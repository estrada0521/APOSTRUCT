"""Command-line interface for APOSTRUCT."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from APOSTRUCT.Backend.pipeline import (
    build_state_from_cif_info,
    build_state_from_parent_state,
)
from APOSTRUCT.Backend.modes.mode_detail_text import render_mode_detail_text
from APOSTRUCT.Backend.parent import (
    build_generic_parent_state,
    build_parent_state_from_cif_info,
    read_cif_summary,
)
from APOSTRUCT.Backend.reciprocal.k_points import select_k_points
from APOSTRUCT.case_input import (
    APOSTRUCTCase,
    case_from_cli,
    case_from_json,
    load_case,
    pipeline_request,
    pipeline_request_from_cif_info,
)
from APOSTRUCT.cli_output import compact_opd_result, compact_result
from APOSTRUCT.cli_parser import parser


def _exact_vector(text: str, name: str) -> tuple[Fraction, Fraction, Fraction]:
    parts = str(text).strip().strip("()").split(",")
    if len(parts) != 3 or any(not part.strip() for part in parts):
        raise ValueError(f"{name} must be an exact X,Y,Z vector")
    try:
        return tuple(Fraction(part.strip()) for part in parts)  # type: ignore[return-value]
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must contain exact rational values") from exc


def _exact_basis(text: str) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    rows = str(text).strip().split(";")
    if len(rows) != 3:
        raise ValueError("--basis must contain three X,Y,Z rows separated by semicolons")
    return tuple(
        _exact_vector(row, f"--basis row {index}")
        for index, row in enumerate(rows, start=1)
    )


def _direct_parent_state(args: argparse.Namespace) -> dict[str, Any] | None:
    sg = getattr(args, "sg", None)
    wyckoff = getattr(args, "wyckoff", ())
    if sg is None:
        if wyckoff:
            raise ValueError("--wyckoff requires --sg")
        return None
    if getattr(args, "structure", None) is not None:
        raise ValueError("structure and --sg are mutually exclusive")
    site_modes = bool(
        getattr(args, "displacive", ()) or getattr(args, "magnetic", ())
    )
    site_free_strain = (
        args.command in {"irreps", "opds", "modes", "invariants"}
        and bool(getattr(args, "strain", False))
        and not site_modes
    )
    if (
        args.command in {"irreps", "opds", "modes", "invariants"}
        and not wyckoff
        and not site_free_strain
    ):
        raise ValueError(f"{args.command} with --sg requires --wyckoff")
    return build_generic_parent_state(int(sg), list(wyckoff))


def _load_cli_case(path: Path) -> APOSTRUCTCase:
    if str(path) != "-":
        return load_case(path)
    try:
        payload = json.load(sys.stdin)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid case JSON from stdin: {exc}") from exc
    return case_from_json(payload, base=Path.cwd())


def _direction_mode_selection(
    args: argparse.Namespace,
) -> tuple[dict[str, Any], int] | None:
    path = getattr(args, "from_directions", None)
    row_number = getattr(args, "direction_row", None)
    if path is None and row_number is None:
        return None
    if path is None or row_number is None:
        raise ValueError("--from-directions and --direction-row must be used together")
    if args.command != "modes":
        raise ValueError("--from-directions is available only for modes")
    if args.case is not None or args.k or args.irrep or args.opd is not None:
        raise ValueError(
            "--from-directions cannot be combined with --case, --k, --irrep, or --opd"
        )
    try:
        if str(path) == "-":
            payload = json.load(sys.stdin)
        else:
            payload = json.loads(path.expanduser().resolve().read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise ValueError(f"could not load directions JSON: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema") != "APOSTRUCT.cli.directions":
        raise ValueError("--from-directions requires a compact directions payload")
    rows = payload.get("directions")
    if not isinstance(rows, list) or not 1 <= int(row_number) <= len(rows):
        raise ValueError(
            f"direction row must be an integer from 1 through {len(rows) if isinstance(rows, list) else 0}"
        )
    row = rows[int(row_number) - 1]
    if not isinstance(row, dict) or row.get("role") != "primary":
        raise ValueError(f"directions row {row_number} is not primary")
    opd = str(row.get("opd") or "").strip()
    if not opd:
        raise ValueError(f"directions row {row_number} has no Source OPD handle")
    parameters = row.get("k_parameters") or {}
    if not isinstance(parameters, dict):
        raise ValueError(f"directions row {row_number} has invalid K parameters")
    args.k = [[str(row["k_label"]), *[f"{name}={value}" for name, value in parameters.items()]]]
    args.irrep = [str(row["irrep"])]
    args.opd = opd
    return payload, int(row_number)


def _cli_case(args: argparse.Namespace) -> APOSTRUCTCase:
    if args.command in {"modes", "invariants", "run"} and args.case is not None:
        supplied = any((
            getattr(args, "structure", None) is not None,
            getattr(args, "sg", None) is not None,
            bool(getattr(args, "wyckoff", ())),
            bool(getattr(args, "k", ())),
            bool(getattr(args, "irrep", ())),
            bool(getattr(args, "displacive", ())),
            bool(getattr(args, "magnetic", ())),
            bool(getattr(args, "strain", False)),
            getattr(args, "opd", None) is not None,
        ))
        if supplied:
            raise ValueError("--case cannot be combined with direct mode-selection arguments")
        case = _load_cli_case(args.case)
        stage = args.upto if args.command == "run" else args.command
        pure_strain_irreps = (
            stage == "irreps"
            and case.strain
            and not any(displacive or magnetic for _, displacive, magnetic in case.sites)
        )
        if pure_strain_irreps and any(
            slot.label.casefold() != "gm" for slot in case.k
        ):
            raise ValueError("strain-only irreps require k label GM or no k")
        if not case.k and stage != "kpoints" and not pure_strain_irreps:
            raise ValueError(f"{stage} case requires one to four ordered k slots")
        if stage in {"modes", "invariants"} and not case.opd:
            raise ValueError(f"{stage} case requires an OPD")
        if (
            stage == "modes"
            and case.space_group is not None
            and not case.wyckoff
            and not (case.strain and not case.sites)
        ):
            raise ValueError("modes with an sg case requires wyckoff")
        return case
    if args.structure is None and args.sg is None:
        if args.command in {"modes", "invariants"}:
            raise ValueError(f"{args.command} requires structure, --sg, or --case")
        raise ValueError(f"{args.command} requires structure or --sg")
    if args.command in {"modes", "invariants"} and (
        not args.k or not args.irrep or not args.opd
    ):
        raise ValueError(
            f"direct {args.command} requires --k, --irrep, and --opd"
        )
    pure_strain = bool(getattr(args, "strain", False)) and not (
        getattr(args, "displacive", ()) or getattr(args, "magnetic", ())
    )
    if args.command == "irreps" and not args.k and not pure_strain:
        raise ValueError("irreps requires --k unless --strain is the only distortion")
    if args.command == "irreps" and pure_strain and args.k and any(
        str(slot[0]).casefold() != "gm" for slot in args.k if slot
    ):
        raise ValueError("strain-only irreps require --k GM or no --k")
    if args.command in {"irreps", "opds", "modes", "invariants"} and not any((
        getattr(args, "displacive", ()),
        getattr(args, "magnetic", ()),
        getattr(args, "strain", False),
    )):
        raise ValueError(
            "select at least one distortion with --displacive SITE, "
            "--magnetic SITE, or --strain; SITE is an info.sites[].label or type value"
        )
    return case_from_cli(
        args.structure if args.structure is not None else Path(f"SG{args.sg}"),
        k=getattr(args, "k", ()) or (),
        irreps=getattr(args, "irrep", ()) or (),
        displacive=getattr(args, "displacive", ()) or (),
        magnetic=getattr(args, "magnetic", ()) or (),
        strain=bool(getattr(args, "strain", False)),
        opd=getattr(args, "opd", None),
    )


def execute_case(
    case: APOSTRUCTCase,
    *,
    stage: str,
    cif_info: dict[str, Any] | None = None,
    parent_state: dict[str, Any] | None = None,
    direction_selection: tuple[dict[str, Any], int] | None = None,
) -> dict[str, Any]:
    require_irreps = stage in {"opds", "modes", "invariants"}
    request_kwargs = {
        "require_irreps": require_irreps,
        "default_displacive": False,
        "allow_empty_distortion": stage == "kpoints",
    }
    if parent_state is not None:
        cif_info = parent_state["input"]
    request = (
        pipeline_request(case, **request_kwargs)
        if cif_info is None
        else pipeline_request_from_cif_info(case, cif_info, **request_kwargs)
    )
    cif_info = request.pop("cif_info")
    parent_state = request.pop("parent_state", parent_state)
    if stage == "kpoints":
        return select_k_points(
            parent_state or build_parent_state_from_cif_info(cif_info),
            k_label=request["k_label"],
            k_index=request["k_index"],
            k_params=request["k_params"],
            display_k_params=request["display_k_params"],
            selections=request["selections"],
        )
    builder = (
        build_state_from_parent_state
        if parent_state is not None
        else build_state_from_cif_info
    )
    return builder(
        parent_state if parent_state is not None else cif_info,
        **request,
        include_opd=stage in {"opds", "modes", "invariants"},
        include_mode_details=stage == "modes",
        selected_opd_only=stage in {"modes", "invariants"},
        **(
            {"direction_selection": direction_selection}
            if direction_selection is not None
            else {}
        ),
    )


def _invariant_result(
    state: dict[str, Any],
    *,
    minimum_degree: int,
    maximum_degree: int,
) -> dict[str, Any]:
    selected = state.get("selected") or {}
    source_factors = selected.get("opd_factors") or []
    if not source_factors:
        raise ValueError("selected OPD has no primary invariant factors")
    orderparam = selected.get("orderparam") or {}
    selection_label = str((orderparam.get("direction") or {}).get("label") or "")
    if not selection_label:
        raise ValueError("selected OPD lost its public label")
    return _invariant_result_from_factors(
        space_group=int(state["space_group"]["number"]),
        source_factors=source_factors,
        selection_label=selection_label,
        distortion_selection=dict(state.get("distortion_selection") or {}),
        opd=compact_opd_result(dict(orderparam)),
        minimum_degree=minimum_degree,
        maximum_degree=maximum_degree,
    )


def _invariant_result_from_factors(
    *,
    space_group: int,
    source_factors: Sequence[dict[str, Any]],
    selection_label: str,
    distortion_selection: dict[str, Any],
    opd: dict[str, Any],
    minimum_degree: int,
    maximum_degree: int,
) -> dict[str, Any]:
    from APOSTRUCT import _invariant_service, handle_api_post

    factors = [
        {
            "gid": factor["gid"],
            "magnetic": bool(factor.get("magnetic")),
            "k_parameters": list(
                factor.get("source_k_parameters")
                or factor.get("k_parameters")
                or []
            ),
            "opd": factor["opd"],
            "domain": factor["domain"],
            **(
                {"parameter_offset": factor["parameter_offset"]}
                if "parameter_offset" in factor
                else {}
            ),
            **(
                {"direction_matrix": factor["direction_matrix"]}
                if "direction_matrix" in factor
                else {}
            ),
        }
        for factor in source_factors
    ]
    request = {
            "space_group": space_group,
            "factors": factors,
            "minimum_degree": minimum_degree,
            "maximum_degree": maximum_degree,
        }
    if any("direction_matrix" in factor for factor in factors):
        result = _invariant_service().compute({
            **request,
            "factors": [
                {
                    **{
                        key: value
                        for key, value in factor.items()
                        if key not in {"opd", "direction_matrix"}
                    },
                    "direction": factor["opd"],
                    "selected_direction_matrix": factor.get("direction_matrix"),
                }
                for factor in factors
            ],
        })
    else:
        result = handle_api_post("invariants", {}, request)
    result.pop("elapsed_ms", None)
    degrees_by_number = {
        int(row["degree"]): row for row in result.get("degrees") or []
    }
    result["degrees"] = [
        degrees_by_number.get(
            degree,
            {"degree": degree, "count": 0, "invariants": [], "polynomials": []},
        )
        for degree in range(minimum_degree, maximum_degree + 1)
    ]
    public_factors = result.get("factors") or []
    if len(public_factors) != len(source_factors):
        raise ValueError("invariant result lost a selected factor")
    for source, public in zip(source_factors, public_factors, strict=True):
        public["slot"] = int(source["slot"])
        public["role"] = str(source.get("role") or "primary")
        public["k_parameters"] = list(source.get("k_parameters") or [])
    return {
        "schema": "APOSTRUCT.cli.invariants",
        "selection": {"opd": selection_label},
        "distortion_selection": distortion_selection,
        "opd": opd,
        **result,
    }


def _state_with_invariant_factors(state: dict[str, Any]) -> dict[str, Any]:
    """Add GUI-equivalent invariant factor records without mutating mode state."""

    from APOSTRUCT import _complete_invariant_factors

    selected = state.get("selected") or {}
    mode_details = selected.get("mode_details") or {}
    if mode_details.get("status") != "ok":
        return state
    factor_state = {
        "space_group": state.get("space_group"),
        "input": state.get("input"),
        "selected": {
            "opd_factors": deepcopy(selected.get("opd_factors") or []),
            "mode_details": {
                "secondary_opd_factors": deepcopy(
                    mode_details.get("secondary_opd_factors") or []
                )
            },
        },
    }
    _complete_invariant_factors(factor_state)
    completed = factor_state["selected"]
    factors = [
        {**factor, "role": "primary"}
        for factor in completed.get("opd_factors") or []
    ]
    factors.extend(completed["mode_details"].get("secondary_opd_factors") or [])
    public_state = dict(state)
    public_state["selected"] = {**selected, "invariant_factors": factors}
    return public_state


def _invariant_result_from_modes(
    payload: dict[str, Any],
    *,
    secondary_slots: Sequence[int],
    minimum_degree: int,
    maximum_degree: int,
) -> dict[str, Any]:
    factors = payload.get("invariant_factors")
    if not isinstance(factors, list) or not factors:
        raise ValueError("compact modes payload has no invariant factor records")
    requested = list(secondary_slots)
    if len(set(requested)) != len(requested):
        raise ValueError("secondary factor slots must not be repeated")
    requested_set = set(requested)
    primary = [factor for factor in factors if factor.get("role") == "primary"]
    secondary = [factor for factor in factors if factor.get("role") == "secondary"]
    available = {int(factor["slot"]): factor for factor in secondary}
    missing = [slot for slot in requested if slot not in available]
    if missing:
        labels = ", ".join(str(slot) for slot in available)
        raise ValueError(
            f"secondary factor slot {missing[0]} not found; available slots: [{labels}]"
        )
    selected = primary + [
        factor for factor in secondary if int(factor["slot"]) in requested_set
    ]
    if not primary:
        raise ValueError("compact modes payload has no primary invariant factors")
    space_group = payload.get("space_group") or {}
    opd = payload.get("opd") or {}
    selection_label = str(opd.get("label") or "")
    if not selection_label:
        raise ValueError("compact modes payload has no selected OPD label")
    return _invariant_result_from_factors(
        space_group=int(space_group["number"]),
        source_factors=selected,
        selection_label=selection_label,
        distortion_selection=dict(payload.get("distortion_selection") or {}),
        opd=dict(opd),
        minimum_degree=minimum_degree,
        maximum_degree=maximum_degree,
    )


def format_result(value: dict[str, Any], *, output_format: str, indent: int) -> str:
    if output_format == "text":
        return render_mode_detail_text(value).rstrip() + "\n"
    if output_format != "json":
        raise ValueError(f"unknown output format: {output_format}")
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, indent=indent, allow_nan=False
    ) + "\n"


def _write_text(text: str, *, output: Path | None) -> None:
    if output is None:
        sys.stdout.write(text)
        return
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")


def _output_format(args: argparse.Namespace, *, stage: str) -> str:
    if getattr(args, "format", None):
        if args.format == "text" and stage != "modes":
            raise ValueError("text output is available only for modes")
        return args.format
    return "json"


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "serve":
            from APOSTRUCT.server import serve

            serve(host=args.host, port=args.port, open_browser=args.open_browser)
            return 0
        if args.command == "show":
            from APOSTRUCT.gui_launch import (
                gui_launch_from_case,
                send_gui_launch,
            )

            result = send_gui_launch(
                gui_launch_from_case(
                    _load_cli_case(args.case),
                    raw_amplitudes=args.amplitude,
                ),
                server=args.server,
            )
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        if args.command == "settings":
            from APOSTRUCT.Backend.settings import international_settings

            _write_text(
                format_result(
                    international_settings(args.sg),
                    output_format="json",
                    indent=args.indent,
                ),
                output=args.output,
            )
            return 0
        if args.command == "directions":
            from APOSTRUCT.Backend.directions import (
                compatible_directions,
                compatible_magnetic_directions,
            )

            basis = _exact_basis(args.basis)
            origin = _exact_vector(args.origin, "--origin")
            result = (
                compatible_magnetic_directions(
                    args.parent_sg,
                    args.subgroup_msg,
                    basis,
                    origin,
                    parent_setting_id=args.parent_setting,
                    subgroup_setting_id=args.subgroup_setting,
                )
                if args.subgroup_msg is not None
                else compatible_directions(
                    args.parent_sg,
                    args.subgroup_sg,
                    basis,
                    origin,
                    parent_setting_id=args.parent_setting,
                    subgroup_setting_id=args.subgroup_setting,
                )
            )
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        if args.command == "combine-modes":
            from APOSTRUCT.cli_mode_combination import (
                combine_mode_payload,
                load_compact_modes,
            )

            result = combine_mode_payload(
                load_compact_modes(args.modes_json),
                raw_weights=args.weight or args.amplitude,
                apply_normfactor=args.amplitude is not None,
            )
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        if args.command == "invariants" and args.from_modes is not None:
            supplied = any((
                args.structure is not None,
                args.sg is not None,
                bool(args.wyckoff),
                args.case is not None,
                bool(args.k),
                bool(args.irrep),
                bool(args.displacive),
                bool(args.magnetic),
                bool(args.strain),
                args.opd is not None,
            ))
            if supplied:
                raise ValueError(
                    "--from-modes cannot be combined with direct mode-selection arguments"
                )
            from APOSTRUCT.cli_mode_combination import load_compact_modes

            result = _invariant_result_from_modes(
                load_compact_modes(args.from_modes),
                secondary_slots=args.secondary,
                minimum_degree=args.minimum_degree,
                maximum_degree=args.maximum_degree,
            )
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        if args.command == "invariants" and args.secondary:
            raise ValueError("--secondary requires --from-modes")
        stage = args.upto if args.command == "run" else args.command
        if stage == "invariants" and getattr(args, "full_state", False):
            raise ValueError("invariants output is compact; --full-state is unavailable")
        if stage != "invariants" and args.command == "run" and (
            args.minimum_degree is not None or args.maximum_degree is not None
        ):
            raise ValueError("degree arguments require --upto invariants")
        direction_selection = _direction_mode_selection(args)
        parent_state = _direct_parent_state(args)
        if args.command == "info":
            if args.structure is None and parent_state is None:
                raise ValueError("info requires structure or --sg")
            result = (
                parent_state
                if parent_state is not None
                else build_parent_state_from_cif_info(
                    read_cif_summary(args.structure.expanduser().resolve())
                )
            )
            if not args.full_state:
                result = compact_result(result, stage="info")
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        case = _cli_case(args)
        result = execute_case(
            case,
            stage=stage,
            parent_state=parent_state,
            direction_selection=direction_selection,
        )
        output_format = _output_format(args, stage=stage)
        if stage == "invariants":
            result = _invariant_result(
                result,
                minimum_degree=(
                    2 if args.minimum_degree is None else args.minimum_degree
                ),
                maximum_degree=(
                    4 if args.maximum_degree is None else args.maximum_degree
                ),
            )
        elif output_format == "json" and not getattr(args, "full_state", False):
            if stage == "modes":
                result = _state_with_invariant_factors(result)
            result = compact_result(result, stage=stage)
        text = format_result(
            result, output_format=output_format, indent=args.indent
        )
        _write_text(text, output=args.output)
        if stage == "modes":
            details = result.get("mode_details")
            if not isinstance(details, dict):
                details = (result.get("selected") or {}).get("mode_details")
            if isinstance(details, dict) and details.get("status") == "error":
                return 2
        return 0
    except (KeyError, ValueError, IndexError, OSError) as exc:
        args._error_parser.error(str(exc))
    return 2


__all__ = ["compact_result", "execute_case", "format_result", "main", "parser"]
