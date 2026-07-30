"""Command-line entrypoint for the local ISODISTORT implementation."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import sys
from typing import Any, Sequence

from ISODISTORT.Assembled.Backend.pipeline import build_state_from_cif_info
from ISODISTORT.Assembled.Backend.modes.mode_detail_text import render_mode_detail_text
from ISODISTORT.Assembled.Backend.parent import (
    build_parent_state_from_cif_info,
    read_cif_summary,
)
from ISODISTORT.Assembled.Backend.reciprocal.k_points import select_k_points
from ISODISTORT.Assembled.Backend.source.tables import source_tables
from ISODISTORT.Assembled.case_input import (
    AssembledCase,
    case_from_cli,
    load_case,
    pipeline_request,
    pipeline_request_from_cif_info,
)


def _add_output_arguments(
    parser: argparse.ArgumentParser, *, selectable_format: bool = False
) -> None:
    parser.add_argument("--output", "-o", type=Path)
    if selectable_format:
        parser.add_argument("--format", choices=("json", "text"))
    parser.add_argument("--indent", type=int, default=2)


def _add_k_arguments(parser: argparse.ArgumentParser, *, require_irrep: bool) -> None:
    parser.add_argument(
        "--k",
        action="append",
        nargs="+",
        required=True,
        metavar=("LABEL", "PARAM"),
        help="k label followed by its exact positional parameters",
    )
    parser.add_argument(
        "--irrep",
        action="append",
        default=[],
        required=require_irrep,
        help="irrep paired by order with each --k",
    )


def _add_distortion_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--displacive", nargs="+", default=[], metavar="SITE")
    parser.add_argument("--magnetic", nargs="+", default=[], metavar="SITE")
    parser.add_argument("--strain", action="store_true")


def _add_structure_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", nargs="?", type=Path)
    parser.add_argument("--sg", type=int, help="parent space-group number")
    parser.add_argument(
        "--wyckoff",
        nargs="+",
        default=[],
        metavar="POSITION",
        help="Source Wyckoff positions used to construct a generic parent structure",
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="isodistort", description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    info = commands.add_parser("info", help="show the normalized parent structure")
    _add_structure_arguments(info)
    _add_output_arguments(info)

    kpoints = commands.add_parser("kpoints", help="enumerate parent k points")
    _add_structure_arguments(kpoints)
    _add_output_arguments(kpoints)

    irreps = commands.add_parser("irreps", help="enumerate irreps at selected k points")
    _add_structure_arguments(irreps)
    _add_k_arguments(irreps, require_irrep=False)
    _add_output_arguments(irreps)

    opds = commands.add_parser("opds", help="enumerate order-parameter directions")
    _add_structure_arguments(opds)
    _add_k_arguments(opds, require_irrep=True)
    _add_distortion_arguments(opds)
    _add_output_arguments(opds)

    modes = commands.add_parser("modes", help="compute complete mode details")
    _add_structure_arguments(modes)
    modes.add_argument("--case", type=Path)
    modes.add_argument("--k", action="append", nargs="+", default=[])
    modes.add_argument("--irrep", action="append", default=[])
    _add_distortion_arguments(modes)
    modes.add_argument("--opd")
    _add_output_arguments(modes, selectable_format=True)

    serve = commands.add_parser("serve", help="serve the standalone frontend")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8300)
    serve.add_argument("--open-browser", action="store_true")
    for command in (info, kpoints, irreps, opds, modes, serve):
        command.set_defaults(_error_parser=command)
    return root


def _generic_cell(sg: int) -> dict[str, str]:
    if sg <= 2:
        values = ("4", "5", "6", "70", "80", "75")
    elif sg <= 15:
        values = ("4", "5", "6", "90", "105", "90")
    elif sg <= 74:
        values = ("4", "5", "6", "90", "90", "90")
    elif sg <= 142:
        values = ("4", "4", "6", "90", "90", "90")
    elif sg <= 194:
        values = ("4", "4", "6", "90", "90", "120")
    else:
        values = ("4", "4", "4", "90", "90", "90")
    return dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), values, strict=True))


def _virtual_cif_info(sg: int, wyckoff: Sequence[str]) -> dict[str, Any]:
    if not 1 <= int(sg) <= 230:
        raise ValueError("--sg must be an integer from 1 through 230")
    data = source_tables()
    parameter_values = (Fraction(1, 7), Fraction(2, 11), Fraction(3, 13))
    sites: list[dict[str, Any]] = []
    used_labels: set[str] = set()
    for ordinal, raw_position in enumerate(wyckoff, 1):
        match = re.fullmatch(r"(?:(\d+))?([A-Za-z]+)", str(raw_position).strip())
        if match is None:
            raise ValueError(f"invalid Wyckoff position {raw_position!r}")
        supplied_multiplicity, raw_label = match.groups()
        label = raw_label.casefold()
        if label in used_labels:
            raise ValueError(f"duplicate Wyckoff position {raw_position!r}")
        used_labels.add(label)
        try:
            row = data.wyckoff_row_by_label(int(sg), label)
        except KeyError as exc:
            raise ValueError(str(exc)) from exc
        multiplicity = data.wyckoff_multiplicity(int(sg), row)
        if supplied_multiplicity is not None and int(supplied_multiplicity) != multiplicity:
            raise ValueError(
                f"Wyckoff position {raw_position!r} has multiplicity {multiplicity} in SG{sg}"
            )
        vectors = data.wyckoff_fraction_vectors(row)
        point = tuple(
            (
                vectors[0][axis]
                + sum(
                    parameter_values[index] * vectors[index + 1][axis]
                    for index in range(3)
                )
            )
            % 1
            for axis in range(3)
        )
        site = {
            "type": label,
            "label": f"{label}{ordinal}",
            "multiplicity": str(multiplicity),
            "fract": [str(float(value)) for value in point],
            "occupancy": "1",
        }
        identified = data.match_wyckoff_site(int(sg), site)
        if identified is None or int(identified["row_id"]) != int(row.row_id):
            raise ValueError(
                f"could not construct a generic representative for SG{sg} Wyckoff {label}"
            )
        sites.append(site)
    return {
        "path": None,
        "block": f"space_group_{sg}",
        "formula": None,
        "parent": {"number": int(sg), "symbol": data.space_symbol(int(sg))},
        "explicit_parent_preference": None,
        "lattice": _generic_cell(int(sg)),
        "atom_sites": sites,
        "symmetry_operations": [],
    }


def _direct_cif_info(args: argparse.Namespace) -> dict[str, Any] | None:
    if args.sg is None:
        if args.wyckoff:
            raise ValueError("--wyckoff requires --sg")
        return None
    if args.structure is not None:
        raise ValueError("structure and --sg are mutually exclusive")
    if args.command in {"irreps", "opds", "modes"} and not args.wyckoff:
        raise ValueError(f"{args.command} with --sg requires --wyckoff")
    return _virtual_cif_info(int(args.sg), args.wyckoff)


def _cli_case(args: argparse.Namespace) -> AssembledCase:
    if args.command == "modes" and args.case is not None:
        supplied = any((
            args.structure is not None,
            args.sg is not None,
            bool(args.wyckoff),
            bool(args.k),
            bool(args.irrep),
            bool(args.displacive),
            bool(args.magnetic),
            bool(args.strain),
            args.opd is not None,
        ))
        if supplied:
            raise ValueError("--case cannot be combined with direct mode-selection arguments")
        case = load_case(args.case)
        if not case.opd:
            raise ValueError("modes case requires an OPD")
        return case
    if args.structure is None and args.sg is None:
        if args.command == "modes":
            raise ValueError("modes requires structure, --sg, or --case")
        raise ValueError(f"{args.command} requires structure or --sg")
    if args.command == "modes" and (not args.k or not args.irrep or not args.opd):
        raise ValueError("direct modes requires --k, --irrep, and --opd")
    return case_from_cli(
        args.structure if args.structure is not None else Path(f"SG{args.sg}"),
        k=getattr(args, "k", ()) or (),
        irreps=getattr(args, "irrep", ()) or (),
        displacive=getattr(args, "displacive", ()) or (),
        magnetic=getattr(args, "magnetic", ()) or (),
        strain=bool(getattr(args, "strain", args.command in {"kpoints", "irreps"})),
        opd=getattr(args, "opd", None),
    )


def execute_case(
    case: AssembledCase,
    *,
    stage: str,
    cif_info: dict[str, Any] | None = None,
) -> dict[str, Any]:
    require_irreps = stage in {"opds", "modes"}
    request_kwargs = {
        "require_irreps": require_irreps,
        "default_displacive": stage in {"kpoints", "irreps"},
        "allow_empty_distortion": stage == "kpoints",
    }
    request = (
        pipeline_request(case, **request_kwargs)
        if cif_info is None
        else pipeline_request_from_cif_info(case, cif_info, **request_kwargs)
    )
    cif_info = request.pop("cif_info")
    if stage == "kpoints":
        return select_k_points(
            build_parent_state_from_cif_info(cif_info),
            k_label=request["k_label"],
            k_index=request["k_index"],
            k_params=request["k_params"],
            display_k_params=request["display_k_params"],
            selections=request["selections"],
        )
    return build_state_from_cif_info(
        cif_info,
        **request,
        include_opd=stage in {"opds", "modes"},
        include_mode_details=stage == "modes",
        selected_opd_only=stage == "modes",
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


def _output_format(args: argparse.Namespace) -> str:
    if getattr(args, "format", None):
        if args.format == "text" and args.command != "modes":
            raise ValueError("text output is available only for modes")
        return args.format
    if args.command != "modes":
        return "json"
    if args.case is not None and args.case.suffix.casefold() == ".json":
        return "json"
    return "text"


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "serve":
            from ISODISTORT.Assembled.server import serve

            serve(host=args.host, port=args.port, open_browser=args.open_browser)
            return 0
        cif_info = _direct_cif_info(args)
        if args.command == "info":
            if args.structure is None and cif_info is None:
                raise ValueError("info requires structure or --sg")
            source = (
                cif_info
                if cif_info is not None
                else read_cif_summary(args.structure.expanduser().resolve())
            )
            result = build_parent_state_from_cif_info(source)
            _write_text(
                format_result(result, output_format="json", indent=args.indent),
                output=args.output,
            )
            return 0
        case = _cli_case(args)
        result = execute_case(case, stage=args.command, cif_info=cif_info)
        text = format_result(
            result, output_format=_output_format(args), indent=args.indent
        )
        _write_text(text, output=args.output)
        return 0
    except (KeyError, ValueError, IndexError, OSError) as exc:
        args._error_parser.error(str(exc))
    return 2


__all__ = ["execute_case", "format_result", "main", "parser"]
