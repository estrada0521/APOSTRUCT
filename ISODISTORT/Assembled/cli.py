"""Command-line entrypoint for the local ISODISTORT implementation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from ISODISTORT.Assembled.Backend.pipeline import build_state_from_cif_info
from ISODISTORT.Assembled.Backend.modes.mode_detail_text import render_mode_detail_text
from ISODISTORT.Assembled.Backend.parent import build_parent_state_from_cif_info
from ISODISTORT.Assembled.Backend.reciprocal.k_points import select_k_points
from ISODISTORT.Assembled.case_input import (
    AssembledCase,
    case_from_cli,
    load_case,
    pipeline_request,
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


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="isodistort", description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)

    kpoints = commands.add_parser("kpoints", help="enumerate parent k points")
    kpoints.add_argument("structure", type=Path)
    _add_output_arguments(kpoints)

    irreps = commands.add_parser("irreps", help="enumerate irreps at selected k points")
    irreps.add_argument("structure", type=Path)
    _add_k_arguments(irreps, require_irrep=False)
    _add_output_arguments(irreps)

    opds = commands.add_parser("opds", help="enumerate order-parameter directions")
    opds.add_argument("structure", type=Path)
    _add_k_arguments(opds, require_irrep=True)
    _add_distortion_arguments(opds)
    _add_output_arguments(opds)

    modes = commands.add_parser("modes", help="compute complete mode details")
    modes.add_argument("structure", nargs="?", type=Path)
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
    return root


def _cli_case(args: argparse.Namespace) -> AssembledCase:
    if args.command == "modes" and args.case is not None:
        supplied = any((
            args.structure is not None,
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
    if args.command == "modes" and args.structure is None:
        raise ValueError("modes requires structure or --case")
    if args.command == "modes" and (not args.k or not args.irrep or not args.opd):
        raise ValueError("direct modes requires --k, --irrep, and --opd")
    return case_from_cli(
        args.structure,
        k=getattr(args, "k", ()) or (),
        irreps=getattr(args, "irrep", ()) or (),
        displacive=getattr(args, "displacive", ()) or (),
        magnetic=getattr(args, "magnetic", ()) or (),
        strain=bool(getattr(args, "strain", args.command in {"kpoints", "irreps"})),
        opd=getattr(args, "opd", None),
    )


def execute_case(case: AssembledCase, *, stage: str) -> dict[str, Any]:
    require_irreps = stage in {"opds", "modes"}
    request = pipeline_request(
        case,
        require_irreps=require_irreps,
        default_displacive=stage in {"kpoints", "irreps"},
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
    output.write_text(text)


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
        case = _cli_case(args)
        result = execute_case(case, stage=args.command)
        text = format_result(
            result, output_format=_output_format(args), indent=args.indent
        )
        _write_text(text, output=args.output)
        return 0
    except (KeyError, ValueError, IndexError) as exc:
        parser().error(str(exc))
    return 2


__all__ = ["execute_case", "format_result", "main", "parser"]
