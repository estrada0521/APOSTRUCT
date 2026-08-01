"""Argument parser for the Distortropy command line."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any


def _add_output_arguments(
    parser: argparse.ArgumentParser,
    *,
    selectable_format: bool = False,
    full_state: bool = False,
) -> None:
    parser.add_argument("--output", "-o", type=Path, help="write output to PATH")
    if selectable_format:
        parser.add_argument(
            "--format",
            choices=("json", "text"),
            help="output format; JSON is the default",
        )
    if full_state:
        parser.add_argument(
            "--full-state",
            action="store_true",
            help="emit the complete internal pipeline state instead of the compact result",
        )
    parser.add_argument("--indent", type=int, default=2, help="JSON indentation")


def _add_k_arguments(
    parser: argparse.ArgumentParser,
    *,
    require_irrep: bool,
    required: bool = True,
) -> None:
    parser.add_argument(
        "--k",
        action="append",
        nargs="+",
        required=required,
        metavar=("LABEL", "PARAM"),
        help=(
            "repeatable k-point selection; exact positional values follow "
            "kpoints[].parameter_names order in kpoints[].kvector, or use "
            "NAME=VALUE entries"
        ),
    )
    if require_irrep:
        parser.add_argument(
            "--irrep",
            action="append",
            required=True,
            help="repeatable irreps row symbol paired by order with each --k",
        )


def _add_distortion_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--displacive",
        nargs="+",
        default=[],
        metavar="SITE",
        help=(
            "enable displacement modes for info.sites[].label or type values; "
            "a type selects every matching site"
        ),
    )
    parser.add_argument(
        "--magnetic",
        nargs="+",
        default=[],
        metavar="SITE",
        help=(
            "enable magnetic modes for info.sites[].label or type values; "
            "a type selects every matching site"
        ),
    )
    parser.add_argument("--strain", action="store_true", help="include strain modes")


def _add_structure_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", nargs="?", type=Path, help="input CIF path")
    parser.add_argument("--sg", type=int, help="parent space-group number")
    parser.add_argument(
        "--wyckoff",
        nargs="+",
        default=[],
        metavar="POSITION",
        help=(
            "ordered symbolic Wyckoff sites used with --sg; letters or "
            "multiplicity+letter forms are accepted, and mode details for a "
            "free site require values such as i:x=1/7,y=2/11,z=3/13"
        ),
    )


def _add_selected_opd_arguments(parser: argparse.ArgumentParser) -> None:
    _add_structure_arguments(parser)
    parser.add_argument("--case", type=Path, help="JSON case path or - for stdin")
    parser.add_argument(
        "--k",
        action="append",
        nargs="+",
        default=[],
        metavar=("LABEL", "PARAM"),
        help=(
            "repeatable k-point selection; exact positional values follow "
            "kpoints[].parameter_names order in kpoints[].kvector, or use "
            "NAME=VALUE entries"
        ),
    )
    parser.add_argument(
        "--irrep",
        action="append",
        default=[],
        help="repeatable irreps row symbol paired by order with each --k",
    )
    _add_distortion_arguments(parser)
    parser.add_argument(
        "--opd",
        help="exact opds[].label to compute; unknown labels are rejected",
    )


def _add_invariant_degree_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_minimum: int | None = 2,
    default_maximum: int | None = 4,
) -> None:
    parser.add_argument(
        "--minimum-degree",
        type=int,
        default=default_minimum,
        help="lowest polynomial degree, from 1 through 12",
    )
    parser.add_argument(
        "--maximum-degree",
        type=int,
        default=default_maximum,
        help="highest polynomial degree, from minimum through 12",
    )


def _command_parser(
    commands: Any,
    name: str,
    *,
    help_text: str,
    description: str,
    example: str,
) -> argparse.ArgumentParser:
    return commands.add_parser(
        name,
        help=help_text,
        description=description,
        epilog=f"example:\n  {example}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="distortropy",
        description=(
            "Compute symmetry-adapted distortion modes from a parent CIF or "
            "a symbolic space-group/Wyckoff parent."
        ),
        epilog="""workflow:
  distortropy settings --sg 14
  distortropy directions --parent-sg 221 --subgroup-sg 140 --basis=... (three rows)
  distortropy info structure.cif
  distortropy kpoints structure.cif
  distortropy irreps structure.cif --k R --displacive O
  distortropy opds structure.cif --k R --irrep R4- --displacive O
  distortropy modes structure.cif --k R --irrep R4- --displacive O --opd P1
  distortropy invariants structure.cif --k R --irrep R4- --displacive O --opd P1
  distortropy invariants --from-modes modes.json --secondary 2
  distortropy combine-modes modes.json --weight displacive-1=1

Use each compact JSON result as the vocabulary for the next command.
Use --full-state only when complete pipeline provenance is required.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    commands = root.add_subparsers(dest="command", required=True)

    settings = _command_parser(
        commands,
        "settings",
        help_text="list International settings for a space group",
        description=(
            "List the International setting IDs accepted by directions.\n"
            "The choices match the setting-level parent and subgroup catalogs used by COPL."
        ),
        example="distortropy settings --sg 14",
    )
    settings.add_argument("--sg", type=int, required=True, help="space-group number")
    _add_output_arguments(settings)

    directions = _command_parser(
        commands,
        "directions",
        help_text="find OP directions compatible with a subgroup embedding",
        description=(
            "Find ordinary or time-odd commensurate order-parameter directions "
            "allowed by an exact parent-to-subgroup embedding.\n"
            "Basis and origin use the selected International settings of both groups; "
            "defaults are listed by settings. A BNS number selects a magnetic subgroup. "
            "No CIF or Wyckoff sites are involved."
        ),
        example=(
            "distortropy directions --parent-sg 221 --subgroup-sg 140 "
            "--basis=-1,1,0 --basis=-1,-1,0 --basis=0,0,2"
        ),
    )
    directions.add_argument("--parent-sg", type=int, required=True)
    directions.add_argument(
        "--parent-setting",
        type=int,
        help="International setting ID from distortropy settings --sg PARENT_SG",
    )
    subgroup = directions.add_mutually_exclusive_group(required=True)
    subgroup.add_argument("--subgroup-sg", type=int)
    subgroup.add_argument(
        "--subgroup-msg",
        metavar="BNS",
        help="BNS magnetic group number, for example 62.448",
    )
    directions.add_argument(
        "--subgroup-setting",
        type=int,
        help="International setting ID for the ordinary subgroup number",
    )
    directions.add_argument(
        "--basis",
        action="append",
        required=True,
        metavar="X,Y,Z",
        help="repeat exactly three times for the subgroup basis in parent coordinates",
    )
    directions.add_argument(
        "--origin",
        default="0,0,0",
        metavar="X,Y,Z",
        help="subgroup origin in parent coordinates; default: 0,0,0",
    )
    _add_output_arguments(directions)

    info = _command_parser(
        commands,
        "info",
        help_text="show the normalized parent structure",
        description=(
            "Show the space group, cell, formula, and selectable site types.\n"
            "With --sg, Wyckoff labels and coordinates use the default International "
            "Tables setting; free coordinates are arbitrary generic placeholders."
        ),
        example="distortropy info --sg 221 --wyckoff 1a 1b 3c",
    )
    _add_structure_arguments(info)
    _add_output_arguments(info, full_state=True)

    kpoints = _command_parser(
        commands,
        "kpoints",
        help_text="enumerate parent k points",
        description="List k-point labels and the exact parameters required to select them.",
        example="distortropy kpoints structure.cif",
    )
    _add_structure_arguments(kpoints)
    _add_output_arguments(kpoints, full_state=True)

    irreps = _command_parser(
        commands,
        "irreps",
        help_text="enumerate irreps at selected k points",
        description=(
            "List selectable irrep symbols for one to four ordered k points.\n"
            "For homogeneous strain alone, --k and --wyckoff may be omitted; "
            "the k point is GM."
        ),
        example="distortropy irreps structure.cif --k DT b=1/3 --displacive O",
    )
    _add_structure_arguments(irreps)
    _add_k_arguments(irreps, require_irrep=False, required=False)
    _add_distortion_arguments(irreps)
    _add_output_arguments(irreps, full_state=True)

    opds = _command_parser(
        commands,
        "opds",
        help_text="enumerate order-parameter directions",
        description=(
            "List OPD labels, subgroups, bases, origins, and active k vectors.\n"
            "Homogeneous strain selections may use --sg without --wyckoff.\n"
            "In JSON, index is group index i and cell_index is supercell factor s.\n"
            "ferroic_properties uses proper/improper ferroelectric, "
            "proper/improper ferroelastic, and ferromagnetic; "
            "ferroic_classified distinguishes a classified empty result."
        ),
        example=(
            "distortropy opds structure.cif --k R --irrep R4- --displacive O"
        ),
    )
    _add_structure_arguments(opds)
    _add_k_arguments(opds, require_irrep=True)
    _add_distortion_arguments(opds)
    _add_output_arguments(opds, full_state=True)

    modes = _command_parser(
        commands,
        "modes",
        help_text="compute complete mode details",
        description=(
            "Compute mode definitions for one exact label returned by opds[].label.\n"
            "Shell-quote coupled labels such as 'P3(1)P3(1)'.\n"
            "In JSON definitions, mode identifies the kind, K vector, irrep, site, "
            "and site irrep. Role is primary or secondary; an unambiguous primary "
            "includes factor with its slot, gid, and global parameter names."
        ),
        example=(
            "distortropy modes structure.cif --k R --irrep R4- "
            "--displacive O --opd P1"
        ),
    )
    _add_selected_opd_arguments(modes)
    modes.add_argument(
        "--from-directions",
        type=Path,
        help=(
            "directions JSON whose exact embedding and primary direction are "
            "reused instead of a catalog OPD representative"
        ),
    )
    modes.add_argument(
        "--direction-row",
        type=int,
        help="1-based directions[] row in --from-directions input",
    )
    _add_output_arguments(modes, selectable_format=True, full_state=True)

    invariants = _command_parser(
        commands,
        "invariants",
        help_text="compute invariants of a selected OPD domain",
        description=(
            "Compute primary-order-parameter invariants for one exact OPD selection.\n"
            "Alternatively, reuse compact modes JSON and include selected secondary "
            "factor slots without rerunning modes.\n"
            "Homogeneous strain selections may use --sg without --wyckoff. "
            "Shell-quote coupled labels such as 'P3(1)P3(1)'.\n"
            "In JSON polynomial terms, coefficient is an exact expression string and "
            "exponents follow variables order. ferroic_classified distinguishes a "
            "classified empty ferroic_properties result."
        ),
        example=(
            "distortropy invariants structure.cif --k R --irrep R4- "
            "--displacive O --opd P2 --minimum-degree 2 --maximum-degree 4"
        ),
    )
    _add_selected_opd_arguments(invariants)
    invariants.add_argument(
        "--from-modes",
        type=Path,
        help="compact modes JSON path or - for stdin; primary factors are included",
    )
    invariants.add_argument(
        "--secondary",
        action="append",
        type=int,
        default=[],
        metavar="SLOT",
        help="repeatable secondary factor slot from --from-modes input",
    )
    _add_invariant_degree_arguments(invariants)
    _add_output_arguments(invariants)

    combine_modes = _command_parser(
        commands,
        "combine-modes",
        help_text="combine definitions from compact modes JSON",
        description=(
            "Apply explicit weights to displacive or magnetic definition rows from one "
            "compact modes payload.\n"
            "Definition IDs are local to that payload. --weight multiplies published "
            "rows directly; --amplitude applies each definition's published "
            "normfactor and does not establish a common physical unit across unrelated "
            "definitions. "
            "Returned vectors use the child "
            "crystallographic dxyz component basis; net_magnetic_vector is their sum "
            "over the returned atoms."
        ),
        example=(
            "distortropy combine-modes modes.json "
            "--weight magnetic-1=1 --weight magnetic-2=3/4"
        ),
    )
    combine_modes.add_argument(
        "modes_json", type=Path, help="compact modes JSON path or - for stdin"
    )
    coefficients = combine_modes.add_mutually_exclusive_group(required=True)
    coefficients.add_argument(
        "--weight",
        action="append",
        metavar="DEFINITION_ID=VALUE",
        help="repeatable coefficient applied directly to a published definition row",
    )
    coefficients.add_argument(
        "--amplitude",
        action="append",
        metavar="DEFINITION_ID=VALUE",
        help="repeatable coefficient in each definition's published normfactor convention",
    )
    _add_output_arguments(combine_modes)

    run = _command_parser(
        commands,
        "run",
        help_text="execute a JSON case through one pipeline stage",
        description=(
            "Read a JSON case from a path or stdin and stop after the requested stage.\n"
            "Use exactly one structure source: structure, cif, or sg with optional "
            "wyckoff. Wyckoff items may be labels or objects with wyckoff and "
            "parameters; free parameters may remain absent through invariants but "
            "are required for modes. Relative structure paths resolve from the case "
            "file directory.\n"
            "Other case fields are sites, strain, ordered k, and opd. k may be omitted "
            "for kpoints and strain-only irreps. Each k item uses label, optional params, "
            "and ir for opds or later. Degree arguments require --upto invariants, where "
            "--full-state is unavailable."
        ),
        example=(
            "printf '%s' '{\"sg\":123,\"wyckoff\":[\"1a\"],"
            "\"sites\":{\"a\":{\"displacive\":true,\"magnetic\":false}},"
            "\"k\":[{\"label\":\"M\"}]}' | "
            "distortropy run --case - --upto irreps"
        ),
    )
    run.add_argument(
        "--case", type=Path, required=True, help="JSON case path or - for stdin"
    )
    run.add_argument(
        "--upto",
        choices=("kpoints", "irreps", "opds", "modes", "invariants"),
        default="modes",
        help="last pipeline stage to execute",
    )
    _add_invariant_degree_arguments(
        run,
        default_minimum=None,
        default_maximum=None,
    )
    _add_output_arguments(run, selectable_format=True, full_state=True)

    serve = _command_parser(
        commands,
        "serve",
        help_text="serve the standalone frontend",
        description="Serve the local graphical interface.",
        example="distortropy serve --host 127.0.0.1 --port 8300",
    )
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8300)
    serve.add_argument("--open-browser", action="store_true")
    for command in (
        settings,
        directions,
        info,
        kpoints,
        irreps,
        opds,
        modes,
        invariants,
        combine_modes,
        run,
        serve,
    ):
        command.set_defaults(_error_parser=command)
    return root


__all__ = ["parser"]
