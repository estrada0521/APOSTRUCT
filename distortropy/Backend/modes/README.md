# Mode Domain

This directory owns local mode construction, structure expansion, and output
presentation for Distortropy.

The code under `engine/` evaluates Source tables and projection routines. This
runtime is self-contained: it must not import from sibling Disassembled
workspaces or invoke external executables.

Boundary rules:

- runtime calls stay inside `distortropy/Backend`;
- oracle generation and implementation-analysis evidence stay under
  `ReverseEngineering/` or `Validation/`;
- the public boundary consumes selected reciprocal-space and isotropy state and
  returns Distortropy mode-detail payloads.

Runtime module boundaries:

- `common.py` holds shared numeric and decoded-data primitives;
- `site_transport.py` and `request_context.py` normalize site and k-vector input;
- `structure_runtime.py` builds child-orbit and Undistorted structure layouts;
- `print_layout.py` and `print_intertwiner.py` construct Source print columns;
- `definition_presentation.py` formats displacive, magnetic, and strain definitions;
- `subduction_specs.py` enumerates single and coupled render specifications;
- `engine/` contains Source-table decoding, subgroup expansion, and projection;
  and
- `assembly.py` is the public domain orchestrator.

Dependencies must point toward the earlier, narrower layers.  Leaf modules
must never import `assembly.py`.
