"""Ported ``project_vector_`` and adjacent MAIN__ bridge loops."""

from __future__ import annotations

from typing import Iterable


class ProjectVectorMixin:
    def project_vector_from_boundary(
        self,
        *,
        site_pg: int,
        vector_basis_id: int,
        target_vector_rep: int,
        atom_count: int,
        vector_dim: int,
        vector_setting: int,
        point_op: int,
        project_count: int,
        basis_function: Iterable[float],
        output_before: Iterable[float] = (),
        output_length: int = 720,
        sign: float = 1.0,
    ) -> tuple[int, list[float]]:
        """Port ``project_vector_`` for a captured routine-boundary state.

        This helper intentionally takes the same scalar values that the binary
        receives at ``project_vector_`` entry.  It is for routine-boundary
        verification first; later the upstream Python pipeline should provide
        these values without GDB injection.
        """

        basis = list(float(x) for x in basis_function)
        out = list(float(x) for x in output_before)
        if len(out) < output_length:
            out.extend([0.0] * (output_length - len(out)))
        else:
            out = out[:output_length]

        reps = self.iso.wyckoff["iwyckoff_pg_vector_reps"]
        basis_codes = self.iso.wyckoff["iwyckoff_pg_vector_basis"]
        point = tuple(int(x) for x in self.iso.space["ipoint_op"][9 * (int(point_op) - 1):9 * int(point_op)])

        emitted = 0
        component_start = 1
        local_94 = 0
        while component_start < 4:
            local_80 = local_94
            rep_index = (int(site_pg) - 1) * 6 + (int(vector_setting) - 1) * 3 + (component_start - 1)
            if int(target_vector_rep) == int(reps[rep_index]):
                if int(project_count) == 0:
                    local_80 = local_94 + 1
                elif int(project_count) > 0:
                    family = 0
                    local_88 = (local_94 * 9 + 9) * 16
                    while family < int(project_count):
                        local_80 += 1
                        for atom in range(max(0, int(atom_count))):
                            local_vec = [0.0, 0.0, 0.0]
                            for basis_col in range(max(0, int(vector_dim))):
                                input_index = family * 144 + atom * 3 + basis_col
                                input_value = basis[input_index] if input_index < len(basis) else 0.0
                                code_base = (
                                    (int(vector_basis_id) - 1) * 18
                                    + (int(vector_setting) - 1) * 9
                                    + (component_start - 1) * 3
                                    + basis_col * 3
                                )
                                for axis in range(3):
                                    code = int(basis_codes[code_base + axis])
                                    local_vec[axis] += float(self.iso.const[code]) * input_value
                            for axis in range(3):
                                # the mode kernel walks point-op columns here:
                                # indices axis, axis+3, axis+6 in row-major storage.
                                value = sum(float(point[inner * 3 + axis]) * sign * local_vec[inner] for inner in range(3))
                                output_index = local_88 - 145 + (axis + 1) + atom * 3
                                if 0 <= output_index < len(out):
                                    out[output_index] = value
                        family += 1
                        local_88 += 144
            local_94 = local_80
            component_start += max(1, int(vector_dim))
        emitted = local_94
        return emitted, out

    @staticmethod
    def project_vector_bridge_source_e8_from_active_values(
        active_values: Iterable[Iterable[object]],
        *,
        family_stride: int = 144,
        output_length: int = 768,
    ) -> list[dict[str, object]]:
        """Synthesize MAIN__ bridge source buffers from a `project_` basis.

        `project_` stores one selected real block per family with the local
        layout `parent_row * 3 + site_col` and a 144-double family stride.
        The first MAIN__ bridge loop does not read that block from index zero;
        it points at a window whose active entries begin at
        `64 - 3 * (bridge_project_count + 1)`.  `bridge_project_count` is the
        highest parent-row index in the current family, not the number of
        selected project families.
        """

        grouped: dict[int, list[tuple[int, float]]] = {}
        for item in active_values:
            values = list(item)
            if len(values) < 2:
                continue
            index = int(values[0])
            value = float(values[1])
            family = index // int(family_stride)
            local_index = index % int(family_stride)
            grouped.setdefault(family, []).append((local_index, value))

        out: list[dict[str, object]] = []
        for family in sorted(grouped):
            entries = grouped[family]
            if not entries:
                continue
            bridge_project_count = max(local_index // 3 for local_index, _value in entries)
            vector_loop_count = max(local_index % 3 for local_index, _value in entries) + 1
            offset = 64 - 3 * (bridge_project_count + 1)
            source = [0.0] * int(output_length)
            for local_index, value in entries:
                target = offset + local_index
                if 0 <= target < len(source):
                    source[target] = value
            out.append({
                "family": family,
                "bridge_project_count": int(bridge_project_count),
                "vector_loop_count": int(vector_loop_count),
                "source_e8_minus512": source,
            })
        return out

    def project_vector_bridge_second_stage(
        *,
        vector_loop_count: int,
        atom_count: int,
        project_count: int,
        source_minus256: Iterable[float],
        coefficients: Iterable[float],
        output_length: int = 360,
    ) -> list[float]:
        """Port the MAIN__ bridge loop immediately before ``project_vector_``.

        This corresponds to the checked address range ``0x4148ea..0x41497f``.
        The routine takes the source matrix pointed to by stack slot ``0x178``
        dumped from ``ptr-0x200`` and the coefficient buffer beginning at
        ``rsp+0x8640``.  It writes the stack buffer beginning at ``rsp+0x81c0``,
        which is later passed as ``project_vector_`` argument 9.
        """

        source = [float(x) for x in source_minus256]
        coeff = [float(x) for x in coefficients]
        out = [0.0] * int(output_length)
        term_count = int(project_count) + 1
        for vector_index in range(max(0, int(vector_loop_count))):
            for atom_index in range(max(0, int(atom_count))):
                value = 0.0
                for term_index in range(max(0, term_count)):
                    source_index = atom_index * 48 + (16 - term_count) + term_index
                    coeff_index = vector_index + term_index * 3
                    source_value = source[source_index] if 0 <= source_index < len(source) else 0.0
                    coeff_value = coeff[coeff_index] if 0 <= coeff_index < len(coeff) else 0.0
                    value += source_value * coeff_value
                output_index = vector_index + atom_index * 3
                if 0 <= output_index < len(out):
                    out[output_index] = value
            for term_index in range(max(0, term_count)):
                coeff_index = vector_index + term_index * 3
                output_index = 144 + term_index * 3 + vector_index
                if 0 <= coeff_index < len(coeff) and 0 <= output_index < len(out):
                    out[output_index] = coeff[coeff_index]
        return out
    def project_vector_bridge_first_stage(
        *,
        vector_loop_count: int,
        supercell_atom_count: int,
        project_count: int,
        source_e8_minus512: Iterable[float],
        stack_0490: Iterable[float],
        output_length: int = 181,
    ) -> list[float]:
        """Port the first MAIN__ bridge loop that builds the coefficient buffer.

        This mirrors ``0x4147c2..0x414839``.  It fills the coefficient buffer
        beginning at ``rsp+0x8640`` from the source pointer stored at stack slot
        ``0xe8`` and the stack work array beginning near ``rsp+0x490``.
        """

        source = [float(x) for x in source_e8_minus512]
        weights = [float(x) for x in stack_0490]
        out = [0.0] * int(output_length)
        term_count = int(project_count) + 1
        for atom_index in range(max(0, int(supercell_atom_count))):
            for vector_index in range(max(0, int(vector_loop_count))):
                value = 0.0
                for term_index in range(max(0, term_count)):
                    source_index = 64 + vector_index - 3 * (term_count - term_index)
                    weight_index = atom_index + term_index * 48
                    source_value = source[source_index] if 0 <= source_index < len(source) else 0.0
                    weight_value = weights[weight_index] if 0 <= weight_index < len(weights) else 0.0
                    value += source_value * weight_value
                output_index = atom_index * 3 + vector_index
                if 0 <= output_index < len(out):
                    out[output_index] = value
        return out
