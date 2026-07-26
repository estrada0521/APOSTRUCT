"""Project-vector expansion and its coefficient-building stages."""

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
        """Evaluate ``project_vector_`` from a materialized boundary state.

        The explicit scalar arguments form the projection boundary contract;
        callers must derive them from Source-backed pipeline state.
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

        component_start = 1
        emitted_count = 0
        while component_start < 4:
            next_emitted_count = emitted_count
            rep_index = (int(site_pg) - 1) * 6 + (int(vector_setting) - 1) * 3 + (component_start - 1)
            if int(target_vector_rep) == int(reps[rep_index]):
                if int(project_count) == 0:
                    next_emitted_count = emitted_count + 1
                elif int(project_count) > 0:
                    family = 0
                    output_family_offset = emitted_count * 144
                    while family < int(project_count):
                        next_emitted_count += 1
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
                                output_index = output_family_offset + axis + atom * 3
                                if 0 <= output_index < len(out):
                                    out[output_index] = value
                        family += 1
                        output_family_offset += 144
            emitted_count = next_emitted_count
            component_start += max(1, int(vector_dim))
        return emitted_count, out

    @staticmethod
    def project_vector_bridge_source_e8_from_active_values(
        active_values: Iterable[Iterable[object]],
        *,
        family_stride: int = 144,
        output_length: int = 768,
    ) -> list[dict[str, object]]:
        """Synthesize bridge source buffers from a projected basis.

        `project_` stores one selected real block per family with the local
        layout `parent_row * 3 + site_col` and a 144-double family stride.
        The coefficient stage reads a shifted window whose active entries begin at
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
        """Build project-vector basis functions from source rows and coefficients.

        Source rows use a 48-value stride and coefficients use three values per
        term.  The leading output block contains one vector per atom; the block
        at offset 144 preserves the coefficient rows consumed by the next stage.
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
        """Build the coefficient buffer from projected basis rows and weights.

        Source families are stored in three-value rows, while operation weights
        use a 48-value column stride.  The output is atom-major with three vector
        components per atom.
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
