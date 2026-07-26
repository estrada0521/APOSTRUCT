"""Decoders for static Source tables used by isotropy and mode kernels."""

from __future__ import annotations

from fractions import Fraction
import math
from numbers import Integral
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from ISODISTORT.Assembled.Backend.source.representation import decode_little_sparse_matrix, real_phase_operator

from ISODISTORT.Assembled.Backend.exactmath import (
    fraction_identity3,
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
)
from ISODISTORT.Assembled.Backend.source import magnetic as magnetic_data
from ISODISTORT.Assembled.Backend.source.iso_data import ISOData
from ISODISTORT.Assembled.Backend.source.tables import (
    SOURCE,
    SPACE_SETTINGS_CINTER_BASE,
    source_tables,
)
from ISODISTORT.Assembled.Backend.isotropy.engine.lattice import _reduce_fraction_vector, get_new_fractionals
from ISODISTORT.Assembled.Backend.isotropy.engine.orderparam import orderparam_check


class SourceData:
    """Shared table decoder for isotropy and subgroup construction."""

    _IDENTITY_AFFINE_RECORD = (
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
    )

    def __init__(self, source: str | Path = SOURCE, *, tables: Any | None = None):
        self.source = Path(source)
        if tables is None:
            tables = (
                source_tables().iso
                if self.source.resolve() == SOURCE.resolve()
                else ISOData(self.source)
            )
        table_dir_value = getattr(tables, "data_dir", None)
        if table_dir_value is None:
            table_dir_value = getattr(tables, "source", None)
        if table_dir_value is None:
            raise ValueError("Source tables must declare their data directory")
        table_dir = Path(table_dir_value)
        if table_dir.resolve() != self.source.resolve():
            raise ValueError(
                f"Source table directory mismatch: {table_dir} != {self.source}"
            )
        self.space = tables.space
        self.little = tables.little
        self.const = tables.const

    def _sequence_key(
        self,
        sequence: Sequence[tuple[int, ...]],
        *,
        cache_name: str,
    ) -> tuple[tuple[int, ...], ...]:
        """Return a stable value key without relying on reusable object ids.

        Dynamic-row generation calls pair multiplication many times with the
        same operation/fraction sequences.  A raw ``id(sequence)`` key is fast
        but unsafe across long comparisons because Python can reuse ids after
        a previous case is freed.  This cache keeps the live sequence object
        alongside its value key, preventing that reuse while avoiding repeated
        tuple materialization in hot loops.
        """

        cache = getattr(self, cache_name, None)
        if cache is None:
            cache = {}
            setattr(self, cache_name, cache)
        ident = id(sequence)
        cached = cache.get(ident)
        if cached is not None and cached[0] is sequence:
            return cached[1]
        key = tuple(tuple(int(item) for item in record) for record in sequence)
        cache[ident] = (sequence, key)
        return key

    def generate_space_group_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return the operation records assembled by ``generate_space_group_``."""

        sg = int(sg)
        point_group = int(self.space["ispace_point_group"][sg - 1])
        count = int(self.space["ipoint_group_order"][point_group - 1])
        pointer = int(self.space["ispace_elements_pointer"][sg - 1])
        raw = self.space["ispace_elements"][(pointer - 1) * 5 : (pointer - 1 + count) * 5]
        return tuple(
            tuple(int(value) for value in raw[offset : offset + 5])
            for offset in range(0, len(raw), 5)
        )

    def magnetic_parent_group_for_space_group(self, sg: int) -> int:
        """Return the first magnetic group used by ``SETTING MAGNETIC``."""

        table = magnetic_data.data().table
        return int(table["mag_space_group_pointer"][int(sg) - 1])

    def magnetic_orderparam_group_for_irrep(self, gid: int) -> int:
        """Return the magnetic group used by ``orderparam_to_subgroup_magnetic_``.

        The ISO routine indexes the parent space group through the selected
        irrep and then advances one slot from ``mag_space_group_pointer``.  For
        grey magnetic irrep calculations this includes the time-reversal point
        operations that are absent from the first ordinary magnetic setting.
        """

        table = magnetic_data.data().table
        sg = int(self.little["little_irr_space_group"][int(gid) - 1])
        return int(table["mag_space_group_pointer"][sg - 1]) + 1

    def generate_magnetic_space_group_records(
        self,
        magnetic_group: int,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Generate the closure used by ``generate_spacegroup_magnetic_``.

        Magnetic records keep magnetic point-operation numbers.  Translation
        rotation uses the corresponding nonmagnetic point operation, matching
        ``gmlt_magnetic_``.
        """

        magnetic_group = int(magnetic_group)
        cache = getattr(self, "_generate_magnetic_space_group_records_cache", None)
        if cache is None:
            cache = {}
            setattr(self, "_generate_magnetic_space_group_records_cache", cache)
        cached = cache.get(magnetic_group)
        if cached is not None:
            return cached

        table = magnetic_data.data().table
        sg = int(table["mag_space_group"][magnetic_group - 1])
        generator_count = int(table["mag_gen_count"][magnetic_group - 1])
        raw = table["mag_gen"][(magnetic_group - 1) * 20 : magnetic_group * 20]
        records: list[tuple[int, int, int, int, int]] = [
            self._reduce_operation_record(tuple(int(value) for value in raw[index * 5 : (index + 1) * 5]))  # type: ignore[arg-type]
            for index in range(generator_count)
        ]
        by_point_op = {int(record[4]): position for position, record in enumerate(records, start=1)}
        changed = True
        while changed:
            changed = False
            scan_count = len(records)
            for left in records[:scan_count]:
                for right in records[:scan_count]:
                    product = self._compose_magnetic_operation_records(left, right, sg=sg)
                    if int(product[4]) in by_point_op:
                        continue
                    by_point_op[int(product[4])] = len(records) + 1
                    records.append(product)
                    changed = True
        result = tuple(sorted(records, key=lambda record: int(record[4])))
        cache[magnetic_group] = result
        return result

    def _compose_magnetic_operation_records(
        self,
        left_record: tuple[int, int, int, int, int],
        right_record: tuple[int, int, int, int, int],
        *,
        sg: int,
    ) -> tuple[int, int, int, int, int]:
        table = magnetic_data.data().table
        nonmag_left_op = int(table["mag_point_op_mag2nonmag"][int(left_record[4]) - 1])
        fake_left = (int(left_record[0]), int(left_record[1]), int(left_record[2]), int(left_record[3]), nonmag_left_op)
        lattice = int(self.space["ispace_lattice"][int(sg) - 1])
        rotated_right = self._rotate_kernel_fraction_by_space_operation(
            (int(right_record[0]), int(right_record[1]), int(right_record[2]), int(right_record[3])),
            fake_left,
            lattice=lattice,
        )
        left_denominator = int(left_record[3])
        right_denominator = int(rotated_right[3])
        if left_denominator == 0 or right_denominator == 0:
            raise ZeroDivisionError("operation record has zero translation denominator")
        denominator = math.lcm(left_denominator, right_denominator)
        left_scale = denominator // left_denominator
        right_scale = denominator // right_denominator
        numerators = [
            int(left_record[axis]) * left_scale + int(rotated_right[axis]) * right_scale
            for axis in range(3)
        ]
        x, y, z, den = _reduce_fraction_vector((numerators[0], numerators[1], numerators[2], denominator))
        # ``mag_point_op_mlt`` is right-major: right operation, then left.
        point_op = int(
            table["mag_point_op_mlt"][(int(right_record[4]) - 1) * 144 + (int(left_record[4]) - 1)]
        )
        return self._reduce_operation_record((x, y, z, den, point_op))

    def inter_setting_basis_transform(self, setting_index: int, half: int = 0) -> tuple[int, ...]:
        """Return the 3x3 basis block from ``ispace_settings_inter``.

        ``data_space:ispace_settings_inter`` stores each setting as two 4x4
        affine integer transforms.  The basis transform used by
        ``id_subgroup_`` is the upper-left 3x3 block of one half.
        ``setting_index`` is the 1-based Source setting id; ``half`` selects
        the first or second affine transform in its 32-value table block.
        """

        if half not in (0, 1):
            raise ValueError(f"invalid inter setting half: {half}")
        start = (int(setting_index) - 1) * 32 + half * 16
        raw = [int(value) for value in self.space["ispace_settings_inter"][start : start + 16]]
        if len(raw) != 16:
            raise IndexError(f"inter setting index out of range: {setting_index}")
        return tuple(raw[row * 4 + col] for row in range(3) for col in range(3))

    def matching_inter_setting_transforms(
        self,
        transform: tuple[int, ...] | list[int],
        *,
        child_sg: int | None = None,
    ) -> tuple[dict[str, Any], ...]:
        """Find ``data_space`` inter-setting records with this 3x3 transform."""

        target = tuple(int(value) for value in transform[:9])
        out: list[dict[str, Any]] = []
        count = len(self.space["ispace_inter_number"])
        for index in range(1, count + 1):
            sg = int(self.space["ispace_inter_number"][index - 1])
            if child_sg is not None and sg != int(child_sg):
                continue
            for half in (0, 1):
                if self.inter_setting_basis_transform(index, half) != target:
                    continue
                out.append({
                    "setting_index": index,
                    "half": half,
                    "space_group": sg,
                    "axis": str(self.space["space_inter_axis"][index - 1]),
                    "abc": str(self.space["space_inter_abc"][index - 1]),
                    "label": str(self.space["space_inter_label_full"][index - 1]).strip(),
                })
        return tuple(out)

    def point_operation_matrix(self, point_op: int) -> tuple[int, ...]:
        """Return a 3x3 point-operation matrix from ``ipoint_op``."""

        start = (int(point_op) - 1) * 9
        raw = self.space["ipoint_op"][start : start + 9]
        if len(raw) != 9:
            raise IndexError(f"point operation out of range: {point_op}")
        return tuple(int(value) for value in raw)

    def ml_lattice_setting_record(self, lattice: int) -> tuple[int, ...]:
        """Return the 36-int ML lattice setting record used by ``id_subgroup_``."""

        record_index = int(lattice)
        raw = self.space["lattice_ml"][(record_index - 1) * 36 : record_index * 36]
        if len(raw) != 36:
            raise IndexError(f"ML lattice setting record out of range: {lattice}")
        return tuple(int(value) for value in raw)

    def ml_lattice_setting_denominators(self, lattice: int) -> tuple[int, int]:
        """Return the two denominators paired with ``ml_lattice_setting_record``."""

        record_index = int(lattice)
        denom = self.space["lattice_ml_denom"][(record_index - 1) * 2 : record_index * 2]
        if len(denom) != 2:
            raise IndexError(f"ML lattice setting denominator out of range: {lattice}")
        return (int(denom[0]), int(denom[1]))

    @staticmethod
    def _reduce_integer_matrix_denominator(
        matrix: tuple[int, ...] | list[int],
        denominator: int,
    ) -> tuple[tuple[int, ...], int]:
        gcd = abs(int(denominator))
        for value in matrix:
            gcd = math.gcd(gcd, abs(int(value)))
        if gcd <= 0:
            gcd = 1
        reduced_denominator = int(denominator) // gcd
        reduced = tuple(int(value) // gcd for value in matrix)
        if reduced_denominator < 0:
            reduced_denominator = -reduced_denominator
            reduced = tuple(-value for value in reduced)
        return reduced, reduced_denominator

    @staticmethod
    def _matmlt_iso_storage(
        left: tuple[int, ...] | list[int],
        right: tuple[int, ...] | list[int],
    ) -> tuple[int, ...]:
        """Apply ``matmlt_`` storage semantics to flat 3x3 integer buffers."""

        out = [0] * 9
        for outer in range(3):
            right_offset = 0
            for group in range(3):
                total = 0
                for k in range(3):
                    total += int(left[outer + 3 * k]) * int(right[right_offset + k])
                out[outer + 3 * group] = total
                right_offset += 3
        return tuple(out)

    def id_subgroup_candidate_generator_matrix(
        self,
        candidate_subgroup: int,
        point_op: int,
    ) -> tuple[tuple[int, int, int], ...]:
        """Build the candidate generator matrix used by ``id_subgroup_``.

        The routine conjugates the canonical point-operation matrix by a
        point-group convention transform selected from Source lattice tables:

        ``S0 * ipoint_op(point_op) * S1``

        followed by ``reduc3_`` with the first setting denominator.  This is
        the 3x3 buffer later used to build the origin equations.
        """

        lattice = int(self.space["ispace_lattice"][int(candidate_subgroup) - 1])
        setting = self.ml_lattice_setting_record(lattice)
        denom0, _denom1 = self.ml_lattice_setting_denominators(lattice)
        left = setting[:9]
        right = setting[9:18]
        raw = self._matmlt_iso_storage(
            self._matmlt_iso_storage(left, self.point_operation_matrix(point_op)),
            right,
        )
        reduced, denominator = self._reduce_integer_matrix_denominator(raw, denom0)
        if denominator != 1:
            raise ValueError(
                f"id_subgroup candidate matrix did not reduce to integer: "
                f"sg={candidate_subgroup} op={point_op} denominator={denominator}"
            )
        return tuple(
            tuple(int(reduced[row * 3 + col]) for col in range(3))
            for row in range(3)
        )  # type: ignore[return-value]

    @staticmethod
    def _normalize_setting_label(setting: str, length: int | None = None) -> str:
        text = str(setting)
        if length is not None:
            text = text[: int(length)]
        return text.rstrip()

    @staticmethod
    def _reduce_operation_record(record: tuple[int, int, int, int, int]) -> tuple[int, int, int, int, int]:
        x, y, z, den, op = (int(value) for value in record)
        if den == 0:
            return record
        x %= den
        y %= den
        z %= den
        rx, ry, rz, rden = _reduce_fraction_vector((x, y, z, den))
        return (rx % rden, ry % rden, rz % rden, rden, op)

    def _pml_to_cml_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        lattice = int(self.space["ispace_lattice"][int(sg) - 1])
        raw = self.space["lattice_ml"][(lattice - 1) * 36 + 9 : (lattice - 1) * 36 + 18]
        den = int(self.space["lattice_ml_denom"][(lattice - 1) * 2])
        if den == 0:
            raise ValueError(f"zero lattice_ml pml->cml denominator for lattice {lattice}")
        return tuple(
            tuple(Fraction(int(raw[3 * row + col]), den) for col in range(3))
            for row in range(3)
        )  # type: ignore[return-value]

    def _cml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        choice = int(self.space["ispace_inter_choice"][int(sg) - 1])
        raw = tuple(
            int(value)
            for value in self.space["ispace_settings_inter"][(choice - 1) * 32 : (choice - 1) * 32 + 16]
        )
        if len(raw) != 16:
            raise IndexError(f"inter setting choice out of range: {choice}")
        den = int(raw[15])
        if den == 0:
            raise ValueError(f"zero inter-setting denominator for SG{sg}")
        return tuple(
            tuple(Fraction(raw[4 * row + col], den) for col in range(3))
            for row in range(3)
        )  # type: ignore[return-value]

    def _pml_to_cinter_matrix(self, sg: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        return fraction_matrix_multiply3(
            self._pml_to_cml_matrix(sg),
            self._cml_to_cinter_matrix(sg),
        )

    def setting_change_matrix(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        identity = fraction_identity3()
        to_cinter = {
            "cinter": identity,
            "cml": self._cml_to_cinter_matrix(sg),
            "pml": self._pml_to_cinter_matrix(sg),
        }
        source = str(from_setting).strip().lower()
        target = str(to_setting).strip().lower()
        return fraction_matrix_multiply3(
            to_cinter[source],
            fraction_matrix_inverse3(to_cinter[target]),
        )

    @staticmethod
    def _fraction_record_from_values(values: tuple[Fraction, Fraction, Fraction]) -> tuple[int, int, int, int]:
        den = 1
        for value in values:
            den = math.lcm(den, value.denominator)
        nums = [int(value * den) for value in values]
        return _reduce_fraction_vector((nums[0], nums[1], nums[2], den))

    def vector_change_setting_record(
        self,
        sg: int,
        from_setting: str,
        to_setting: str,
        record: tuple[int, int, int, int],
    ) -> tuple[int, int, int, int]:
        x, y, z, den = (int(value) for value in record)
        values = (Fraction(x, den), Fraction(y, den), Fraction(z, den))
        matrix = self.setting_change_matrix(sg, from_setting, to_setting)
        out = tuple(
            sum(values[row] * matrix[row][col] for row in range(3))
            for col in range(3)
        )
        return self._fraction_record_from_values(out)  # type: ignore[arg-type]

    @staticmethod
    def _fraction_vector_to_record(values: tuple[Fraction, Fraction, Fraction], op: int) -> tuple[int, int, int, int, int]:
        den = 1
        for value in values:
            den = den * value.denominator // math.gcd(den, value.denominator)
        nums = [int(value * den) for value in values]
        x, y, z, rden = _reduce_fraction_vector((nums[0], nums[1], nums[2], den))
        return (x % rden, y % rden, z % rden, rden, int(op))

    def cml_to_pml_operation_record(self, sg: int, record: tuple[int, int, int, int, int]) -> tuple[int, int, int, int, int]:
        """Apply the cml->pml translation part used after ``op_change_setting_``.

        The initialized generator block already supplies pml point-operation
        numbers.  This helper therefore applies only the affine translation
        part of ``op_change_setting_(cml,pml)``.
        """

        x, y, z, den, op = (int(value) for value in record)
        cml_to_pml = fraction_matrix_inverse3(self._pml_to_cml_matrix(sg))
        values = (Fraction(x, den), Fraction(y, den), Fraction(z, den))
        out = tuple(
            sum(values[row] * cml_to_pml[row][col] for row in range(3))
            for col in range(3)
        )
        return self._fraction_vector_to_record(out, op)

    def get_generators_records_from_initialized_block(
        self,
        sg: int,
        setting: str,
        generators_conv: list[int] | tuple[int, ...],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Reproduce ``get_generators_`` from an initialized common-space block."""

        sg = int(sg)
        requested = self._normalize_setting_label(setting)
        if not requested:
            raise ValueError("empty generator setting")
        table_setting = "c" + requested[1:]
        labels = [self._normalize_setting_label(value) for value in self.space["space_gen_setting_label"]]
        try:
            setting_slot = labels.index(table_setting) + 1
        except ValueError as exc:
            raise KeyError(f"unknown generator setting: {table_setting!r}") from exc
        point_group = int(self.space["ispace_point_group"][sg - 1])
        count = int(self.space["ipoint_group_gen_count"][point_group - 1])
        pointer = int(self.space["ispace_gen_pointer"][(setting_slot - 1) * 230 + (sg - 1)])
        records = [
            tuple(int(value) for value in generators_conv[(pointer + index) * 5 : (pointer + index + 1) * 5])
            for index in range(count)
        ]
        if requested == table_setting:
            return tuple(records)
        if requested == "pml" and table_setting == "cml":
            return tuple(self.cml_to_pml_operation_record(sg, record) for record in records)
        raise NotImplementedError(f"generator setting transform {table_setting!r}->{requested!r}")

    @staticmethod
    def _source_table_integer(value: Any, *, field: str) -> int:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError(f"{field} must be an exact integer")
        return int(value)

    def _space_generator_setting_is_identity(self, sg: int) -> bool:
        start = SPACE_SETTINGS_CINTER_BASE + int(sg) * 32
        raw = self.space["ispace_settings"][start : start + 32]
        if len(raw) != 32:
            raise IndexError(f"generator setting block out of range for SG{sg}")
        values = tuple(
            self._source_table_integer(value, field="ispace_settings")
            for value in raw
        )
        identity = self._IDENTITY_AFFINE_RECORD
        return values[:16] == identity and values[16:] == identity

    def _space_element_generator_records(
        self,
        sg: int,
        *,
        setting_slot: int,
        count: int,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        setting_slot = self._source_table_integer(
            setting_slot, field="generator setting slot"
        )
        count = self._source_table_integer(count, field="generator count")
        pointer_index = (setting_slot - 1) * 230 + (int(sg) - 1)
        pointers = self.space["ispace_gen_pointer"]
        if not 0 <= pointer_index < len(pointers):
            raise IndexError(f"generator pointer out of range for SG{sg}")
        pointer = self._source_table_integer(
            pointers[pointer_index], field="ispace_gen_pointer"
        )
        if pointer <= 0 or count <= 0:
            raise ValueError(f"invalid generator pointer/count for SG{sg}: {pointer}/{count}")

        start = (pointer - 1) * 5
        raw = self.space["ispace_elements"][start : start + count * 5]
        if len(raw) != count * 5:
            raise IndexError(f"generator records out of range for SG{sg}")

        generated = {
            self._reduce_operation_record(tuple(record))
            for record in self.generate_space_group_records(sg)
        }
        records: list[tuple[int, int, int, int, int]] = []
        point_operation_count = len(self.space["ipoint_op_code"])
        for index in range(count):
            record = tuple(
                self._source_table_integer(value, field="ispace_elements")
                for value in raw[index * 5 : (index + 1) * 5]
            )
            if len(record) != 5:
                raise IndexError(f"incomplete generator record for SG{sg}")
            if record[3] <= 0:
                raise ValueError(f"invalid generator denominator for SG{sg}: {record[3]}")
            if not 1 <= record[4] <= point_operation_count:
                raise ValueError(f"invalid generator point operation for SG{sg}: {record[4]}")
            reduced = self._reduce_operation_record(record)
            if reduced not in generated:
                raise ValueError(
                    f"generator record is not a generated SG{sg} element: {reduced}"
                )
            records.append(reduced)
        if len(set(records)) != len(records):
            raise ValueError(f"duplicate generator records for SG{sg}")
        return tuple(records)

    def get_generators_records(self, sg: int, setting: str = "pml") -> tuple[tuple[int, int, int, int, int], ...]:
        """Read generator records with the ``get_generators_`` table contract.

        The routine maps requested primitive settings to the corresponding
        conventional generator table by replacing the first setting character
        with ``c``.  When Source declares the conventional-to-primitive affine
        setting as exact identity, the 1-based generator pointers address
        records directly in the shared space-element table.
        """

        sg = self._source_table_integer(sg, field="space group")
        if not 1 <= sg <= 230:
            raise ValueError(f"space group out of range: {sg}")
        requested = self._normalize_setting_label(setting)
        if not requested:
            raise ValueError("empty generator setting")
        table_setting = "c" + requested[1:]
        labels = [self._normalize_setting_label(value) for value in self.space["space_gen_setting_label"]]
        try:
            setting_slot = labels.index(table_setting) + 1
        except ValueError as exc:
            raise KeyError(f"unknown generator setting: {table_setting!r}") from exc
        point_group = self._source_table_integer(
            self.space["ispace_point_group"][sg - 1], field="ispace_point_group"
        )
        if not 1 <= point_group <= len(self.space["ipoint_group_gen_count"]):
            raise ValueError(f"point group out of range for SG{sg}: {point_group}")
        count = self._source_table_integer(
            self.space["ipoint_group_gen_count"][point_group - 1],
            field="ipoint_group_gen_count",
        )
        if (
            requested == "pml"
            and table_setting == "cml"
            and self._space_generator_setting_is_identity(sg)
        ):
            return self._space_element_generator_records(
                sg,
                setting_slot=setting_slot,
                count=count,
            )
        if requested.startswith("p"):
            point_ops = [
                int(value)
                for value in self.space["ipoint_group_gen"][(point_group - 1) * 5 : (point_group - 1) * 5 + count]
            ]
            generated = [tuple(int(value) for value in record) for record in self.generate_space_group_records(sg)]
            by_op = {int(record[4]): record for record in generated}
            out: list[tuple[int, int, int, int, int]] = []
            used_ops: set[int] = set()
            for op in point_ops:
                record = by_op.get(op)
                if record is None:
                    target_code = int(self.space["ipoint_op_code"][op - 1])
                    candidates = [
                        candidate
                        for candidate in generated
                        if int(candidate[4]) not in used_ops
                        and int(self.space["ipoint_op_code"][int(candidate[4]) - 1]) == target_code
                    ]
                    if candidates:
                        # Multiple generated operations can share a point-op
                        # code. This mapper selects the last representative,
                        # except PG18/op58, where it selects the first
                        # (55/57/59 -> 55). Keep this policy local to generator
                        # setting conversion.
                        if point_group == 18 and int(op) == 58:
                            record = candidates[0]
                        else:
                            record = candidates[-1]
                if record is None:
                    raise KeyError(op)
                used_ops.add(int(record[4]))
                out.append(self._reduce_operation_record(tuple(record)))
            return tuple(out)
        pointer = int(self.space["ispace_gen_pointer"][(setting_slot - 1) * 230 + (sg - 1)])
        raw = self.space["ispace_generators_conv"][pointer * 5 : (pointer + count) * 5]
        records = [
            tuple(int(value) for value in raw[offset : offset + 5])
            for offset in range(0, len(raw), 5)
        ]
        return tuple(records)

    def _space_group_record_count(self, sg: int) -> int:
        point_group = int(self.space["ispace_point_group"][sg - 1])
        return int(self.space["ipoint_group_order"][point_group - 1])

    def generated_space_operation_position(self, sg: int, record: tuple[int, int, int, int, int]) -> int:
        """Return the 1-based generated operation position used by ``get_irrep4_``."""

        request_tau = tuple(Fraction(int(record[i]), int(record[3])) % 1 for i in range(3))
        point_op = int(record[4])
        pointer = int(self.space["ispace_elements_pointer"][sg - 1])
        count = self._space_group_record_count(sg)
        for index in range(count):
            raw = self.space["ispace_elements"][(pointer - 1 + index) * 5 : (pointer + index) * 5]
            if len(raw) != 5:
                continue
            x, y, z, den, op = (int(value) for value in raw)
            if op != point_op:
                continue
            tau = (Fraction(x, den) % 1, Fraction(y, den) % 1, Fraction(z, den) % 1)
            if tau == request_tau:
                return index + 1
        raise KeyError(f"operation record {record} not found in generated SG{sg} records")

    def generated_space_tau(self, sg: int, point_op: int) -> tuple[Fraction, Fraction, Fraction]:
        """Return the generated representative translation for a point operation."""

        for x, y, z, den, op in self.generate_space_group_records(sg):
            if int(op) == int(point_op):
                return (Fraction(int(x), int(den)), Fraction(int(y), int(den)), Fraction(int(z), int(den)))
        raise KeyError(f"point op {point_op} not found in generated SG{sg} records")

    def _little_sparse_matrix_by_gid_record(self, gid: int, record: tuple[int, int, int, int, int]) -> np.ndarray:
        gid = int(gid)
        dim = int(self.little["little_irr_full_dim"][gid - 1])
        sg = int(self.little["little_irr_space_group"][gid - 1])
        position = self.generated_space_operation_position(sg, record)
        return decode_little_sparse_matrix(
            self.little, self.const, gid=gid, dim=dim, position=position
        )

    @staticmethod
    def _phase_operator(dim: int, phases: tuple[Fraction, ...]) -> np.ndarray:
        return real_phase_operator(dim, phases)

    @staticmethod
    def _kparam_fractions(kparam: tuple[int, ...] | list[int]) -> tuple[Fraction, Fraction, Fraction]:
        if len(kparam) < 4:
            return (Fraction(0), Fraction(0), Fraction(0))
        den = int(kparam[3])
        if den == 0:
            return (Fraction(0), Fraction(0), Fraction(0))
        return (
            Fraction(int(kparam[0]), den),
            Fraction(int(kparam[1]), den),
            Fraction(int(kparam[2]), den),
        )

    def _little_k_star_vectors(self, gid: int, kparam: tuple[int, ...] | list[int]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        gid = int(gid)
        sg = int(self.little["little_irr_space_group"][gid - 1])
        kslot = int(self.little["little_irr_k"][gid - 1])
        slot = (sg - 1) * 27 + kslot - 1
        count = int(self.little["little_k_star_count"][slot])
        pointer = int(self.little["little_k_star_ml_pointer"][slot])
        params = self._kparam_fractions(kparam)
        lattice = int(self.space["ispace_lattice"][sg - 1])
        lattice_slot = (lattice - 1) * 27 + kslot - 1
        if int(self.little["little_k_dim"][lattice_slot]) <= 0:
            raw_special = self.little["little_k_star_special"]
            vectors: list[tuple[Fraction, Fraction, Fraction]] = []
            for arm in range(count):
                start = 4 * (pointer - 1 + arm)
                x, y, z, den = (int(value) for value in raw_special[start:start + 4])
                denominator = den or 1
                vectors.append((Fraction(x, denominator), Fraction(y, denominator), Fraction(z, denominator)))
            return tuple(vectors)
        raw = self.little["little_k_star"]
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for arm in range(count):
            start = 16 * (pointer - 1 + arm)
            record = tuple(int(value) for value in raw[start : start + 16])
            pieces: list[tuple[Fraction, Fraction, Fraction]] = []
            for offset in range(0, 16, 4):
                den = int(record[offset + 3])
                if den == 0:
                    pieces.append((Fraction(0), Fraction(0), Fraction(0)))
                else:
                    pieces.append(tuple(Fraction(record[offset + axis], den) for axis in range(3)))  # type: ignore[arg-type]
            out = [pieces[0][axis] for axis in range(3)]
            for param_index, param in enumerate(params):
                for axis in range(3):
                    out[axis] += param * pieces[param_index + 1][axis]
            vectors.append(tuple(out))  # type: ignore[arg-type]
        return tuple(vectors)

    def operation_record_phases(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: tuple[int, ...] | list[int],
    ) -> tuple[Fraction, ...]:
        gid = int(gid)
        sg = int(self.little["little_irr_space_group"][gid - 1])
        kslot = int(self.little["little_irr_k"][gid - 1])
        lattice = int(self.space["ispace_lattice"][sg - 1])
        slot = (lattice - 1) * 27 + kslot - 1
        request_tau = tuple(Fraction(int(record[i]), int(record[3])) for i in range(3))
        if int(self.little["little_k_dim"][slot]) > 0:
            delta = request_tau
        else:
            generated_tau = self.generated_space_tau(sg, int(record[4]))
            delta = tuple(request_tau[i] - generated_tau[i] for i in range(3))
        return tuple(sum(kvec[i] * delta[i] for i in range(3)) for kvec in self._little_k_star_vectors(gid, kparam))

    def get_irreps_matrix_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: tuple[int, ...] | list[int],
    ) -> np.ndarray:
        """Evaluate a nonmagnetic operation in the ``get_irreps_`` matrix frame."""

        key = (int(gid), tuple(int(value) for value in record), tuple(int(value) for value in kparam))
        cache = getattr(self, "_get_irreps_matrix_for_record_cache", None)
        if cache is None:
            cache = {}
            setattr(self, "_get_irreps_matrix_for_record_cache", cache)
        cached = cache.get(key)
        if cached is not None:
            return cached

        dim = int(self.little["little_irr_full_dim"][int(gid) - 1])
        phases = self.operation_record_phases(gid, record, kparam)
        matrix = self._phase_operator(dim, phases) @ self._little_sparse_matrix_by_gid_record(gid, record)
        cache[key] = matrix
        return matrix

    def get_irreps_matrix_for_magnetic_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        kparam: tuple[int, ...] | list[int],
    ) -> np.ndarray:
        """Evaluate a magnetic operation in the ``get_irreps_`` matrix frame.

        Convert the magnetic point operation to its nonmagnetic partner before
        calling ``get_irrep4_``. If the magnetic
        operation carries time reversal, it multiplies the resulting matrix by
        ``-1``.
        """

        table = magnetic_data.data().table
        key = (int(gid), tuple(int(value) for value in record), tuple(int(value) for value in kparam))
        cache = getattr(self, "_get_irreps_matrix_for_magnetic_record_cache", None)
        if cache is None:
            cache = {}
            setattr(self, "_get_irreps_matrix_for_magnetic_record_cache", cache)
        cached = cache.get(key)
        if cached is not None:
            return cached
        mag_op = int(record[4])
        nonmag_op = int(table["mag_point_op_mag2nonmag"][mag_op - 1])
        nonmag_record = (
            int(record[0]),
            int(record[1]),
            int(record[2]),
            int(record[3]),
            nonmag_op,
        )
        matrix = self.get_irreps_matrix_for_record(gid, nonmag_record, kparam)
        if bool(table["mag_point_op_r"][mag_op - 1]):
            matrix = -matrix
        cache[key] = matrix
        return matrix

    @staticmethod
    def _matrix_to_stride48(matrix: np.ndarray) -> tuple[float, ...]:
        out = [0.0] * (48 * 48)
        rows, cols = matrix.shape
        for row in range(rows):
            for col in range(cols):
                out[row + col * 48] = float(matrix[row, col])
        return tuple(out)

    @staticmethod
    def _reduce_basis_rows(rows: list[list[int]]) -> tuple[int, ...]:
        for _ in range(3):
            for row_index in range(3):
                current = rows[row_index]
                current_size = max(abs(value) for value in current)
                for other_index in range(3):
                    if other_index == row_index:
                        continue
                    for sign in (-1, 1):
                        candidate = [
                            current[axis] + sign * rows[other_index][axis]
                            for axis in range(3)
                        ]
                        candidate_size = max(abs(value) for value in candidate)
                        if candidate_size < current_size:
                            rows[row_index] = candidate
                            current = candidate
                            current_size = candidate_size
        return tuple(value for row in rows for value in row)

    def _identity_point_op(self, sg: int) -> int:
        return int(self.generate_space_group_records(sg)[0][4])

    def _invariant_translation(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        vector: tuple[int, int, int],
        orderparam: tuple[float, ...] | list[float],
    ) -> bool:
        sg = int(self.little["little_irr_space_group"][int(gid) - 1])
        record = (int(vector[0]), int(vector[1]), int(vector[2]), 1, self._identity_point_op(sg))
        matrix = self.get_irreps_matrix_for_record(gid, record, kparam)
        dim = int(matrix.shape[0])
        return orderparam_check(dim, 1, self._matrix_to_stride48(matrix), orderparam) == 1

    @staticmethod
    def _inverse_fraction_matrix(matrix: tuple[int, ...] | list[int]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        rows = tuple(
            tuple(Fraction(int(matrix[3 * row + column])) for column in range(3))
            for row in range(3)
        )
        try:
            return fraction_matrix_inverse3(rows)
        except ValueError:
            raise ValueError("singular basis")

    @staticmethod
    def _vector_in_basis(vector: tuple[int, int, int], basis: tuple[int, ...] | list[int]) -> bool:
        inv = SourceData._inverse_fraction_matrix(basis)
        coords = [
            sum(Fraction(vector[row]) * inv[row][col] for row in range(3))
            for col in range(3)
        ]
        return all(value.denominator == 1 for value in coords)

    def _search_newlat_basis(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        orderparam: tuple[float, ...] | list[float],
        *,
        contains: tuple[int, ...] | list[int] | None = None,
        max_iter: int = 256,
    ) -> tuple[int, ...]:
        def allowed(vector: tuple[int, int, int]) -> bool:
            if contains is not None and not self._vector_in_basis(vector, contains):
                return False
            return self._invariant_translation(gid, kparam, vector, orderparam)

        periods = [0, 0, 0]
        for axis in range(3):
            trial = [0, 0, 0]
            for value in range(1, max_iter + 1):
                trial[axis] = value
                if allowed(tuple(trial)):  # type: ignore[arg-type]
                    periods[axis] = value
                    break
            if periods[axis] == 0:
                raise ValueError(f"newlat_order_: failed to find period for axis {axis}")

        rows: list[list[int]] = []
        for axis in range(3):
            trial = [0, 0, 0]
            trial[axis] = 1
            for _ in range(max_iter**3):
                if allowed(tuple(trial)):  # type: ignore[arg-type]
                    rows.append(list(trial))
                    break
                if axis < 0:
                    raise AssertionError("unreachable")
                trial[0] += 1
                if trial[0] != periods[0] or axis == 0:
                    continue
                trial[0] = 0
                trial[1] += 1
                if trial[1] != periods[1] or axis == 1:
                    continue
                trial[1] = 0
                trial[2] += 1
                if trial[2] != periods[2] or axis == 2:
                    continue
                raise ValueError(f"newlat_order_: failed to find basis row {axis}")
            else:
                raise ValueError(f"newlat_order_: exceeded search guard for basis row {axis}")
        return self._reduce_basis_rows(rows)

    def newlat_order(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        orderparam: tuple[float, ...] | list[float],
    ) -> tuple[int, ...]:
        """Select a lattice with ``newlat_order_`` rules for one OPD row."""

        return self._search_newlat_basis(gid, kparam, orderparam)

    def newlat_order3(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        orderparam: tuple[float, ...] | list[float],
        previous_basis: tuple[int, ...] | list[int],
    ) -> tuple[int, ...]:
        """Select a constrained lattice with ``newlat_order3_`` rules."""

        return self._search_newlat_basis(gid, kparam, orderparam, contains=previous_basis)

    @staticmethod
    def _matmul3(left: tuple[int, ...] | list[int], right: tuple[int, ...] | list[int]) -> tuple[int, ...]:
        out: list[int] = []
        for row in range(3):
            for col in range(3):
                out.append(sum(int(left[row * 3 + k]) * int(right[k * 3 + col]) for k in range(3)))
        return tuple(out)

    @staticmethod
    def _basis_transform_matrix(old_basis: tuple[int, ...] | list[int], new_basis: tuple[int, ...] | list[int]) -> tuple[int, ...]:
        inv = SourceData._inverse_fraction_matrix(old_basis)
        out: list[int] = []
        for row in range(3):
            for col in range(3):
                value = sum(Fraction(int(new_basis[row * 3 + k])) * inv[k][col] for k in range(3))
                if value.denominator != 1:
                    raise ValueError("new basis is not an integer transform of old basis")
                out.append(int(value))
        return tuple(out)

    @staticmethod
    def _vadd_fraction_operation(
        fraction: tuple[int, int, int, int],
        operation: tuple[int, int, int, int, int],
    ) -> tuple[int, int, int, int, int]:
        values = [
            Fraction(int(fraction[axis]), int(fraction[3]))
            + Fraction(int(operation[axis]), int(operation[3]))
            for axis in range(3)
        ]
        den = 1
        for value in values:
            den = math.lcm(den, value.denominator)
        return (
            int(values[0] * den),
            int(values[1] * den),
            int(values[2] * den),
            den,
            int(operation[4]),
        )

    @staticmethod
    def _operation_fraction_part(record: tuple[int, int, int, int, int]) -> tuple[int, int, int, int]:
        return (int(record[0]), int(record[1]), int(record[2]), int(record[3]))

    @staticmethod
    def _fraction_values(fraction: tuple[int, int, int, int]) -> tuple[Fraction, Fraction, Fraction]:
        return tuple(Fraction(int(fraction[axis]), int(fraction[3])) for axis in range(3))  # type: ignore[return-value]

    @staticmethod
    def _fraction_record_from_values(values: tuple[Fraction, Fraction, Fraction]) -> tuple[int, int, int, int]:
        denominator = 1
        for value in values:
            denominator = math.lcm(denominator, value.denominator)
        numerators = [int(value * denominator) for value in values]
        return _reduce_fraction_vector((numerators[0], numerators[1], numerators[2], denominator))

    @staticmethod
    def _fraction_add(
        left: tuple[int, int, int, int],
        right: tuple[int, int, int, int],
    ) -> tuple[int, int, int, int]:
        left_values = SourceData._fraction_values(left)
        right_values = SourceData._fraction_values(right)
        return SourceData._fraction_record_from_values(tuple(left_values[i] + right_values[i] for i in range(3)))

    def _kernel_fraction_index_for_delta(
        self,
        *,
        kernel_basis: tuple[int, ...] | list[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        delta: tuple[Fraction, Fraction, Fraction],
    ) -> int:
        """Return the 1-based kernel-fraction coset index used by ``iso_find_op_``."""

        basis_key = tuple(int(value) for value in kernel_basis)
        fraction_key = self._sequence_key(kernel_fractions, cache_name="_kernel_fraction_sequence_key_cache")
        delta_key = tuple(Fraction(value) for value in delta)
        cache_key = (basis_key, fraction_key, delta_key)
        cache = getattr(self, "_kernel_fraction_index_cache", None)
        if cache is None:
            cache = {}
            setattr(self, "_kernel_fraction_index_cache", cache)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached

        inverse_cache = getattr(self, "_inverse_fraction_matrix_cache", None)
        if inverse_cache is None:
            inverse_cache = {}
            setattr(self, "_inverse_fraction_matrix_cache", inverse_cache)
        inverse = inverse_cache.get(basis_key)
        if inverse is None:
            inverse = self._inverse_fraction_matrix(kernel_basis)
            inverse_cache[basis_key] = inverse

        def residue_key(vector: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction]:
            coords = [
                sum(vector[row] * inverse[row][col] for row in range(3))
                for col in range(3)
            ]
            return tuple(value % 1 for value in coords)  # type: ignore[return-value]

        residue_cache = getattr(self, "_kernel_fraction_residue_index_cache", None)
        if residue_cache is None:
            residue_cache = {}
            setattr(self, "_kernel_fraction_residue_index_cache", residue_cache)
        residue_cache_key = (basis_key, fraction_key)
        residue_index = residue_cache.get(residue_cache_key)
        if residue_index is None:
            residue_index = {}
            duplicate = False
            for fraction_index, fraction in enumerate(kernel_fractions, start=1):
                key = residue_key(self._fraction_values(fraction))
                if key in residue_index:
                    duplicate = True
                    break
                residue_index[key] = fraction_index
            residue_cache[residue_cache_key] = None if duplicate else residue_index
        if residue_index is not None:
            direct = residue_index.get(residue_key(delta))
            if direct is not None:
                cache[cache_key] = direct
                return direct

        values_cache = getattr(self, "_kernel_fraction_value_index_cache", None)
        if values_cache is None:
            values_cache = {}
            setattr(self, "_kernel_fraction_value_index_cache", values_cache)
        values_by_fraction = values_cache.get(fraction_key)
        if values_by_fraction is None:
            values_by_fraction = {
                self._fraction_values(fraction): fraction_index
                for fraction_index, fraction in enumerate(kernel_fractions, start=1)
            }
            values_cache[fraction_key] = values_by_fraction

        for candidate, fraction_index in values_by_fraction.items():
            residual = tuple(delta[axis] - candidate[axis] for axis in range(3))
            coords = [
                sum(residual[row] * inverse[row][col] for row in range(3))
                for col in range(3)
            ]
            if all(value.denominator == 1 for value in coords):
                cache[cache_key] = fraction_index
                return fraction_index
        raise KeyError(f"kernel fraction coset not found: {delta}")

    def _operation_pair_for_record(
        self,
        *,
        space_operations: Sequence[tuple[int, int, int, int, int]],
        kernel_basis: tuple[int, ...] | list[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        record: tuple[int, int, int, int, int],
    ) -> tuple[int, int]:
        """Return 1-based ``(ispace, kernel_fraction)`` matching ``iso_find_op_``."""

        point_op = int(record[4])
        target = tuple(Fraction(int(record[axis]), int(record[3])) for axis in range(3))
        point_cache_key = self._sequence_key(space_operations, cache_name="_space_operation_sequence_key_cache")
        point_cache = getattr(self, "_space_index_by_point_op_cache", None)
        if point_cache is None:
            point_cache = {}
            setattr(self, "_space_index_by_point_op_cache", point_cache)
        by_point = point_cache.get(point_cache_key)
        if by_point is None:
            by_point = {}
            for space_index, space_operation in enumerate(space_operations, start=1):
                by_point.setdefault(int(space_operation[4]), []).append((space_index, space_operation))
            point_cache[point_cache_key] = by_point
        items = by_point.get(point_op)
        if not items:
            raise KeyError(f"operation record not represented by ispace/kernel pair: {record}")
        for space_index, space_operation in items:
            base = tuple(Fraction(int(space_operation[axis]), int(space_operation[3])) for axis in range(3))
            delta = tuple(target[axis] - base[axis] for axis in range(3))
            try:
                fraction_index = self._kernel_fraction_index_for_delta(
                    kernel_basis=kernel_basis,
                    kernel_fractions=kernel_fractions,
                    delta=delta,
                )
            except KeyError:
                continue
            return space_index, fraction_index
        raise KeyError(f"operation record not represented by ispace/kernel pair: {record}")

    def _rotate_kernel_fraction_by_space_operation(
        self,
        fraction: tuple[int, int, int, int],
        space_operation: tuple[int, int, int, int, int],
        *,
        lattice: int,
    ) -> tuple[int, int, int, int]:
        """Apply the nonmagnetic ``vrot_(setting, point_op, fraction)`` contract."""

        fraction_values = tuple(int(value) for value in fraction)
        fraction_denominator = int(fraction_values[3])
        if fraction_denominator == 0:
            raise ZeroDivisionError("fraction has zero denominator")
        setting = self.ml_lattice_setting_record(int(lattice))
        setting_denominator, _ = self.ml_lattice_setting_denominators(int(lattice))
        backward = setting[:9]
        forward = setting[9:18]
        matrix = self.point_operation_matrix(int(space_operation[4]))
        first = tuple(
            sum(fraction_values[row] * int(forward[col + 3 * row]) for row in range(3))
            for col in range(3)
        )
        second = tuple(
            sum(first[row] * int(matrix[col + 3 * row]) for row in range(3))
            for col in range(3)
        )
        rotated = tuple(
            sum(second[row] * int(backward[col + 3 * row]) for row in range(3))
            for col in range(3)
        )
        return _reduce_fraction_vector(
            (
                int(rotated[0]),
                int(rotated[1]),
                int(rotated[2]),
                fraction_denominator * int(setting_denominator),
            )
        )

    def _compose_operation_records(
        self,
        left_record: tuple[int, int, int, int, int],
        right_record: tuple[int, int, int, int, int],
        space_operations: Sequence[tuple[int, int, int, int, int]],
        *,
        lattice: int,
    ) -> tuple[int, int, int, int, int]:
        left_t = tuple(Fraction(int(left_record[axis]), int(left_record[3])) for axis in range(3))
        rotated_right = self._fraction_values(
            self._rotate_kernel_fraction_by_space_operation(
                (int(right_record[0]), int(right_record[1]), int(right_record[2]), int(right_record[3])),
                left_record,
                lattice=int(lattice),
            )
        )
        product_t = tuple(
            (
                left_t[col]
                + rotated_right[col]
            )
            for col in range(3)
        )
        left_matrix = self.point_operation_matrix(int(left_record[4]))
        right_matrix = self.point_operation_matrix(int(right_record[4]))
        product_matrix = tuple(
            sum(int(right_matrix[row * 3 + k]) * int(left_matrix[k * 3 + col]) for k in range(3))
            for row in range(3)
            for col in range(3)
        )
        point_matrix_cache_key = self._sequence_key(
            space_operations,
            cache_name="_space_operation_sequence_key_cache",
        )
        point_matrix_cache = getattr(self, "_point_op_by_matrix_cache", None)
        if point_matrix_cache is None:
            point_matrix_cache = {}
            setattr(self, "_point_op_by_matrix_cache", point_matrix_cache)
        point_by_matrix = point_matrix_cache.get(point_matrix_cache_key)
        if point_by_matrix is None:
            point_by_matrix = {
                self.point_operation_matrix(int(record[4])): int(record[4])
                for record in space_operations
            }
            point_matrix_cache[point_matrix_cache_key] = point_by_matrix
        point_op = point_by_matrix.get(product_matrix)
        if point_op is None:
            raise KeyError(f"point-operation product not found: {left_record} * {right_record}")
        denominator = 1
        for value in product_t:
            denominator = math.lcm(denominator, value.denominator)
        numerators = [int(value * denominator) for value in product_t]
        x, y, z, den = _reduce_fraction_vector((numerators[0], numerators[1], numerators[2], denominator))
        return (x, y, z, den, point_op)

    def operation_pair_multiply(
        self,
        *,
        space_operations: Sequence[tuple[int, int, int, int, int]],
        kernel_basis: tuple[int, ...] | list[int],
        kernel_fractions: Sequence[tuple[int, int, int, int]],
        lattice: int = 1,
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> tuple[int, int]:
        """Return the ``iso_mlt_ops_`` membership key for one pair product.

        ``get_isotropy_`` keeps operation membership in a table keyed by
        ``(ispace_index, kernel_fraction_index)``.  Matrix evaluation still
        uses ``vadd(kernel_fraction, ispace_op)``; this helper is only for the
        closure/seen table.
        """

        cache_key = (
            self._sequence_key(space_operations, cache_name="_space_operation_sequence_key_cache"),
            tuple(int(value) for value in kernel_basis),
            self._sequence_key(kernel_fractions, cache_name="_kernel_fraction_sequence_key_cache"),
            int(lattice),
            (int(left[0]), int(left[1])),
            (int(right[0]), int(right[1])),
        )
        cache = getattr(self, "_operation_pair_multiply_cache", None)
        if cache is None:
            cache = {}
            setattr(self, "_operation_pair_multiply_cache", cache)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached

        left_space = int(left[0])
        left_fraction = int(left[1])
        right_space = int(right[0])
        right_fraction = int(right[1])

        left_record = self._vadd_fraction_operation(
            kernel_fractions[left_fraction - 1],
            space_operations[left_space - 1],
        )
        right_record = self._vadd_fraction_operation(
            kernel_fractions[right_fraction - 1],
            space_operations[right_space - 1],
        )
        product_record = self._compose_operation_records(
            left_record,
            right_record,
            space_operations,
            lattice=int(lattice),
        )
        result = self._operation_pair_for_record(
            space_operations=space_operations,
            kernel_basis=kernel_basis,
            kernel_fractions=kernel_fractions,
            record=product_record,
        )
        cache[cache_key] = result
        return result

    def orderparam_to_subgroup(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        orderparam: tuple[float, ...] | list[float],
        row_count: int,
    ) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
        """Evaluate the nonmagnetic ``orderparam_to_subgroup_`` loop."""

        gid = int(gid)
        dim = int(self.little["little_irr_full_dim"][gid - 1])
        sg = int(self.little["little_irr_space_group"][gid - 1])
        candidates = list(self.generate_space_group_records(sg))
        basis: tuple[int, ...] | None = None
        fractions: tuple[tuple[int, int, int, int], ...] = ()

        for row_index in range(int(row_count)):
            row = tuple(float(orderparam[row_index * 48 + col]) for col in range(dim))
            if row_index == 0:
                basis = self.newlat_order(gid, kparam, row)
                fractions = get_new_fractionals(basis)
            else:
                if basis is None:
                    raise ValueError("missing previous basis")
                next_basis = self.newlat_order3(gid, kparam, row, basis)
                transform = self._basis_transform_matrix(basis, next_basis)
                fractions = tuple(
                    self._reduce_fraction_parent(basis, fraction)
                    for fraction in get_new_fractionals(transform)
                )
                basis = next_basis

            next_candidates: list[tuple[int, int, int, int, int]] = []
            for operation in candidates:
                for fraction in fractions:
                    shifted = self._vadd_fraction_operation(fraction, operation)
                    matrix = self.get_irreps_matrix_for_record(gid, shifted, kparam)
                    if orderparam_check(dim, 1, self._matrix_to_stride48(matrix), row):
                        next_candidates.append(shifted)
                        break
            candidates = next_candidates
        if basis is None:
            basis = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        return basis, tuple(candidates)

    def orderparam_to_subgroup_magnetic(
        self,
        gid: int,
        kparam: tuple[int, ...] | list[int],
        orderparam: tuple[float, ...] | list[float],
        row_count: int,
    ) -> tuple[tuple[int, ...], tuple[tuple[int, int, int, int, int], ...]]:
        """Evaluate the magnetic ``orderparam_to_subgroup_magnetic_`` loop."""

        gid = int(gid)
        dim = int(self.little["little_irr_full_dim"][gid - 1])
        magnetic_group = self.magnetic_orderparam_group_for_irrep(gid)
        candidates = list(self.generate_magnetic_space_group_records(magnetic_group))
        basis: tuple[int, ...] | None = None
        fractions: tuple[tuple[int, int, int, int], ...] = ()

        for row_index in range(int(row_count)):
            row = tuple(float(orderparam[row_index * 48 + col]) for col in range(dim))
            if row_index == 0:
                basis = self.newlat_order(gid, kparam, row)
                fractions = get_new_fractionals(basis)
            else:
                if basis is None:
                    raise ValueError("missing previous basis")
                next_basis = self.newlat_order3(gid, kparam, row, basis)
                transform = self._basis_transform_matrix(basis, next_basis)
                fractions = tuple(
                    self._reduce_fraction_parent(basis, fraction)
                    for fraction in get_new_fractionals(transform)
                )
                basis = next_basis

            next_candidates: list[tuple[int, int, int, int, int]] = []
            for operation in candidates:
                for fraction in fractions:
                    shifted = self._vadd_fraction_operation(fraction, operation)
                    matrix = self.get_irreps_matrix_for_magnetic_record(gid, shifted, kparam)
                    if orderparam_check(dim, 1, self._matrix_to_stride48(matrix), row):
                        next_candidates.append(shifted)
                        break
            candidates = next_candidates
        if basis is None:
            basis = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        return basis, tuple(candidates)

    @staticmethod
    def _reduce_fraction_parent(
        basis: tuple[int, ...] | list[int],
        fraction: tuple[int, int, int, int],
    ) -> tuple[int, int, int, int]:
        x, y, z, den = (int(value) for value in fraction)
        values = [
            x * int(basis[0 + axis]) + y * int(basis[3 + axis]) + z * int(basis[6 + axis])
            for axis in range(3)
        ]
        factor = 0
        for value in values + [den]:
            factor = math.gcd(factor, abs(value))
        if factor == 0:
            factor = 1
        return (values[0] // factor, values[1] // factor, values[2] // factor, den // factor)
