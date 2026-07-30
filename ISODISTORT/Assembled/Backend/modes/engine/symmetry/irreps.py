from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
import math
import operator
from typing import Iterable, overload

import numpy as np

from ISODISTORT.Assembled.Backend.source.representation import decode_little_sparse_matrix, real_phase_operator
from ISODISTORT.Assembled.Backend.modes.engine.input import Case
from ISODISTORT.Assembled.Backend.modes.engine.records import (
    LittleKVectorRecord,
)


@dataclass(frozen=True, slots=True)
class BridgeWeightView(Sequence[float]):
    """Immutable view of a small matrix in the bridge's 48-stride layout."""

    _matrix: tuple[tuple[float, ...], ...]
    _output_length: int
    _dense_fallback: tuple[float, ...] | None

    @classmethod
    def from_matrix(
        cls,
        matrix: np.ndarray,
        *,
        output_length: int,
    ) -> BridgeWeightView:
        values = np.asarray(matrix)
        if values.ndim != 2:
            raise ValueError("bridge weight matrix must be two-dimensional")
        if np.max(np.abs(values.imag)) > 1e-10:
            raise ValueError("bridge weight matrix has non-real entries")
        frozen = tuple(
            tuple(float(value) for value in row)
            for row in values.real
        )
        length = max(0, int(output_length))
        if len(frozen) > 48:
            dense = [0.0] * length
            for row, row_values in enumerate(frozen):
                for column, value in enumerate(row_values):
                    index = row + column * 48
                    if index < length:
                        dense[index] = value
            return cls((), length, tuple(dense))
        return cls(frozen, length, None)

    def __len__(self) -> int:
        return self._output_length

    @overload
    def __getitem__(self, index: int) -> float: ...

    @overload
    def __getitem__(self, index: slice) -> list[float]: ...

    def __getitem__(self, index: int | slice) -> float | list[float]:
        if isinstance(index, slice):
            return [self[offset] for offset in range(*index.indices(len(self)))]
        offset = operator.index(index)
        if offset < 0:
            offset += len(self)
        if offset < 0 or offset >= len(self):
            raise IndexError("bridge weight index out of range")
        if self._dense_fallback is not None:
            return self._dense_fallback[offset]
        row = offset % 48
        column = offset // 48
        if row >= len(self._matrix):
            return 0.0
        values = self._matrix[row]
        return values[column] if column < len(values) else 0.0

class IrrepMixin:
    def generated_space_point_op_positions(self, sg: int) -> dict[int, int]:
        """Return the point-op -> internal position map used by `get_irrep4_ssgn_`.

        `get_irrep4_ssgn_` calls `generate_space_group_` and builds a `mapop`
        array from the generated operation records.  The sparse
        `little_irr_full_matrices_*` tables are indexed by that internal
        position, not by the user-facing operation order used by
        `irrep_matrices_with_ops`.
        """

        pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
        count = len(self.iso.space_ops(sg))
        raw = self.iso.space["ispace_elements"][(pointer - 1) * 5:(pointer - 1 + count) * 5]
        out: dict[int, int] = {}
        for index in range(count):
            record = raw[5 * index:5 * index + 5]
            if len(record) == 5:
                # The Fortran map is overwritten if a point op appears more
                # than once; mirror that behavior.
                out[int(record[4])] = index + 1
        return out

    def generated_space_operation_position(self, sg: int, record: tuple[int, int, int, int, int]) -> int:
        """Return the generated operation position matching a full 5-int record.

        Parent space groups can have several operations with the same point-op
        label but different centering/screw/glide translations.  The parent
        `get_irrep4_ssgn_` matrix lookup is keyed by the generated operation
        position, so point-op-only lookup is insufficient for centered groups.
        """

        tau = tuple(Fraction(int(record[i]), int(record[3])) % 1 for i in range(3))
        point_op = int(record[4])
        pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
        count = len(self.iso.space_ops(sg))
        for index in range(count):
            x, y, z, den, op = self.iso.space["ispace_elements"][
                5 * (pointer - 1 + index):5 * (pointer - 1 + index + 1)
            ]
            if int(op) != point_op:
                continue
            candidate = (Fraction(int(x), int(den)) % 1, Fraction(int(y), int(den)) % 1, Fraction(int(z), int(den)) % 1)
            if candidate == tau:
                return index + 1
        raise KeyError(f"operation record {record} not found in generated SG{sg} records")

    def space_operation_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return generated space-operation records as native 5-int tuples."""

        pointer = int(self.iso.space["ispace_elements_pointer"][sg - 1])
        count = len(self.iso.space_ops(sg))
        out: list[tuple[int, int, int, int, int]] = []
        for index in range(count):
            raw = self.iso.space["ispace_elements"][
                5 * (pointer - 1 + index):5 * (pointer - 1 + index + 1)
            ]
            if len(raw) == 5:
                out.append(tuple(int(x) for x in raw))  # type: ignore[arg-type]
        return tuple(out)

    def generated_space_tau(self, sg: int, point_op: int) -> tuple[Fraction, Fraction, Fraction]:
        """Return the generated operation translation for a point operation.

        This mirrors the helper used by the parent `get_irrep4_` phase
        verifier.  If a point operation appears more than once, the generated
        table's first matching operation is the one used for the sparse matrix
        representative; extra translation is carried as a Bloch phase.
        """

        for x, y, z, den, op in self.space_operation_records(sg):
            if int(op) == int(point_op):
                return (Fraction(x, den), Fraction(y, den), Fraction(z, den))
        raise KeyError(f"point op {point_op} not found in generated SG{sg} records")

    def operation_record_phases(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
    ) -> tuple[Fraction, ...]:
        """Return the k-star phase vector used by `get_irrep4_` for a record."""

        k_record = self.little_k_vectors_by_gid(gid)
        request_tau = tuple(Fraction(int(record[i]), int(record[3])) for i in range(3))
        generated_tau = self.generated_space_tau(k_record.sg, int(record[4]))
        delta = tuple(request_tau[i] - generated_tau[i] for i in range(3))
        return tuple(sum(kvec[i] * delta[i] for i in range(3)) for kvec in k_record.vectors)

    def little_sparse_matrix_by_gid_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
    ) -> np.ndarray:
        """Read the sparse matrix table for a full parent operation record."""

        dim = self.little_record_by_gid(gid).full_dim
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        position = self.generated_space_operation_position(sg, record)
        return decode_little_sparse_matrix(
            self.iso.little, self.iso.const, gid=gid, dim=dim, position=position
        )

    def little_sparse_matrix_by_gid_full_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
    ) -> np.ndarray:
        """Resolve a sparse matrix while retaining a noncanonical full record.

        Generated Source records remain the primary matrix authority because
        centered groups may contain several records with the same point-op.
        A representative transport can legitimately produce a full
        translation absent from that table; only then is the canonical matrix
        selected by point-op while the caller retains the full translation for
        its Bloch phase.
        """

        try:
            return self.little_sparse_matrix_by_gid_record(int(gid), record)
        except KeyError:
            sg = int(self.iso.little["little_irr_space_group"][int(gid) - 1])
            canonical = next(
                (
                    candidate
                    for candidate in self.generate_space_group_records(sg)
                    if int(candidate[4]) == int(record[4])
                ),
                None,
            )
            if canonical is None:
                raise KeyError(
                    f"point op {record[4]} not found in generated SG{sg} records"
                ) from None
            return self.little_sparse_matrix_by_gid_record(int(gid), canonical)

    def little_sparse_matrix_by_gid_op(self, gid: int, point_op: int) -> np.ndarray:
        """Read the exact sparse matrix table used by `get_irrep4_ssgn_`."""

        dim = self.little_record_by_gid(gid).full_dim
        # The internal positions are keyed by the irrep's actual space group,
        # not by the native little row id.
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        positions = self.generated_space_point_op_positions(sg)
        position = positions[int(point_op)]
        return decode_little_sparse_matrix(
            self.iso.little, self.iso.const, gid=gid, dim=dim, position=position
        )

    def phase_operator(self, gid: int, phases: Iterable[Fraction]) -> np.ndarray:
        """Return the real row-side phase operator used by ``get_irrep4_ssgn_``."""

        phase_values = tuple(phases)
        dim = self.little_record_by_gid(gid).full_dim
        return real_phase_operator(dim, phase_values, gid=gid)

    def little_phase_matrix_by_gid_op(
        self,
        gid: int,
        point_op: int,
        phases: Iterable[Fraction],
    ) -> np.ndarray:
        """Return the phase-adjusted parent matrix from ``get_irrep4_ssgn_``."""

        return self.phase_operator(gid, phases) @ self.little_sparse_matrix_by_gid_op(gid, point_op)

    def little_phase_matrix_by_gid_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        phases: Iterable[Fraction],
    ) -> np.ndarray:
        """Return the phase-adjusted parent matrix for a full operation record."""

        return self.phase_operator(gid, phases) @ self.little_sparse_matrix_by_gid_full_record(gid, record)

    def site_little_sparse_matrices(
        self,
        site_pg: int,
        pg_irrep: int,
    ) -> dict[int, np.ndarray]:
        """Return site matrices from the same `data_little` path as mode kernel."""

        old_id = self.site_pg_irrep_old_id(site_pg, pg_irrep)
        gid = self.little_gid_for_old_id(old_id)
        return {
            int(op): self.little_sparse_matrix_by_gid_op(gid, int(op))
            for op in self.site_pg_element_settings(site_pg)[0]
        }

    def site_irrep_operation_order(self, site_pg: int, pg_irrep: int) -> tuple[int, ...]:
        """Return the site-irrep operation order used by `get_irreps_`."""

        old_id = self.site_pg_irrep_old_id(site_pg, pg_irrep)
        sg = int(self.iso.irreps["irrep_space_group"][old_id - 1])
        label = self.iso.irreps["irrep_label"][old_id - 1].strip()
        elements = self.iso.irrep_matrices_with_ops(sg, label)
        return tuple(int(element.op_index) for element in elements[3:] if element.op_index is not None)

    def site_pg_element_settings(self, site_pg: int) -> list[tuple[int, ...]]:
        """Return candidate site point-group operation settings.

        Each setting is a row from ``iwyckoff_pg_elements`` with trailing zeros
        removed.  In ``project_`` the selected row is converted into 5-int
        operation records ``(0,0,0,1, point_op)`` before calling
        ``get_irreps_`` for the site point-group irrep.
        """

        return list(self.iso.site_pg_element_settings(int(site_pg)))

    def site_pg_project_records(self, site_pg: int, setting_index0: int = 0) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return 5-int records for one site point-group operation setting."""

        settings = self.site_pg_element_settings(site_pg)
        point_ops = settings[setting_index0]
        return tuple((0, 0, 0, 1, op) for op in point_ops)

    def little_parent_point_ops(self, sg: int, label: str) -> tuple[int, ...]:
        """Return the point operations in the little group for an irrep label."""

        little = self.little_record(sg, label)
        return self.little_parent_point_ops_by_gid(sg, little.gid)

    def little_parent_point_ops_by_gid(self, sg: int, gid: int) -> tuple[int, ...]:
        """Return little-group point operations using the native little row id."""

        little = self.little_record_by_gid(gid)
        slot0 = int(self.iso.little["little_irr_k"][little.gid - 1]) - 1
        slot = (sg - 1) * 27 + slot0
        count = int(self.iso.little["little_ops_count"][slot])
        return tuple(int(op) for op in self.iso.little["little_ops"][slot * 48:slot * 48 + count])

    def little_k_vectors_by_gid(self, gid: int) -> LittleKVectorRecord:
        """Return the k-vector records copied by ``get_irrep4_``.

        For fixed special-k rows, the relevant static table is ``little_k_star_special``, not the later
        ``little_kvec`` section.  ``read_little`` loads it as 2100 records of
        four integers ``(h, k, l, denominator)``.  ``get_irrep4_`` addresses it
        with the 1-based ``little_k_star_ml_pointer`` and copies ``nmod``
        consecutive records into the first column of the runtime ``kvec.1``
        common block.

        Parametric rows require runtime ``Case.k_params`` and are handled by
        ``little_k_star_vectors_for_case`` using the 16-int ``little_k_star``
        records.
        """

        gid = int(gid)
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        slot = (sg - 1) * 27 + kslot - 1
        nmod = int(self.iso.little["little_k_star_count"][slot])
        pointer = int(self.iso.little["little_k_star_ml_pointer"][slot])
        if nmod < 0:
            raise ValueError(f"negative k-star count for gid={gid}: {nmod}")
        if nmod and pointer <= 0:
            raise ValueError(f"missing little_k_star_ml_pointer for gid={gid} SG{sg} kslot={kslot}")

        raw = self.iso.little["little_k_star_special"]
        records: list[tuple[int, int, int, int]] = []
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for offset in range(nmod):
            start = 4 * (pointer - 1 + offset)
            item = tuple(int(x) for x in raw[start:start + 4])
            if len(item) != 4:
                raise ValueError(f"k-vector pointer out of range for gid={gid}: pointer={pointer} nmod={nmod}")
            den = item[3]
            if den == 0:
                raise ValueError(f"zero k-vector denominator for gid={gid}: pointer={pointer} offset={offset}")
            records.append(item)  # type: ignore[arg-type]
            vectors.append((Fraction(item[0], den), Fraction(item[1], den), Fraction(item[2], den)))

        return LittleKVectorRecord(
            gid=gid,
            sg=sg,
            kslot=kslot,
            nmod=nmod,
            pointer=pointer,
            records=tuple(records),
            vectors=tuple(vectors),
        )

    def little_k_star_vectors_for_case(
        self,
        gid: int,
        case: Case,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Evaluate runtime k-star arms for a possibly parametric case.

        For parametric little irreps (`little_irr_old == 0`), the mode kernel uses
        `little_k_star` as 16-int records:
        `base, direction_1, direction_2, direction_3`.  Runtime k parameters
        are substituted per arm before phase construction.
        """

        gid = int(gid)
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        if sg != int(case.sg):
            raise ValueError(f"case SG{case.sg} does not match gid={gid} SG{sg}")
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        slot = (sg - 1) * 27 + kslot - 1
        nmod = int(self.iso.little["little_k_star_count"][slot])
        pointer = int(self.iso.little["little_k_star_ml_pointer"][slot])
        if nmod < 0:
            raise ValueError(f"negative k-star count for gid={gid}: {nmod}")
        if nmod and pointer <= 0:
            raise ValueError(f"missing little_k_star_ml_pointer for gid={gid} SG{sg} kslot={kslot}")

        raw = self.iso.little["little_k_star"]
        params = tuple(Fraction(value) for value in case.k_params)
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for arm in range(nmod):
            start = 16 * (pointer - 1 + arm)
            record = tuple(int(x) for x in raw[start:start + 16])
            if len(record) != 16:
                raise ValueError(f"k-star pointer out of range for gid={gid}: pointer={pointer} nmod={nmod}")
            pieces: list[tuple[Fraction, Fraction, Fraction]] = []
            for offset in range(0, 16, 4):
                den = record[offset + 3]
                if den == 0:
                    pieces.append((Fraction(0), Fraction(0), Fraction(0)))
                else:
                    pieces.append(tuple(Fraction(record[offset + axis], den) for axis in range(3)))  # type: ignore[arg-type]
            out = [pieces[0][axis] for axis in range(3)]
            for param_index, param in enumerate(params[:3]):
                for axis in range(3):
                    out[axis] += param * pieces[param_index + 1][axis]
            vectors.append(tuple(out))  # type: ignore[arg-type]
        return tuple(vectors)

    def operation_record_phases_for_case(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        case: Case,
    ) -> tuple[Fraction, ...]:
        """Return `get_irrep4_` phase vector using runtime case parameters."""

        k_vectors = self.little_k_star_vectors_for_case(gid, case)
        sg = int(self.iso.little["little_irr_space_group"][gid - 1])
        request_tau = tuple(Fraction(int(record[i]), int(record[3])) for i in range(3))
        kslot = int(self.iso.little["little_irr_k"][gid - 1])
        lattice = int(self.iso.space["ispace_lattice"][sg - 1])
        slot = (lattice - 1) * 27 + kslot - 1
        if int(self.iso.little["little_k_dim"][slot]) > 0:
            # Parametric-k rows follow the `get_irrep4_` branch that feeds the
            # request record directly into the phase calculation.  Fixed
            # special-k rows subtract the generated space-group representative.
            delta = request_tau
        else:
            generated_tau = self.generated_space_tau(sg, int(record[4]))
            delta = tuple(request_tau[i] - generated_tau[i] for i in range(3))
        return tuple(sum(kvec[i] * delta[i] for i in range(3)) for kvec in k_vectors)

    def little_phase_matrix_by_gid_record_for_case(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        case: Case,
    ) -> np.ndarray:
        """Return phase-adjusted matrix for a runtime mode-kernel case."""

        phases = self.operation_record_phases_for_case(gid, record, case)
        return self.phase_operator(gid, phases) @ self.little_sparse_matrix_by_gid_full_record(gid, record)

    @staticmethod
    def add_translation_to_operation_record(
        translation: tuple[int, int, int, int],
        operation: tuple[int, int, int, int, int],
    ) -> tuple[int, int, int, int, int]:
        """Add a fractional translation to a 5-int operation record.

        The supercell expansion supplies a 4-int fractional translation. The
        sum keeps the operation record's point-operation label for subsequent
        irrep matrix evaluation.
        """

        tx, ty, tz, tden = (int(value) for value in translation)
        ox, oy, oz, oden, op = (int(value) for value in operation)
        values = (
            Fraction(tx, tden) + Fraction(ox, oden),
            Fraction(ty, tden) + Fraction(oy, oden),
            Fraction(tz, tden) + Fraction(oz, oden),
        )
        den = 1
        for value in values:
            den = math.lcm(den, value.denominator)
        return tuple(int(value * den) for value in values) + (den, op)  # type: ignore[return-value]

    @lru_cache(maxsize=32768)
    def project_vector_bridge_weight_view_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        case: Case,
        *,
        output_length: int = 2304,
    ) -> BridgeWeightView:
        """Return the exact bridge weights without materializing zero padding."""

        matrix = self._bridge_irrep_matrix_for_record(gid, record, case)
        return BridgeWeightView.from_matrix(matrix, output_length=output_length)

    def _bridge_irrep_matrix_for_record(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        case: Case,
    ) -> np.ndarray:
        """Return the exact parent-irrep matrix before OPD contraction."""

        if case.k_params:
            return self.little_phase_matrix_by_gid_record_for_case(gid, record, case)
        return self.little_phase_matrix_by_gid_record(
            gid,
            record,
            self.operation_record_phases(gid, record),
        )

    def get_irreps_matrix_for_case(
        self,
        gid: int,
        record: tuple[int, int, int, int, int],
        case: Case,
    ) -> np.ndarray:
        """Evaluate the type-2 matrix branch of ``get_irreps_`` for nonmagnetic calls."""

        return self.little_phase_matrix_by_gid_record_for_case(gid, record, case)
