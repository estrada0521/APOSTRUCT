"""K-vector and find_isotropy_subgroup_ input helpers."""

from __future__ import annotations

from fractions import Fraction
import math
from ISODISTORT.Assembled.Backend.modes.engine.input import Case

class IsotropyKVectorMixin:
    def k_slot_for_label(self, sg: int, label: str) -> int:
        """Return the Miller-Love k-slot for a displayed k label."""

        lattice = int(self.iso.space["ispace_lattice"][sg - 1])
        count = int(self.iso.little["little_k_count"][lattice - 1])
        matches = [
            kslot
            for kslot in range(1, count + 1)
            if self.iso.little["little_k_label"][(lattice - 1) * 27 + kslot - 1].strip() == label
        ]
        if len(matches) != 1:
            raise KeyError(f"expected one k label {label!r} for SG{sg}, got {matches}")
        return matches[0]

    def k_vector_from_case(self, case: Case) -> tuple[Fraction, Fraction, Fraction]:
        """Evaluate the conventional fractional k vector for a mode-kernel case.

        `data_little:little_k` stores four rational 3-vectors:
        `base + p1*v1 + p2*v2 + p3*v3`.  Fixed special k labels have
        `little_k_dim == 0` and use only the base row.
        """

        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        dim = int(self.iso.little["little_k_dim"][slot])
        if len(case.k_params) < dim:
            raise ValueError(
                f"k label {case.k_label} for SG{case.sg} expects {dim} parameter(s), got {len(case.k_params)}"
            )
        raw = tuple(int(x) for x in self.iso.little["little_k"][slot * 16:(slot + 1) * 16])
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for start in range(0, 16, 4):
            den = raw[start + 3]
            if den == 0:
                vectors.append((Fraction(0), Fraction(0), Fraction(0)))
            else:
                vectors.append(tuple(Fraction(raw[start + axis], den) for axis in range(3)))  # type: ignore[arg-type]
        out = [vectors[0][axis] for axis in range(3)]
        for param_index in range(dim):
            param = Fraction(case.k_params[param_index])
            for axis in range(3):
                out[axis] += param * vectors[param_index + 1][axis]
        return tuple(value % 1 for value in out)  # type: ignore[return-value]

    def _first_little_gid_for_case(self, case: Case) -> int:
        """Return the first little-irrep row for the case's SG/k slot.

        The final stdout k-vector line is printed once per k selection before
        the irrep loop.  The binary uses the k-star attached to that selection;
        all little-irrep rows for the same SG/k share the same first star arm,
        so the first row is sufficient and avoids depending on a specific irrep.
        """

        kslot = self.k_slot_for_label(case.sg, case.k_label)
        for gid, sg in enumerate(self.iso.little["little_irr_space_group"], start=1):
            if int(sg) == int(case.sg) and int(self.iso.little["little_irr_k"][gid - 1]) == kslot:
                return gid
        raise KeyError(f"no little-irrep row for SG{case.sg} k={case.k_label}")

    def ml_star_lead_k_vector_from_case(self, case: Case) -> tuple[Fraction, Fraction, Fraction]:
        """Return the first Miller-Love star arm used as stdout k-vector source."""

        gid = self._first_little_gid_for_case(case)
        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        if int(self.iso.little["little_k_dim"][slot]) > 0:
            vectors = self.little_k_star_vectors_for_case(gid, case)
        else:
            vectors = self.little_k_vectors_by_gid(gid).vectors
        if not vectors:
            raise ValueError(f"empty k-star for SG{case.sg} k={case.k_label}")
        return tuple(vectors[0])  # type: ignore[return-value]

    def display_k_vector_from_case(self, case: Case) -> tuple[Fraction, Fraction, Fraction]:
        """Return the k-vector coordinates printed by mode-kernel stdout.

        ``mode kernel.txt`` states that input labels follow Miller-Love, while the
        printed coordinates are expressed in the reciprocal basis of the
        conventional International Tables lattice.  The binary's stdout line is
        reproduced by taking the first Miller-Love k-star arm and applying the
        reciprocal-coordinate ``pml -> cinter`` setting transform.  Values are
        intentionally not reduced modulo reciprocal lattice vectors; the binary
        prints entries such as ``1`` and ``8/5``.
        """

        return self.reciprocal_setting_change_vector(
            case.sg,
            "pml",
            "cinter",
            self.ml_star_lead_k_vector_from_case(case),
        )

    def raw_k_vector_from_case(self, case: Case) -> tuple[Fraction, Fraction, Fraction]:
        """Evaluate k without reducing components modulo reciprocal vectors."""

        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        dim = int(self.iso.little["little_k_dim"][slot])
        if len(case.k_params) < dim:
            raise ValueError(
                f"k label {case.k_label} for SG{case.sg} expects {dim} parameter(s), got {len(case.k_params)}"
            )
        vectors = self.little_k_basis_vectors_for_case(case)
        out = [vectors[0][axis] for axis in range(3)]
        for param_index in range(dim):
            param = Fraction(case.k_params[param_index])
            for axis in range(3):
                out[axis] += param * vectors[param_index + 1][axis]
        return tuple(out)  # type: ignore[return-value]

    def little_k_basis_vectors_for_case(
        self,
        case: Case,
    ) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
        """Return the base and parameter-direction k vectors for a case.

        `data_little:little_k` stores four records:
        `base, direction_1, direction_2, direction_3`.  This helper exposes the
        raw rational vectors before substituting the runtime k parameters.
        """

        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        raw = tuple(int(x) for x in self.iso.little["little_k"][slot * 16:(slot + 1) * 16])
        vectors: list[tuple[Fraction, Fraction, Fraction]] = []
        for start in range(0, 16, 4):
            den = raw[start + 3]
            if den == 0:
                vectors.append((Fraction(0), Fraction(0), Fraction(0)))
            else:
                vectors.append(tuple(Fraction(raw[start + axis], den) for axis in range(3)))  # type: ignore[arg-type]
        return tuple(vectors)

    def newlat_input_records_for_case(
        self,
        case: Case,
    ) -> tuple[tuple[int, int, int, int], ...]:
        """Build the four 4-int records passed from `find_isotropy_subgroup_` to `newlat_`.

        For parametric k labels, the first record is the evaluated k vector
        modulo reciprocal lattice vectors.  Following records are the individual
        parameter contributions `p_i * direction_i` without modulo reduction.
        The binary passes exactly four records and normalizes unused zero rows
        as `(0,0,0,1)`.
        """

        vectors = self.little_k_basis_vectors_for_case(case)
        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        dim = int(self.iso.little["little_k_dim"][slot])
        if len(case.k_params) < dim:
            raise ValueError(
                f"k label {case.k_label} for SG{case.sg} expects {dim} parameter(s), got {len(case.k_params)}"
            )
        evaluated = self.raw_k_vector_from_case(case)
        rows: list[tuple[int, int, int, int]] = [self._fraction_vector_to_record(evaluated)]
        for param_index in range(dim):
            contribution = tuple(
                Fraction(case.k_params[param_index]) * vectors[param_index + 1][axis]
                for axis in range(3)
            )
            rows.append(self._fraction_vector_to_record(contribution))
        while len(rows) < 4:
            rows.append((0, 0, 0, 1))
        return tuple(rows[:4])

    def find_isotropy_orderparam_records_for_case(
        self,
        case: Case,
    ) -> tuple[tuple[int, int, int, int], ...]:
        """Build the k-vector records passed to ``find_isotropy_subgroup_``.

        The public ``newlat_input_records_for_case`` helper exposes the compact
        records that ``find_isotropy_subgroup_`` constructs internally before
        calling ``newlat_``.  The routine entry receives a larger
        order-parameter buffer instead.  For parametric k, its first record is
        the sum of the runtime parameter contributions transformed by the
        per-SG/k ``little_k_star_conv2ml`` 4x4 row-vector matrix.  The following
        records retain the evaluated k and individual parameter contributions
        in the Miller-Love basis.  This mirrors the entry buffer observed at
        the binary boundary and keeps the upstream data available for the
        surrounding subgroup-selection port.
        """

        lattice = int(self.iso.space["ispace_lattice"][case.sg - 1])
        kslot = self.k_slot_for_label(case.sg, case.k_label)
        slot = (lattice - 1) * 27 + kslot - 1
        dim = int(self.iso.little["little_k_dim"][slot])
        if dim <= 0:
            return ()
        if len(case.k_params) < dim:
            raise ValueError(
                f"k label {case.k_label} for SG{case.sg} expects {dim} parameter(s), got {len(case.k_params)}"
            )

        sg_slot = (case.sg - 1) * 27 + kslot - 1
        pointer = int(self.iso.little["little_k_star_conv2ml_pointer"][sg_slot])
        if pointer <= 0:
            raise ValueError(f"missing little_k_star_conv2ml pointer for SG{case.sg} {case.k_label}")
        raw = tuple(int(x) for x in self.iso.little["little_k_star"][16 * (pointer - 1):16 * pointer])
        if len(raw) != 16:
            raise ValueError(f"little_k_star_conv2ml pointer out of range: SG{case.sg} {case.k_label}")
        matrix = tuple(
            tuple(Fraction(raw[4 * row + col]) for col in range(4))
            for row in range(4)
        )

        newlat_rows = self.newlat_input_records_for_case(case)
        transformed = [Fraction(0), Fraction(0), Fraction(0), Fraction(0)]
        for param_record in newlat_rows[1:1 + dim]:
            vector = (
                Fraction(param_record[0], param_record[3]),
                Fraction(param_record[1], param_record[3]),
                Fraction(param_record[2], param_record[3]),
                Fraction(0),
            )
            for col in range(4):
                transformed[col] += sum(vector[row] * matrix[row][col] for row in range(4))
        scale = transformed[3] if transformed[3] else Fraction(1)
        first = self._fraction_vector_to_record(transformed[col] / scale for col in range(3))

        zero = (0, 0, 0, 0)
        records: list[tuple[int, int, int, int]] = [
            first,
            newlat_rows[0],
            zero,
            newlat_rows[0],
        ]
        records.extend(newlat_rows[1:1 + dim])
        while len(records) < 8:
            records.append(zero)
        return tuple(records)

    def find_isotropy_initial_basis_for_case(self, case: Case) -> tuple[int, ...]:
        """Return the 3x3 pml basis produced by the initial `find_isotropy_subgroup_` path.

        In the decoded parametric-k branch, `find_isotropy_subgroup_` first
        builds the four `newlat_` records from `little_k` and runtime k
        parameters, calls `newlat_` with mode count 1, then passes the resulting
        basis into `subgroup_to_orderparam_`.  Current GDB assets show the
        returned basis remains this initial `newlat_` basis; later comparisons
        should catch any case where the fallback search replaces it.
        """

        return self.newlat_basis(1, self.newlat_input_records_for_case(case))

    def find_isotropy_candidate_operation_records_for_case(
        self,
        case: Case,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return candidate subgroup operations in ``find_isotropy_subgroup_`` order.

        After the initial parametric-k superlattice basis is built, the binary
        enumerates internal supercell translations from that basis and adds
        each of them to every non-identity generated space operation.  The
        resulting records are the trial operations accumulated in
        ``subgroup_to_orderparam_``'s final argument.
        """

        basis = self.find_isotropy_initial_basis_for_case(case)
        translations = self.get_new_fractionals_from_basis_columns(basis)
        candidates: list[tuple[int, int, int, int, int]] = []
        for space_record in self.generate_space_group_records(case.sg)[1:]:
            sx, sy, sz, sden, point_op = (int(value) for value in space_record)
            tau = (Fraction(sx, sden), Fraction(sy, sden), Fraction(sz, sden))
            for tx, ty, tz, tden in translations:
                shift = (Fraction(tx, tden), Fraction(ty, tden), Fraction(tz, tden))
                vx, vy, vz, den = self._fraction_vector_to_record(
                    tau[axis] + shift[axis] for axis in range(3)
                )
                candidates.append((vx, vy, vz, den, point_op))
        return tuple(candidates)

    def k_supercell_index(self, case: Case) -> int:
        """Return the minimal diagonal-index needed to fold k to Gamma.

        Here we only need the index, not a particular supercell basis:
        it is the least common multiple of denominators of the evaluated
        conventional fractional k-vector components.
        """

        index = 1
        for value in self.k_vector_from_case(case):
            index = math.lcm(index, value.denominator)
        return index
