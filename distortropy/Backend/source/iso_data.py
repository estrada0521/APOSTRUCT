"""
Decoder for selected ISOTROPY/ISO-IR Source data tables.

This module decodes
  data_space.txt, data_irreps.txt, data_little.txt, data_wyckoff.txt, const.dat
directly from the repository Source files.

All data pointers in the ISOTROPY tables are Fortran 1-based and are converted
locally where used.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import math
import re
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np

from distortropy.Backend.exactmath import (
    fraction_identity3,
    fraction_matrix_inverse3,
    fraction_matrix_multiply3,
)

_SECTION = re.compile(r"^[a-z][a-z0-9_]*$")
_NUMBER = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")

Frac3 = Tuple[Fraction, Fraction, Fraction]
SpaceOp = Tuple[int, Frac3]  # 1-based point-op index, fractional translation tau
Matrix3 = tuple[tuple[Fraction, Fraction, Fraction], ...]


def pml_to_cml_matrix_from_table(
    space: Mapping[str, Sequence[Any]],
    sg: int,
) -> Matrix3:
    lattice = int(space["ispace_lattice"][int(sg) - 1])
    raw = space["lattice_ml"][(lattice - 1) * 36 + 9:(lattice - 1) * 36 + 18]
    denominator = int(space["lattice_ml_denom"][(lattice - 1) * 2])
    if denominator == 0:
        raise ValueError(f"zero lattice_ml pml->cml denominator for lattice {lattice}")
    return tuple(
        tuple(Fraction(int(raw[3 * row + column]), denominator) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def setting_to_cinter_affine_from_table(
    space: Mapping[str, Sequence[Any]],
    sg: int,
    setting: str,
    setting_id: int | None = None,
) -> tuple[Matrix3, Frac3]:
    """Decode a Source setting as an affine row transform to cinter."""

    name = str(setting).strip().lower()
    if name == "cinter":
        return fraction_identity3(), (Fraction(0), Fraction(0), Fraction(0))
    if name not in {"cml", "pml"}:
        raise KeyError(f"unsupported setting {name!r}")

    pml_to_cml = pml_to_cml_matrix_from_table(space, sg) if name == "pml" else None
    choice = int(setting_id or space["ispace_inter_choice"][int(sg) - 1])
    raw = tuple(int(value) for value in space["ispace_settings_inter"][(choice - 1) * 32:(choice - 1) * 32 + 16])
    if len(raw) != 16:
        raise IndexError(f"inter setting choice out of range: {choice}")
    denominator = raw[15]
    if denominator == 0:
        raise ValueError(f"zero ispace_settings_inter denominator for SG{sg} choice {choice}")
    matrix: Matrix3 = tuple(
        tuple(Fraction(raw[4 * row + column], denominator) for column in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    origin: Frac3 = tuple(Fraction(raw[12 + axis], denominator) for axis in range(3))  # type: ignore[assignment]
    if pml_to_cml is not None:
        matrix = fraction_matrix_multiply3(pml_to_cml, matrix)
    return matrix, origin


def reciprocal_to_cinter_matrix_from_table(
    space: Mapping[str, Sequence[Any]],
    sg: int,
    setting: str,
    setting_id: int | None = None,
) -> Matrix3:
    """Decode the reciprocal-coordinate transform to cinter."""

    inverse = fraction_matrix_inverse3(
        setting_to_cinter_affine_from_table(space, sg, setting, setting_id)[0]
    )
    return tuple(tuple(inverse[column][row] for column in range(3)) for row in range(3))  # type: ignore[return-value]


def setting_change_matrix_from_table(
    space: Mapping[str, Sequence[Any]],
    sg: int,
    from_setting: str,
    to_setting: str,
) -> Matrix3:
    """Decode the Source direct-space row transform between named settings."""

    cml_to_cinter = setting_to_cinter_affine_from_table(space, sg, "cml")[0]
    to_cinter = {
        "cinter": fraction_identity3(),
        "cml": cml_to_cinter,
        "pml": setting_to_cinter_affine_from_table(space, sg, "pml")[0],
    }
    source = str(from_setting).strip().lower()
    target = str(to_setting).strip().lower()
    if source not in to_cinter or target not in to_cinter:
        raise KeyError(f"unsupported setting transform {source!r}->{target!r}")
    return fraction_matrix_multiply3(
        to_cinter[source],
        fraction_matrix_inverse3(to_cinter[target]),
    )


def generate_space_group_records_from_table(
    space: Mapping[str, Sequence[Any]],
    sg: int,
) -> tuple[tuple[int, int, int, int, int], ...]:
    """Decode the Source operation records for one ordinary space group."""

    sg = int(sg)
    point_group = int(space["ispace_point_group"][sg - 1])
    count = int(space["ipoint_group_order"][point_group - 1])
    pointer = int(space["ispace_elements_pointer"][sg - 1])
    out: list[tuple[int, int, int, int, int]] = []
    for index in range(count):
        raw = space["ispace_elements"][5 * (pointer - 1 + index):5 * (pointer + index)]
        if len(raw) == 5:
            out.append(tuple(int(value) for value in raw))  # type: ignore[arg-type]
    return tuple(out)


def _split_sections(path: Path):
    sections = []
    name = None
    buf: List[str] = []
    for line in path.read_text(errors="replace").splitlines(True):
        s = line.strip()
        if _SECTION.match(s):
            if name is not None:
                sections.append((name, "".join(buf)))
            name = s
            buf = []
        else:
            if name is not None:
                buf.append(line)
    if name is not None:
        sections.append((name, "".join(buf)))
    return sections


def _parse_block(data: str):
    strings = re.findall(r'"([^"]*)"', data)
    rest = re.sub(r'"[^"]*"', " ", data)
    nums = _NUMBER.findall(rest)
    if strings and not nums:
        return strings
    if nums:
        if all("." not in x and "e" not in x.lower() for x in nums):
            return [int(x) for x in nums]
        return [float(x) for x in nums]
    tokens = rest.split()
    if tokens and all(token in {"T", "F"} for token in tokens):
        return [token == "T" for token in tokens]
    return strings


def _load_table(path: Path) -> Dict[str, list]:
    return {name: _parse_block(data) for name, data in _split_sections(path)}


def _load_constants(path: Path) -> Dict[int, float]:
    const: Dict[int, float] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        code = int(parts[0])
        value = float(parts[1])
        const[code] = value
        if code >= 2 and code % 2 == 0:
            const[code + 1] = -value
    return const


@dataclass(frozen=True)
class IrrepElement:
    op_index: int | None          # None for primitive translation generators
    op_label: str
    R: np.ndarray                 # 3x3 integer point operation
    tau: Frac3                    # fractional translation of coset rep
    D: np.ndarray                 # complex representation matrix


class ISOData:
    TABLES = ("space", "irreps", "little", "isotropy", "wyckoff", "images", "const")

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        # Single source of truth: the distributed data_ text (parsed each run;
        # <1s). No pickle cache — data_ is the only authority.
        self.space = _load_table(self._find("data_space.txt", "data_space(1).txt"))
        self.irreps = _load_table(self._find("data_irreps.txt", "data_irreps(1).txt"))
        self.little = _load_table(self._find("data_little.txt"))
        self.isotropy = _load_table(self._find("data_isotropy.txt"))
        self.wyckoff = _load_table(self._find("data_wyckoff.txt", "data_wyckoff(1).txt"))
        self.images = _load_table(self._find("data_images.txt", "data_images(1).txt"))
        self.const = _load_constants(self._find("const.dat", "const(1).dat"))

        self.point_labels = [s.strip() for s in self.space["point_op_label"]]
        self.point_ops = [
            np.array(self.space["ipoint_op"][9 * i:9 * (i + 1)], dtype=int).reshape(3, 3)
            for i in range(72)
        ]
        self._matrix_to_point_indices: Dict[Tuple[int, ...], List[int]] = {}
        for i, R in enumerate(self.point_ops, start=1):
            self._matrix_to_point_indices.setdefault(tuple(R.flatten()), []).append(i)

    def _find(self, *names: str) -> Path:
        for name in names:
            p = self.data_dir / name
            if p.exists():
                return p
        raise FileNotFoundError(f"none of these files exists in {self.data_dir}: {names}")

    # ------------------------------------------------------------------
    # Basic data_irreps path
    # ------------------------------------------------------------------
    def old_irrep_id(self, sg: int, label: str) -> int:
        label = label.strip()
        start = self.irreps["irrep_space_group_pointer"][sg - 1]
        end_excl = self.irreps["irrep_space_group_pointer"][sg]
        for oid in range(start, end_excl):
            if self.irreps["irrep_label"][oid - 1].strip() == label:
                return oid
        raise KeyError(f"irrep label {label!r} not found for SG{sg}")

    def irrep_dim(self, old_id: int) -> int:
        # The reliable dimension table is in data_little, keyed through little_irr_old.
        for i, oid in enumerate(self.little["little_irr_old"]):
            if oid == old_id:
                return int(self.little["little_irr_full_dim"][i])
        raise KeyError(f"dimension for old irrep id {old_id} not found")

    def irrep_matrices(self, sg: int, label: str) -> List[np.ndarray]:
        old_id = self.old_irrep_id(sg, label)
        dim = self.irrep_dim(old_id)
        ptr_start = self.irreps["irrep_elements_pointer"][old_id - 1]
        if old_id < len(self.irreps["irrep_elements_pointer"]):
            ptr_end = self.irreps["irrep_elements_pointer"][old_id]
        else:
            # The final irrep has no `irrep_elements_pointer` sentinel.  The
            # corresponding matrix-pointer table is zero padded, so use the
            # first zero matrix pointer as the exclusive end.
            matrix_pointer_table = self.irreps["irrep_matrices_pointer"]
            ptr_end = ptr_start
            while ptr_end <= len(matrix_pointer_table) and int(matrix_pointer_table[ptr_end - 1]) > 0:
                ptr_end += 1
        matrix_ptrs = self.irreps["irrep_matrices_pointer"][ptr_start - 1:ptr_end - 1]
        mats = []
        for p in matrix_ptrs:
            raw = self.irreps["irrep_matrices"][p - 1:p - 1 + dim * dim]
            mats.append(np.array([self.const[x] for x in raw], dtype=complex).reshape(dim, dim))
        return mats

    # ------------------------------------------------------------------
    # data_space: generate finite coset representatives {R|tau}
    # ------------------------------------------------------------------
    def _point_index_from_matrix(self, R: np.ndarray, sg: int) -> int:
        candidates = self._matrix_to_point_indices[tuple(np.array(R, dtype=int).flatten())]
        lattice = self.space["ispace_lattice"][sg - 1]
        # The first 48 labels are the orthorhombic/tetragonal/cubic set.
        # The final 24 labels are the hexagonal/trigonal set.  Some matrices
        # are duplicated between these two sets, so choose the appropriate family.
        if lattice in (10, 11):  # HR, HP in data_space lattice list
            good = [i for i in candidates if i >= 49]
        else:
            good = [i for i in candidates if i <= 48]
        return good[0] if good else candidates[0]

    def _compose(self, g: SpaceOp, h: SpaceOp, sg: int) -> SpaceOp:
        # Active convention: {Rg|tg} after {Rh|th} = {Rg Rh | Rg th + tg}.
        ig, tg = g
        ih, th = h
        Rg = self.point_ops[ig - 1]
        Rh = self.point_ops[ih - 1]
        R = Rg @ Rh
        idx = self._point_index_from_matrix(R, sg)
        tau = tuple(
            (sum(Fraction(int(Rg[i, j])) * th[j] for j in range(3)) + tg[i]) % 1
            for i in range(3)
        )
        return idx, tau  # type: ignore[return-value]

    def space_generators_conv(self, sg: int) -> List[SpaceOp]:
        choice = self.space["ispace_inter_choice"][sg - 1]  # 1-based inter-setting id
        count = self.space["ispace_generators_conv_count"][choice - 1]
        raw = self.space["ispace_generators_conv"][20 * (choice - 1):20 * choice]
        gens = []
        for k in range(count):
            x, y, z, den, op = raw[5 * k:5 * k + 5]
            gens.append((op, (Fraction(x, den), Fraction(y, den), Fraction(z, den))))
        return gens

    def space_ops(self, sg: int) -> List[SpaceOp]:
        """Return finite coset representatives derived from static ``data_space``.

        The reps come from the generator closure of the conventional
        space-group generators (``_space_ops_from_generators``); no separate
        precomputed coset table is consulted.

        Known limitation: the mod-1 closure over-generates for centered
        lattices (it counts centering translations / high-order screw powers as
        distinct cosets), so for centered space groups ``len(...)`` can exceed
        ``ipoint_group_order``.  The authoritative coset count is the point
        group order in ``data_space`` (see ``generate_space_group_records``).
        """
        ops_cache = getattr(self, "_space_ops_cache", None)
        if ops_cache is None:
            ops_cache = {}
            object.__setattr__(self, "_space_ops_cache", ops_cache)
        if sg not in ops_cache:
            ops_cache[sg] = self._space_ops_from_generators(sg)
        return ops_cache[sg]

    def _space_ops_from_generators(self, sg: int) -> List[SpaceOp]:
        lattice = self.space["ispace_lattice"][sg - 1]
        identity_op = 49 if lattice in (10, 11) else 1
        E: SpaceOp = (identity_op, (Fraction(0), Fraction(0), Fraction(0)))
        gens = self.space_generators_conv(sg)
        seen = {E}
        queue = [E]
        while queue:
            a = queue.pop(0)
            for g in gens:
                for prod in (self._compose(g, a, sg), self._compose(a, g, sg)):
                    if prod not in seen:
                        seen.add(prod)
                        queue.append(prod)
        return sorted(seen, key=lambda x: (x[0], x[1]))

    # ------------------------------------------------------------------
    # Source operation and Wyckoff geometry authority
    # ------------------------------------------------------------------
    @lru_cache(maxsize=None)
    def generate_space_group_records(self, sg: int) -> tuple[tuple[int, int, int, int, int], ...]:
        return generate_space_group_records_from_table(self.space, int(sg))

    @lru_cache(maxsize=None)
    def site_pg_element_settings(self, site_pg: int) -> tuple[tuple[int, ...], ...]:
        site_pg = int(site_pg)
        pointer = int(self.wyckoff["iwyckoff_pg_elements_pointer"][site_pg - 1])
        count = int(self.wyckoff["iwyckoff_pg_settings_count"][site_pg - 1])
        settings: list[tuple[int, ...]] = []
        for index in range(count):
            start = (pointer - 1 + index) * 48
            row = tuple(int(value) for value in self.wyckoff["iwyckoff_pg_elements"][start:start + 48])
            settings.append(tuple(value for value in row if value))
        return tuple(settings)

    @lru_cache(maxsize=None)
    def wyckoff_fraction_vectors(self, row_id: int) -> tuple[Frac3, ...]:
        row_id = int(row_id)
        raw = tuple(
            int(value)
            for value in self.wyckoff["iwyckoff_fract"][(row_id - 1) * 16:row_id * 16]
        )
        vectors: list[Frac3] = []
        for start in range(0, 16, 4):
            denominator = int(raw[start + 3])
            if denominator == 0:
                vectors.append((Fraction(0), Fraction(0), Fraction(0)))
            else:
                vectors.append(
                    tuple(Fraction(raw[start + axis], denominator) for axis in range(3))  # type: ignore[arg-type]
                )
        return tuple(vectors)

    @staticmethod
    def vmlt_fraction(matrix: Iterable[int], vector: Frac3) -> Frac3:
        values = tuple(int(value) for value in matrix)
        return tuple(
            sum(Fraction(values[axis + 3 * component]) * vector[component] for component in range(3))
            for axis in range(3)
        )  # type: ignore[return-value]

    def vrot_fraction(self, sg: int, point_op: int, vector: Frac3) -> Frac3:
        lattice = int(self.space["ispace_lattice"][int(sg) - 1])
        record = tuple(
            int(value)
            for value in self.space["lattice_ml"][(lattice - 1) * 36:lattice * 36]
        )
        back = record[0:9]
        forward = record[9:18]
        denominator = int(self.space["lattice_ml_denom"][(lattice - 1) * 2])
        if denominator == 0:
            raise ValueError(f"zero lattice_ml denominator for lattice {lattice}")
        point = tuple(
            int(value)
            for value in self.space["ipoint_op"][9 * (int(point_op) - 1):9 * int(point_op)]
        )
        rotated = self.vmlt_fraction(
            back,
            self.vmlt_fraction(point, self.vmlt_fraction(forward, vector)),
        )
        return tuple(value / denominator for value in rotated)  # type: ignore[return-value]

    @staticmethod
    def fraction_operation_record(vector: Frac3, point_op: int) -> tuple[int, int, int, int, int]:
        denominator = 1
        for value in vector:
            denominator = math.lcm(denominator, value.denominator)
        numerators = [int(value * denominator) for value in vector]
        divisor = denominator
        for numerator in numerators:
            divisor = math.gcd(divisor, abs(numerator))
        if divisor > 1:
            numerators = [numerator // divisor for numerator in numerators]
            denominator //= divisor
        if denominator < 0:
            numerators = [-numerator for numerator in numerators]
            denominator = -denominator
        return (*numerators, denominator, int(point_op))  # type: ignore[return-value]

    @lru_cache(maxsize=None)
    def wyc_pg_elements_records(
        self,
        sg: int,
        row_id: int,
        site_pg: int,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        base, *parameters = self.wyckoff_fraction_vectors(int(row_id))
        records: list[tuple[int, int, int, int, int]] = []
        for x, y, z, denominator, point_op in self.generate_space_group_records(int(sg)):
            rotated_base = self.vrot_fraction(int(sg), point_op, base)
            shift = tuple(base[axis] - rotated_base[axis] for axis in range(3))
            generated = (
                Fraction(x, denominator),
                Fraction(y, denominator),
                Fraction(z, denominator),
            )
            if any(((shift[axis] - generated[axis]) % 1) != 0 for axis in range(3)):
                continue
            if all(
                self.vrot_fraction(int(sg), point_op, parameter) == parameter
                for parameter in parameters
            ):
                records.append(self.fraction_operation_record(shift, point_op))
        by_op = {record[4]: record for record in records}
        record_ops = set(by_op)
        for setting in self.site_pg_element_settings(int(site_pg)):
            if len(setting) == len(records) and set(setting) == record_ops:
                return tuple(by_op[point_op] for point_op in setting)
        return tuple(records)

    def point_op_product(self, left: int, right: int) -> int:
        return int(self.space["ipoint_op_mlt"][(int(left) - 1) * 72 + int(right) - 1])

    @lru_cache(maxsize=None)
    def wyc_pg_cosets_records(
        self,
        sg: int,
        site_records: tuple[tuple[int, int, int, int, int], ...],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        sg = int(sg)
        generated = self.generate_space_group_records(sg)
        visited = [False] * len(generated)
        index_by_op = {int(record[4]): index for index, record in enumerate(generated)}
        representatives: list[tuple[int, int, int, int, int]] = []
        for index, record in enumerate(generated):
            if visited[index]:
                continue
            representatives.append(record)
            current_op = int(record[4])
            for site in site_records:
                product = self.point_op_product(int(site[4]), current_op)
                if product not in index_by_op:
                    raise KeyError(f"product op {product} not in generated SG{sg} records")
                mark = index_by_op[product]
                if visited[mark]:
                    raise ValueError(f"duplicate coset visit for SG{sg} op={product}")
                visited[mark] = True
        if len(representatives) * len(site_records) != len(generated):
            raise ValueError(
                f"coset decomposition mismatch for SG{sg}: reps={len(representatives)} "
                f"sites={len(site_records)} generated={len(generated)}"
            )
        return tuple(representatives)

    def irrep_matrices_with_ops(self, sg: int, label: str) -> List[IrrepElement]:
        mats = self.irrep_matrices(sg, label)
        ops = self.space_ops(sg)
        if len(mats) != 3 + len(ops):
            raise ValueError(
                f"matrix count mismatch for SG{sg} {label}: "
                f"data_irreps has {len(mats)}, expected 3 + {len(ops)}"
            )
        result: List[IrrepElement] = []
        for ax, D in enumerate(mats[:3]):
            result.append(IrrepElement(None, f"T{ax + 1}", np.eye(3, dtype=int), (Fraction(0),) * 3, D))
        for (op, tau), D in zip(ops, mats[3:]):
            result.append(IrrepElement(op, self.point_labels[op - 1], self.point_ops[op - 1], tau, D))
        return result

    # ------------------------------------------------------------------
