"""Source-only atom transport between equivalent selected child bases."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import math

from APOSTRUCT.Backend.modes.structure.magnetic_wyckoff import (
    MagneticWyckoffIdentification,
    identify_magnetic_wyckoff_branch,
    magnetic_group_setting,
    magnetic_orbit_points,
)
from Verification.comparison.basis_lattice import (
    child_coordinate_origin_shift,
)


FractionPoint = tuple[Fraction, Fraction, Fraction]
IntMatrix3 = tuple[
    tuple[int, int, int],
    tuple[int, int, int],
    tuple[int, int, int],
]


@dataclass(frozen=True)
class BasisEquivalentAtom:
    """One local atom expressed and reidentified in the request child frame."""

    label: str
    site: str
    xyz: FractionPoint
    orbit_points: tuple[FractionPoint, ...]


@dataclass(frozen=True)
class BasisEquivalentAtomFrameTransport:
    """A complete atom sequence transported by one exact ``GL(3,Z)`` map."""

    basis_change: IntMatrix3
    request_origin: FractionPoint
    child_space_group: int
    magnetic_group: int
    atoms: tuple[BasisEquivalentAtom, ...]


@dataclass(frozen=True)
class BasisEquivalentOrdinaryAtomFrameTransport:
    """Complete ordinary atom orbits transported by one exact GL(3,Z) map."""

    basis_change: IntMatrix3
    request_origin: FractionPoint
    child_space_group: int
    atoms: tuple[BasisEquivalentAtom, ...]


def _exact_integer(value: object) -> int:
    if isinstance(value, bool):
        raise ValueError("boolean matrix entry")
    fraction = value if isinstance(value, Fraction) else Fraction(str(value))
    if fraction.denominator != 1:
        raise ValueError("non-integral matrix entry")
    return int(fraction)


def exact_unimodular_matrix3(value: object) -> IntMatrix3 | None:
    """Return one exact integral 3x3 matrix only when its determinant is +/-1."""

    if isinstance(value, (str, bytes, bytearray)):
        return None
    try:
        if not isinstance(value, Sequence) or len(value) != 3:
            return None
        matrix = tuple(
            tuple(_exact_integer(item) for item in row)
            for row in value
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray))
        )
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        return None
    determinant = (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    return matrix if abs(determinant) == 1 else None  # type: ignore[return-value]


def _fraction_vector3(value: object) -> FractionPoint | None:
    if isinstance(value, (str, bytes, bytearray)):
        return None
    try:
        if not isinstance(value, Sequence) or len(value) != 3:
            return None
        out = []
        for item in value:
            if isinstance(item, bool):
                return None
            if isinstance(item, float) and not math.isfinite(item):
                return None
            out.append(item if isinstance(item, Fraction) else Fraction(str(item)))
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    return tuple(out)  # type: ignore[return-value]


def _row_multiply(point: FractionPoint, matrix: IntMatrix3) -> FractionPoint:
    return tuple(
        sum(point[row] * matrix[row][column] for row in range(3)) % 1
        for column in range(3)
    )  # type: ignore[return-value]


def _request_point(
    point: FractionPoint,
    matrix: IntMatrix3,
    shift: FractionPoint,
) -> FractionPoint:
    rotated = _row_multiply(point, matrix)
    return tuple((rotated[index] + shift[index]) % 1 for index in range(3))  # type: ignore[return-value]


def _same_periodic_point(
    left: FractionPoint,
    right: FractionPoint,
    *,
    tol: float = 1.0e-12,
) -> bool:
    return all(
        abs(((float(first - second) + 0.5) % 1.0) - 0.5) <= tol
        for first, second in zip(left, right, strict=True)
    )


def _periodic_points_are_distinct(points: Sequence[FractionPoint]) -> bool:
    return all(
        not _same_periodic_point(points[left], points[right])
        for left in range(len(points))
        for right in range(left)
    )


def _site_multiplicity(value: object) -> int | None:
    if not isinstance(value, str) or not value:
        return None
    prefix = value[: len(value) - len(value.lstrip("0123456789"))]
    if not prefix or prefix.startswith("0"):
        return None
    try:
        multiplicity = int(prefix)
    except ValueError:
        return None
    return multiplicity if multiplicity > 0 else None


def basis_equivalent_ordinary_atom_frame_transport(
    local_atoms: object,
    *,
    basis_change: object,
    request_origin: object,
    local_origin: object,
    child_space_group: object,
    tol: float,
    request_basis: object = None,
) -> BasisEquivalentOrdinaryAtomFrameTransport | None:
    """Move complete ordinary atom-orbit witnesses into the request basis.

    The operation is intentionally narrower than ordinary Wyckoff
    reidentification: every local representative must already carry its full
    presentation orbit, and the orbit cardinality must equal the displayed
    site multiplicity.  The caller remains responsible for proving the
    transformed rows against the complete Web atom set.
    """

    matrix = exact_unimodular_matrix3(basis_change)
    request_origin_value = _fraction_vector3(request_origin)
    local_origin_value = _fraction_vector3(local_origin)
    if request_origin_value is None or local_origin_value is None:
        return None
    if request_origin_value == local_origin_value:
        origin_shift: FractionPoint = (Fraction(0), Fraction(0), Fraction(0))
    elif isinstance(request_basis, Sequence) and not isinstance(
        request_basis, (str, bytes, bytearray)
    ):
        origin_shift = child_coordinate_origin_shift(
            request_basis,
            request_origin_value,
            local_origin_value,
        )
        if origin_shift is None:
            return None
    else:
        return None
    try:
        child = _exact_integer(child_space_group)
    except (OverflowError, TypeError, ValueError, ZeroDivisionError):
        return None
    if (
        matrix is None
        or not 1 <= child <= 230
        or not math.isfinite(tol)
        or tol <= 0.0
        or isinstance(local_atoms, (str, bytes, bytearray))
        or not isinstance(local_atoms, Sequence)
        or not local_atoms
    ):
        return None

    labels: set[str] = set()
    all_orbit_points: list[FractionPoint] = []
    transformed: list[BasisEquivalentAtom] = []
    for row in local_atoms:
        if not isinstance(row, Mapping):
            return None
        label = row.get("label")
        site = row.get("site")
        xyz = _fraction_vector3(row.get("xyz"))
        orbit_value = row.get("_presentation_orbit_points")
        multiplicity = _site_multiplicity(site)
        if (
            not isinstance(label, str)
            or not label
            or label in labels
            or not isinstance(site, str)
            or xyz is None
            or multiplicity is None
            or isinstance(orbit_value, (str, bytes, bytearray))
            or not isinstance(orbit_value, Sequence)
            or not orbit_value
        ):
            return None
        orbit: list[FractionPoint] = []
        for point_value in orbit_value:
            point = _fraction_vector3(point_value)
            if point is None:
                return None
            orbit.append(tuple(value % 1 for value in point))  # type: ignore[arg-type]
        folded_xyz = tuple(value % 1 for value in xyz)
        representative_matches = [
            point
            for point in orbit
            if _same_periodic_point(folded_xyz, point, tol=tol)
        ]
        if (
            len(orbit) != multiplicity
            or not _periodic_points_are_distinct(orbit)
            or len(representative_matches) != 1
            or any(
                _same_periodic_point(point, prior)
                for point in orbit
                for prior in all_orbit_points
            )
        ):
            return None
        folded_xyz = representative_matches[0]
        request_xyz = _request_point(folded_xyz, matrix, origin_shift)  # type: ignore[arg-type]
        request_orbit = tuple(
            _request_point(point, matrix, origin_shift) for point in orbit
        )
        if not _periodic_points_are_distinct(request_orbit) or request_xyz not in request_orbit:
            return None
        labels.add(label)
        all_orbit_points.extend(orbit)
        transformed.append(
            BasisEquivalentAtom(
                label=label,
                site=site,
                xyz=request_xyz,
                orbit_points=request_orbit,
            )
        )
    if len(transformed) != len(local_atoms):
        return None
    return BasisEquivalentOrdinaryAtomFrameTransport(
        basis_change=matrix,
        request_origin=request_origin_value,
        child_space_group=child,
        atoms=tuple(transformed),
    )


def basis_equivalent_atom_frame_transport(
    local_atoms: object,
    *,
    basis_change: object,
    request_origin: object,
    local_origin: object,
    child_space_group: int,
    magnetic_group: int,
    request_basis: object = None,
) -> BasisEquivalentAtomFrameTransport | None:
    """Move local child coordinates into the request basis and reidentify sites.

    Provenance: M.  The coordinate relation follows exactly from
    ``B_local = U B_request``.  Magnetic Wyckoff identification and centered
    orbit multiplicity come from the Source-backed cinter tables.  No Web
    coordinate or label participates in this transform.
    """

    matrix = exact_unimodular_matrix3(basis_change)
    request_origin_value = _fraction_vector3(request_origin)
    local_origin_value = _fraction_vector3(local_origin)
    if request_origin_value is None or local_origin_value is None:
        return None
    if request_origin_value == local_origin_value:
        origin_shift: FractionPoint = (Fraction(0), Fraction(0), Fraction(0))
    elif isinstance(request_basis, Sequence) and not isinstance(
        request_basis, (str, bytes, bytearray)
    ):
        origin_shift = child_coordinate_origin_shift(
            request_basis,
            request_origin_value,
            local_origin_value,
        )
        if origin_shift is None:
            return None
    else:
        return None
    if (
        matrix is None
        or isinstance(local_atoms, (str, bytes, bytearray))
        or not isinstance(local_atoms, Sequence)
        or not local_atoms
    ):
        return None
    try:
        child = _exact_integer(child_space_group)
        group = _exact_integer(magnetic_group)
        setting = magnetic_group_setting(group)
    except (IndexError, KeyError, OverflowError, TypeError, ValueError):
        return None
    if child <= 0 or setting.ordinary_space_group != child:
        return None

    identifications: dict[FractionPoint, MagneticWyckoffIdentification] = {}
    orbits: dict[FractionPoint, tuple[FractionPoint, ...]] = {}
    transformed: list[BasisEquivalentAtom] = []
    for row in local_atoms:
        if not isinstance(row, Mapping):
            return None
        label = row.get("label")
        xyz = _fraction_vector3(row.get("xyz"))
        if not isinstance(label, str) or not label or xyz is None:
            return None
        request_xyz = _request_point(xyz, matrix, origin_shift)
        try:
            identification = identifications.get(request_xyz)
            if identification is None:
                identification = identify_magnetic_wyckoff_branch(
                    group,
                    request_xyz,
                    setting="cinter",
                )
                identifications[request_xyz] = identification
            orbit = orbits.get(request_xyz)
            if orbit is None:
                orbit = magnetic_orbit_points(
                    group,
                    request_xyz,
                    record_setting="cinter",
                )
                if not orbit or len(set(orbit)) != len(orbit) or request_xyz not in orbit:
                    return None
                for point in orbit:
                    orbits[point] = orbit
            source_ordinal = identification.row.source_ordinal
            for point in orbit:
                orbit_identification = identifications.get(point)
                if orbit_identification is None:
                    orbit_identification = identify_magnetic_wyckoff_branch(
                        group,
                        point,
                        setting="cinter",
                    )
                    identifications[point] = orbit_identification
                if orbit_identification.row.source_ordinal != source_ordinal:
                    return None
            letter = str(identification.row.label).strip()
        except (IndexError, KeyError, OverflowError, TypeError, ValueError, ZeroDivisionError):
            return None
        if not letter:
            return None
        transformed.append(
            BasisEquivalentAtom(
                label=label,
                site=f"{len(orbit)}{letter}",
                xyz=request_xyz,
                orbit_points=orbit,
            )
        )
    if len(transformed) != len(local_atoms):
        return None
    return BasisEquivalentAtomFrameTransport(
        basis_change=matrix,
        request_origin=request_origin_value,
        child_space_group=child,
        magnetic_group=group,
        atoms=tuple(transformed),
    )
