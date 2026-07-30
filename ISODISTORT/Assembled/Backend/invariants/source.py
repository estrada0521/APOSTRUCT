"""Source-table adapter for exact crystallographic DISPLAY INVARIANT actions."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd
from typing import Sequence

from sympy import Expr, Rational, expand, nsimplify, sqrt, sqrtdenest

from .algebra import ExactMatrix
from .authority import InvariantSource, invariant_source


@dataclass(frozen=True, eq=False)
class Quadratic3:
    """Canonical element of ``Q(sqrt(2), sqrt(3))``.

    The historical name remains part of the internal API.  Fixed Source
    representations also use ``sqrt(2)`` and ``sqrt(6)``, while parametric
    representations occupy its ``Q(sqrt(3))`` subfield.
    """

    rational: Fraction = Fraction(0)
    radical: Fraction = Fraction(0)
    root2: Fraction = Fraction(0)
    root6: Fraction = Fraction(0)

    @staticmethod
    def _parts(value: object) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        if isinstance(value, Quadratic3):
            return value.rational, value.radical, value.root2, value.root6
        return Fraction(value), Fraction(0), Fraction(0), Fraction(0)  # type: ignore[arg-type]

    @classmethod
    def from_float(cls, value: float, *, tolerance: float = 1e-10) -> Quadratic3:
        scalar = float(value)
        if abs(scalar) <= tolerance:
            return cls()
        root2_symbol = sqrt(2)
        root3_symbol = sqrt(3)
        root6_symbol = sqrt(6)
        expression = expand(nsimplify(
            scalar,
            [root2_symbol, root3_symbol, root6_symbol],
            tolerance=tolerance,
            full=False,
        ))
        root2_part = expression.coeff(root2_symbol)
        root3_part = expression.coeff(root3_symbol)
        root6_part = expression.coeff(root6_symbol)
        rational = (
            expression
            - root2_part * root2_symbol
            - root3_part * root3_symbol
            - root6_part * root6_symbol
        )
        parts = (rational, root3_part, root2_part, root6_part)
        if any(part.is_Rational is not True for part in parts):
            raise ValueError(f"entry is outside Q(sqrt(2), sqrt(3)): {scalar!r}")
        fraction_parts = tuple(
            Fraction(int(part.p), int(part.q))
            for part in parts
        )
        if any(
            abs(part.numerator) > 16 or part.denominator > 96
            for part in fraction_parts
        ):
            raise ValueError(f"entry is outside Q(sqrt(2), sqrt(3)): {scalar!r}")
        out = cls(
            fraction_parts[0],
            fraction_parts[1],
            fraction_parts[2],
            fraction_parts[3],
        )
        if abs(float(out) - scalar) > tolerance:
            raise ValueError(f"entry does not round-trip in Q(sqrt(2), sqrt(3)): {scalar!r}")
        return out

    def __bool__(self) -> bool:
        return bool(self.rational or self.radical or self.root2 or self.root6)

    def __eq__(self, other: object) -> bool:
        try:
            rational, radical, root2, root6 = self._parts(other)
        except (TypeError, ValueError, ZeroDivisionError):
            return False
        return (
            self.rational == rational
            and self.radical == radical
            and self.root2 == root2
            and self.root6 == root6
        )

    def __hash__(self) -> int:
        if not (self.radical or self.root2 or self.root6):
            return hash(self.rational)
        return hash((self.rational, self.radical, self.root2, self.root6, "sqrt2sqrt3"))

    def __float__(self) -> float:
        return (
            float(self.rational)
            + float(self.radical) * (3.0**0.5)
            + float(self.root2) * (2.0**0.5)
            + float(self.root6) * (6.0**0.5)
        )

    @property
    def expression(self) -> Expr:
        return (
            Rational(self.rational.numerator, self.rational.denominator)
            + Rational(self.radical.numerator, self.radical.denominator) * sqrt(3)
            + Rational(self.root2.numerator, self.root2.denominator) * sqrt(2)
            + Rational(self.root6.numerator, self.root6.denominator) * sqrt(6)
        )

    def __neg__(self) -> Quadratic3:
        return Quadratic3(-self.rational, -self.radical, -self.root2, -self.root6)

    def __add__(self, other: object) -> Quadratic3:
        rational, radical, root2, root6 = self._parts(other)
        return Quadratic3(
            self.rational + rational,
            self.radical + radical,
            self.root2 + root2,
            self.root6 + root6,
        )

    def __radd__(self, other: object) -> Quadratic3:
        return self + other

    def __sub__(self, other: object) -> Quadratic3:
        rational, radical, root2, root6 = self._parts(other)
        return Quadratic3(
            self.rational - rational,
            self.radical - radical,
            self.root2 - root2,
            self.root6 - root6,
        )

    def __rsub__(self, other: object) -> Quadratic3:
        return -self + other

    def __mul__(self, other: object) -> Quadratic3:
        rational, radical, root2, root6 = self._parts(other)
        return Quadratic3(
            self.rational * rational
            + 3 * self.radical * radical
            + 2 * self.root2 * root2
            + 6 * self.root6 * root6,
            self.rational * radical
            + self.radical * rational
            + 2 * (self.root2 * root6 + self.root6 * root2),
            self.rational * root2
            + self.root2 * rational
            + 3 * (self.radical * root6 + self.root6 * radical),
            self.rational * root6
            + self.root6 * rational
            + self.radical * root2
            + self.root2 * radical,
        )

    def __rmul__(self, other: object) -> Quadratic3:
        return self * other

    def __truediv__(self, other: object) -> Quadratic3:
        rational, radical, root2, root6 = self._parts(other)
        divisor = Quadratic3(rational, radical, root2, root6)
        if not divisor:
            raise ZeroDivisionError("division by zero in Q(sqrt(2), sqrt(3))")
        if not (radical or root2 or root6):
            return Quadratic3(
                self.rational / rational,
                self.radical / rational,
                self.root2 / rational,
                self.root6 / rational,
            )
        conjugate_root2 = Quadratic3(rational, radical, -root2, -root6)
        norm_root2 = divisor * conjugate_root2
        if norm_root2.root2 or norm_root2.root6:
            raise AssertionError("sqrt(2) conjugation did not eliminate its field component")
        conjugate_root3 = Quadratic3(norm_root2.rational, -norm_root2.radical)
        field_norm = norm_root2 * conjugate_root3
        if field_norm.radical or field_norm.root2 or field_norm.root6 or not field_norm.rational:
            raise AssertionError("biquadratic field norm is not a nonzero rational")
        inverse_numerator = conjugate_root2 * conjugate_root3
        inverse = Quadratic3(
            inverse_numerator.rational / field_norm.rational,
            inverse_numerator.radical / field_norm.rational,
            inverse_numerator.root2 / field_norm.rational,
            inverse_numerator.root6 / field_norm.rational,
        )
        return self * inverse

    def __rtruediv__(self, other: object) -> Quadratic3:
        return Quadratic3(*self._parts(other)) / self

    def as_fraction(self) -> Fraction:
        if self.radical or self.root2 or self.root6:
            raise ValueError("algebraic Reynolds coefficient did not reduce to a rational")
        return self.rational


@dataclass(frozen=True)
class ExactAlgebraicNumber:
    """Sparse exact sum of rational multiples of square-free radicals."""

    parts: tuple[tuple[int, Fraction], ...]
    symbolic: Expr | None = None

    def __post_init__(self) -> None:
        if self.symbolic is not None:
            expression = expand(sqrtdenest(self.symbolic))
            if expression.is_real is not True or expression.is_algebraic is not True:
                raise ValueError("symbolic entry must be a real algebraic number")
            if expression.is_Rational is True:
                object.__setattr__(
                    self,
                    "parts",
                    ((1, Fraction(int(expression.p), int(expression.q))),),
                )
                object.__setattr__(self, "symbolic", None)
                return
            object.__setattr__(self, "parts", ())
            object.__setattr__(self, "symbolic", expression)
            return
        normalized = tuple(sorted(
            (int(radicand), Fraction(coefficient))
            for radicand, coefficient in self.parts
            if coefficient
        ))
        if any(radicand < 1 for radicand, _coefficient in normalized):
            raise ValueError("algebraic radicands must be positive")
        if len({radicand for radicand, _coefficient in normalized}) != len(normalized):
            raise ValueError("algebraic radicands must be unique")
        object.__setattr__(self, "parts", normalized)

    @staticmethod
    def _square_free(value: int) -> tuple[int, int]:
        outside = 1
        inside = 1
        remaining = int(value)
        factor = 2
        while factor * factor <= remaining:
            count = 0
            while remaining % factor == 0:
                remaining //= factor
                count += 1
            outside *= factor ** (count // 2)
            if count % 2:
                inside *= factor
            factor += 1
        if remaining > 1:
            inside *= remaining
        return outside, inside

    @classmethod
    def _from_expression(cls, value: Expr) -> ExactAlgebraicNumber:
        coefficients: dict[int, Fraction] = {}
        expression = expand(sqrtdenest(value))
        terms = expression.as_ordered_terms() if expression.is_Add else (expression,)
        for term in terms:
            coefficient, radical = term.as_coeff_Mul(rational=True)
            if coefficient.is_Rational is not True:
                return cls((), symbolic=expression)
            radicand = 1
            factors = radical.as_ordered_factors() if radical.is_Mul else (radical,)
            for item in factors:
                if item == 1:
                    continue
                if not (
                    item.is_Pow
                    and item.exp == Rational(1, 2)
                    and item.base.is_Integer is True
                    and int(item.base) > 0
                ):
                    return cls((), symbolic=expression)
                radicand *= int(item.base)
            outside, radicand = cls._square_free(radicand)
            fraction = Fraction(int(coefficient.p), int(coefficient.q)) * outside
            coefficients[radicand] = coefficients.get(radicand, Fraction(0)) + fraction
        return cls(tuple(coefficients.items()))

    @staticmethod
    def _coerce(value: object) -> ExactAlgebraicNumber:
        if isinstance(value, ExactAlgebraicNumber):
            return value
        if isinstance(value, Quadratic3):
            return ExactAlgebraicNumber(tuple(
                (radicand, coefficient)
                for radicand, coefficient in (
                    (1, value.rational),
                    (3, value.radical),
                    (2, value.root2),
                    (6, value.root6),
                )
                if coefficient
            ))
        return ExactAlgebraicNumber(((1, Fraction(value)),))  # type: ignore[arg-type]

    @classmethod
    @lru_cache(maxsize=2048)
    def from_float(
        cls,
        value: float,
        *,
        tolerance: float = 1e-10,
    ) -> ExactAlgebraicNumber:
        scalar = float(value)
        if abs(scalar) <= tolerance:
            return cls(())
        expression = nsimplify(scalar, tolerance=tolerance, full=False)
        if expression.is_real is not True or expression.is_algebraic is not True:
            raise ValueError(f"entry is not a real algebraic number: {scalar!r}")
        out = cls._from_expression(expression)
        if any(
            abs(coefficient.numerator) > 4096 or coefficient.denominator > 4096
            for _radicand, coefficient in out.parts
        ):
            raise ValueError(f"algebraic entry has implausible rational content: {scalar!r}")
        if abs(float(out) - scalar) > tolerance:
            raise ValueError(f"entry does not round-trip as an algebraic number: {scalar!r}")
        return out

    @property
    def is_rational(self) -> bool:
        if self.symbolic is not None:
            return self.symbolic.is_Rational is True
        return all(radicand == 1 for radicand, _coefficient in self.parts)

    @property
    def expression(self) -> Expr:
        if self.symbolic is not None:
            return self.symbolic
        return sum(
            (
                Rational(coefficient.numerator, coefficient.denominator)
                * (1 if radicand == 1 else sqrt(radicand))
                for radicand, coefficient in self.parts
            ),
            Rational(0),
        )

    def __bool__(self) -> bool:
        if self.symbolic is not None:
            return self.symbolic != 0
        return bool(self.parts)

    def __eq__(self, other: object) -> bool:
        try:
            converted = self._coerce(other)
            if self.symbolic is not None or converted.symbolic is not None:
                return self.expression == converted.expression
            return self.parts == converted.parts
        except (TypeError, ValueError, ZeroDivisionError):
            return False

    def __hash__(self) -> int:
        if self.symbolic is not None:
            return hash((self.symbolic, "algebraic"))
        if self.is_rational:
            return hash(self.as_fraction())
        return hash((self.parts, "multiquadratic"))

    def __float__(self) -> float:
        if self.symbolic is not None:
            return float(self.symbolic)
        return float(sum(
            float(coefficient) * (radicand ** 0.5)
            for radicand, coefficient in self.parts
        ))

    def __neg__(self) -> ExactAlgebraicNumber:
        if self.symbolic is not None:
            return self._from_expression(-self.symbolic)
        return ExactAlgebraicNumber(tuple(
            (radicand, -coefficient) for radicand, coefficient in self.parts
        ))

    def __add__(self, other: object) -> ExactAlgebraicNumber:
        converted = self._coerce(other)
        if self.symbolic is not None or converted.symbolic is not None:
            return self._from_expression(self.expression + converted.expression)
        coefficients = dict(self.parts)
        for radicand, coefficient in converted.parts:
            coefficients[radicand] = coefficients.get(radicand, Fraction(0)) + coefficient
        return ExactAlgebraicNumber(tuple(coefficients.items()))

    def __radd__(self, other: object) -> ExactAlgebraicNumber:
        return self + other

    def __sub__(self, other: object) -> ExactAlgebraicNumber:
        return self + -self._coerce(other)

    def __rsub__(self, other: object) -> ExactAlgebraicNumber:
        return self._coerce(other) - self

    def __mul__(self, other: object) -> ExactAlgebraicNumber:
        converted = self._coerce(other)
        if self.symbolic is not None or converted.symbolic is not None:
            return self._from_expression(self.expression * converted.expression)
        coefficients: dict[int, Fraction] = {}
        for left_radicand, left_coefficient in self.parts:
            for right_radicand, right_coefficient in converted.parts:
                common = gcd(left_radicand, right_radicand)
                radicand = (left_radicand // common) * (right_radicand // common)
                coefficient = left_coefficient * right_coefficient * common
                coefficients[radicand] = coefficients.get(radicand, Fraction(0)) + coefficient
        return ExactAlgebraicNumber(tuple(coefficients.items()))

    def __rmul__(self, other: object) -> ExactAlgebraicNumber:
        return self * other

    def __truediv__(self, other: object) -> ExactAlgebraicNumber:
        divisor = self._coerce(other)
        if not divisor:
            raise ZeroDivisionError("division by zero in an algebraic number field")
        if self.symbolic is not None or divisor.symbolic is not None:
            return self._from_expression(self.expression / divisor.expression)
        if divisor.is_rational:
            rational = divisor.as_fraction()
            return ExactAlgebraicNumber(tuple(
                (radicand, coefficient / rational)
                for radicand, coefficient in self.parts
            ))
        quotient = nsimplify(self.expression / divisor.expression, full=False)
        return ExactAlgebraicNumber._from_expression(quotient)

    def __rtruediv__(self, other: object) -> ExactAlgebraicNumber:
        return self._coerce(other) / self

    def as_fraction(self) -> Fraction:
        if not self.is_rational:
            raise ValueError("algebraic Reynolds coefficient did not reduce to a rational")
        if self.symbolic is not None:
            return Fraction(int(self.symbolic.p), int(self.symbolic.q))
        return dict(self.parts).get(1, Fraction(0))


def _source_display_entry(value: object) -> Fraction | ExactAlgebraicNumber:
    """Restore an exact scalar from a rounded Source OPD display value."""

    token = str(value).strip()
    scalar = float(token)
    lower = token.lower()
    mantissa, _, exponent_text = lower.partition("e")
    decimals = len(mantissa.rsplit(".", 1)[1]) if "." in mantissa else 0
    exponent = int(exponent_text) if exponent_text else 0
    tolerance = 0.51 * (10.0 ** (exponent - decimals))
    rational = Fraction(scalar).limit_denominator(96)
    if abs(float(rational) - scalar) <= tolerance:
        return rational
    exact = ExactAlgebraicNumber.from_float(scalar, tolerance=tolerance)
    if exact.is_rational:
        return exact.as_fraction()
    return exact


def _source_gid(data: object, space_group: int, irrep: str) -> int:
    matches = [
        gid
        for gid, parent in enumerate(data.little["little_irr_space_group"], start=1)
        if int(parent) == int(space_group)
        and str(data.little["little_irr_full_label"][gid - 1]).strip() == str(irrep)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one Source irrep {irrep!r} in SG{space_group}, found {len(matches)}"
        )
    return matches[0]


def _rational_entry(
    value: complex,
    *,
    maximum_denominator: int = 96,
) -> Fraction | ExactAlgebraicNumber:
    scalar = complex(value)
    if abs(scalar.imag) > 1e-12:
        raise ValueError("real DISPLAY INVARIANT matrix has a complex entry")
    rational = Fraction(float(scalar.real)).limit_denominator(maximum_denominator)
    if abs(float(rational) - scalar.real) <= 1e-12:
        return rational
    exact = ExactAlgebraicNumber.from_float(scalar.real)
    if exact.is_rational:
        return exact.as_fraction()
    return exact


def _multiply(left: ExactMatrix, right: ExactMatrix) -> ExactMatrix:
    size = len(left)
    if (
        len(right) != size
        or any(len(row) != size for row in left)
        or any(len(row) != size for row in right)
    ):
        raise ValueError("representation matrix dimensions differ")
    rows = []
    for row in range(size):
        values: list[object] = [Fraction(0) for _column in range(size)]
        for inner, left_value in enumerate(left[row]):
            if not left_value:
                continue
            for column, right_value in enumerate(right[inner]):
                if right_value:
                    values[column] = values[column] + left_value * right_value  # type: ignore[operator]
        rows.append(tuple(values))
    return tuple(rows)  # type: ignore[return-value]


def _generated_matrix_group(generators: Sequence[ExactMatrix]) -> tuple[ExactMatrix, ...]:
    if not generators:
        raise ValueError("representation generator set is empty")
    size = len(generators[0])
    identity: ExactMatrix = tuple(
        tuple(Fraction(int(row == column)) for column in range(size))
        for row in range(size)
    )
    group = [identity]
    seen = {identity}
    selected: list[ExactMatrix] = []
    for generator in dict.fromkeys(generators):
        if generator in seen:
            continue
        selected.append(generator)
        cursor = 0
        while cursor < len(group):
            current = group[cursor]
            cursor += 1
            for selected_generator in selected:
                product = _multiply(selected_generator, current)
                if product not in seen:
                    seen.add(product)
                    group.append(product)
    return tuple(group)


def _exact_record_matrix(
    data: object,
    gid: int,
    record: tuple[int, int, int, int, int],
) -> ExactMatrix:
    raw = data.get_irreps_matrix_for_record(gid, record, ())
    rows, columns = raw.shape
    if rows != columns:
        raise ValueError("Source irrep matrix is not square")
    return tuple(
        tuple(_rational_entry(raw[row, column]) for column in range(columns))
        for row in range(rows)
    )


def _quadratic3_record_matrix(
    projection: InvariantSource,
    gid: int,
    record: tuple[int, int, int, int, int],
    kparam: Sequence[float | Fraction | int],
) -> ExactMatrix:
    raw = projection.get_irreps_matrix_for_record(gid, record, kparam=kparam)
    rows, columns = raw.shape
    if rows != columns:
        raise ValueError("Source irrep matrix is not square")
    out = []
    for row in range(rows):
        values = []
        for column in range(columns):
            value = complex(raw[row, column])
            if abs(value.imag) > 1e-12:
                raise ValueError("real DISPLAY INVARIANT matrix has a complex entry")
            exact = ExactAlgebraicNumber.from_float(value.real)
            values.append(exact.as_fraction() if exact.is_rational else exact)
        out.append(tuple(values))
    return tuple(out)


def _block_diagonal(blocks: Sequence[ExactMatrix]) -> ExactMatrix:
    size = sum(len(block) for block in blocks)
    rows = [[Fraction(0) for _column in range(size)] for _row in range(size)]
    offset = 0
    for block in blocks:
        width = len(block)
        for row in range(width):
            if len(block[row]) != width:
                raise ValueError("Source irrep matrix is not square")
            for column in range(width):
                rows[offset + row][offset + column] = block[row][column]
        offset += width
    return tuple(tuple(row) for row in rows)


def coupled_fixed_irrep_matrices(
    space_group: int,
    irreps: Sequence[str],
    *,
    source_data: object | None = None,
) -> tuple[ExactMatrix, ...]:
    """Return the synchronized finite image of fixed rational Source irreps.

    Parent point operations alone are insufficient away from Gamma: lattice
    translations can have a nontrivial finite image on the K star.  The binary
    includes that translation kernel before invariant averaging.  Generate the
    same image from every parent operation plus the three unit translations.
    For coupled irreps, each operation is evaluated once in every factor and
    assembled block-diagonally before closure; taking an independent Cartesian
    product of the factor images would describe the wrong group action.
    """

    labels = tuple(str(irrep) for irrep in irreps)
    if not labels:
        raise ValueError("at least one Source irrep is required")
    data = source_data or invariant_source().source_data
    gids = tuple(_source_gid(data, int(space_group), irrep) for irrep in labels)
    lattice = int(data.space["ispace_lattice"][int(space_group) - 1])
    for gid in gids:
        kslot = int(data.little["little_irr_k"][gid - 1])
        slot = (lattice - 1) * 27 + kslot - 1
        if int(data.little["little_k_dim"][slot]) != 0:
            raise ValueError(
                "parametric K requires the later exact-parameter representation frontier"
            )

    space_records = data.generate_space_group_records(int(space_group))
    identity_operation = int(space_records[0][4])
    translation_records = (
        (1, 0, 0, 1, identity_operation),
        (0, 1, 0, 1, identity_operation),
        (0, 0, 1, 1, identity_operation),
    )
    generators: list[ExactMatrix] = []
    for record in (*space_records, *translation_records):
        generators.append(
            _block_diagonal(
                tuple(_exact_record_matrix(data, gid, record) for gid in gids)
            )
        )
    return _generated_matrix_group(generators)


def coupled_parametric_irrep_matrices(
    space_group: int,
    irreps: Sequence[str],
    k_parameters: Sequence[Sequence[float | Fraction | int]],
    *,
    projection_source: InvariantSource | None = None,
) -> tuple[ExactMatrix, ...]:
    """Return the synchronized exact multiquadratic parametric image."""

    labels = tuple(str(irrep) for irrep in irreps)
    parameters = tuple(tuple(values) for values in k_parameters)
    if not labels or len(labels) != len(parameters):
        raise ValueError("irreps and k-parameter rows must have the same nonzero length")
    projection = projection_source or invariant_source()
    data = projection.source_data
    gids = tuple(_source_gid(data, int(space_group), irrep) for irrep in labels)
    source_parameters = []
    for gid, values in zip(gids, parameters, strict=True):
        little = projection.little_record_by_gid(gid)
        if int(little.old_id) > 0:
            if values:
                raise ValueError("fixed Source irrep does not accept explicit k parameters")
            source_parameters.append(None)
        else:
            if not values:
                raise ValueError("parametric Source irrep requires explicit k parameters")
            source_parameters.append(projection.source_kparam_for_gid(gid, values))

    space_records = projection.generate_space_group_records(int(space_group))
    identity_operation = int(space_records[0][4])
    translation_records = (
        (1, 0, 0, 1, identity_operation),
        (0, 1, 0, 1, identity_operation),
        (0, 0, 1, 1, identity_operation),
    )
    generators = []
    for record in (*space_records, *translation_records):
        blocks = []
        for gid, values in zip(gids, source_parameters, strict=True):
            if values is None:
                blocks.append(_exact_record_matrix(data, gid, record))
            else:
                blocks.append(_quadratic3_record_matrix(projection, gid, record, values))
        generators.append(_block_diagonal(tuple(blocks)))
    return _generated_matrix_group(generators)


def fixed_irrep_dimensions(
    space_group: int,
    irreps: Sequence[str],
    *,
    source_data: object | None = None,
) -> tuple[int, ...]:
    data = source_data or invariant_source().source_data
    return tuple(
        int(data.little["little_irr_full_dim"][_source_gid(data, space_group, irrep) - 1])
        for irrep in irreps
    )


def fixed_irrep_matrices(
    space_group: int,
    irrep: str,
    *,
    source_data: object | None = None,
) -> tuple[ExactMatrix, ...]:
    return coupled_fixed_irrep_matrices(
        int(space_group),
        (str(irrep),),
        source_data=source_data,
    )
