"""Exact homogeneous-polynomial operations used by ``coupled_invar_``.

The binary constructs invariants degree by degree.  Products of invariants
already found at lower degrees are admitted first; Reynolds projections of
ordered monomial seeds then complete the invariant subspace.  Keeping these
two phases separate reproduces the observable generator order without making
the printed polynomial basis itself a mathematical uniqueness claim.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from typing import Iterable, Mapping, Protocol, Sequence

from sympy import Expr, Matrix, Rational, sstr


Exponent = tuple[int, ...]


class ExactAlgebraic(Protocol):
    @property
    def expression(self) -> Expr: ...

    def as_fraction(self) -> Fraction: ...


ExactScalar = Fraction | ExactAlgebraic
Polynomial = dict[Exponent, ExactScalar]
ExactMatrix = tuple[tuple[ExactScalar, ...], ...]


def _coefficient(value: object) -> ExactScalar:
    if hasattr(value, "as_fraction"):
        return value  # type: ignore[return-value]
    return Fraction(value)  # type: ignore[arg-type]


def _as_fraction(value: object) -> Fraction:
    if hasattr(value, "as_fraction"):
        return value.as_fraction()  # type: ignore[union-attr,no-any-return]
    return Fraction(value)  # type: ignore[arg-type]


def _sympy_scalar(value: object) -> Expr:
    try:
        fraction = _as_fraction(value)
    except ValueError:
        expression = getattr(value, "expression", None)
        if expression is None:
            raise
        return expression
    return Rational(fraction.numerator, fraction.denominator)


def _divide(left: ExactScalar, right: ExactScalar) -> ExactScalar:
    try:
        return _coefficient(left) / _coefficient(right)  # type: ignore[operator,return-value]
    except TypeError:
        return _coefficient(right).__rtruediv__(_coefficient(left))  # type: ignore[attr-defined,no-any-return]


def _clean(polynomial: Mapping[Exponent, ExactScalar]) -> Polynomial:
    return {
        tuple(int(power) for power in exponent): _coefficient(coefficient)
        for exponent, coefficient in polynomial.items()
        if coefficient
    }


def _constant(variable_count: int) -> Polynomial:
    return {(0,) * int(variable_count): Fraction(1)}


def add(left: Mapping[Exponent, ExactScalar], right: Mapping[Exponent, ExactScalar]) -> Polynomial:
    out = dict(_clean(left))
    for exponent, coefficient in right.items():
        key = tuple(exponent)
        out[key] = out.get(key, Fraction(0)) + _coefficient(coefficient)  # type: ignore[operator]
    return _clean(out)


def multiply(left: Mapping[Exponent, ExactScalar], right: Mapping[Exponent, ExactScalar]) -> Polynomial:
    if not left or not right:
        return {}
    out: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            if len(left_exponent) != len(right_exponent):
                raise ValueError("polynomial variable counts differ")
            exponent = tuple(a + b for a, b in zip(left_exponent, right_exponent, strict=True))
            out[exponent] = out.get(exponent, Fraction(0)) + (
                _coefficient(left_coefficient) * _coefficient(right_coefficient)  # type: ignore[operator]
            )
    return _clean(out)


def power(polynomial: Mapping[Exponent, ExactScalar], exponent: int) -> Polynomial:
    if type(exponent) is not int or exponent < 0:
        raise ValueError("polynomial exponent must be a nonnegative integer")
    if not polynomial:
        return {} if exponent else {(0,): Fraction(1)}
    variable_count = len(next(iter(polynomial)))
    result = _constant(variable_count)
    factor = _clean(polynomial)
    count = exponent
    while count:
        if count & 1:
            result = multiply(result, factor)
        count //= 2
        if count:
            factor = multiply(factor, factor)
    return result


def primitive(polynomial: Mapping[Exponent, ExactScalar]) -> Polynomial:
    """Clear denominators, common content, and the leading sign exactly."""

    clean = _clean(polynomial)
    if not clean:
        return {}
    try:
        rational = {key: _as_fraction(value) for key, value in clean.items()}
    except ValueError:
        leading = max(clean)
        scaled = {
            key: _divide(value, clean[leading])
            for key, value in clean.items()
        }
        try:
            rational = {key: _as_fraction(value) for key, value in scaled.items()}
        except ValueError:
            return _clean(scaled)
    denominator = reduce(lcm, (value.denominator for value in rational.values()), 1)
    integers = {key: int(value * denominator) for key, value in rational.items()}
    content = reduce(gcd, (abs(value) for value in integers.values()))
    integers = {key: value // content for key, value in integers.items()}
    leading = max(integers)
    if integers[leading] < 0:
        integers = {key: -value for key, value in integers.items()}
    return {key: Fraction(value) for key, value in integers.items() if value}


def homogeneous_exponents(variable_count: int, degree: int) -> tuple[Exponent, ...]:
    if type(variable_count) is not int or variable_count < 1:
        raise ValueError("variable_count must be positive")
    if type(degree) is not int or degree < 0:
        raise ValueError("degree must be nonnegative")
    if variable_count == 1:
        return ((degree,),)
    return tuple(
        (first, *tail)
        for first in range(degree, -1, -1)
        for tail in homogeneous_exponents(variable_count - 1, degree - first)
    )


def _linear_form(row: Sequence[ExactScalar], variable: int) -> Polynomial:
    count = len(row)
    out: Polynomial = {}
    for column, coefficient in enumerate(row):
        if coefficient:
            exponent = [0] * count
            exponent[column] = int(variable)
            out[tuple(exponent)] = _coefficient(coefficient)
    return out


def transform_monomial(exponent: Exponent, matrix: ExactMatrix) -> Polynomial:
    count = len(exponent)
    if len(matrix) != count or any(len(row) != count for row in matrix):
        raise ValueError("representation matrix and monomial dimensions differ")
    result = _constant(count)
    for row, row_power in zip(matrix, exponent, strict=True):
        if row_power:
            result = multiply(result, power(_linear_form(row, 1), row_power))
    return result


def reynolds_seed(exponent: Exponent, matrices: Sequence[ExactMatrix]) -> Polynomial:
    if not matrices:
        raise ValueError("a nonempty finite representation group is required")
    projected: Polynomial = {}
    for matrix in matrices:
        projected = add(projected, transform_monomial(exponent, matrix))
    return primitive(projected)


def _vector(polynomial: Mapping[Exponent, ExactScalar], monomials: Sequence[Exponent]) -> list[Expr]:
    return [_sympy_scalar(polynomial.get(exponent, 0)) for exponent in monomials]


def _independent(
    polynomial: Mapping[Exponent, ExactScalar],
    accepted: Sequence[Mapping[Exponent, ExactScalar]],
    monomials: Sequence[Exponent],
) -> bool:
    if not polynomial:
        return False
    columns = [_vector(item, monomials) for item in accepted]
    candidate = _vector(polynomial, monomials)
    if not columns:
        return any(candidate)
    before = Matrix.hstack(*(Matrix(column) for column in columns)).rank()
    after = Matrix.hstack(*(Matrix(column) for column in (*columns, candidate))).rank()
    return after > before


def _decomposable_candidates(
    degree: int,
    bases: Mapping[int, Sequence[Polynomial]],
    *,
    block_dimensions: Sequence[int],
    target_multidegree: Exponent,
) -> Iterable[Polynomial]:
    seen: set[tuple[tuple[Exponent, Fraction], ...]] = set()
    for left_degree in range(1, degree):
        right_degree = degree - left_degree
        if left_degree > right_degree:
            break
        for left_index, left in enumerate(bases.get(left_degree, ())):
            for right_index, right in enumerate(bases.get(right_degree, ())):
                if left_degree == right_degree and right_index < left_index:
                    continue
                candidate = primitive(multiply(left, right))
                if _polynomial_multidegree(candidate, block_dimensions) != target_multidegree:
                    continue
                key = tuple(sorted(candidate.items()))
                if candidate and key not in seen:
                    seen.add(key)
                    yield candidate


def _scalar_decomposable_candidates(
    degree: int,
    bases: Mapping[int, Sequence[Polynomial]],
) -> Iterable[Polynomial]:
    """Yield ``invar_`` products for a list of one-dimensional OPD factors.

    The binary scans the chronological invariant list with the later factor
    outside the earlier factor.  This differs from the per-irrep multidegree
    scan used when any selected OPD has more than one free parameter.
    """

    prior = tuple(
        (prior_degree, polynomial)
        for prior_degree in range(1, degree)
        for polynomial in bases.get(prior_degree, ())
    )
    seen: set[tuple[tuple[Exponent, Fraction], ...]] = set()
    for outer_index, (outer_degree, outer) in enumerate(prior):
        for inner_degree, inner in prior[: outer_index + 1]:
            if outer_degree + inner_degree != degree:
                continue
            candidate = primitive(multiply(outer, inner))
            key = tuple(sorted(candidate.items()))
            if candidate and key not in seen:
                seen.add(key)
                yield candidate


def _exponent_multidegree(
    exponent: Exponent,
    block_dimensions: Sequence[int],
) -> Exponent:
    out = []
    offset = 0
    for width in block_dimensions:
        out.append(sum(exponent[offset : offset + width]))
        offset += width
    if offset != len(exponent):
        raise ValueError("block dimensions do not cover the polynomial variables")
    return tuple(out)


def _polynomial_multidegree(
    polynomial: Mapping[Exponent, Fraction],
    block_dimensions: Sequence[int],
) -> Exponent | None:
    degrees = {
        _exponent_multidegree(exponent, block_dimensions)
        for exponent, coefficient in polynomial.items()
        if coefficient
    }
    if not degrees:
        return None
    if len(degrees) != 1:
        raise ValueError("polynomial is not homogeneous in the irrep blocks")
    return next(iter(degrees))


def invariant_basis(
    matrices: Sequence[ExactMatrix],
    *,
    minimum_degree: int,
    maximum_degree: int,
    block_dimensions: Sequence[int] | None = None,
) -> dict[int, tuple[Polynomial, ...]]:
    """Construct the binary-ordered homogeneous invariant basis."""

    if not matrices:
        raise ValueError("a nonempty finite representation group is required")
    variable_count = len(matrices[0])
    if any(
        len(matrix) != variable_count or any(len(row) != variable_count for row in matrix)
        for matrix in matrices
    ):
        raise ValueError("representation matrices must be square and equally sized")
    if type(minimum_degree) is not int or type(maximum_degree) is not int:
        raise TypeError("degree bounds must be integers")
    if minimum_degree < 1 or maximum_degree < minimum_degree:
        raise ValueError("invalid degree range")
    blocks = tuple(int(width) for width in (block_dimensions or (variable_count,)))
    if not blocks or any(width < 1 for width in blocks) or sum(blocks) != variable_count:
        raise ValueError("block dimensions must be positive and cover every variable")

    bases: dict[int, tuple[Polynomial, ...]] = {}
    for degree in range(1, maximum_degree + 1):
        all_monomials = homogeneous_exponents(variable_count, degree)
        if all(width == 1 for width in blocks):
            accepted: list[Polynomial] = []
            for candidate in _scalar_decomposable_candidates(degree, bases):
                if _independent(candidate, accepted, all_monomials):
                    accepted.append(candidate)
            for exponent in all_monomials:
                candidate = reynolds_seed(exponent, matrices)
                if _independent(candidate, accepted, all_monomials):
                    accepted.append(candidate)
            bases[degree] = tuple(accepted)
            continue
        degree_basis: list[Polynomial] = []
        for multidegree in homogeneous_exponents(len(blocks), degree):
            monomials = tuple(
                exponent
                for exponent in all_monomials
                if _exponent_multidegree(exponent, blocks) == multidegree
            )
            accepted: list[Polynomial] = []
            for candidate in _decomposable_candidates(
                degree,
                bases,
                block_dimensions=blocks,
                target_multidegree=multidegree,
            ):
                if _independent(candidate, accepted, monomials):
                    accepted.append(candidate)
            for exponent in monomials:
                candidate = reynolds_seed(exponent, matrices)
                if _independent(candidate, accepted, monomials):
                    accepted.append(candidate)
            degree_basis.extend(accepted)
        bases[degree] = tuple(degree_basis)
    return {
        degree: bases[degree]
        for degree in range(minimum_degree, maximum_degree + 1)
        if bases[degree]
    }


def polynomial_text(polynomial: Mapping[Exponent, ExactScalar]) -> str:
    """Render one primitive polynomial in the binary's compact variable syntax."""

    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = _sympy_scalar(polynomial[exponent])
        if coefficient == 0:
            continue
        variables = "".join(
            f"n{index}" + (f"^{value}" if value != 1 else "")
            for index, value in enumerate(exponent, start=1)
            if value
        )
        negative = coefficient.is_negative is True
        magnitude = -coefficient if negative else coefficient
        if variables and magnitude == 1:
            body = variables
        else:
            number = sstr(magnitude).replace("**", "^")
            if variables and magnitude.is_Add:
                number = f"({number})"
            body = number + variables
        if not terms:
            terms.append(("-" if negative else "") + body)
        else:
            terms.append((" -" if negative else " +") + body)
    return "".join(terms) or "0"
