"""Exact fixed-K gradient invariants for ``DISPLAY INVARIANT``.

The binary combines each selected irrep with the symmetric Cartesian tensor
representations of derivative orders ``0..gradient``.  Invariants are then
selected at one polynomial degree and one total derivative order.  Keeping
the synchronized ``(irrep, Cartesian)`` group small before expanding the jet
representation avoids closing a much larger block matrix group.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from itertools import combinations_with_replacement
from math import gcd, lcm
from typing import Mapping, Sequence

from .algebra import (
    ExactMatrix,
    ExactScalar,
    Exponent,
    Polynomial,
    _as_fraction,
    _independent,
    add,
    homogeneous_exponents,
    primitive,
    transform_monomial,
)
from .authority import InvariantSource, invariant_source
from APOSTRUCT.Backend.modes.engine.tensor import _cartesian_point_matrix
from .source import (
    ExactAlgebraicNumber,
    _block_diagonal,
    _exact_record_matrix,
    _multiply,
    _quadratic3_record_matrix,
    _source_gid,
)


DerivativeWord = tuple[int, ...]
GradientVariable = tuple[int, DerivativeWord]
_ROOT3 = ExactAlgebraicNumber.from_float(3.0**0.5)
_HEX2: ExactMatrix = (
    (_ROOT3 * Fraction(2, 3), _ROOT3 * Fraction(1, 3), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)
_INVERSE_HEX2: ExactMatrix = (
    (_ROOT3 * Fraction(1, 2), Fraction(-1, 2), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)


def _transpose3(value: ExactMatrix) -> ExactMatrix:
    return tuple(tuple(value[column][row] for column in range(3)) for row in range(3))


def _cartesian_action(
    projection: InvariantSource,
    space_group: int,
    point_operation: int,
) -> ExactMatrix:
    """Return the exact Cartesian carrier initialized by binary ``get_trep_``.

    For operations 1..48, the caller's setting bridge is the same one used by
    ``tensor_basis_``.  ``get_trep_`` initializes entries 49..72 separately by
    conjugating the canonical Source table with its ``hex2`` matrix; routing
    those entries through the tensor helper would apply the setting conversion
    twice.  The final transpose converts Source row-vector storage to the
    column action used by the polynomial substitution below.
    """

    if not 49 <= int(point_operation) <= 72:
        raw = _cartesian_point_matrix(projection, int(space_group), int(point_operation))
        return tuple(
            tuple(
                exact.as_fraction() if exact.is_rational else exact
                for exact in (
                    ExactAlgebraicNumber.from_float(
                        float(raw[row, column]), tolerance=1e-10
                    )
                    for column in range(3)
                )
            )
            for row in range(3)
        )

    flat = projection.source_data.point_operation_matrix(int(point_operation))
    matrix: ExactMatrix = tuple(
        tuple(Fraction(flat[column * 3 + row]) for column in range(3))
        for row in range(3)
    )

    return _multiply(_multiply(_transpose3(_INVERSE_HEX2), matrix), _transpose3(_HEX2))


def _paired_group(
    generators: Sequence[tuple[ExactMatrix, ExactMatrix]],
) -> tuple[tuple[ExactMatrix, ExactMatrix], ...]:
    if not generators:
        raise ValueError("gradient representation generator set is empty")
    irrep_size = len(generators[0][0])
    identity_irrep: ExactMatrix = tuple(
        tuple(Fraction(int(row == column)) for column in range(irrep_size))
        for row in range(irrep_size)
    )
    identity_cartesian: ExactMatrix = tuple(
        tuple(Fraction(int(row == column)) for column in range(3))
        for row in range(3)
    )
    ordered_generators = tuple(dict.fromkeys(generators))
    group = [(identity_irrep, identity_cartesian)]
    seen = set(group)
    cursor = 0
    while cursor < len(group):
        current_irrep, current_cartesian = group[cursor]
        cursor += 1
        for generator_irrep, generator_cartesian in ordered_generators:
            product = (
                _multiply(generator_irrep, current_irrep),
                _multiply(generator_cartesian, current_cartesian),
            )
            if product not in seen:
                seen.add(product)
                group.append(product)
    return tuple(group)


def _derivative_words(order: int) -> tuple[DerivativeWord, ...]:
    if type(order) is not int or order < 0:
        raise ValueError("derivative order must be a nonnegative integer")
    # ``get_trep_`` advances the last axis first and carries to the left.
    # Sorting by the reversed word reproduces xx,yx,yy,zx,zy,zz at order 2.
    return tuple(
        sorted(combinations_with_replacement(range(3), order), key=lambda word: word[::-1])
    )


def _fixed_record_matrix(
    projection: InvariantSource,
    gid: int,
    record: tuple[int, int, int, int, int],
) -> ExactMatrix:
    try:
        return _exact_record_matrix(projection.source_data, gid, record)
    except ValueError as rational_error:
        try:
            return _quadratic3_record_matrix(projection, gid, record, ())
        except ValueError:
            raise rational_error


def _symmetric_derivative_action(matrix: ExactMatrix, order: int) -> ExactMatrix:
    """Return the binary ``get_trep_`` symmetric Cartesian action."""

    words = _derivative_words(order)
    word_index = {word: index for index, word in enumerate(words)}
    rows: list[list[ExactScalar]] = [
        [Fraction(0) for _word in words]
        for _word in words
    ]
    for source_index, source_word in enumerate(words):
        terms: dict[DerivativeWord, ExactScalar] = {(): Fraction(1)}
        for source_axis in source_word:
            expanded: dict[DerivativeWord, ExactScalar] = {}
            for word, coefficient in terms.items():
                for target_axis in range(3):
                    target_word = tuple(sorted((*word, target_axis)))
                    expanded[target_word] = (
                        expanded.get(target_word, Fraction(0))
                        + coefficient * matrix[source_axis][target_axis]  # type: ignore[operator]
                    )
            terms = expanded
        for target_word, coefficient in terms.items():
            rows[source_index][word_index[target_word]] = coefficient
    return tuple(tuple(row) for row in rows)


def _tensor_product(left: ExactMatrix, right: ExactMatrix) -> ExactMatrix:
    return tuple(
        tuple(
            left[row // len(right)][column // len(right)]
            * right[row % len(right)][column % len(right)]  # type: ignore[operator]
            for column in range(len(left) * len(right))
        )
        for row in range(len(left) * len(right))
    )


def _jet_action_from_paired_group(
    paired: Sequence[tuple[ExactMatrix, ExactMatrix]],
    irrep_dimensions: Sequence[int],
    gradient_order: int,
) -> tuple[tuple[ExactMatrix, ...], tuple[GradientVariable, ...]]:
    derivative_cache: dict[tuple[ExactMatrix, int], ExactMatrix] = {}
    matrices = []
    for irrep_matrix, cartesian_matrix in paired:
        blocks = []
        offset = 0
        for width in irrep_dimensions:
            irrep_block = tuple(
                tuple(irrep_matrix[offset + row][offset + column] for column in range(width))
                for row in range(width)
            )
            offset += width
            for order in range(gradient_order + 1):
                key = (cartesian_matrix, order)
                derivative = derivative_cache.get(key)
                if derivative is None:
                    derivative = _symmetric_derivative_action(cartesian_matrix, order)
                    derivative_cache[key] = derivative
                blocks.append(_tensor_product(irrep_block, derivative))
        matrices.append(_block_diagonal(tuple(blocks)))

    variables_list: list[GradientVariable] = []
    component_offset = 0
    for width in irrep_dimensions:
        for order in range(gradient_order + 1):
            for component in range(component_offset, component_offset + width):
                variables_list.extend((component, word) for word in _derivative_words(order))
        component_offset += width
    variables = tuple(variables_list)
    if len(variables) != len(matrices[0]):
        raise AssertionError("gradient variable identities do not cover the jet action")
    return tuple(matrices), variables


def coupled_fixed_gradient_action(
    space_group: int,
    irreps: Sequence[str],
    gradient_order: int,
    *,
    projection_source: InvariantSource | None = None,
) -> tuple[tuple[ExactMatrix, ...], tuple[GradientVariable, ...]]:
    """Return the synchronized fixed-K jet action and printed variable identities."""

    labels = tuple(str(irrep) for irrep in irreps)
    if not labels:
        raise ValueError("at least one Source irrep is required")
    if type(gradient_order) is not int or gradient_order < 1:
        raise ValueError("gradient_order must be a positive integer")
    projection = projection_source or invariant_source()
    data = projection.source_data
    gids = tuple(_source_gid(data, int(space_group), irrep) for irrep in labels)
    lattice = int(data.space["ispace_lattice"][int(space_group) - 1])
    for gid in gids:
        kslot = int(data.little["little_irr_k"][gid - 1])
        slot = (lattice - 1) * 27 + kslot - 1
        if int(data.little["little_k_dim"][slot]) != 0:
            raise ValueError("parametric K gradient invariants are not yet supported")

    records = data.generate_space_group_records(int(space_group))
    identity_operation = int(records[0][4])
    translation_records = (
        (1, 0, 0, 1, identity_operation),
        (0, 1, 0, 1, identity_operation),
        (0, 0, 1, 1, identity_operation),
    )
    generators = []
    for record in (*records, *translation_records):
        generators.append(
            (
                _block_diagonal(
                    tuple(_fixed_record_matrix(projection, gid, record) for gid in gids)
                ),
                _cartesian_action(projection, int(space_group), int(record[4])),
            )
        )
    paired = _paired_group(generators)
    irrep_dimensions = tuple(
        len(_fixed_record_matrix(projection, gid, records[0]))
        for gid in gids
    )
    return _jet_action_from_paired_group(paired, irrep_dimensions, gradient_order)


def coupled_parametric_gradient_action(
    space_group: int,
    irreps: Sequence[str],
    k_parameters: Sequence[Sequence[float | Fraction | int]],
    gradient_order: int,
    *,
    projection_source: InvariantSource | None = None,
) -> tuple[tuple[ExactMatrix, ...], tuple[GradientVariable, ...]]:
    """Return the synchronized exact multiquadratic parametric-K jet action."""

    labels = tuple(str(irrep) for irrep in irreps)
    parameter_rows = tuple(tuple(values) for values in k_parameters)
    if not labels or len(labels) != len(parameter_rows):
        raise ValueError("irreps and k-parameter rows must have the same nonzero length")
    if type(gradient_order) is not int or gradient_order < 1:
        raise ValueError("gradient_order must be a positive integer")
    projection = projection_source or invariant_source()
    data = projection.source_data
    gids = tuple(_source_gid(data, int(space_group), irrep) for irrep in labels)
    source_parameters = []
    for gid, values in zip(gids, parameter_rows, strict=True):
        little = projection.little_record_by_gid(gid)
        if int(little.old_id) > 0:
            if values:
                raise ValueError("fixed Source irrep does not accept explicit k parameters")
            source_parameters.append(None)
        else:
            if not values:
                raise ValueError("parametric Source irrep requires explicit k parameters")
            source_parameters.append(projection.source_kparam_for_gid(gid, values))

    records = projection.generate_space_group_records(int(space_group))
    identity_operation = int(records[0][4])
    translation_records = (
        (1, 0, 0, 1, identity_operation),
        (0, 1, 0, 1, identity_operation),
        (0, 0, 1, 1, identity_operation),
    )
    generators = []
    for record in (*records, *translation_records):
        blocks = []
        for gid, values in zip(gids, source_parameters, strict=True):
            if values is None:
                blocks.append(_fixed_record_matrix(projection, gid, record))
            else:
                blocks.append(_quadratic3_record_matrix(projection, gid, record, values))
        generators.append(
            (
                _block_diagonal(tuple(blocks)),
                _cartesian_action(projection, int(space_group), int(record[4])),
            )
        )
    paired = _paired_group(generators)
    irrep_dimensions = tuple(
        len(
            _fixed_record_matrix(projection, gid, records[0])
            if values is None
            else _quadratic3_record_matrix(projection, gid, records[0], values)
        )
        for gid, values in zip(gids, source_parameters, strict=True)
    )
    return _jet_action_from_paired_group(paired, irrep_dimensions, gradient_order)


def gradient_invariant_basis(
    matrices: Sequence[ExactMatrix],
    variables: Sequence[GradientVariable],
    *,
    degree: int,
    gradient_order: int,
) -> tuple[Polynomial, ...]:
    """Construct the binary-ordered basis at one degree and derivative order."""

    if not matrices or len(variables) != len(matrices[0]):
        raise ValueError("gradient action and variables are empty or inconsistent")
    if type(degree) is not int or degree < 1:
        raise ValueError("degree must be a positive integer")
    if type(gradient_order) is not int or gradient_order < 1:
        raise ValueError("gradient_order must be a positive integer")
    weights = tuple(len(word) for _component, word in variables)
    block_dimensions: list[int] = []
    block_weights: list[int] = []
    for weight in weights:
        if block_weights and block_weights[-1] == weight:
            block_dimensions[-1] += 1
        else:
            block_weights.append(weight)
            block_dimensions.append(1)

    def block_degree(exponent: Exponent) -> Exponent:
        degrees = []
        offset = 0
        for width in block_dimensions:
            degrees.append(sum(exponent[offset : offset + width]))
            offset += width
        if offset != len(exponent):
            raise AssertionError("gradient blocks do not cover the jet variables")
        return tuple(degrees)

    all_monomials = homogeneous_exponents(len(variables), degree)
    block_multidegrees = tuple(
        multidegree
        for multidegree in homogeneous_exponents(len(block_dimensions), degree)
        if sum(
            power * weight
            for power, weight in zip(multidegree, block_weights, strict=True)
        )
        == gradient_order
    )
    monomial_groups = tuple(
        tuple(exponent for exponent in all_monomials if block_degree(exponent) == multidegree)
        for multidegree in block_multidegrees
    )
    monomials = tuple(exponent for group in monomial_groups for exponent in group)
    accepted: list[Polynomial] = []
    for group in monomial_groups:
        for exponent in group:
            projected: Polynomial = {}
            for matrix in matrices:
                projected = add(projected, transform_monomial(exponent, matrix))
            candidate = primitive(projected)
            # ``coupled_invar_grad_`` retains the orientation of the Reynolds
            # seed selected by its monomial scan.  Unlike the ordinary invariant
            # printer, a later lexicographic leading term may therefore be
            # negative (observable for hexagonal rank-two gradients).
            if candidate.get(exponent, Fraction(0)) < 0:  # type: ignore[operator]
                candidate = {key: -value for key, value in candidate.items()}
            if _independent(candidate, accepted, monomials):
                accepted.append(candidate)
    return tuple(accepted)


def gradient_polynomial_text(
    polynomial: Mapping[Exponent, ExactScalar],
    variables: Sequence[GradientVariable],
) -> str:
    """Render the compact derivative suffix convention used by the binary."""

    clean = {tuple(exponent): _as_fraction(value) for exponent, value in polynomial.items() if value}
    if not clean:
        return "0"
    denominator = reduce(lcm, (value.denominator for value in clean.values()), 1)
    integers = {key: int(value * denominator) for key, value in clean.items()}
    content = reduce(gcd, (abs(value) for value in integers.values()))
    integers = {key: value // content for key, value in integers.items()}
    terms = []
    for exponent in sorted(integers, reverse=True):
        coefficient = integers[exponent]
        if not coefficient:
            continue
        factors = []
        for position, power in enumerate(exponent):
            if not power:
                continue
            component, derivative = variables[position]
            # ``get_trep_`` stores a symmetric word in ascending axis order;
            # the printer emits its coordinate suffix in reverse order.
            suffix = "".join("xyz"[axis] for axis in reversed(derivative))
            factors.append(
                f"n{component + 1}{suffix}" + (f"^{power}" if power != 1 else "")
            )
        magnitude = abs(coefficient)
        body = "".join(factors)
        if magnitude != 1 or not body:
            body = f"{magnitude}{body}"
        if not terms:
            terms.append(("-" if coefficient < 0 else "") + body)
        else:
            terms.append((" -" if coefficient < 0 else " +") + body)
    return "".join(terms)
