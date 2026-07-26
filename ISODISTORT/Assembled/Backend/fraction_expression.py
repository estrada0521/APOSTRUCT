"""Restricted exact arithmetic for Source display expressions."""

from __future__ import annotations

import ast
from collections.abc import Mapping
from fractions import Fraction


def split_coordinate_expression3(value: object) -> tuple[str, str, str] | None:
    text = str(value or "").strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    parts = tuple(part.strip() for part in text.split(","))
    return parts if len(parts) == 3 else None


def evaluate_fraction_expression(
    expression: str,
    parameters: Mapping[str, object],
    unknown_name_message: str,
    unsupported_message: str,
    *,
    decimal_names: bool = True,
) -> Fraction:
    tree = ast.parse(expression, mode="eval")

    def visit(node: ast.AST) -> Fraction:
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return Fraction(str(node.value))
        if isinstance(node, ast.Name):
            if node.id not in parameters:
                raise ValueError(unknown_name_message.replace("{name}", repr(node.id)))
            value = parameters[node.id]
            return Fraction(str(value)) if decimal_names else value  # type: ignore[return-value]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return visit(node.operand)
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            right = visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        raise ValueError(unsupported_message)

    return visit(tree)
