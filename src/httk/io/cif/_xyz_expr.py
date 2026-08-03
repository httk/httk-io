"""Private mCIF-only xyz evaluator.

This small duplicate is intentional: the IO layer cannot depend on atomistic's public
operation model, and mCIF needs matrices for moment arithmetic rather than affine objects.
"""

import re
from decimal import Decimal
from fractions import Fraction
from typing import Any

_NUMBER = r"(?:\d+/\d+|\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def _number(text: str, use_fractions: bool) -> Fraction | float:
    if use_fractions:
        return Fraction(text) if "/" in text else Fraction(Decimal(text))
    return float(Fraction(text)) if "/" in text else float(text)


def _parse_linear_expr(expr: str, *, use_fractions: bool = False) -> tuple[tuple[Any, ...], Any]:
    return _parse(expr, ("x", "y", "z"), use_fractions=use_fractions, integer_coefficients=False)


def _parse_linear_expr_algebraic(
    expr: str, *, allowed_vars: tuple[str, ...] = ("x1", "x2", "x3"), use_fractions: bool = False
) -> tuple[tuple[int, ...], Any]:
    row, constant = _parse(expr, allowed_vars, use_fractions=use_fractions, integer_coefficients=True)
    return tuple(int(value) for value in row), constant


def _parse(
    expr: str,
    variables: tuple[str, ...],
    *,
    use_fractions: bool,
    integer_coefficients: bool,
) -> tuple[tuple[Any, ...], Any]:
    text = expr.replace(" ", "")
    if not text:
        raise ValueError("Empty expression")
    if text[0] not in "+-":
        text = "+" + text
    names = "|".join(re.escape(name) for name in sorted(variables, key=len, reverse=True))
    token = re.compile(rf"([+-])(?:(?:(?P<coefficient>{_NUMBER})\*?)?(?P<variable>{names})|(?P<constant>{_NUMBER}))")
    rows: dict[str, Fraction] = {name: Fraction(0) for name in variables}
    constant: Fraction | float = Fraction(0) if use_fractions else 0.0
    position = 0
    for match in token.finditer(text):
        if match.start() != position:
            raise ValueError(f"Unparsed tail in {expr!r} near {text[position:]!r}")
        position = match.end()
        sign = 1 if match.group(1) == "+" else -1
        variable = match.group("variable")
        if variable is None:
            constant += sign * _number(match.group("constant"), use_fractions)
            continue
        coefficient_text = match.group("coefficient")
        if coefficient_text is None:
            coefficient = Fraction(1)
        elif "/" in coefficient_text:
            coefficient = Fraction(coefficient_text)
        else:
            coefficient = Fraction(Decimal(coefficient_text))
        if integer_coefficients and coefficient.denominator != 1:
            raise ValueError(f"Non-integer coefficient on {variable} in {expr!r}")
        rows[variable] += sign * coefficient
    if position != len(text):
        raise ValueError(f"Unparsed tail in {expr!r} near {text[position:]!r}")
    return tuple(rows[name] for name in variables), constant
