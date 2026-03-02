"""Unit-aware dimension input parser (CT-7a through CT-7j).

Parses user-entered dimension strings and returns a value in feet.
Supports feet, inches, centimetres, millimetres, and metres, with
flexible whitespace and up to three tokens per input line.
"""

import re

# Conversion factors: unit → feet
_CONVERSIONS = {
    "ft": 1.0,
    "feet": 1.0,
    "in": 1.0 / 12.0,
    "inches": 1.0 / 12.0,
    "cm": 1.0 / 30.48,
    "centimeters": 1.0 / 30.48,
    "mm": 1.0 / 304.8,
    "millimeters": 1.0 / 304.8,
    "m": 1.0 / 0.3048,
    "meters": 1.0 / 0.3048,
}

# Pattern: number (int or float), optional whitespace, optional unit suffix.
# The foot-mark (') and inch-mark (") are handled as single-char units.
# Word units must not be followed by more word chars (to avoid partial match).
_TOKEN_RE = re.compile(
    r"""
    (-?(?:\d+\.?\d*|\.\d+))   # group 1: number (leading minus allowed)
    \s*                         # optional whitespace
    (?:                         # optional unit group
      (['\u2032])               # group 2: foot-mark (' or ′)
      |(["\u2033])              # group 3: inch-mark (" or ″)
      |(feet|ft|inches|in|centimeters|cm|millimeters|mm|meters|m)  # group 4: word unit
      (?![a-z])                 # negative lookahead: no trailing letters
    )?
    """,
    re.VERBOSE | re.IGNORECASE,
)


def parse_dimension(text: str) -> float | None:
    """Parse a dimension string and return value in feet, or None on failure.

    Rules (CT-7a–CT-7j):
    - A bare number is feet.
    - A number + unit suffix is converted per that unit.
    - Two bare numbers → feet and inches.
    - Up to 3 tokens are summed; 4th+ tokens are ignored.
    - Fractions (e.g. "1/3") are evaluated as division (CT-8).
    """
    if text is None:
        return None
    text = text.strip()
    if not text:
        return None

    # CT-8: fraction shortcut — "1/3", "3/4", etc.
    if "/" in text and not any(c in text for c in "' \" \t"):
        parts = text.split("/")
        if len(parts) == 2:
            try:
                return float(parts[0]) / float(parts[1])
            except (ValueError, ZeroDivisionError):
                return None

    tokens = _TOKEN_RE.findall(text)
    if not tokens:
        return None

    # Limit to 3 tokens (CT-7i)
    tokens = tokens[:3]

    result = 0.0
    bare_count = 0  # count of numbers with no unit

    for num_str, foot_mark, inch_mark, word_unit in tokens:
        try:
            value = float(num_str)
        except ValueError:
            return None

        if foot_mark:
            result += value  # feet
        elif inch_mark:
            result += value / 12.0  # inches → feet
        elif word_unit:
            factor = _CONVERSIONS.get(word_unit.lower())
            if factor is None:
                return None
            result += value * factor
        else:
            # Bare number — first is feet, second is inches (CT-7g)
            bare_count += 1
            if bare_count == 1:
                result += value  # feet
            else:
                result += value / 12.0  # inches

    return result
