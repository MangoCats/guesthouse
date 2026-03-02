"""Tests for app.parse_input — CT-7a through CT-7j, CT-8.

Verifies the unit-aware dimension parser converts user input strings
to feet values correctly.
"""

import pytest
from app.parse_input import parse_dimension


# ── CT-7a  Raw Number → feet ────────────────────────────────────────

class TestCT7a_RawNumber:
    def test_integer(self):
        assert parse_dimension("6") == 6.0

    def test_float(self):
        assert parse_dimension("6.5") == 6.5

    def test_leading_decimal(self):
        assert parse_dimension(".5") == 0.5

    def test_zero(self):
        assert parse_dimension("0") == 0.0

    def test_negative(self):
        assert parse_dimension("-3") == -3.0


# ── CT-7b  Feet Suffixes ────────────────────────────────────────────

class TestCT7b_FeetSuffixes:
    def test_single_quote(self):
        assert parse_dimension("6.5'") == 6.5

    def test_ft_suffix(self):
        assert parse_dimension("6.5 ft") == 6.5

    def test_ft_no_space(self):
        assert parse_dimension("6.5ft") == 6.5

    def test_feet_suffix(self):
        assert parse_dimension("6.5feet") == 6.5

    def test_feet_with_space(self):
        assert parse_dimension("6.5 feet") == 6.5


# ── CT-7c  Inch Suffixes ────────────────────────────────────────────

class TestCT7c_InchSuffixes:
    def test_double_quote(self):
        assert parse_dimension('80"') == pytest.approx(80 / 12.0)

    def test_in_suffix(self):
        assert parse_dimension("6.5 in") == pytest.approx(6.5 / 12.0)

    def test_in_no_space(self):
        assert parse_dimension("6.5in") == pytest.approx(6.5 / 12.0)

    def test_inches_suffix(self):
        assert parse_dimension("78inches") == pytest.approx(78 / 12.0)

    def test_inches_with_space(self):
        assert parse_dimension("78 inches") == pytest.approx(78 / 12.0)


# ── CT-7d  Centimetre Suffixes ──────────────────────────────────────

class TestCT7d_CentimetreSuffixes:
    def test_cm_suffix(self):
        assert parse_dimension("30.48cm") == pytest.approx(1.0)

    def test_cm_with_space(self):
        assert parse_dimension("30.48 cm") == pytest.approx(1.0)

    def test_centimeters_suffix(self):
        assert parse_dimension("60.96centimeters") == pytest.approx(2.0)


# ── CT-7e  Millimetre Suffixes ──────────────────────────────────────

class TestCT7e_MillimetreSuffixes:
    def test_mm_suffix(self):
        assert parse_dimension("304.8mm") == pytest.approx(1.0)

    def test_mm_with_space(self):
        assert parse_dimension("304.8 mm") == pytest.approx(1.0)

    def test_millimeters_suffix(self):
        assert parse_dimension("609.6millimeters") == pytest.approx(2.0)


# ── CT-7f  Metre Suffixes ───────────────────────────────────────────

class TestCT7f_MetreSuffixes:
    def test_m_suffix(self):
        assert parse_dimension("0.3048m") == pytest.approx(1.0)

    def test_m_with_space(self):
        assert parse_dimension("0.3048 m") == pytest.approx(1.0)

    def test_meters_suffix(self):
        assert parse_dimension("0.6096meters") == pytest.approx(2.0)


# ── CT-7g  Two Bare Numbers as Feet-Inches ──────────────────────────

class TestCT7g_TwoBareNumbers:
    def test_six_six(self):
        assert parse_dimension("6 6") == pytest.approx(6.5)

    def test_five_zero(self):
        assert parse_dimension("5 0") == pytest.approx(5.0)

    def test_five_six_inches(self):
        # 5 feet 6 inches
        assert parse_dimension("5 6") == pytest.approx(5.5)

    def test_zero_eight(self):
        # 0 feet 8 inches = 8/12
        assert parse_dimension("0 8") == pytest.approx(8 / 12.0)


# ── CT-7h  Multi-Token Summation ────────────────────────────────────

class TestCT7h_MultiTokenSummation:
    def test_feet_and_inches(self):
        assert parse_dimension("6' 6\"") == pytest.approx(6.5)

    def test_mixed_units(self):
        # 6' + 6.7" + 2cm
        expected = 6.0 + 6.7 / 12.0 + 2.0 / 30.48
        assert parse_dimension("6.7in 6' 2cm") == pytest.approx(expected)

    def test_three_tokens(self):
        # 1' + 2" + 3cm
        expected = 1.0 + 2.0 / 12.0 + 3.0 / 30.48
        assert parse_dimension("1' 2\" 3cm") == pytest.approx(expected)

    def test_any_order(self):
        # Same tokens in different order should give same result
        expected = 1.0 + 2.0 / 12.0 + 3.0 / 30.48
        assert parse_dimension("3cm 1' 2\"") == pytest.approx(expected)
        assert parse_dimension("2\" 3cm 1'") == pytest.approx(expected)


# ── CT-7i  Fourth Token Ignored ─────────────────────────────────────

class TestCT7i_FourthTokenIgnored:
    def test_four_tokens(self):
        # 1' + 2" + 3cm, ignore 999ft
        expected = 1.0 + 2.0 / 12.0 + 3.0 / 30.48
        assert parse_dimension("1' 2\" 3cm 999ft") == pytest.approx(expected)

    def test_five_tokens(self):
        # Only first 3 count
        expected = 1.0 + 2.0 + 3.0
        assert parse_dimension("1 2 3 4 5") == pytest.approx(
            1.0 + 2.0 / 12.0 + 3.0 / 12.0  # bare: 1=ft, 2=in, 3=in
        )


# ── CT-7j  Whitespace Flexibility ───────────────────────────────────

class TestCT7j_WhitespaceFlexibility:
    def test_no_space_feet_inches(self):
        # 6'6" — no space between tokens
        assert parse_dimension("6'6\"") == pytest.approx(6.5)

    def test_space_before_unit(self):
        assert parse_dimension("80 \"") == pytest.approx(80 / 12.0)

    def test_leading_trailing_whitespace(self):
        assert parse_dimension("  6.5  ") == 6.5

    def test_multiple_spaces(self):
        assert parse_dimension("6   6") == pytest.approx(6.5)


# ── CT-8  Fractions ─────────────────────────────────────────────────

class TestCT8_Fractions:
    def test_one_third(self):
        assert parse_dimension("1/3") == pytest.approx(1 / 3)

    def test_three_quarters(self):
        assert parse_dimension("3/4") == pytest.approx(0.75)

    def test_seven_eighths(self):
        assert parse_dimension("7/8") == pytest.approx(7 / 8)


# ── Edge Cases ──────────────────────────────────────────────────────

class TestEdgeCases:
    def test_empty_string(self):
        assert parse_dimension("") is None

    def test_none(self):
        assert parse_dimension(None) is None

    def test_whitespace_only(self):
        assert parse_dimension("   ") is None

    def test_non_numeric(self):
        assert parse_dimension("abc") is None

    def test_unit_only(self):
        assert parse_dimension("ft") is None

    def test_case_insensitive_units(self):
        assert parse_dimension("6.5 FT") == 6.5
        assert parse_dimension("80 IN") == pytest.approx(80 / 12.0)
        assert parse_dimension("30.48 CM") == pytest.approx(1.0)

    def test_m_not_matched_in_mm(self):
        # "mm" should be matched as millimetres, not "m" + leftover "m"
        assert parse_dimension("304.8mm") == pytest.approx(1.0)

    def test_fraction_divide_by_zero(self):
        assert parse_dimension("1/0") is None
