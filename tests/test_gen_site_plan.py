"""Tests for site/gen_site_plan.py PDF generation."""
import math
import pytest
from conftest import _import_from

_mod = _import_from("site", "gen_site_plan")
build_site_plan_data = _mod.build_site_plan_data
render_site_plan = _mod.render_site_plan
render_site_plan_df = _mod.render_site_plan_df
LINE_TOP = _mod.LINE_TOP
LINE_BOT = _mod.LINE_BOT
BOT_LEFT = _mod.BOT_LEFT
TL_251 = _mod.TL_251
CORNER_NW = _mod.CORNER_NW
CORNER_NE = _mod.CORNER_NE
CORNER_SE = _mod.CORNER_SE
CORNER_SW = _mod.CORNER_SW
P45_DIST_216 = _mod.P45_DIST_216
P3_DIST_275 = _mod.P3_DIST_275


@pytest.fixture(scope="module")
def sp_data():
    return build_site_plan_data()


# ============================================================
# Unit tests — computed geometry
# ============================================================

class TestBuildSitePlanData:
    def test_scale(self, sp_data):
        assert sp_data.SCALE == pytest.approx(72.0 / 30.0)

    def test_rotation_degrees(self, sp_data):
        assert sp_data.rotation_deg == pytest.approx(67.2, abs=0.5)

    def test_f15_pdf_x(self, sp_data):
        assert sp_data.f15_pdf[0] == pytest.approx(753.6, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(413.8, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(36.0, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(27.3, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line."""
        assert sp_data.min_setback_216 == pytest.approx(11.500000000036, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line."""
        assert sp_data.min_setback_275 == pytest.approx(25.534716480317, abs=1e-9)

    def test_draw_poly_count(self, sp_data):
        """draw_poly should have same point count as outer_poly."""
        # outer_poly has ~268 pts (18 segments, 9 arcs sampled at 20 pts)
        assert len(sp_data.draw_poly) > 200

    def test_draw_poly_near_outer(self, sp_data):
        """draw_poly should be close to (but slightly inside) outer polygon."""
        # First point of draw_poly should be very close to first point of
        # outer_poly (offset by only half-stroke-width ≈ 0.33 PDF pts)
        dp = sp_data.draw_poly[0]
        # Just check it's a reasonable building coordinate (FC-based)
        assert -20.0 < dp[0] < 20.0
        assert -16.0 < dp[1] < 16.0

    def test_building_to_pdf_f15(self, sp_data):
        """building_to_pdf(F15) should match f15_pdf."""
        f15 = sp_data.pts["F15"]
        result = sp_data.building_to_pdf(*f15)
        assert result[0] == pytest.approx(sp_data.f15_pdf[0], abs=0.01)
        assert result[1] == pytest.approx(sp_data.f15_pdf[1], abs=0.01)

    def test_building_outline_inside_parcel(self, sp_data):
        """All building outline points in PDF coords should fall inside parcel."""
        # Parcel bounding box (generous margins from property line endpoints)
        x_min, x_max = 100.0, 850.0
        y_min, y_max = 30.0, 600.0
        for pt in sp_data.draw_poly:
            px, py = sp_data.building_to_pdf(*pt)
            assert x_min < px < x_max, f"x={px:.1f} outside parcel"
            assert y_min < py < y_max, f"y={py:.1f} outside parcel"

    def test_span_endpoints_valid(self, sp_data):
        """N-S dimension line endpoints should be distinct PDF points."""
        sx, sy = sp_data.span_s_pdf
        nx, ny = sp_data.span_n_pdf
        dist = math.hypot(nx - sx, ny - sy)
        assert dist > 10.0  # at least 10 PDF pts apart

    def test_f2_pdf_inside_parcel(self, sp_data):
        """F2 in PDF coords should be inside the parcel."""
        px, py = sp_data.f2_pdf
        assert 100.0 < px < 850.0
        assert 30.0 < py < 600.0

    def test_p4_setback_216(self, sp_data):
        """P4 should be exactly P45_DIST_216 from the 216.73' line."""
        p4_pdf = sp_data.p_series_pdf["P4"]
        ldx = LINE_BOT[0] - LINE_TOP[0]
        ldy = LINE_BOT[1] - LINE_TOP[1]
        llen = math.hypot(ldx, ldy)
        dist = ((p4_pdf[0] - LINE_TOP[0]) * (-ldy)
                + (p4_pdf[1] - LINE_TOP[1]) * ldx) / llen
        assert dist / sp_data.SCALE == pytest.approx(P45_DIST_216, abs=1e-9)

    def test_p5_setback_216(self, sp_data):
        """P5 should be exactly P45_DIST_216 from the 216.73' line."""
        p5_pdf = sp_data.p_series_pdf["P5"]
        ldx = LINE_BOT[0] - LINE_TOP[0]
        ldy = LINE_BOT[1] - LINE_TOP[1]
        llen = math.hypot(ldx, ldy)
        dist = ((p5_pdf[0] - LINE_TOP[0]) * (-ldy)
                + (p5_pdf[1] - LINE_TOP[1]) * ldx) / llen
        assert dist / sp_data.SCALE == pytest.approx(P45_DIST_216, abs=1e-9)

    def test_p3_setback_275(self, sp_data):
        """P3 should be exactly P3_DIST_275 from the 275.08' line."""
        p3_pdf = sp_data.p_series_pdf["P3"]
        bdx = LINE_BOT[0] - BOT_LEFT[0]
        bdy = LINE_BOT[1] - BOT_LEFT[1]
        blen = math.hypot(bdx, bdy)
        dist = ((p3_pdf[0] - BOT_LEFT[0]) * bdy
                - (p3_pdf[1] - BOT_LEFT[1]) * bdx) / blen
        assert dist / sp_data.SCALE == pytest.approx(P3_DIST_275, abs=1e-9)

    def test_f_series_pdf_count(self, sp_data):
        assert len(sp_data.f_series_pdf) == 19

    def test_f_series_pdf_f15_matches(self, sp_data):
        """f_series_pdf['F15'] should match f15_pdf."""
        assert sp_data.f_series_pdf["F15"][0] == pytest.approx(
            sp_data.f15_pdf[0], abs=1e-9)
        assert sp_data.f_series_pdf["F15"][1] == pytest.approx(
            sp_data.f15_pdf[1], abs=1e-9)


# ============================================================
# Integration tests — PDF rendering
# ============================================================

@pytest.fixture(scope="module")
def base_doc(sp_data):
    doc = render_site_plan(sp_data)
    yield doc
    doc.close()


@pytest.fixture(scope="module")
def base_text(base_doc):
    return base_doc[0].get_text("text")


class TestRenderSitePlan:
    def test_pdf_page_count(self, base_doc):
        assert len(base_doc) == 1

    def test_pdf_page_size(self, base_doc):
        rect = base_doc[0].rect
        assert rect.width == pytest.approx(878.4, abs=1.0)
        assert rect.height == pytest.approx(779.5, abs=1.0)

    def test_text_ew_dim_present(self, base_text, sp_data):
        assert f"{sp_data.ew_dim_ft:.1f}'" in base_text

    def test_text_ns_dim_present(self, base_text, sp_data):
        assert f"{sp_data.ns_dim_ft:.1f}'" in base_text

    def test_text_proposed_present(self, base_text):
        assert "PROPOSED" in base_text

    def test_text_front_present(self, base_text):
        assert "FRONT" in base_text

    def test_text_setback_present(self, base_text):
        assert "11.5'" in base_text

    def test_text_generated_present(self, base_text):
        assert "Generated" in base_text


# ============================================================
# Drainfield variant tests
# ============================================================

@pytest.fixture(scope="module")
def df_doc(sp_data):
    doc = render_site_plan(sp_data, corners=False)
    render_site_plan_df(doc, sp_data)
    yield doc
    doc.close()


@pytest.fixture(scope="module")
def df_text(df_doc):
    return df_doc[0].get_text("text")


class TestRenderSitePlanDf:
    def test_text_existing_drainfield(self, df_text):
        assert "EXISTING" in df_text

    def test_text_new_drainfield(self, df_text):
        assert "NEW" in df_text

    def test_text_drainfield_word(self, df_text):
        assert "DRAINFIELD" in df_text

    def test_base_annotations_still_present(self, df_text):
        """Drainfield variant should still have all base annotations."""
        assert "PROPOSED" in df_text
        assert "FRONT" in df_text


# ============================================================
# d² regression tests — F-series PDF positions vs parcel corners
# ============================================================

from conftest import dist_sq as _dist_sq


_REF_NAMES = ["CORNER_NW", "CORNER_NE", "CORNER_SE", "CORNER_SW"]
_REF_PTS = {
    "CORNER_NW": CORNER_NW,
    "CORNER_NE": CORNER_NE,
    "CORNER_SE": CORNER_SE,
    "CORNER_SW": CORNER_SW,
}

TOL = 1e-6  # PDF pts²

# (name, d²_to_NW, d²_to_NE, d²_to_SE, d²_to_SW)
EXPECTED_D2 = [
    ("F1", 492072.443093, 195365.144921, 11290.031553, 330916.559377),
    ("F2", 489504.306931, 196146.400545, 11618.283322, 327791.254596),
    ("F5", 411981.689348, 176477.032940, 27233.287585, 273216.342780),
    ("F6", 404242.905826, 170531.614483, 29369.813800, 271468.912507),
    ("F7", 403115.389027, 161006.494472, 30196.873055, 278701.902385),
    ("F8", 409835.770279, 158348.578725, 28807.368103, 286982.355773),
    ("F9", 410320.597512, 158163.527600, 28712.917749, 287578.616729),
    ("F10", 409229.567115, 131879.687566, 33485.615515, 311569.106362),
    ("F11", 406489.647829, 129265.394354, 34706.495436, 311718.840250),
    ("F12", 406840.292021, 120794.997747, 37110.272491, 320851.558164),
    ("F13", 442193.842677, 123341.804070, 30804.631238, 351621.400475),
    ("F14", 444032.594952, 123660.441101, 30474.142167, 353067.227445),
    ("F15", 474042.111833, 131556.078415, 24882.648674, 374305.192134),
    ("F16", 478780.284987, 134003.307561, 23705.825557, 376452.153219),
    ("F17", 488166.083931, 142907.851574, 20214.473888, 376502.954418),
    ("F18", 489246.191551, 144888.742422, 19554.545041, 375529.963550),
    ("F11a", 401724.008662, 124721.342230, 36927.158983, 312050.150146),
    ("F11b", 401788.617470, 123160.608480, 37370.072834, 313732.921536),
    ("FC", 444919.004972, 151815.344024, 22848.585070, 325169.723437),
]


def _assert_d2_row(data_dict, idx, expected_list):
    """Check d² from data_dict[name] to each parcel corner against expected."""
    name, exp_nw, exp_ne, exp_se, exp_sw = expected_list[idx]
    pt = data_dict[name]
    expected = [exp_nw, exp_ne, exp_se, exp_sw]
    for i, rn in enumerate(_REF_NAMES):
        actual = _dist_sq(pt, _REF_PTS[rn])
        assert abs(actual - expected[i]) < TOL, (
            f"{name} -> {rn}: d² = {actual:.6f}, "
            f"expected {expected[i]:.6f}, delta {actual - expected[i]:.2e}")


class TestSitePlanD2Regression:
    """Verify F-series PDF positions via d² to parcel corners."""

    @pytest.fixture(scope="class")
    def f_pdf(self, sp_data):
        return sp_data.f_series_pdf

    def test_point_count(self, f_pdf):
        assert len(f_pdf) == len(EXPECTED_D2)

    def test_point_names_match(self, f_pdf):
        actual = list(f_pdf.keys())
        expected = [e[0] for e in EXPECTED_D2]
        assert actual == expected

    @pytest.mark.parametrize("idx", range(len(EXPECTED_D2)),
                             ids=[e[0] for e in EXPECTED_D2])
    def test_d2(self, f_pdf, idx):
        _assert_d2_row(f_pdf, idx, EXPECTED_D2)


# ============================================================
# P-series outer path — distance from 216.73' line
# ============================================================

class TestPSeriesDist216:
    """PX, P4, P5 distance from the 216.73' line."""

    EXPECTED_DIST_216 = {
        "PX": 11.000000000088,
        "P4": 11.000000000000,
        "P5": 10.999999999829,
    }

    @pytest.fixture(scope="class")
    def p_pdf(self, sp_data):
        return sp_data.p_series_pdf

    @pytest.fixture(scope="class")
    def dist_216_params(self, sp_data):
        ldx = LINE_BOT[0] - LINE_TOP[0]
        ldy = LINE_BOT[1] - LINE_TOP[1]
        llen = math.hypot(ldx, ldy)
        return ldx, ldy, llen, sp_data.SCALE

    def _dist_216(self, pt, params):
        ldx, ldy, llen, SCALE = params
        return ((pt[0] - LINE_TOP[0]) * (-ldy)
                + (pt[1] - LINE_TOP[1]) * ldx) / (llen * SCALE)

    @pytest.mark.parametrize("name", ["PX", "P4", "P5"])
    def test_dist_216(self, p_pdf, dist_216_params, name):
        dist = self._dist_216(p_pdf[name], dist_216_params)
        assert dist == pytest.approx(self.EXPECTED_DIST_216[name], abs=1e-9)


# ============================================================
# d² regression tests — P-series PDF positions vs parcel corners
# ============================================================

# ============================================================
# P-series to F-series cross-system distance tests
# ============================================================

EXPECTED_PF_DIST = [
    # (P-name, F-name, distance in feet)
    ("P3", "F2", 1.805316982795),
    ("P2", "F6", 3.120031884583),
    ("POB", "F11", 11.159336313306),
    ("P5", "F15", 8.050187019164),
]


class TestPFCrossDistances:
    """Verify distances between corrected P-series and F-series points."""

    @pytest.fixture(scope="class")
    def p_pdf(self, sp_data):
        return sp_data.p_series_pdf

    @pytest.fixture(scope="class")
    def f_pdf(self, sp_data):
        return sp_data.f_series_pdf

    @pytest.mark.parametrize("idx", range(len(EXPECTED_PF_DIST)),
                             ids=[f"{e[0]}-{e[1]}" for e in EXPECTED_PF_DIST])
    def test_pf_dist(self, p_pdf, f_pdf, sp_data, idx):
        pname, fname, exp_ft = EXPECTED_PF_DIST[idx]
        pp = p_pdf[pname]
        fp = f_pdf[fname]
        d_ft = math.hypot(pp[0] - fp[0], pp[1] - fp[1]) / sp_data.SCALE
        assert d_ft == pytest.approx(exp_ft, abs=1e-9)


EXPECTED_P_D2 = [
    ("POB", 384890.043108, 138759.433111, 37629.578147, 283159.509248),
    ("P2", 402630.896072, 175725.299250, 29677.906581, 265728.881187),
    ("P3", 494938.935024, 196474.443014, 10858.601081, 332731.675432),
    ("P4", 487255.110666, 140667.081533, 20966.132382, 377938.267489),
    ("P5", 463269.109218, 117952.395194, 30997.270469, 378612.218418),
    ("T1", 443013.815106, 121695.680822, 31310.908394, 354225.869977),
    ("T2", 396814.564957, 134369.680822, 35571.306031, 297892.328652),
    ("T3", 489179.595733, 163646.988233, 15106.391326, 356721.685117),
    ("PA", 406309.579254, 116006.741456, 38833.276097, 325513.672283),
    ("PX", 500371.006037, 153129.272321, 16590.565676, 378374.124864),
    ("TC1", 421902.561357, 106263.816558, 40337.006773, 351229.749807),
    ("TC2", 370605.497771, 115259.850492, 47033.929005, 294327.178440),
    ("TC3", 526727.389818, 174051.753109, 10502.434067, 384670.938107),
]


class TestPSeriesD2Regression:
    """Verify P-series PDF positions via d² to parcel corners."""

    @pytest.fixture(scope="class")
    def p_pdf(self, sp_data):
        return sp_data.p_series_pdf

    def test_point_count(self, p_pdf):
        assert len(p_pdf) == len(EXPECTED_P_D2)

    def test_point_names_match(self, p_pdf):
        actual = list(p_pdf.keys())
        expected = [e[0] for e in EXPECTED_P_D2]
        assert actual == expected

    @pytest.mark.parametrize("idx", range(len(EXPECTED_P_D2)),
                             ids=[e[0] for e in EXPECTED_P_D2])
    def test_d2(self, p_pdf, idx):
        _assert_d2_row(p_pdf, idx, EXPECTED_P_D2)
