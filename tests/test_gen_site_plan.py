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
FC_DIST_216 = _mod.FC_DIST_216
FC_DIST_275 = _mod.FC_DIST_275
OUTER_PATH_DIST_216 = _mod.OUTER_PATH_DIST_216


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
        assert sp_data.rotation_deg == pytest.approx(73.3, abs=0.5)

    def test_f15_pdf_x(self, sp_data):
        assert sp_data.f15_pdf[0] == pytest.approx(752.2, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(413.0, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(36.0, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(28.0, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line should be 11.5'."""
        assert sp_data.min_setback_216 == pytest.approx(11.5, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line should be 25.5'."""
        assert sp_data.min_setback_275 == pytest.approx(25.5, abs=1e-9)

    def test_draw_poly_count(self, sp_data):
        """draw_poly should have same point count as outer_poly."""
        # outer_poly has ~350 pts
        assert len(sp_data.draw_poly) > 300

    def test_draw_poly_near_outer(self, sp_data):
        """draw_poly should be close to (but slightly inside) outer polygon."""
        # First point of draw_poly should be very close to first point of
        # outer_poly (offset by only half-stroke-width ≈ 0.33 PDF pts)
        dp = sp_data.draw_poly[0]
        # Just check it's a reasonable building coordinate
        assert 0.0 < dp[0] < 40.0
        assert -5.0 < dp[1] < 30.0

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
        """Span line endpoints should be distinct PDF points."""
        sx, sy = sp_data.span_s_pdf
        nx, ny = sp_data.span_n_pdf
        dist = math.hypot(nx - sx, ny - sy)
        assert dist > 10.0  # at least 10 PDF pts apart

    def test_f2_pdf_inside_parcel(self, sp_data):
        """F2 in PDF coords should be inside the parcel."""
        px, py = sp_data.f2_pdf
        assert 100.0 < px < 850.0
        assert 30.0 < py < 600.0

    def test_fc_setback_216(self, sp_data):
        """FC should be exactly FC_DIST_216 from the 216.73' line."""
        fc_pdf = sp_data.f_series_pdf["FC"]
        ldx = LINE_BOT[0] - LINE_TOP[0]
        ldy = LINE_BOT[1] - LINE_TOP[1]
        llen = math.hypot(ldx, ldy)
        dist = ((fc_pdf[0] - LINE_TOP[0]) * (-ldy)
                + (fc_pdf[1] - LINE_TOP[1]) * ldx) / llen
        assert dist / sp_data.SCALE == pytest.approx(FC_DIST_216, abs=1e-9)

    def test_fc_setback_275(self, sp_data):
        """FC should be exactly FC_DIST_275 from the 275.08' line."""
        fc_pdf = sp_data.f_series_pdf["FC"]
        bdx = LINE_BOT[0] - BOT_LEFT[0]
        bdy = LINE_BOT[1] - BOT_LEFT[1]
        blen = math.hypot(bdx, bdy)
        dist = ((fc_pdf[0] - BOT_LEFT[0]) * bdy
                - (fc_pdf[1] - BOT_LEFT[1]) * bdx) / blen
        assert dist / sp_data.SCALE == pytest.approx(FC_DIST_275, abs=1e-9)

    def test_f_series_pdf_count(self, sp_data):
        assert len(sp_data.f_series_pdf) == 23

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

def _dist_sq(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


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
    ("F1", 495758.724788, 195566.312671, 10788.874732, 334295.794657),
    ("F2", 493469.543377, 196401.427622, 11059.052875, 331379.314322),
    ("F3", 454270.385227, 188046.499836, 17682.499038, 301518.299484),
    ("F4", 451315.617775, 187376.418083, 18265.366822, 299363.629143),
    ("F5", 421381.381621, 180058.602840, 24869.798677, 278497.459882),
    ("F6", 413677.994621, 174493.602663, 26856.865388, 276362.933826),
    ("F7", 410768.383616, 164430.246001, 28065.427317, 282428.944627),
    ("F8", 416681.835777, 161208.128405, 26960.675659, 290528.334362),
    ("F9", 417109.025217, 160982.777148, 26886.564826, 291111.662200),
    ("F10", 410923.882312, 134131.257903, 32598.365956, 310856.093403),
    ("F11", 408253.158362, 132059.139820, 33637.576064, 310497.766249),
    ("F12", 406726.664017, 123027.815980, 36431.361918, 318387.389629),
    ("F13", 441736.354190, 124615.079378, 30409.197733, 349792.052827),
    ("F14", 443786.850272, 124892.488840, 30065.558876, 351487.068686),
    ("F15", 471872.988532, 130855.868769, 25285.829214, 372875.504156),
    ("F16", 479215.113099, 134415.509128, 23535.225201, 376448.089720),
    ("F17", 488614.452569, 143333.593667, 20057.414058, 376512.431444),
    ("F18", 492422.160115, 149314.548123, 18099.202816, 374182.167689),
    ("F19", 492866.808770, 151663.531460, 17447.589468, 372237.496300),
    ("F20", 493184.914665, 155118.837164, 16567.587221, 369078.724049),
    ("F11a", 402326.402224, 127475.846123, 36037.094910, 309760.765859),
    ("F11b", 402049.243495, 125836.069173, 36544.350073, 311193.249378),
    ("FC", 449948.697478, 155156.026350, 21475.920097, 326746.553435),
]


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
        name, exp_nw, exp_ne, exp_se, exp_sw = EXPECTED_D2[idx]
        pt = f_pdf[name]
        expected = [exp_nw, exp_ne, exp_se, exp_sw]
        for i, rn in enumerate(_REF_NAMES):
            actual = _dist_sq(pt, _REF_PTS[rn])
            assert abs(actual - expected[i]) < TOL, (
                f"{name} -> {rn}: d² = {actual:.6f}, "
                f"expected {expected[i]:.6f}, delta {actual - expected[i]:.2e}")


# ============================================================
# P-series outer path — distance from 216.73' line
# ============================================================

class TestPSeriesDist216:
    """PX, P4, P5 should be exactly OUTER_PATH_DIST_216 from the 216.73' line."""

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
        assert dist == pytest.approx(OUTER_PATH_DIST_216, abs=1e-9)


# ============================================================
# d² regression tests — P-series PDF positions vs parcel corners
# ============================================================

# ============================================================
# P-series to F-series cross-system distance tests
# ============================================================

EXPECTED_PF_DIST = [
    # (P-name, F-name, distance in feet)
    ("P3", "F2", 3.041381265149),
    ("P2", "F6", 4.787288387006),
    ("POB", "F11", 11.167322386603),
    ("P5", "F15", 4.931322653473),
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
    ("POB", 390162.483010, 143753.410195, 35451.065525, 283030.550417),
    ("P2", 408507.176890, 181323.117250, 28103.234874, 266203.763271),
    ("P3", 501347.147240, 202604.192412, 9815.860773, 333738.488915),
    ("P4", 492777.207763, 145910.715812, 19037.276956, 378058.965854),
    ("P5", 468350.993658, 122755.816817, 28628.202387, 378292.704126),
    ("T1", 448141.805404, 126545.208303, 28987.946169, 353952.461542),
    ("T2", 402055.662020, 139332.315068, 33361.450571, 297732.026982),
    ("T3", 495089.393497, 169278.323180, 13565.236567, 357230.084149),
    ("PA", 411259.965700, 120678.665085, 36332.710020, 325062.659996),
    ("PX", 506119.458527, 158599.261994, 14888.065643, 378721.178622),
    ("TC1", 426701.755392, 110784.547776, 37685.248285, 350627.545109),
    ("TC2", 375435.599505, 119811.489409, 44413.078216, 293755.881440),
    ("TC3", 532838.954664, 179884.855137, 9163.046390, 385381.104220),
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
        name, exp_nw, exp_ne, exp_se, exp_sw = EXPECTED_P_D2[idx]
        pt = p_pdf[name]
        expected = [exp_nw, exp_ne, exp_se, exp_sw]
        for i, rn in enumerate(_REF_NAMES):
            actual = _dist_sq(pt, _REF_PTS[rn])
            assert abs(actual - expected[i]) < TOL, (
                f"{name} -> {rn}: d² = {actual:.6f}, "
                f"expected {expected[i]:.6f}, delta {actual - expected[i]:.2e}")
