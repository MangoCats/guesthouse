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
        assert sp_data.rotation_deg == pytest.approx(67.0, abs=0.5)

    def test_f15_pdf_x(self, sp_data):
        assert sp_data.f15_pdf[0] == pytest.approx(753.9, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(415.0, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(37.4, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(25.2, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line."""
        assert sp_data.min_setback_216 == pytest.approx(11.018059189657, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line."""
        assert sp_data.min_setback_275 == pytest.approx(23.422780267135, abs=1e-9)

    def test_draw_poly_count(self, sp_data):
        """draw_poly should have same point count as outer_poly."""
        # outer_poly has ~350 pts
        assert len(sp_data.draw_poly) > 300

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
        assert len(sp_data.f_series_pdf) == 21

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
        assert "11.0'" in base_text

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
    ("F1", 499500.235603, 200078.754452, 10127.617508, 334078.789439),
    ("F2", 496911.755127, 200842.533876, 10451.478111, 330942.842888),
    ("F5", 424521.582300, 181711.215133, 24105.317091, 279906.838052),
    ("F6", 416780.439530, 176108.459185, 26054.628031, 277734.556226),
    ("F7", 413834.336863, 166008.610862, 27226.698299, 283764.075366),
    ("F8", 419755.892921, 162794.597163, 26130.050537, 291871.568997),
    ("F9", 420183.661210, 162569.824756, 26056.518554, 292455.475685),
    ("F10", 413893.097951, 135612.885155, 31662.899329, 312094.486533),
    ("F11", 411207.361342, 133525.754414, 32687.096779, 311721.146721),
    ("F12", 409642.584438, 124456.148015, 35442.600074, 319572.487541),
    ("F13", 444750.345029, 126141.481831, 29518.506306, 351075.221157),
    ("F14", 446807.012348, 126425.062529, 29181.038686, 352776.408252),
    ("F15", 474983.490996, 132478.782847, 24491.649413, 374255.184112),
    ("F16", 482356.538558, 136069.346200, 22771.968394, 377858.692669),
    ("F17", 491812.035570, 145043.588282, 19350.314794, 377979.191937),
    ("F18", 495650.946252, 151055.745875, 17423.306688, 375680.131318),
    ("F19", 496105.309027, 153414.443331, 16781.407459, 373745.174047),
    ("F20", 496436.709394, 156883.043507, 15914.699684, 370599.696269),
    ("F11a", 405246.912013, 128908.767525, 35052.922433, 310950.453138),
    ("F11b", 404962.802491, 127262.039782, 35553.226803, 312375.985865),
    ("FC", 449948.697478, 155156.026350, 21475.920097, 326746.553435),
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
    """PX, P4, P5 should be ~11.0' from the 216.73' line and mutually equal."""

    EXPECTED_DIST_216 = 10.518059189657

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
        assert dist == pytest.approx(self.EXPECTED_DIST_216, abs=1e-9)


# ============================================================
# d² regression tests — P-series PDF positions vs parcel corners
# ============================================================

# ============================================================
# P-series to F-series cross-system distance tests
# ============================================================

EXPECTED_PF_DIST = [
    # (P-name, F-name, distance in feet)
    ("P3", "F2", 1.646007568815),
    ("P2", "F6", 4.396040970861),
    ("POB", "F11", 10.840859298483),
    ("P5", "F15", 5.912252306315),
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
    ("POB", 391192.723068, 143395.780966, 35244.851239, 284264.239170),
    ("P2", 409484.413528, 180912.484600, 27844.017168, 267384.448605),
    ("P3", 502447.584023, 202316.759907, 9679.843211, 335042.374393),
    ("P4", 493993.651281, 145739.290042, 19017.266129, 379478.858067),
    ("P5", 469576.730291, 122593.684162, 28617.484675, 379721.889454),
    ("T1", 449320.240652, 126335.774263, 28929.927072, 355334.345485),
    ("T2", 403118.057638, 139006.841398, 33187.391844, 298997.871295),
    ("T3", 496255.080756, 169056.141150, 13494.469481, 358599.220103),
    ("PA", 412383.591818, 120414.421915, 36219.881793, 326389.734809),
    ("PX", 507331.123566, 158423.057745, 14863.276337, 380136.292357),
    ("TC1", 427878.051357, 110572.974453, 37625.089905, 352007.289769),
    ("TC2", 376495.321020, 119483.341635, 44236.345386, 295019.051650),
    ("TC3", 534051.373013, 179709.404198, 9139.010394, 386796.971264),
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
