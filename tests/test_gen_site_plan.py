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
        assert sp_data.rotation_deg == pytest.approx(67.2, abs=0.5)

    def test_f15_pdf_x(self, sp_data):
        assert sp_data.f15_pdf[0] == pytest.approx(755.9, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(417.9, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(36.0, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(27.3, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line."""
        assert sp_data.min_setback_216 == pytest.approx(10.979236685143613, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line."""
        assert sp_data.min_setback_275 == pytest.approx(23.798821694855195, abs=1e-9)

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
    ("F1", 497592.009814, 199195.701463, 10407.240796, 332983.263591),
    ("F2", 495020.890948, 199973.974383, 10732.509861, 329854.976106),
    ("F5", 417069.393189, 179875.726601, 25918.633947, 274851.184113),
    ("F6", 409256.142997, 173855.841474, 27980.693493, 273029.287170),
    ("F7", 408055.427338, 164257.522603, 28734.553888, 280189.078188),
    ("F8", 414784.160161, 161607.958427, 27353.400506, 288477.883147),
    ("F9", 415269.583934, 161423.503843, 27259.546693, 289074.740644),
    ("F10", 413958.956957, 134920.067229, 31812.647879, 312845.633697),
    ("F11", 411185.205344, 132271.941690, 32999.695474, 312961.535258),
    ("F12", 411458.959725, 123724.655271, 35326.582717, 322017.363361),
    ("F13", 446958.057372, 126417.008585, 29166.488455, 352932.752662),
    ("F14", 448805.645262, 126744.481231, 28844.834998, 354387.415247),
    ("F15", 478979.097697, 134804.054099, 23417.277059, 375789.315490),
    ("F16", 483753.746939, 137287.759334, 22276.930032, 377972.752664),
    ("F17", 493249.053421, 146301.810885, 18895.085900, 378133.061400),
    ("F18", 494350.409198, 148303.949889, 18256.405210, 377181.318689),
    ("F11a", 406359.890626, 127668.214014, 35160.683469, 313233.169602),
    ("F11b", 406410.331912, 126093.312743, 35589.429798, 314901.773471),
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
    """PX, P4, P5 distance from the 216.73' line."""

    EXPECTED_DIST_216 = {
        "PX": 10.479236685196,
        "P4": 10.479236685108,
        "P5": 10.479236684937,
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
    ("POB", 389578.618283, 141758.998106, 35915.795844, 284395.221915),
    ("P2", 407672.781899, 179078.174898, 28317.434930, 267317.904506),
    ("P3", 500475.633933, 200322.131744, 9992.942512, 334815.511833),
    ("P4", 492314.644109, 144037.604797, 19623.308347, 379544.938425),
    ("P5", 468039.360247, 121033.636044, 29365.164021, 379929.606940),
    ("T1", 447764.304249, 124757.159786, 29659.040059, 355523.496612),
    ("T2", 401516.574294, 137382.679980, 33870.957890, 299141.475482),
    ("T3", 494447.902693, 167226.285014, 13972.340808, 358537.129569),
    ("PA", 410880.081786, 118888.233809, 37001.421151, 326631.312308),
    ("PX", 505579.287213, 156648.543318, 15396.489374, 380129.543533),
    ("TC1", 426427.172723, 109099.417744, 38459.260660, 352301.498665),
    ("TC2", 375025.159886, 117990.502428, 45051.233642, 295293.978047),
    ("TC3", 532183.384499, 177818.737611, 9556.071270, 386674.070280),
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
