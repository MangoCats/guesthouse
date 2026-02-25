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
        assert sp_data.rotation_deg == pytest.approx(73.3, abs=0.5)

    def test_f15_pdf_x(self, sp_data):
        assert sp_data.f15_pdf[0] == pytest.approx(752.2, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(415.6, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(36.0, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(27.3, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line."""
        assert sp_data.min_setback_216 == pytest.approx(11.483624293815556, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line."""
        assert sp_data.min_setback_275 == pytest.approx(24.387148809491602, abs=1e-9)

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
    ("F1", 503255.090953, 198074.891240, 9729.425159, 339416.307067),
    ("F2", 500989.543752, 199072.047524, 9970.379323, 336370.067876),
    ("F5", 423670.027452, 183440.991953, 24296.071650, 277726.838198),
    ("F6", 415123.214396, 177670.863626, 26428.060579, 274972.418064),
    ("F7", 412163.067692, 167670.522784, 27520.718986, 280845.422979),
    ("F8", 418036.199854, 164408.085189, 26375.647327, 288904.492714),
    ("F9", 418460.509294, 164179.853932, 26298.656495, 289484.940552),
    ("F10", 411871.589182, 136470.351406, 31868.151715, 309395.475298),
    ("F11", 408694.130551, 133891.777776, 33105.287491, 309062.599642),
    ("F12", 407119.736998, 124922.197982, 35788.056170, 316766.757224),
    ("F13", 441058.095038, 125259.000587, 30270.103823, 348437.346101),
    ("F14", 442859.117614, 125473.522765, 29974.901246, 349955.029063),
    ("F15", 474256.448040, 132156.076938, 24655.626991, 373869.918857),
    ("F16", 479636.522641, 134770.010831, 23386.165128, 376493.319612),
    ("F17", 489047.561612, 143699.794871, 19920.053486, 376569.360838),
    ("F18", 490112.127129, 145378.101519, 19356.280549, 375900.517497),
    ("F11a", 403156.732415, 129410.730061, 35354.624121, 308542.652964),
    ("F11b", 402866.638660, 127758.018086, 35848.944259, 309962.201458),
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
        "PX": 10.983624519962,
        "P4": 10.983624671598,
        "P5": 10.983624966496,
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
    ("P3", "F2", 1.883462155565),
    ("P2", "F6", 3.653256062376),
    ("POB", "F11", 10.361849087670),
    ("P5", "F15", 5.555344567743),
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
    ("POB", 390236.484601, 143778.212249, 35428.385724, 283071.833122),
    ("P2", 408583.809361, 181350.548998, 28083.184366, 266247.676628),
    ("P3", 501431.868527, 202639.711997, 9803.899443, 333790.492164),
    ("P4", 492859.366873, 145943.674975, 19022.756174, 378108.407546),
    ("P5", 468430.237864, 122785.861928, 28610.767553, 378339.230859),
    ("T1", 448219.781133, 126573.984830, 28969.242289, 353997.719457),
    ("T2", 402130.525924, 139357.979507, 33339.633469, 297774.172232),
    ("T3", 495172.673633, 169312.402600, 13551.835618, 357280.646595),
    ("PA", 411334.776172, 120704.276676, 36310.840620, 325104.752192),
    ("PX", 506203.116468, 158633.719549, 14875.043253, 378772.119174),
    ("TC1", 426777.245509, 110810.839324, 37664.059336, 350670.317302),
    ("TC2", 375507.356394, 119834.047625, 44388.154778, 293794.919539),
    ("TC3", 532925.302972, 179922.002358, 9152.713750, 385434.735247),
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
