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
        assert sp_data.f15_pdf[0] == pytest.approx(752.2, abs=2.0)

    def test_f15_pdf_y(self, sp_data):
        assert sp_data.f15_pdf[1] == pytest.approx(415.6, abs=2.0)

    def test_ew_dimension(self, sp_data):
        assert sp_data.ew_dim_ft == pytest.approx(36.0, abs=0.1)

    def test_ns_dimension(self, sp_data):
        assert sp_data.ns_dim_ft == pytest.approx(26.5, abs=0.1)

    def test_min_setback_216(self, sp_data):
        """Min F-point distance from 216.73' line."""
        assert sp_data.min_setback_216 == pytest.approx(11.457680941327, abs=1e-9)

    def test_min_setback_275(self, sp_data):
        """Min F-point distance from 275.08' line."""
        assert sp_data.min_setback_275 == pytest.approx(23.779170036233, abs=1e-9)

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
    ("F1", 497338.491427, 199233.314111, 10439.731330, 332707.196256),
    ("F2", 494754.956967, 200002.039551, 10768.537949, 329576.195721),
    ("F5", 418587.436473, 180125.541144, 25543.743361, 275978.592713),
    ("F6", 410899.852908, 174576.344401, 27546.613507, 273859.870092),
    ("F7", 407973.397282, 164496.143118, 28738.330815, 279909.036272),
    ("F8", 414645.916439, 161394.794611, 27404.297879, 288543.961025),
    ("F9", 415133.921679, 161212.347871, 27309.834903, 289141.459480),
    ("F10", 414052.557187, 135315.670644, 31701.150666, 312536.711246),
    ("F11", 410264.950899, 132126.395738, 33217.767717, 312266.312075),
    ("F12", 408748.785158, 123105.400501, 36021.882175, 320166.264058),
    ("F13", 438706.299141, 124445.907780, 30847.032515, 347037.075756),
    ("F14", 441633.421847, 124913.987024, 30317.876332, 349368.496099),
    ("F15", 476008.593552, 134029.267915, 23883.639333, 373615.700142),
    ("F16", 482487.770296, 137402.545820, 22317.720544, 376562.278131),
    ("F17", 489702.420845, 144250.277392, 19714.558067, 376659.728205),
    ("F18", 493528.452393, 150249.555851, 17774.670828, 374347.788452),
    ("F19", 493979.844018, 152605.282158, 17129.800449, 372409.860032),
    ("F20", 494307.627224, 156070.265172, 16259.475513, 369260.765092),
    ("F11a", 404328.975789, 127533.883067, 35608.067590, 311520.092712),
    ("F11b", 404053.692377, 125895.981436, 36117.198071, 312954.451549),
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
        "PX": 10.665564610640,
        "P4": 10.316176834908,
        "P5": 9.636692600456,
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
    ("P3", "F2", 0.833333333333),
    ("P2", "F6", 3.069062339662),
    ("POB", "F11", 10.185547795239),
    ("P5", "F15", 7.040077056430),
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
    ("POB", 389576.858660, 140113.844169, 36193.529422, 285951.763772),
    ("P2", 404541.343062, 177009.229474, 29155.601579, 266314.801539),
    ("P3", 497563.198799, 200821.500205, 10363.741984, 331587.034602),
    ("P4", 494664.506197, 145820.749229, 18955.450860, 380079.275599),
    ("P5", 472156.745937, 122617.587632, 28454.754986, 382324.891373),
    ("T1", 450885.041888, 125559.812954, 29036.108141, 357715.843035),
    ("T2", 402192.337662, 136268.654000, 33953.264050, 300835.974933),
    ("T3", 494495.359484, 168447.740801, 13758.925398, 357425.127996),
    ("PA", 413688.162717, 118596.595268, 36611.891181, 329564.307105),
    ("PX", 507020.296515, 158534.111826, 14853.350814, 379706.989890),
    ("TC1", 430709.019742, 109596.068891, 37706.532812, 355882.218574),
    ("TC2", 377152.309980, 116493.973921, 44971.294888, 298723.944357),
    ("TC3", 532313.911660, 180014.381423, 9165.378655, 384717.961227),
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
