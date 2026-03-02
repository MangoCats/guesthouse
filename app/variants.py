"""Compute furniture/appliance items for each floorplan variant.

Replicates positioning math from floorplan/gen_floorplan.py for the
interactive canvas, producing JSON-serialisable item dicts instead of SVG.
This duplication is intentional per NF-4 (see ARCHITECTURE.md § NF-4).
"""
import math

from shared.geometry import seg_vecs, offset_pt, line_isect
from app.apputil import point_to_list, bbox_from_poly

# ---------------------------------------------------------------------------
# Variant registry
# ---------------------------------------------------------------------------
VARIANTS = {
    "standard": {"label": "Standard", "flags": {}},
    "minik":    {"label": "Small Kitchen", "flags": {"minik": True}},
    "daybed":   {"label": "Daybed", "flags": {"db": True}},
    "bare":     {"label": "Room Dimensions", "flags": {"bare": True}},
    "sf":       {"label": "Square Footage", "flags": {"sf": True}},
}

# ---------------------------------------------------------------------------
# Item dimensions hardcoded in gen_floorplan.py (not in floorplan/constants.py)
# ---------------------------------------------------------------------------
HAMPER_W = 31.5 / 12.0
HAMPER_D = 19.0 / 12.0
MINIK_APPL_W = 32.0 / 12.0
MINIK_APPL_D = 27.0 / 12.0
MICROWAVE_W = 19.5 / 12.0
MICROWAVE_D = 16.625 / 12.0
COFFEE_W = 7.2 / 12.0
COFFEE_D = 9.2 / 12.0
COOKTOP_W = 13.4 / 12.0
COOKTOP_D = 20.5 / 12.0
TOASTER_W = 13.7 / 12.0
TOASTER_D = 12.5 / 12.0
DINING_CHAIR_W = 18.0 / 12.0
DINING_CHAIR_D = 21.0 / 12.0
DINING_TBL_BASE = 31.5 / 12.0
DINING_TBL_H = 35.25 / 12.0
DAYBED_W = 86.0 / 12.0
DAYBED_D = 43.0 / 12.0
WORK_CTR_W = 60.0 / 12.0
WORK_CTR_D = 18.0 / 12.0
STD_FRIDGE_W = 32.75 / 12.0
STD_FRIDGE_D = 35.0 / 12.0
SOFA_FULL_W = 80.75 / 12.0
SOFA_FULL_D = 34.625 / 12.0

# Realistic toilet plan-view polygon (from gen_floorplan.py _TOILET_SVG).
# Coordinates in SVG units; dx = across width (centered on 0), dy = depth from wall.
_TOILET_SVG = [
    (-1.905, 0), (-1.905, 2.032), (-0.841, 2.032),
    (-1.078, 2.224), (-1.267, 2.455), (-1.408, 2.719),
    (-1.495, 3.005), (-1.524, 3.302),
    (-1.732, 5.461), (-1.699, 5.799), (-1.600, 6.124),
    (-1.440, 6.423), (-1.225, 6.686), (-0.962, 6.901),
    (-0.663, 7.061), (-0.338, 7.160), (0, 7.193),
    (0.338, 7.160), (0.663, 7.061), (0.962, 6.901),
    (1.225, 6.686), (1.440, 6.423), (1.600, 6.124),
    (1.699, 5.799), (1.732, 5.461),
    (1.524, 3.302), (1.495, 3.005), (1.408, 2.719),
    (1.267, 2.455), (1.078, 2.224), (0.847, 2.035),
    (0.841, 2.032), (1.905, 2.032), (1.905, 0),
]
_SVG_TO_FT = 10.0 / 30.48


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _toilet_poly(center, facing, width):
    """Transform the toilet plan-view polygon to building coordinates.

    center: (e, n) point on wall face at toilet center-line
    facing: unit vector from wall toward bowl
    width: unit vector along toilet width (perpendicular to facing)
    """
    return [(center[0] + dx * _SVG_TO_FT * width[0] + dy * _SVG_TO_FT * facing[0],
             center[1] + dx * _SVG_TO_FT * width[1] + dy * _SVG_TO_FT * facing[1])
            for dx, dy in _TOILET_SVG]


def _item(name, item_type, poly, label=None, shape="rect"):
    """Build a standard item dict."""
    return {
        "name": name,
        "type": item_type,
        "poly": [point_to_list(p) for p in poly],
        "bbox": bbox_from_poly(poly),
        "label": label or name,
        "shape": shape,
    }


def _circle_item(name, item_type, center, radius, label=None):
    """Build a circle item dict with bounding poly."""
    n_pts = 24
    poly = []
    for i in range(n_pts):
        a = 2 * math.pi * i / n_pts
        poly.append((center[0] + radius * math.cos(a),
                      center[1] + radius * math.sin(a)))
    return {
        "name": name,
        "type": item_type,
        "poly": [point_to_list(p) for p in poly],
        "bbox": bbox_from_poly(poly),
        "label": label or name,
        "shape": "circle",
        "center": point_to_list(center),
        "radius": round(radius, 6),
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_variant_items(pts, inner_poly, layout, radii, variant="standard"):
    """Compute all furniture/appliance items for a variant.

    Returns dict of items keyed by name.
    """
    info = VARIANTS.get(variant, VARIANTS["standard"])
    flags = info["flags"]
    minik = flags.get("minik", False)
    db = flags.get("db", False)
    bare = flags.get("bare", False)
    sf = flags.get("sf", False)

    if bare or sf:
        return {}

    import floorplan.constants as C

    items = {}

    # Common direction vectors
    w2w5_al, w2w5_in = seg_vecs(pts["W2"], pts["W5"])
    w9w10_al, w9w10_in = seg_vecs(pts["W9"], pts["W10"])
    w12w13_al, _ = seg_vecs(pts["W12"], pts["W13"])
    w11w12_al, w11w12_in = seg_vecs(pts["W11"], pts["W12"])

    _iw2_e_al, _iw2_e_out = seg_vecs(layout.iw2s.poly[1], layout.iw2s.poly[2])
    _iw1_n_al, _iw1_n_cw = seg_vecs(layout.iw1.poly[3], layout.iw1.poly[2])
    _iw1_n_out = (-_iw1_n_cw[0], -_iw1_n_cw[1])
    _iw1_w = (-_iw1_n_al[0], -_iw1_n_al[1])

    # IW8 direction vectors
    _iw8_al, _iw8_in = seg_vecs(layout.iw8.poly[0], layout.iw8.poly[1])
    _iw8_out = (-_iw8_in[0], -_iw8_in[1])
    _iw8_s_ref = layout.iw8.poly[0]

    # North wall anchor and helper
    _nw_anchor = pts["W9"]

    def _nwp(d_along, d_inward=0):
        return offset_pt(offset_pt(_nw_anchor, d_along, w9w10_al),
                         d_inward, w9w10_in)

    # IW2s distance along north wall
    _iw2s_ne = layout.iw2s.poly[2]
    _iw2_d = ((_iw2s_ne[0] - _nw_anchor[0]) * w9w10_al[0] +
              (_iw2s_ne[1] - _nw_anchor[1]) * w9w10_al[1])

    # IW2 east face lower along
    _iw2_lo_e_al, _ = seg_vecs(layout.iw2.poly[1], layout.iw2.poly[2])
    _iw12_corner = line_isect(layout.iw2.poly[1], _iw2_lo_e_al,
                              layout.iw1.poly[3], _iw1_n_al)

    def _iwp(d_e, d_n=0):
        return offset_pt(offset_pt(_iw12_corner, d_e, _iw1_n_al),
                         d_n, _iw1_n_out)

    # IW4 / IW1 corner
    _iw4_w_al, _ = seg_vecs(layout.iw4.poly[3], layout.iw4.poly[0])
    _iw41_corner = line_isect(layout.iw4.poly[3], _iw4_w_al,
                              layout.iw1.poly[3], _iw1_n_al)

    def _lwp(d_w=0, d_n=0):
        return offset_pt(offset_pt(_iw41_corner, d_w, _iw1_w),
                         d_n, _iw1_n_out)

    _small_wd = minik or db

    # ===================================================================
    # APPLIANCES (from _render_appliances)
    # ===================================================================

    # --- Dryer & Washer ---
    _shift = (4.0 / 12.0 * w2w5_in[0] + (-2.0 / 12.0) * w2w5_al[0],
              4.0 / 12.0 * w2w5_in[1] + (-2.0 / 12.0) * w2w5_al[1])
    _dryer_nw_mk = None

    for label, wall_obj in [("dryer", layout.dryer), ("washer", layout.washer)]:
        poly = list(wall_obj.poly)
        if not _small_wd:
            poly = [(p[0] + _shift[0], p[1] + _shift[1]) for p in poly]
        if _small_wd:
            if label == "dryer":
                _sw = poly[0]
                _se = offset_pt(_sw, MINIK_APPL_W, w2w5_in)
                _nw = offset_pt(_sw, MINIK_APPL_D, w2w5_al)
                _ne = offset_pt(_se, MINIK_APPL_D, w2w5_al)
                poly = [_sw, _se, _ne, _nw]
                _dryer_nw_mk = _nw
            else:
                _sw = offset_pt(_dryer_nw_mk, 1.0 / 12.0, w2w5_al)
                _se = offset_pt(_sw, MINIK_APPL_W, w2w5_in)
                _nw = offset_pt(_sw, MINIK_APPL_D, w2w5_al)
                _ne = offset_pt(_se, MINIK_APPL_D, w2w5_al)
                poly = [_sw, _se, _ne, _nw]
        items[label] = _item(label, "appliance", poly, label.upper())

    washer_poly = items["washer"]["poly"]  # [[e,n], ...]

    # --- Hamper (standard only) ---
    if not _small_wd:
        washer_nw_raw = layout.washer.poly[3]
        washer_nw_shifted = (washer_nw_raw[0] + _shift[0],
                             washer_nw_raw[1] + _shift[1])
        _washer_nw_d = ((washer_nw_shifted[0] - pts["W2"][0]) * w2w5_al[0] +
                        (washer_nw_shifted[1] - pts["W2"][1]) * w2w5_al[1])
        _hm_sw = offset_pt(offset_pt(pts["W2"], _washer_nw_d + 2.0 / 12.0, w2w5_al),
                            2.0 / 12.0, w2w5_in)
        _hm_se = offset_pt(_hm_sw, HAMPER_W, w2w5_in)
        _hm_nw = offset_pt(_hm_sw, HAMPER_D, w2w5_al)
        _hm_ne = offset_pt(_hm_se, HAMPER_D, w2w5_al)
        items["hamper"] = _item("hamper", "appliance",
                                [_hm_sw, _hm_se, _hm_ne, _hm_nw], "HAMPER")

    # --- Counter (standard only) ---
    if not _small_wd and layout.ctr_clip:
        items["counter"] = _item("counter", "appliance",
                                 list(layout.ctr_clip), "COUNTER")

    # --- Water heater ---
    _wh_ref = offset_pt(layout.iw2s.poly[2], C.WH_RADIUS, _iw2_e_out)
    wh_tangent_r = (radii["R_a7"] - C.WALL_OUTER) - C.WH_RADIUS
    _c7 = pts["C7"]
    _wh_d = (_wh_ref[0] - _c7[0], _wh_ref[1] - _c7[1])
    _wh_d_al = _wh_d[0] * _iw2_e_al[0] + _wh_d[1] * _iw2_e_al[1]
    _wh_d2 = _wh_d[0]**2 + _wh_d[1]**2
    _wh_disc = wh_tangent_r**2 - _wh_d2 + _wh_d_al**2
    if _wh_disc >= 0:
        _wh_t = -_wh_d_al + math.sqrt(_wh_disc)
        wh_center = offset_pt(_wh_ref, _wh_t, _iw2_e_al)
        items["water_heater"] = _circle_item("water_heater", "appliance",
                                             wh_center, C.WH_RADIUS, "WH")

    # --- Toilets ---
    _dryer_cx = sum(p[0] for p in layout.dryer.poly) / 4
    _dryer_cy = sum(p[1] for p in layout.dryer.poly) / 4
    _d_dryer_al = ((_dryer_cx - _iw8_s_ref[0]) * _iw8_al[0] +
                   (_dryer_cy - _iw8_s_ref[1]) * _iw8_al[1])

    # South toilet: realistic polygon on south face of IW8
    _toilet_s_center = offset_pt(_iw8_s_ref, _d_dryer_al - 4.0 / 12.0, _iw8_al)
    _toilet_s_poly = _toilet_poly(_toilet_s_center, _iw8_in, _iw8_al)
    items["toilet_s"] = _item("toilet_s", "fixture", _toilet_s_poly, "TOILET")

    # North toilet: realistic polygon on north face of IW8
    _iw8_n_ref = layout.iw8.poly[3]
    _d_toilet_n_al = ((_dryer_cx - _iw8_n_ref[0]) * _iw8_al[0] +
                      (_dryer_cy - _iw8_n_ref[1]) * _iw8_al[1])
    _toilet_n_center = offset_pt(_iw8_n_ref, _d_toilet_n_al - 4.0 / 12.0, _iw8_al)
    _toilet_n_poly = _toilet_poly(_toilet_n_center, _iw8_out, _iw8_al)
    items["toilet_n"] = _item("toilet_n", "fixture", _toilet_n_poly, "TOILET")

    # --- Utility sink (south face of IW8) ---
    _ctr_cx = sum(p[0] for p in layout.ctr.poly) / 4
    _ctr_cy = sum(p[1] for p in layout.ctr.poly) / 4
    _sink_mid = ((_dryer_cx + _ctr_cx) / 2, (_dryer_cy + _ctr_cy) / 2)
    _d_sink_al = ((_sink_mid[0] - _iw8_s_ref[0]) * _iw8_al[0] +
                  (_sink_mid[1] - _iw8_s_ref[1]) * _iw8_al[1])
    _sink_s_pt = offset_pt(_iw8_s_ref, _d_sink_al, _iw8_al)
    _sk_s = offset_pt(_sink_s_pt, C.SINK_RY, _iw8_in)
    # Bounding rect for utility sink ellipse
    _us_sw = offset_pt(offset_pt(_sk_s, -C.SINK_RX, _iw8_al), -C.SINK_RY, _iw8_in)
    _us_se = offset_pt(_us_sw, 2 * C.SINK_RX, _iw8_al)
    _us_nw = offset_pt(_us_sw, 2 * C.SINK_RY, _iw8_in)
    _us_ne = offset_pt(_us_se, 2 * C.SINK_RY, _iw8_in)
    items["util_sink"] = _item("util_sink", "fixture",
                               [_us_sw, _us_se, _us_ne, _us_nw], "SINK")

    # --- Bath sink (north face of IW8) ---
    _iw2_w = layout.iw2.poly[0]
    _d_iw2_al = ((_iw2_w[0] - _iw8_n_ref[0]) * _iw8_al[0] +
                 (_iw2_w[1] - _iw8_n_ref[1]) * _iw8_al[1])
    _bath_sink_east_d = _d_iw2_al - 9.0 / 12.0
    _bath_sink_ctr_d = _bath_sink_east_d - C.BATH_SINK_LENGTH / 2
    _bath_sink_anchor = offset_pt(_iw8_n_ref, _bath_sink_ctr_d, _iw8_al)
    _bs_sw = offset_pt(offset_pt(_bath_sink_anchor, -C.BATH_SINK_LENGTH / 2, _iw8_al),
                        0, _iw8_out)
    _bs_se = offset_pt(_bs_sw, C.BATH_SINK_LENGTH, _iw8_al)
    _bs_nw = offset_pt(_bs_sw, C.BATH_SINK_DEPTH, _iw8_out)
    _bs_ne = offset_pt(_bs_se, C.BATH_SINK_DEPTH, _iw8_out)
    items["bath_sink"] = _item("bath_sink", "fixture",
                               [_bs_sw, _bs_se, _bs_ne, _bs_nw], "BATH SINK")

    # ===================================================================
    # KITCHEN (from _render_kitchen)
    # ===================================================================

    # Kitchen appliance chain along north wall
    _st_d = _iw2_d + C.NORTH_CTR_LENGTH + C.KITCHEN_APPL_GAP
    _ks_d = _st_d + C.STOVE_WIDTH + C.KITCHEN_APPL_GAP + 2.0 / 12.0
    _dw_d = _ks_d + C.KITCHEN_SINK_WIDTH + C.KITCHEN_APPL_GAP

    # Kitchen appliances on north wall: stove, sink, D/W
    if not minik:
        # Stove
        _st_c = [_nwp(_st_d, C.KITCHEN_APPL_GAP + C.STOVE_DEPTH),
                 _nwp(_st_d + C.STOVE_WIDTH, C.KITCHEN_APPL_GAP + C.STOVE_DEPTH),
                 _nwp(_st_d + C.STOVE_WIDTH, C.KITCHEN_APPL_GAP),
                 _nwp(_st_d, C.KITCHEN_APPL_GAP)]
        items["stove"] = _item("stove", "appliance", _st_c, "STOVE")

        # D/W
        _dw_c = [_nwp(_dw_d, C.DW_DEPTH), _nwp(_dw_d + C.DW_WIDTH, C.DW_DEPTH),
                 _nwp(_dw_d + C.DW_WIDTH, 0), _nwp(_dw_d, 0)]
        items["dishwasher"] = _item("dishwasher", "appliance", _dw_c, "D/W")

    # Kitchen sink (all variants including minik)
    _ks_c = [_nwp(_ks_d, C.KITCHEN_SINK_DEPTH),
             _nwp(_ks_d + C.KITCHEN_SINK_WIDTH, C.KITCHEN_SINK_DEPTH),
             _nwp(_ks_d + C.KITCHEN_SINK_WIDTH, 0), _nwp(_ks_d, 0)]
    items["kitchen_sink"] = _item("kitchen_sink", "fixture", _ks_c, "SINK")

    # Fridge
    if minik:
        _fr_mk_d = _ks_d + C.KITCHEN_SINK_WIDTH + 3.0 / 12.0
        _fr_mk_i = 3.0 / 12.0
        fr_sw = _nwp(_fr_mk_d, _fr_mk_i + C.MINIK_FRIDGE_D)
        fr_se = _nwp(_fr_mk_d + C.MINIK_FRIDGE_W, _fr_mk_i + C.MINIK_FRIDGE_D)
        fr_ne = _nwp(_fr_mk_d + C.MINIK_FRIDGE_W, _fr_mk_i)
        fr_nw = _nwp(_fr_mk_d, _fr_mk_i)
        items["fridge"] = _item("fridge", "appliance",
                                [fr_sw, fr_se, fr_ne, fr_nw], "FRIDGE")
    else:
        fr_nw = _iwp(C.KITCHEN_APPL_GAP, C.KITCHEN_APPL_GAP + STD_FRIDGE_D)
        fr_se = _iwp(C.KITCHEN_APPL_GAP + STD_FRIDGE_W, C.KITCHEN_APPL_GAP)
        fr_ne = _iwp(C.KITCHEN_APPL_GAP + STD_FRIDGE_W, C.KITCHEN_APPL_GAP + STD_FRIDGE_D)
        fr_sw = _iwp(C.KITCHEN_APPL_GAP, C.KITCHEN_APPL_GAP)
        items["fridge"] = _item("fridge", "appliance",
                                [fr_sw, fr_se, fr_ne, fr_nw], "FRIDGE")

    # ICE
    if db:
        _ice_d = _dw_d + C.DW_WIDTH + 2.0 / 12.0
        _ice_i = C.KITCHEN_APPL_GAP
    elif minik:
        _fr_mk_d_ice = _ks_d + C.KITCHEN_SINK_WIDTH + 3.0 / 12.0
        _ice_d = _fr_mk_d_ice + C.MINIK_FRIDGE_W + 3.0 / 12.0
        _ice_i = 3.0 / 12.0
    else:
        _ice_d = _dw_d + C.DW_WIDTH + 6.0 / 12.0
        _ice_i = C.KITCHEN_APPL_GAP
    _ice_c = [_nwp(_ice_d, _ice_i + C.ICE_DEPTH),
              _nwp(_ice_d + C.ICE_WIDTH, _ice_i + C.ICE_DEPTH),
              _nwp(_ice_d + C.ICE_WIDTH, _ice_i), _nwp(_ice_d, _ice_i)]
    items["ice_maker"] = _item("ice_maker", "appliance", _ice_c, "ICE")

    # Work counter (standard/db only)
    if not minik:
        _wc_d2 = C.KITCHEN_APPL_GAP + STD_FRIDGE_W + C.KITCHEN_APPL_GAP
        _wc_c = [_iwp(_wc_d2, 0), _iwp(_wc_d2 + WORK_CTR_W, 0),
                 _iwp(_wc_d2 + WORK_CTR_W, WORK_CTR_D), _iwp(_wc_d2, WORK_CTR_D)]
        items["work_counter"] = _item("work_counter", "appliance", _wc_c, "COUNTER")

    # Microwave
    if not minik:
        _wc_d2_mw = C.KITCHEN_APPL_GAP + STD_FRIDGE_W + C.KITCHEN_APPL_GAP
        _mw_d2 = _wc_d2_mw + 2.0 / 12.0
        _mw_d1 = 2.0 / 12.0
        _mw_c = [_iwp(_mw_d2, _mw_d1), _iwp(_mw_d2 + MICROWAVE_W, _mw_d1),
                 _iwp(_mw_d2 + MICROWAVE_W, _mw_d1 + MICROWAVE_D),
                 _iwp(_mw_d2, _mw_d1 + MICROWAVE_D)]
        items["microwave"] = _item("microwave", "appliance", _mw_c, "MICRO")
    else:
        _mw_mk_d = _iw2_d + 2.0 / 12.0
        _mw_mk_i = 3.0 / 12.0
        _mw_mk_c = [_nwp(_mw_mk_d, _mw_mk_i + MICROWAVE_D),
                     _nwp(_mw_mk_d + MICROWAVE_W, _mw_mk_i + MICROWAVE_D),
                     _nwp(_mw_mk_d + MICROWAVE_W, _mw_mk_i),
                     _nwp(_mw_mk_d, _mw_mk_i)]
        items["microwave"] = _item("microwave", "appliance", _mw_mk_c, "MICRO")

    # Kitchen counter (minik) / North counter (standard)
    if minik:
        _kc_c = [_nwp(_iw2_d, C.KITCHEN_CTR_DEPTH),
                 _nwp(_iw2_d + C.KITCHEN_CTR_LENGTH, C.KITCHEN_CTR_DEPTH),
                 _nwp(_iw2_d + C.KITCHEN_CTR_LENGTH, 0), _nwp(_iw2_d, 0)]
        items["kitchen_counter"] = _item("kitchen_counter", "appliance",
                                         _kc_c, "COUNTER")
    else:
        _nc_c = [_nwp(_iw2_d, C.NORTH_CTR_DEPTH),
                 _nwp(_iw2_d + C.NORTH_CTR_LENGTH, C.NORTH_CTR_DEPTH),
                 _nwp(_iw2_d + C.NORTH_CTR_LENGTH, 0), _nwp(_iw2_d, 0)]
        items["north_counter"] = _item("north_counter", "appliance",
                                        _nc_c, "COUNTER")

    # Coffee maker
    if minik:
        _cm_d = _iw2_d + 2.0 / 12.0 + MICROWAVE_W + 3.0 / 12.0
        _cm_i = 3.0 / 12.0
        _cm_c = [_nwp(_cm_d, _cm_i + COFFEE_D),
                 _nwp(_cm_d + COFFEE_W, _cm_i + COFFEE_D),
                 _nwp(_cm_d + COFFEE_W, _cm_i), _nwp(_cm_d, _cm_i)]
        items["coffee_maker"] = _item("coffee_maker", "appliance",
                                       _cm_c, "C")
    else:
        _cm_d = _iw2_d + C.NORTH_CTR_LENGTH - 2.0 / 12.0 - COFFEE_W
        _cm_i = 2.0 / 12.0
        _cm_c = [_nwp(_cm_d, _cm_i + COFFEE_D),
                 _nwp(_cm_d + COFFEE_W, _cm_i + COFFEE_D),
                 _nwp(_cm_d + COFFEE_W, _cm_i), _nwp(_cm_d, _cm_i)]
        items["coffee_maker"] = _item("coffee_maker", "appliance",
                                       _cm_c, "C")

    # Cooktop (minik only)
    if minik:
        _cp_d = _iw2_d + 2.0 / 12.0 + MICROWAVE_W + 3.0 / 12.0 + COFFEE_W + 3.0 / 12.0
        _cp_i_far = C.KITCHEN_CTR_DEPTH - 2.0 / 12.0
        _cp_i_near = _cp_i_far - COOKTOP_D
        _cp_c = [_nwp(_cp_d, _cp_i_far), _nwp(_cp_d + COOKTOP_W, _cp_i_far),
                 _nwp(_cp_d + COOKTOP_W, _cp_i_near), _nwp(_cp_d, _cp_i_near)]
        items["cooktop"] = _item("cooktop", "appliance", _cp_c, "COOKTOP")

    # Toaster (minik only)
    if minik:
        _ts_d = _cp_d + COOKTOP_W + 3.0 / 12.0
        _ts_i = 3.0 / 12.0
        _ts_c = [_nwp(_ts_d, _ts_i + TOASTER_D),
                 _nwp(_ts_d + TOASTER_W, _ts_i + TOASTER_D),
                 _nwp(_ts_d + TOASTER_W, _ts_i), _nwp(_ts_d, _ts_i)]
        items["toaster"] = _item("toaster", "appliance", _ts_c, "TOASTER")

    # --- Dining set ---
    _space_ne_ref = _iwp(C.RO1_OFFSET_FROM_IW2)
    _space_n_ref = _nwp(0, C.KITCHEN_CTR_DEPTH) if minik else _nwp(0, 0)

    if minik:
        _tbl_ref = _nwp(_st_d + C.STOVE_WIDTH + C.KITCHEN_APPL_GAP)
    else:
        _tbl_ref = _nwp(_ks_d + C.KITCHEN_SINK_WIDTH / 2)
    _tbl_d_al = ((_tbl_ref[0] - _iw12_corner[0]) * _iw1_n_al[0] +
                 (_tbl_ref[1] - _iw12_corner[1]) * _iw1_n_al[1])
    _space_n_d_out = ((_space_n_ref[0] - _iw12_corner[0]) * _iw1_n_out[0] +
                      (_space_n_ref[1] - _iw12_corner[1]) * _iw1_n_out[1])
    _tbl_n_offset = 30.0 / 12.0 + (28.0 / 12.0 if not minik else 0)
    _tbl_n_d_out = _space_n_d_out - _tbl_n_offset
    tbl_bc = offset_pt(offset_pt(_iw12_corner, _tbl_d_al, _iw1_n_al),
                       _tbl_n_d_out, _iw1_n_out)
    _to_apex = (-_iw1_n_out[0], -_iw1_n_out[1])

    # Simplified table as bounding rectangle
    tbl_ne = offset_pt(tbl_bc, DINING_TBL_BASE / 2, _iw1_n_al)
    tbl_nw = offset_pt(tbl_bc, -DINING_TBL_BASE / 2, _iw1_n_al)
    tbl_se = offset_pt(tbl_ne, DINING_TBL_H, _to_apex)
    tbl_sw = offset_pt(tbl_nw, DINING_TBL_H, _to_apex)
    items["dining_table"] = _item("dining_table", "furniture",
                                   [tbl_nw, tbl_ne, tbl_se, tbl_sw], "TABLE")

    # Dining chairs (simplified as rects along the two tangent sides)
    # Compute tangent points for chair placement
    apex_r = 12.0 / 12.0
    _tbl_apex = offset_pt(tbl_bc, DINING_TBL_H, _to_apex)
    arc_c = offset_pt(_tbl_apex, apex_r, _iw1_n_out)

    def _sym(p):
        vx = p[0] - tbl_bc[0]
        vy = p[1] - tbl_bc[1]
        v_dot = vx * _to_apex[0] + vy * _to_apex[1]
        return (tbl_bc[0] + 2 * v_dot * _to_apex[0] - vx,
                tbl_bc[1] + 2 * v_dot * _to_apex[1] - vy)

    fillet_r = 6.0 / 12.0
    d_base_ne = (-_iw1_n_al[0], -_iw1_n_al[1])
    dx_r = tbl_ne[0] - arc_c[0]
    dn_r = tbl_ne[1] - arc_c[1]
    dist_r = math.sqrt(dx_r**2 + dn_r**2)
    angle_cp = math.atan2(dn_r, dx_r)
    delta = math.acos(max(-1, min(1, apex_r / dist_r)))
    alpha_r = angle_cp - delta
    t_right = (arc_c[0] + apex_r * math.cos(alpha_r),
               arc_c[1] + apex_r * math.sin(alpha_r))
    t_left = _sym(t_right)

    dtr = (t_right[0] - tbl_ne[0], t_right[1] - tbl_ne[1])
    dtr_len = math.sqrt(dtr[0]**2 + dtr[1]**2)
    d_tang_ne = (dtr[0] / dtr_len, dtr[1] / dtr_len)
    cos_th = d_base_ne[0] * d_tang_ne[0] + d_base_ne[1] * d_tang_ne[1]
    half_angle = math.acos(max(-1, min(1, cos_th))) / 2
    fillet_dist = fillet_r / math.sin(half_angle)
    bis_ne = (d_base_ne[0] + d_tang_ne[0], d_base_ne[1] + d_tang_ne[1])
    bis_ne_len = math.sqrt(bis_ne[0]**2 + bis_ne[1]**2)
    bis_ne = (bis_ne[0] / bis_ne_len, bis_ne[1] / bis_ne_len)
    _fc_ne = (tbl_ne[0] + fillet_dist * bis_ne[0], tbl_ne[1] + fillet_dist * bis_ne[1])
    v_ne = (_fc_ne[0] - tbl_ne[0], _fc_ne[1] - tbl_ne[1])
    t_proj = v_ne[0] * d_tang_ne[0] + v_ne[1] * d_tang_ne[1]
    f_ne_tang = (tbl_ne[0] + t_proj * d_tang_ne[0], tbl_ne[1] + t_proj * d_tang_ne[1])
    f_nw_tang = _sym(f_ne_tang)

    ch_short = DINING_CHAIR_W
    ch_long = DINING_CHAIR_D
    chair_gap = 2.0 / 12.0

    for idx, (side_start, side_end) in enumerate([(f_ne_tang, t_right), (t_left, f_nw_tang)]):
        mid_e = (side_start[0] + side_end[0]) / 2
        mid_n = (side_start[1] + side_end[1]) / 2
        se_d = (side_end[0] - side_start[0], side_end[1] - side_start[1])
        sl = math.sqrt(se_d[0]**2 + se_d[1]**2)
        su = (se_d[0] / sl, se_d[1] / sl)
        sn = (-su[1], su[0])
        to_ctr = (tbl_bc[0] - mid_e, tbl_bc[1] - mid_n)
        if sn[0] * to_ctr[0] + sn[1] * to_ctr[1] > 0:
            sn = (-sn[0], -sn[1])
        cc_e = mid_e + sn[0] * (ch_long / 2 + chair_gap)
        cc_n = mid_n + sn[1] * (ch_long / 2 + chair_gap)
        corners = []
        for ds, dn in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            ce = cc_e + su[0] * ds * ch_short / 2 + sn[0] * dn * ch_long / 2
            cn = cc_n + su[1] * ds * ch_short / 2 + sn[1] * dn * ch_long / 2
            corners.append((ce, cn))
        name = f"dining_chair_{idx + 1}"
        items[name] = _item(name, "furniture", corners, "CHAIR")

    # ===================================================================
    # FURNITURE (from _render_furniture)
    # ===================================================================

    # --- Bed ---
    items["bed"] = _item("bed", "furniture", list(layout.bed.poly), "KING BED")

    # --- Dresser ---
    items["dresser"] = _item("dresser", "furniture",
                             list(layout.dresser.poly), "DRESSER")

    # --- Shelves ---
    items["shelves"] = _item("shelves", "furniture",
                             list(layout.shelves.poly), "SHELVES")

    # --- Variant-specific seating ---
    if minik:
        # SOFA
        _cx_d = -(6.0 / 12.0 + (C.SOFA_WIDTH - 24.0 / 12.0) / 2)
        sofa_nw = _lwp(_cx_d + SOFA_FULL_W / 2, 2.0 / 12.0 + SOFA_FULL_D)
        sofa_ne = _lwp(_cx_d - SOFA_FULL_W / 2, 2.0 / 12.0 + SOFA_FULL_D)
        sofa_se = _lwp(_cx_d - SOFA_FULL_W / 2, 2.0 / 12.0)
        sofa_sw = _lwp(_cx_d + SOFA_FULL_W / 2, 2.0 / 12.0)
        items["sofa"] = _item("sofa", "furniture",
                               [sofa_sw, sofa_se, sofa_ne, sofa_nw], "SOFA")

        # ROCKER (minik): midpoint between ICE SE and SOFA NW
        _ice_se = _nwp(_ice_d + C.ICE_WIDTH, 3.0 / 12.0 + C.ICE_DEPTH)
        _rk_mid = ((sofa_nw[0] + _ice_se[0]) / 2,
                    (sofa_nw[1] + _ice_se[1]) / 2)
        rk_center = offset_pt(_rk_mid, 18.0 / 12.0, w9w10_in)
        _rk_cr = (-w12w13_al[1], w12w13_al[0])
        rk_hw = C.ROCKER_DEPTH / 2
        rk_hh = C.ROCKER_WIDTH / 2
        rk_cx, rk_cy = rk_center
        rk_poly = [
            (rk_cx - rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
            (rk_cx - rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
        ]
        items["rocker"] = _item("rocker", "furniture", rk_poly, "ROCKER")

    elif db:
        # Shelves2 (KALLAX east of IW9)
        _iw9_e_al, _iw9_e_cw = seg_vecs(layout.iw9.poly[1], layout.iw9.poly[2])
        _sh2_nw = line_isect(
            offset_pt(layout.iw1.poly[0], C.SHELVES_GAP_IW1, w9w10_in), w9w10_al,
            offset_pt(layout.iw9.poly[2], C.SHELVES_GAP_IW9, _iw9_e_cw), _iw9_e_al)
        _sh2_ne = offset_pt(_sh2_nw, C.SHELVES_LENGTH, _iw9_e_cw)
        _sh2_sw = offset_pt(_sh2_nw, C.SHELVES_DEPTH, w9w10_in)
        _sh2_se = offset_pt(_sh2_ne, C.SHELVES_DEPTH, w9w10_in)
        items["shelves2"] = _item("shelves2", "furniture",
                                  [_sh2_sw, _sh2_se, _sh2_ne, _sh2_nw], "SHELVES")

        # ET (east end table, db variant)
        et_r = (C.ET_RADIUS_CM / 2.54) / 12.0
        _et_from_iw1 = offset_pt(layout.iw1.poly[3], C.STD_GAP + et_r, _iw1_n_out)
        _et_t_max = 0
        for _i in range(len(inner_poly)):
            _j = (_i + 1) % len(inner_poly)
            _dx = inner_poly[_j][0] - inner_poly[_i][0]
            _dy = inner_poly[_j][1] - inner_poly[_i][1]
            _det = _iw1_n_al[0] * _dy - _iw1_n_al[1] * _dx
            if abs(_det) < 1e-12:
                continue
            _ox = inner_poly[_i][0] - _et_from_iw1[0]
            _oy = inner_poly[_i][1] - _et_from_iw1[1]
            _t = (_ox * _dy - _oy * _dx) / _det
            _s = (_ox * _iw1_n_al[1] - _oy * _iw1_n_al[0]) / _det
            if 0 <= _s <= 1 and _t > 0 and _t > _et_t_max:
                _et_t_max = _t
        et_cx, et_cy = offset_pt(_et_from_iw1, _et_t_max - C.STD_GAP - et_r, _iw1_n_al)
        items["et_east"] = _circle_item("et_east", "furniture",
                                        (et_cx, et_cy), et_r, "ET")

        # DAYBED
        _db_s_ref = offset_pt(layout.iw1.poly[3], C.STD_GAP, _iw1_n_out)
        _neg_al = (-_iw1_n_al[0], -_iw1_n_al[1])
        db_se = offset_pt((et_cx, et_cy), et_r + 3.0 / 12.0, _neg_al)
        _db_se_proj_out = ((db_se[0] - _db_s_ref[0]) * _iw1_n_out[0] +
                           (db_se[1] - _db_s_ref[1]) * _iw1_n_out[1])
        db_se = offset_pt(db_se, -_db_se_proj_out, _iw1_n_out)
        db_sw = offset_pt(db_se, -DAYBED_W, _iw1_n_al)
        db_ne = offset_pt(db_se, DAYBED_D, _iw1_n_out)
        db_nw = offset_pt(db_sw, DAYBED_D, _iw1_n_out)
        items["daybed"] = _item("daybed", "furniture",
                                [db_sw, db_se, db_ne, db_nw], "DAYBED")

        # ET west (db variant)
        _neg_out = (-_iw1_n_out[0], -_iw1_n_out[1])
        et2_center = offset_pt(offset_pt(db_nw, -(6.0 / 12.0 + et_r), _iw1_n_al),
                               et_r, _neg_out)
        items["et_west"] = _circle_item("et_west", "furniture",
                                        et2_center, et_r, "ET")

        # ROCKER (db variant)
        _ro1_e_pt = _nwp(_iw2_d + C.RO1_OFFSET_FROM_IW2 + C.IW1_RO_WIDTH, 0)
        _ref = layout.iw1.poly[3]
        _db_sw_d_al = (db_sw[0] - _ref[0]) * _iw1_n_al[0] + (db_sw[1] - _ref[1]) * _iw1_n_al[1]
        _ro1_d_al = (_ro1_e_pt[0] - _ref[0]) * _iw1_n_al[0] + (_ro1_e_pt[1] - _ref[1]) * _iw1_n_al[1]
        _rk_d_al = (_db_sw_d_al + _ro1_d_al) / 2 - 8.0 / 12.0
        _fr_s_pt = _nwp(0, 3.0 / 12.0 + STD_FRIDGE_D)
        _fr_door_s_pt = offset_pt(_fr_s_pt, STD_FRIDGE_W, w9w10_in)
        _fr_d_out = (_fr_door_s_pt[0] - _ref[0]) * _iw1_n_out[0] + (_fr_door_s_pt[1] - _ref[1]) * _iw1_n_out[1]
        _rk_d_out = _fr_d_out / 2 + 26.0 / 12.0
        rk_center = offset_pt(offset_pt(_ref, _rk_d_al, _iw1_n_al), _rk_d_out, _iw1_n_out)
        rk_cx, rk_cy = rk_center
        rk_hw = C.ROCKER_DEPTH / 2
        rk_hh = C.ROCKER_WIDTH / 2
        _rk_cr = (-w12w13_al[1], w12w13_al[0])
        rk_poly = [
            (rk_cx - rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
            (rk_cx - rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy - rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] + rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] + rk_hw * _rk_cr[1]),
            (rk_cx + rk_hh * w12w13_al[0] - rk_hw * _rk_cr[0],
             rk_cy + rk_hh * w12w13_al[1] - rk_hw * _rk_cr[1]),
        ]
        items["rocker"] = _item("rocker", "furniture", rk_poly, "ROCKER")

    else:
        # Standard: loveseat, ET, loveseat2
        lv_width = C.LOVESEAT_WIDTH
        lv_height = C.LOVESEAT_LENGTH
        lv_nw = _lwp(C.LOVESEAT_OFFSET_IW4, C.LOVESEAT_OFFSET_IW1)
        _lv_perp = (-w12w13_al[1], w12w13_al[0])
        lv_sw_pt = offset_pt(lv_nw, lv_height, w12w13_al)
        lv_ne_pt = offset_pt(lv_nw, lv_width, _lv_perp)
        lv_se_pt = offset_pt(lv_sw_pt, lv_width, _lv_perp)
        items["loveseat"] = _item("loveseat", "furniture",
                                   [lv_sw_pt, lv_se_pt, lv_ne_pt, lv_nw], "LOVESEAT")

        # ET (standard)
        et_r = (C.ET_RADIUS_CM / 2.54) / 12.0
        et_gap = et_r + C.STD_GAP
        _et_from_iw1 = offset_pt(layout.iw1.poly[3], C.STD_GAP + et_r, _iw1_n_out)
        _et_dx = _et_from_iw1[0] - lv_se_pt[0]
        _et_dy = _et_from_iw1[1] - lv_se_pt[1]
        _et_b = 2 * (_et_dx * _iw1_n_al[0] + _et_dy * _iw1_n_al[1])
        _et_c = _et_dx**2 + _et_dy**2 - et_gap**2
        _et_disc = _et_b**2 - 4 * _et_c
        if _et_disc >= 0:
            _et_t = (-_et_b + math.sqrt(_et_disc)) / 2
            et_cx = _et_from_iw1[0] + _et_t * _iw1_n_al[0]
            et_cy = _et_from_iw1[1] + _et_t * _iw1_n_al[1]
            items["et"] = _circle_item("et", "furniture",
                                       (et_cx, et_cy), et_r, "ET")

            # LOVESEAT2
            _lv2_sw = offset_pt(
                offset_pt((et_cx, et_cy), et_r + C.STD_GAP, _iw1_n_al),
                -et_r, _iw1_n_out)
            _lv2_nw = offset_pt(_lv2_sw, lv_width, _iw1_n_out)
            _lv2_se = offset_pt(_lv2_sw, lv_height, _iw1_n_al)
            _lv2_ne = offset_pt(_lv2_nw, lv_height, _iw1_n_al)
            items["loveseat2"] = _item("loveseat2", "furniture",
                                        [_lv2_sw, _lv2_se, _lv2_ne, _lv2_nw], "LOVESEAT")

    # --- Chair + Ottoman (all non-bare/sf variants) ---
    _ch_theta = math.atan2(w12w13_al[1], w12w13_al[0]) - math.pi / 4
    _ch_along = (math.cos(_ch_theta), math.sin(_ch_theta))
    _ch_cross = (-math.sin(_ch_theta), math.cos(_ch_theta))
    _ch_mid = ((pts["W11"][0] + pts["W12"][0]) / 2,
               (pts["W11"][1] + pts["W12"][1]) / 2)
    _ch_base = offset_pt(offset_pt(_ch_mid, -1.0 / 12.0, w11w12_al),
                         8.0 / 12.0, w11w12_in)
    ch_cx, ch_cy = offset_pt(_ch_base, 4.0 / 12.0, _ch_along)
    _ch_hw = C.CHAIR_WIDTH / 2
    _ch_hd = C.CHAIR_DEPTH / 2
    ch_poly = [
        (ch_cx - _ch_hw * _ch_cross[0] - _ch_hd * _ch_along[0],
         ch_cy - _ch_hw * _ch_cross[1] - _ch_hd * _ch_along[1]),
        (ch_cx + _ch_hw * _ch_cross[0] - _ch_hd * _ch_along[0],
         ch_cy + _ch_hw * _ch_cross[1] - _ch_hd * _ch_along[1]),
        (ch_cx + _ch_hw * _ch_cross[0] + _ch_hd * _ch_along[0],
         ch_cy + _ch_hw * _ch_cross[1] + _ch_hd * _ch_along[1]),
        (ch_cx - _ch_hw * _ch_cross[0] + _ch_hd * _ch_along[0],
         ch_cy - _ch_hw * _ch_cross[1] + _ch_hd * _ch_along[1]),
    ]
    items["chair"] = _item("chair", "furniture", ch_poly, "CHAIR")

    # Ottoman
    ot_dist = 39.0 / 12.0
    ot_cx, ot_cy = offset_pt((ch_cx, ch_cy), ot_dist, _ch_along)
    _ot_hs = C.OTTOMAN_SIZE / 2
    ot_poly = [
        (ot_cx - _ot_hs * _ch_cross[0] - _ot_hs * _ch_along[0],
         ot_cy - _ot_hs * _ch_cross[1] - _ot_hs * _ch_along[1]),
        (ot_cx + _ot_hs * _ch_cross[0] - _ot_hs * _ch_along[0],
         ot_cy + _ot_hs * _ch_cross[1] - _ot_hs * _ch_along[1]),
        (ot_cx + _ot_hs * _ch_cross[0] + _ot_hs * _ch_along[0],
         ot_cy + _ot_hs * _ch_cross[1] + _ot_hs * _ch_along[1]),
        (ot_cx - _ot_hs * _ch_cross[0] + _ot_hs * _ch_along[0],
         ot_cy - _ot_hs * _ch_cross[1] + _ot_hs * _ch_along[1]),
    ]
    items["ottoman"] = _item("ottoman", "furniture", ot_poly, "OTTO")

    # --- Desk + Desk chair ---
    w16w17_al, w16w17_in = seg_vecs(pts["W16"], pts["W17"])
    _neg_w16w17_al = (-w16w17_al[0], -w16w17_al[1])
    dk_sw_pt = pts["W17"]
    dk_se_pt = offset_pt(dk_sw_pt, C.DESK_WIDTH, _neg_w16w17_al)
    dk_nw_pt = offset_pt(dk_sw_pt, C.DESK_DEPTH, w16w17_in)
    dk_ne_pt = offset_pt(dk_se_pt, C.DESK_DEPTH, w16w17_in)
    items["desk"] = _item("desk", "furniture",
                           [dk_sw_pt, dk_se_pt, dk_ne_pt, dk_nw_pt], "DESK")

    dc_sw_pt = offset_pt(
        offset_pt(dk_sw_pt, C.DESK_WIDTH / 2 - C.DESK_CHAIR_WIDTH / 2, _neg_w16w17_al),
        C.DESK_DEPTH + C.DESK_CHAIR_GAP, w16w17_in)
    dc_se_pt = offset_pt(dc_sw_pt, C.DESK_CHAIR_WIDTH, _neg_w16w17_al)
    dc_nw_pt = offset_pt(dc_sw_pt, C.DESK_CHAIR_DEPTH, w16w17_in)
    dc_ne_pt = offset_pt(dc_se_pt, C.DESK_CHAIR_DEPTH, w16w17_in)
    items["desk_chair"] = _item("desk_chair", "furniture",
                                 [dc_sw_pt, dc_se_pt, dc_ne_pt, dc_nw_pt], "CHAIR")

    return items
