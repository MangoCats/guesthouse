"""Per-element style management: defaults, validation, and resolution.

Style properties are stored in the element ``properties`` JSON column.
All keys are optional — absent means "use CSS class default".

Keys: fill_color, stroke_color, stroke_width, stroke_style, opacity.

Per-view overrides are stored under ``view_overrides`` as a dict keyed by
variant name, each value a partial style dict that overrides the base.

Product URLs are stored as ``product_url`` (string or None).
"""

import re

STYLE_KEYS = ("fill_color", "stroke_color", "stroke_width", "stroke_style", "opacity")

VALID_STROKE_STYLES = ("solid", "dashed", "dotted")

# Match current CSS class defaults exactly
TYPE_DEFAULTS = {
    "appliance": {
        "fill_color": "#2a3a4a",
        "stroke_color": "#4682B4",
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "furniture": {
        "fill_color": "#2a3a2a",
        "stroke_color": "#5a8a5a",
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "fixture": {
        "fill_color": "#3a2a3a",
        "stroke_color": "#8a5a8a",
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "wall": {
        "fill_color": "#334",
        "stroke_color": "#556",
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "opening": {
        "fill_color": "#2a4a6a",
        "stroke_color": "#4488cc",
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "dimension": {
        "fill_color": None,
        "stroke_color": None,
        "stroke_width": 0.02,
        "stroke_style": "solid",
        "opacity": 100,
    },
    "label": {
        "fill_color": None,
        "stroke_color": None,
        "stroke_width": None,
        "stroke_style": None,
        "opacity": 100,
    },
}

_HEX_RE = re.compile(r"^#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$")


def validate_color(value):
    """Return True if *value* is None (use default) or a valid CSS hex colour."""
    if value is None:
        return True
    if not isinstance(value, str):
        return False
    return bool(_HEX_RE.match(value))


def validate_opacity(value):
    """Return True if *value* is a number in [0, 100]."""
    if not isinstance(value, (int, float)):
        return False
    return 0 <= value <= 100


def validate_stroke_style(value):
    """Return True if *value* is a recognised stroke style string."""
    return value in VALID_STROKE_STYLES


def validate_stroke_width(value):
    """Return True if *value* is a non-negative number."""
    if not isinstance(value, (int, float)):
        return False
    return value >= 0


def validate_style_props(props):
    """Validate all style-related keys in *props*.

    Returns ``(True, None)`` if valid, ``(False, error_message)`` otherwise.
    Only checks keys that are present — absent keys are fine.
    """
    if not isinstance(props, dict):
        return False, "properties must be a dict"

    fc = props.get("fill_color")
    if "fill_color" in props and not validate_color(fc):
        return False, f"invalid fill_color: {fc!r}"

    sc = props.get("stroke_color")
    if "stroke_color" in props and not validate_color(sc):
        return False, f"invalid stroke_color: {sc!r}"

    sw = props.get("stroke_width")
    if "stroke_width" in props and sw is not None and not validate_stroke_width(sw):
        return False, f"invalid stroke_width: {sw!r}"

    ss = props.get("stroke_style")
    if "stroke_style" in props and ss is not None and not validate_stroke_style(ss):
        return False, f"invalid stroke_style: {ss!r}"

    op = props.get("opacity")
    if "opacity" in props and op is not None and not validate_opacity(op):
        return False, f"invalid opacity: {op!r}"

    return True, None


def get_defaults(element_type):
    """Return a copy of the default style dict for *element_type*."""
    defaults = TYPE_DEFAULTS.get(element_type)
    if defaults is None:
        return {
            "fill_color": None,
            "stroke_color": None,
            "stroke_width": 0.02,
            "stroke_style": "solid",
            "opacity": 100,
        }
    return dict(defaults)


def resolve_style(element_type, props, variant=None):
    """Merge style from three layers: type defaults → base props → view overrides.

    *props* is the element's ``properties`` dict (or None).
    *variant* is the current variant name (or None for base style only).

    Returns a dict with all five STYLE_KEYS populated.
    """
    result = get_defaults(element_type)

    if props:
        for key in STYLE_KEYS:
            if key in props and props[key] is not None:
                result[key] = props[key]

        if variant:
            view_ov = props.get("view_overrides")
            if isinstance(view_ov, dict):
                variant_ov = view_ov.get(variant)
                if isinstance(variant_ov, dict):
                    for key in STYLE_KEYS:
                        if key in variant_ov and variant_ov[key] is not None:
                            result[key] = variant_ov[key]

    return result
