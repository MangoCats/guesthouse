"""SVG transform factory and page constants."""

W, H = 792, 612
_s = (368.79 - 151.26) / 18.66

def make_svg_transform(p3_trav):
    """Create to_svg closure from P3 traverse position."""
    px = 368.79 + p3_trav[0] * _s
    py = 124.12 - p3_trav[1] * _s
    def to_svg(e, n):
        return (px + e * _s, py - n * _s)
    return to_svg
