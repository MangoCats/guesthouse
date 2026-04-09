"""
Least-squares adjustment of pentagon (POB, P2, P3, P4, P5) to fit field measurements.

Theoretical distances from distances.md are internally consistent (verified to 0.03").
Field measurements of 4 pairs differ; remaining 6 matched within ~1".
This script finds adjusted coordinates minimizing total distance residuals.
"""

import math
import numpy as np
from scipy.optimize import least_squares


def dms_to_deg(d, m, s):
    return d + m / 60.0 + s / 3600.0


def ft_in(ft, inches):
    return ft + inches / 12.0


# --- Theoretical distances and bearings from distances.md ---
theoretical = {
    ("POB", "P2"): (ft_in(19, 6.2), dms_to_deg(258, 25, 11.4)),
    ("POB", "P3"): (ft_in(38, 0.8), dms_to_deg(210, 8, 49.5)),
    ("POB", "P4"): (ft_in(35, 3.5), dms_to_deg(158, 51, 47.3)),
    ("POB", "P5"): (ft_in(34, 11.5), dms_to_deg(137, 11, 43.6)),
    ("P2", "P3"):  (ft_in(29, 0.0), dms_to_deg(180, 0, 0)),
    ("P2", "P4"):  (ft_in(43, 0.8), dms_to_deg(132, 19, 26.3)),
    ("P2", "P5"):  (ft_in(48, 0.8), dms_to_deg(116, 52, 45.5)),
    ("P3", "P4"):  (ft_in(31, 10.1), dms_to_deg(90, 0, 0)),
    ("P3", "P5"):  (ft_in(43, 5.8), dms_to_deg(80, 22, 35.5)),
    ("P4", "P5"):  (ft_in(13, 2.5), dms_to_deg(56, 36, 31.2)),
}

# --- Field measurements (where they differ from theoretical) ---
field_overrides = {
    ("P2", "P3"): ft_in(28, 11.0),
    ("P2", "P4"): ft_in(42, 11.0),
    ("P3", "P5"): ft_in(43, 8.5),
    ("P4", "P5"): ft_in(13, 5.0),
}

# Build observed distances: field where available, theoretical otherwise
pairs = list(theoretical.keys())
observed = {}
for pair in pairs:
    if pair in field_overrides:
        observed[pair] = field_overrides[pair]
    else:
        observed[pair] = theoretical[pair][0]

# --- Theoretical coordinates (P3 at origin, from earlier verification) ---
# Derived by placing P3=(0,0) and using bearings from distances.md
def bearing_to_en(dist, bearing_deg):
    """Convert distance + bearing (CW from N) to (Easting, Northing)."""
    rad = math.radians(bearing_deg)
    return (dist * math.sin(rad), dist * math.cos(rad))


P3_theo = np.array([0.0, 0.0])

d, b = theoretical[("P2", "P3")]
P2_theo = np.array(bearing_to_en(d, (b + 180) % 360))  # reverse bearing

d, b = theoretical[("P3", "P4")]
P4_theo = np.array(bearing_to_en(d, b))

d, b = theoretical[("P3", "P5")]
P5_theo = np.array(bearing_to_en(d, b))

d, b = theoretical[("POB", "P3")]
POB_theo = np.array(bearing_to_en(d, (b + 180) % 360))  # reverse bearing

theo_coords = {"POB": POB_theo, "P2": P2_theo, "P3": P3_theo, "P4": P4_theo, "P5": P5_theo}

print("=== Theoretical Coordinates (E, N) feet ===")
for name, coord in theo_coords.items():
    print(f"  {name}: E={coord[0]:.4f}, N={coord[1]:.4f}")

# --- Least-squares optimization ---
# Free parameters: POB(E,N), P2(E,N), P4(E,N), P5(E,N) — 8 values
# P3 is fixed at origin.
# Also fix the orientation by constraining P4 to lie on the E axis (N=0),
# reducing to 7 free parameters.

point_names = ["POB", "P2", "P3", "P4", "P5"]


def pack(coords):
    """Pack free coordinates into parameter vector.
    Free: POB(E,N), P2(E,N), P4(E), P5(E,N) = 7 params.
    P3 fixed at (0,0). P4 N-coordinate fixed at 0."""
    return np.array([
        coords["POB"][0], coords["POB"][1],
        coords["P2"][0], coords["P2"][1],
        coords["P4"][0],  # P4 N fixed at 0
        coords["P5"][0], coords["P5"][1],
    ])


def unpack(x):
    """Unpack parameter vector to coordinate dict."""
    return {
        "POB": np.array([x[0], x[1]]),
        "P2":  np.array([x[2], x[3]]),
        "P3":  np.array([0.0, 0.0]),
        "P4":  np.array([x[4], 0.0]),
        "P5":  np.array([x[5], x[6]]),
    }


def residuals(x):
    """Compute distance residuals: computed - observed."""
    coords = unpack(x)
    res = []
    for pair in pairs:
        p1, p2 = pair
        d_computed = np.linalg.norm(coords[p2] - coords[p1])
        d_observed = observed[pair]
        res.append(d_computed - d_observed)
    return np.array(res)


x0 = pack(theo_coords)
result = least_squares(residuals, x0, method="lm")
adj_coords = unpack(result.x)

print("\n=== Adjusted Coordinates (E, N) feet ===")
for name in point_names:
    e, n = adj_coords[name]
    print(f"  {name}: E={e:.4f}, N={n:.4f}")

# --- Point displacements ---
print("\n=== Point Displacements (theoretical -> adjusted) ===")
for name in point_names:
    de = adj_coords[name][0] - theo_coords[name][0]
    dn = adj_coords[name][1] - theo_coords[name][1]
    dist = math.sqrt(de**2 + dn**2)
    print(f"  {name}: dE={de*12:+.2f}\", dN={dn*12:+.2f}\", total={dist*12:.2f}\"")

# --- Distance & bearing table ---
def en_to_bearing(de, dn):
    """Easting/Northing delta to bearing (CW from N) in degrees."""
    return math.degrees(math.atan2(de, dn)) % 360


def deg_to_dms(deg):
    """Convert decimal degrees to (d, m, s) tuple."""
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return d, m, s


def fmt_dms(deg):
    d, m, s = deg_to_dms(deg)
    return f"{d}deg {m:02d}' {s:04.1f}\""


def fmt_ft_in(feet):
    ft = int(feet)
    inches = (feet - ft) * 12
    # handle rounding to avoid e.g. 12.0"
    if inches >= 11.95:
        ft += 1
        inches = 0.0
    return f"{ft}' {inches:.1f}\""


print("\n=== Adjusted Distance & Bearing Table ===")
print(f"{'Pair':<12} {'Distance':>12} {'Fwd Bearing':>18} {'Rev Bearing':>18}")
print("-" * 62)
for pair in pairs:
    p1, p2 = pair
    de = adj_coords[p2][0] - adj_coords[p1][0]
    dn = adj_coords[p2][1] - adj_coords[p1][1]
    dist = math.sqrt(de**2 + dn**2)
    fwd = en_to_bearing(de, dn)
    rev = (fwd + 180) % 360
    print(f"  {p1} -{p2:<5} {fmt_ft_in(dist):>12} {fmt_dms(fwd):>18} {fmt_dms(rev):>18}")

# --- Residuals ---
print("\n=== Residuals (adjusted distance - observed) ===")
final_res = residuals(result.x)
print(f"{'Pair':<12} {'Observed':>12} {'Adjusted':>12} {'Residual':>10} {'Source':>10}")
print("-" * 58)
for i, pair in enumerate(pairs):
    p1, p2 = pair
    d_obs = observed[pair]
    d_adj = d_obs + final_res[i]
    source = "field" if pair in field_overrides else "theo"
    print(f"  {p1}-{p2:<5} {fmt_ft_in(d_obs):>12} {fmt_ft_in(d_adj):>12} {final_res[i]*12:>+8.2f}\"   {source:>6}")

rms = math.sqrt(np.mean(final_res**2)) * 12
print(f"\n  RMS residual: {rms:.2f}\"")
print(f"  Max residual: {max(abs(final_res))*12:.2f}\"")
