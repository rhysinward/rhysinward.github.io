#!/usr/bin/env python3
"""
Render the home-page animation: DENV-2 genotype II (Cosmopolitan) moving
through South and Southeast Asia and southern China, from a discrete-trait
(DTA) phylogeography MCC tree.

The approach follows Gytis Dudas's baltic animation examples (curonia,
https://github.com/evogytis/baltic): lineages that change location on the tree
travel along Bezier paths between location centres, circles track how many
lineages are present in each place, and the time-scaled tree is drawn in step
with the map. This is a fresh implementation styled for rhysinward.github.io:

  - an equal-area (Lambert azimuthal) projection instead of raw lon/lat
  - the site's palette: five regional hues, each location a shade of its hue
  - paths taper from tail to head and blend from origin to destination colour,
    with a paper-coloured halo (as in Gytis's work) so they read over land
  - fainter paths where the destination location is less certain (the MCC
    tree's location.prob)
  - circle area, not radius, proportional to the number of lineages
  - the tree band fills with colour as time passes; the loop ends on a
    crossfade back to the first frame

Usage (from anywhere):
  python3 render.py --data /path/to/Beast_Vietnam/visulisation --out ./out

Then encode with encode.sh. Needs baltic, numpy, scipy, matplotlib,
geopandas, shapely and pyproj, plus ffmpeg for encoding.
"""

import argparse
import json
import math
import os
import sys
import types
from pathlib import Path

import numpy as np

# Homebrew's Python 3.14 can ship with a pyexpat that fails to load, which
# breaks matplotlib's macOS font discovery (it parses a plist). Stub plistlib
# just long enough for matplotlib to fall back to scanning font folders.
class _InvalidPlist(Exception):
    pass

_stub = types.ModuleType("plistlib")
_stub.InvalidFileException = _InvalidPlist
_stub.loads = lambda *a, **k: (_ for _ in ()).throw(_InvalidPlist())
sys.modules["plistlib"] = _stub
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.font_manager  # noqa: E402,F401
del sys.modules["plistlib"]

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects as pe  # noqa: E402
from matplotlib.collections import LineCollection, PatchCollection  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Polygon as MplPolygon  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

import baltic as bt  # noqa: E402
import geopandas as gpd  # noqa: E402
from shapely.geometry import LineString, MultiPolygon, Polygon, box  # noqa: E402

# --------------------------------------------------------------------------
# Look

PAPER = "#FBFAF7"
INK = "#1D2930"
INK2 = "#56636B"
LAND = "#ECEAE3"        # countries outside the analysis
COAST = "#CFCCC2"
GRATICULE = "#EEECE6"

HUES = {                # site tokens: --G, --A, --T, --C, --P
    "South Asia": "#8F5F12",
    "Mainland Southeast Asia": "#2C7652",
    "Vietnam": "#B2433A",
    "Maritime Southeast Asia": "#2B5C8A",
    "China and Taiwan": "#7A4D80",
}
GROUPS = {
    "India": "South Asia",
    "Myanmar": "Mainland Southeast Asia", "Thailand": "Mainland Southeast Asia",
    "Laos": "Mainland Southeast Asia", "Cambodia": "Mainland Southeast Asia",
    "Northern_Vietnam": "Vietnam", "Central_Vietnam": "Vietnam", "Southern_Vietnam": "Vietnam",
    "Singapore": "Maritime Southeast Asia", "Malaysia": "Maritime Southeast Asia",
    "Indonesia": "Maritime Southeast Asia", "Philippines": "Maritime Southeast Asia",
    "Taiwan": "China and Taiwan", "Guangdong": "China and Taiwan", "Yunnan": "China and Taiwan",
    "Zhejiang": "China and Taiwan", "Hunan": "China and Taiwan", "Shandong": "China and Taiwan",
    "Shanghai": "China and Taiwan",
}
# Shade of the group hue for each location: <0 lighter, >0 darker.
SHADE = {
    "Myanmar": -0.30, "Thailand": 0.0, "Laos": -0.45, "Cambodia": 0.25,
    "Northern_Vietnam": 0.25, "Central_Vietnam": -0.35, "Southern_Vietnam": 0.0,
    "Singapore": 0.0, "Malaysia": -0.30, "Indonesia": 0.28, "Philippines": -0.50,
    "Taiwan": -0.35, "Guangdong": 0.0, "Yunnan": 0.28, "Zhejiang": -0.55,
    "Hunan": -0.20, "Shandong": 0.40, "Shanghai": -0.70,
}
LABELS = {  # map label: (anchor location, dx km, dy km, ha); shown at about half size on the page
    "India": ("India", 120, -420, "center"),
    "Myanmar": ("Myanmar", -120, 420, "center"),
    "Thailand": ("Thailand", -260, 120, "right"),
    "Vietnam": ("Central_Vietnam", 210, 30, "left"),
    "Malaysia": ("Malaysia", -170, 230, "right"),
    "Singapore": ("Singapore", -120, -170, "right"),
    "Indonesia": ("Indonesia", 380, -150, "left"),
    "Philippines": ("Philippines", 230, 160, "left"),
    "Taiwan": ("Taiwan", 160, 30, "left"),
    "China": ("Hunan", -150, 260, "center"),
}
SERIF = "Iowan Old Style"
SANS = "Seravek"

PROJ = "+proj=laea +lat_0=14 +lon_0=99 +datum=WGS84 +units=km +no_defs"


def mix(a, b, t):
    a, b = np.array(to_rgb(a)), np.array(to_rgb(b))
    return tuple(a + (b - a) * t)


def loc_colour(loc):
    base = HUES[GROUPS[loc]]
    s = SHADE.get(loc, 0.0)
    return mix(base, PAPER, -s) if s < 0 else mix(base, INK, s * 0.6)


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


# --------------------------------------------------------------------------
# Data

def load_tree(path, rename):
    tree = bt.loadNexus(str(path), absoluteTime=True)
    tree.sortBranches(descending=True)
    for k in tree.Objects:
        loc = k.traits.get("location")
        if loc in rename:
            k.traits["location"] = rename[loc]
    return tree


def load_centres(path):
    centres = {}
    for line in open(path):
        parts = line.strip().split("\t")
        if len(parts) == 3:
            name, lat, lon = parts
            centres[name] = (float(lon), float(lat))
    return centres


def project_points(lonlat):
    pts = gpd.GeoSeries(gpd.points_from_xy([p[0] for p in lonlat], [p[1] for p in lonlat]), crs=4326).to_crs(PROJ)
    return [(p.x, p.y) for p in pts]


def polygons_to_patches(geom):
    """Shapely (Multi)Polygon -> list of matplotlib polygons (exteriors only)."""
    if geom is None or geom.is_empty:
        return []
    parts = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    return [MplPolygon(np.asarray(p.exterior.coords), closed=True) for p in parts if p.area > 40]


def bezier(a, b, n, bend):
    """Quadratic Bezier from a to b, bowed to the left of the direction of travel,
    so A->B and B->A curve apart (as in Gytis's control-point construction)."""
    a, b = np.asarray(a), np.asarray(b)
    d = np.linalg.norm(b - a)
    mid = (a + b) / 2
    perp = np.array([-(b - a)[1], (b - a)[0]]) / max(d, 1e-9)
    c = mid + perp * (bend * d + 90.0)
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * a + 2 * (1 - t) * t * c + t ** 2 * b


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Beast_Vietnam/visulisation folder")
    ap.add_argument("--tree", default="Dengue_2_AL_DENV_2II_dta_MCC.tre")
    ap.add_argument("--out", required=True)
    ap.add_argument("--size", type=int, default=1080)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--start", type=float, default=1998.0)
    ap.add_argument("--end", type=float, default=None, help="default: last tip")
    ap.add_argument("--sec-per-year", type=float, default=1.0)
    ap.add_argument("--hold", type=float, default=1.4, help="seconds to hold the last frame")
    ap.add_argument("--fade", type=float, default=0.8, help="seconds of crossfade back to the start")
    ap.add_argument("--only", type=str, default=None, help="comma-separated decimal years to render as stills")
    args = ap.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tree = load_tree(data / args.tree, {"North_Vietnam": "Northern_Vietnam"})
    leaves = [k for k in tree.Objects if k.branchType == "leaf"]
    last_tip = max(k.absoluteTime for k in leaves)
    end = args.end or last_tip + 0.08
    locs = sorted({k.traits["location"] for k in tree.Objects})
    missing = [l for l in locs if l not in GROUPS]
    assert not missing, f"no colour group for {missing}"

    # --- map geometry
    lon0, lon1, lat0, lat1 = 66.0, 128.0, -11.5, 39.0
    frame_ll = gpd.GeoSeries([box(lon0, lat0, lon1, lat1)], crs=4326).to_crs(PROJ)
    fx0, fy0, fx1, fy1 = frame_ll.total_bounds
    study = gpd.read_file(data / "combined_asia_with_all_china.geojson").to_crs(PROJ)
    study["geometry"] = study.geometry.buffer(0).simplify(3.0)

    # location centres: capitals from the tsv, else a point inside the polygon
    centres_ll = load_centres(data / "countries_sea_lat_long.tsv")
    have = [l for l in locs if l in centres_ll]
    centres = dict(zip(have, project_points([centres_ll[l] for l in have])))
    for l in locs:
        if l not in centres:
            g = study[study["location"] == l].geometry.union_all().representative_point()
            centres[l] = (g.x, g.y)
            print(f"no centre for {l} in the tsv; using a point inside its polygon")
    context = gpd.read_file(data / "asia.geojson").to_crs(PROJ)
    context["geometry"] = context.geometry.buffer(0).simplify(3.0)

    W = H = args.size
    dpi = 100
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi, facecolor=PAPER)
    map_ax = fig.add_axes([0, 0.25, 1, 0.75])
    tree_ax = fig.add_axes([0.075, 0.045, 0.87, 0.165])
    for ax in (map_ax, tree_ax):
        ax.set_facecolor(PAPER)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])

    # fit the frame around the places that must show, then pad to the axes' aspect
    keep = project_points([(68.0, 23.5), (77.5, 7.5), (80.5, 6.0), (106.0, -8.6), (114.0, -8.8),
                           (126.6, 7.0), (122.2, 25.5), (122.5, 37.6), (97.0, 28.5)])
    kx, ky = np.array(keep).T
    pad = 0.035
    x0, x1 = kx.min(), kx.max()
    y0, y1 = ky.min(), ky.max()
    x0, x1 = x0 - (x1 - x0) * pad, x1 + (x1 - x0) * pad
    y0, y1 = y0 - (y1 - y0) * pad, y1 + (y1 - y0) * pad
    aspect = (0.75 * H) / W
    if (y1 - y0) / (x1 - x0) < aspect:
        grow = (aspect * (x1 - x0) - (y1 - y0)) / 2
        y0, y1 = y0 - grow, y1 + grow
    else:
        grow = ((y1 - y0) / aspect - (x1 - x0)) / 2
        x0, x1 = x0 - grow, x1 + grow
    map_ax.set_xlim(x0, x1)
    map_ax.set_ylim(y0, y1)
    map_ax.set_aspect("equal", adjustable="box")
    half_w = (x1 - x0) / 2

    # graticule every 10 degrees, clipped to the frame
    frame_box = box(x0, y0, x1, y1)
    lines = [LineString([(lon, y) for y in np.linspace(-20, 50, 60)]) for lon in range(60, 140, 10)]
    lines += [LineString([(x, lat) for x in np.linspace(55, 140, 80)]) for lat in range(-10, 50, 10)]
    for g in gpd.GeoSeries(lines, crs=4326).to_crs(PROJ):
        g = g.intersection(frame_box)
        for part in getattr(g, "geoms", [g]):
            if not part.is_empty:
                map_ax.plot(*part.xy, color=GRATICULE, lw=0.8, zorder=0)
    px = (2 * half_w) / W  # km per pixel

    # context land
    ctx_patches = [p for geom in context.geometry for p in polygons_to_patches(geom)]
    map_ax.add_collection(PatchCollection(ctx_patches, facecolor=LAND, edgecolor=COAST, lw=0.45, zorder=1))

    # study polygons: one collection per location so they can tint over time
    study_patch = {}
    for name, geom in zip(study["location"], study.geometry):
        patches = polygons_to_patches(geom)
        if not patches:
            continue
        if name in GROUPS:
            base = mix(PAPER, HUES[GROUPS[name]], 0.06)
        else:
            base = mix(LAND, PAPER, 0.25)
        pc = PatchCollection(patches, facecolor=base, edgecolor=PAPER, lw=0.9, zorder=2)
        map_ax.add_collection(pc)
        if name in GROUPS:
            study_patch.setdefault(name, []).append(pc)
    # crisp outer coastline over the white borders
    map_ax.add_collection(PatchCollection([p for geom in study.geometry for p in polygons_to_patches(geom)],
                                          facecolor="none", edgecolor=COAST, lw=0.45, zorder=3))

    halo = [pe.withStroke(linewidth=5.5, foreground=PAPER)]
    for label, (anchor, dx, dy, ha) in LABELS.items():
        x, y = centres[anchor]
        map_ax.text(x + dx, y + dy, label, fontsize=15, fontfamily=SANS, color=INK2, ha=ha, va="center",
                    zorder=40, path_effects=halo)

    year_text = map_ax.text(0.04, 0.07, "", transform=map_ax.transAxes, fontsize=62, fontfamily=SERIF,
                            fontweight="bold", color=INK, zorder=60, path_effects=[pe.withStroke(linewidth=8, foreground=PAPER)])

    # --- lineage presence per location (smoothed counts on the frame grid)
    n_main = int(round((end - args.start) * args.sec_per_year * args.fps)) + 1
    times = np.linspace(args.start, end, n_main)
    branches = [k for k in tree.Objects if k.parent is not None and k.parent.absoluteTime is not None]
    starts = np.array([k.parent.absoluteTime for k in branches])
    ends = np.array([k.absoluteTime for k in branches])
    bl = np.array([k.traits["location"] for k in branches])
    counts = {}
    for loc in locs:
        sel = bl == loc
        c = ((starts[sel][None, :] <= times[:, None]) & (times[:, None] < ends[sel][None, :])).sum(axis=1).astype(float)
        counts[loc] = gaussian_filter1d(c, sigma=0.35 * args.fps * args.sec_per_year)
    cmax = max(c.max() for c in counts.values())

    circles = map_ax.scatter([centres[l][0] for l in locs], [centres[l][1] for l in locs], s=np.zeros(len(locs)),
                             facecolors=[(*loc_colour(l), 0.30) for l in locs],
                             edgecolors=[(*loc_colour(l), 0.95) for l in locs], linewidths=1.3, zorder=30)
    cores = map_ax.scatter([centres[l][0] for l in locs], [centres[l][1] for l in locs], s=np.zeros(len(locs)),
                           c=[loc_colour(l) for l in locs], linewidths=0, zorder=31)

    # --- travelling lineages
    travel = []
    for k in branches:
        if k.parent is tree.root and "location" not in k.parent.traits:
            continue
        a, b = k.parent.traits.get("location"), k.traits.get("location")
        if a and b and a != b:
            T = k.parent.absoluteTime + 0.5 * k.length
            prob = k.traits.get("location.prob", 1.0)
            if isinstance(prob, (list, tuple)):
                prob = max(prob)
            path = bezier(centres[a], centres[b], 64, 0.16)
            travel.append((T, a, b, float(prob), path))
    print(f"{len(travel)} movements, {len(leaves)} tips, root {tree.root.absoluteTime:.2f}, last tip {last_tip:.2f}")
    D = 0.75  # years a movement takes on screen

    arc_halo = LineCollection([], linewidths=[], colors=[], capstyle="round", zorder=34)
    arc_line = LineCollection([], linewidths=[], colors=[], capstyle="round", zorder=35)
    map_ax.add_collection(arc_halo)
    map_ax.add_collection(arc_line)
    heads = map_ax.scatter([], [], s=[], c=[], linewidths=0, zorder=36)

    # --- tree band
    t_lo, t_hi = 1990.0, math.ceil(end)
    ys = {k: k.y for k in tree.Objects}
    yspan = tree.ySpan
    hseg, hcol, vseg, vcol, vx = [], [], [], [], []
    for k in tree.Objects:
        if k.parent is None or k.parent.absoluteTime is None:
            continue
        x0, x1 = max(k.parent.absoluteTime, t_lo), k.absoluteTime
        if x1 < t_lo:
            continue
        hseg.append(((x0, k.y), (x1, k.y)))
        hcol.append(loc_colour(k.traits["location"]))
        if k.branchType == "node":
            yl, yr = k.children[0].y, k.children[-1].y
            vseg.append(((x1, min(yl, yr)), (x1, max(yl, yr))))
            vcol.append(loc_colour(k.traits["location"]))
            vx.append(x1)
    hseg = np.array(hseg)
    hcol = np.array(hcol)
    vseg = np.array(vseg)
    vcol = np.array(vcol)
    vx = np.array(vx)
    tree_ax.add_collection(LineCollection(np.concatenate([hseg, vseg]), colors="#E4E2DB", linewidths=0.55, zorder=1))
    past_h = LineCollection([], linewidths=0.75, zorder=2)
    past_v = LineCollection([], linewidths=0.75, zorder=2)
    tree_ax.add_collection(past_h)
    tree_ax.add_collection(past_v)
    cursor = tree_ax.axvline(args.start, color=INK, lw=1.4, zorder=5)
    tree_ax.set_xlim(t_lo, t_hi)
    tree_ax.set_ylim(-yspan * 0.02, yspan * 1.02)
    for yr in range(int(t_lo), int(t_hi) + 1, 5):
        if yr % 10 == 0:
            tree_ax.text(yr, -yspan * 0.09, str(yr), fontsize=14, fontfamily=SANS, color=INK2, ha="center", va="top")
        tree_ax.plot([yr, yr], [-yspan * 0.035, -yspan * 0.01], color=COAST, lw=0.8, clip_on=False)
    tree_ax.set_xticks([])
    tree_ax.set_yticks([])

    def draw(fi, t):
        # year
        year_text.set_text(f"{int(math.floor(t))}")

        # presence: circles and polygon tint
        idx = min(int(round((t - args.start) / (end - args.start) * (n_main - 1))), n_main - 1) if t >= args.start else 0
        vals = np.array([counts[l][idx] for l in locs])
        area = (vals / cmax) * 2600.0  # points^2 at full presence
        circles.set_sizes(area)
        cores.set_sizes(np.where(vals > 0.05, 9.0, 0.0))
        for l, v in zip(locs, vals):
            if l in study_patch:
                tint = 0.06 + 0.12 * math.sqrt(min(1.0, v / (0.35 * cmax)))
                for pc in study_patch[l]:
                    pc.set_facecolor(mix(PAPER, HUES[GROUPS[l]], tint))

        # movements
        segs, widths, cols, hsegs, hwidths, hcols = [], [], [], [], [], []
        hx, hy, hs, hc = [], [], [], []
        for T, a, b, prob, path in travel:
            p = (t - (T - D)) / D
            if p <= 0 or p >= 1.75:
                continue
            head = smoothstep(min(p, 1.0))
            tail = smoothstep((p - 0.55) / 1.2)
            if head - tail < 0.01:
                continue
            n = 22
            fr = np.linspace(tail, head, n)
            ii = fr * (len(path) - 1)
            pts = np.column_stack([np.interp(ii, np.arange(len(path)), path[:, 0]),
                                   np.interp(ii, np.arange(len(path)), path[:, 1])])
            alpha = 0.35 + 0.65 * max(0.0, min(1.0, (prob - 0.3) / 0.6))
            ca, cb = np.array(loc_colour(a)), np.array(loc_colour(b))
            for q in range(n - 1):
                f = (q + 1) / (n - 1)
                segs.append(pts[q:q + 2])
                widths.append(0.5 + 2.8 * f ** 1.4)
                cols.append((*(ca + (cb - ca) * fr[q + 1]), alpha))
                hsegs.append(pts[q:q + 2])
                hwidths.append(2.8 + 2.8 * f ** 1.4)
                hcols.append((*to_rgb(PAPER), 0.85 * alpha))
            if p < 1.0:
                hx.append(pts[-1, 0]); hy.append(pts[-1, 1]); hs.append(22.0); hc.append((*cb, alpha))
        arc_line.set_segments(segs); arc_line.set_linewidths(widths); arc_line.set_colors(cols)
        arc_halo.set_segments(hsegs); arc_halo.set_linewidths(hwidths); arc_halo.set_colors(hcols)
        heads.set_offsets(np.column_stack([hx, hy]) if hx else np.empty((0, 2)))
        heads.set_sizes(hs); heads.set_facecolors(hc if hc else [])

        # tree band
        on = hseg[:, 0, 0] < t
        h = hseg[on].copy()
        h[:, 1, 0] = np.minimum(h[:, 1, 0], t)
        past_h.set_segments(h); past_h.set_colors(hcol[on])
        von = vx <= t
        past_v.set_segments(vseg[von]); past_v.set_colors(vcol[von])
        cursor.set_xdata([t, t])

    frames = out / "frames"
    frames.mkdir(exist_ok=True)
    if args.only:
        for y in [float(v) for v in args.only.split(",")]:
            draw(0, y)
            fig.savefig(out / f"still_{y:.2f}.png", dpi=dpi, facecolor=PAPER)
            print("still", y)
        return

    for fi, t in enumerate(times):
        draw(fi, t)
        fig.savefig(frames / f"f_{fi:05d}.png", dpi=dpi, facecolor=PAPER)
        if fi % 60 == 0:
            print(f"frame {fi}/{len(times)}  {t:.2f}", flush=True)

    # hold, then crossfade back to the first frame so the loop is seamless
    from PIL import Image
    last = Image.open(frames / f"f_{len(times) - 1:05d}.png").convert("RGB")
    first = Image.open(frames / "f_00000.png").convert("RGB")
    fi = len(times)
    for _ in range(int(args.hold * args.fps)):
        last.save(frames / f"f_{fi:05d}.png"); fi += 1
    nf = int(args.fade * args.fps)
    for j in range(1, nf + 1):
        Image.blend(last, first, smoothstep(j / (nf + 1))).save(frames / f"f_{fi:05d}.png"); fi += 1
    print(f"done: {fi} frames in {frames}")


if __name__ == "__main__":
    main()
