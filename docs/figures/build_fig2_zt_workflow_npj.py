#!/usr/bin/env python3
"""Fig. 2 of the AutoZT paper (npj Comput. Mater.): the zT workflow, 180 mm wide vector SVG.

Layout:
    a-c  phonon highlights (kl-dft-cpu: random displacements, GPU force-constant fit, GPU phonon BTE)
    centre: POSCAR -> phonon chain and electron chain -> zT(n, T)
    d-f  electron highlights (ke-dft-cpu: scattering inputs, deformation potential, carrier transport)
Units: 1 user unit = 0.1 mm, so the 1800-unit width prints at 180 mm (double column).
Type sizes: 1 pt = 3.528 units; nothing but sub/superscripts is set below 5 pt (17.6 units).
PowerPoint: insert the SVG and use Convert to Shape; the primitives are the ones of the Fig. 1 script
(build_fig1_overview_npj.py), chosen so that sub/superscripts, arrows and rounded boxes survive it.
Usage: python3 build_fig2_zt_workflow_npj.py [out.svg]   (reads fig1_logos.json next to this script; needs Pillow)
"""
import bisect
import json
import math
import os
import re
import sys
from xml.sax.saxutils import escape

from PIL import ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'fig1_logos.json')) as _fh:
    LOGOS = json.load(_fh)

# ---------------------------------------------------------------- tokens
FONT = 'Arial'          # a single family name, which PowerPoint keeps; fontconfig maps it to Liberation Sans
PT = 3.528
FS_LETTER = 9 * PT      # panel letters
FS_TITLE = 8 * PT       # panel titles
FS_H = 7 * PT           # headers inside panels
FS_B = 6 * PT           # body
FS_S = 5.4 * PT         # small
FS_XS = 5.02 * PT       # minimum (5 pt after rounding)

INK = '#1D2733'
MUTED = '#5B6673'
FAINT = '#8A96A3'
CODE = '#2B2F36'
PANEL_BG = '#F7F8FA'
PANEL_LINE = '#D6DCE3'
PURPLE = '#5B3FA0'
LAV = '#EFEBF9'
LAV2 = '#F5F3FB'
LAV_LINE = '#D8D2EA'
SLATE = '#2F4A5E'
GREEN = '#238266'
GREEN_BG = '#EEF7F2'
GREEN_SOFT = '#8FB9A5'
BLUE = '#2F5FA7'
BLUE_TX = '#1F4E9A'
BLUE_BG = '#EDF3FC'
ORANGE = '#B8660F'
ORANGE_TX = '#A3570A'
ORANGE_BG = '#FDF3E7'
CRIMSON = '#B8244A'
CRIMSON_BG = '#FCF0F3'
RED = '#A3203F'
PEACH = '#F8E0D1'
PEACH_TX = '#7A3410'
USER_BG = '#F6E1E8'
GATE = '#3A4654'

W = 1800
H = 1300

# ---------------------------------------------------------------- text metrics (Arial, or Liberation Sans = same metrics)
_FONTS = {}
_FONT_FILES = {
    'Regular': ['/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf', '/mnt/c/Windows/Fonts/arial.ttf',
                'C:/Windows/Fonts/arial.ttf', '/System/Library/Fonts/Supplemental/Arial.ttf'],
    'Bold': ['/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf', '/mnt/c/Windows/Fonts/arialbd.ttf',
             'C:/Windows/Fonts/arialbd.ttf', '/System/Library/Fonts/Supplemental/Arial Bold.ttf'],
    'Italic': ['/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf', '/mnt/c/Windows/Fonts/ariali.ttf',
               'C:/Windows/Fonts/ariali.ttf', '/System/Library/Fonts/Supplemental/Arial Italic.ttf'],
    'BoldItalic': ['/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf', '/mnt/c/Windows/Fonts/arialbi.ttf',
                   'C:/Windows/Fonts/arialbi.ttf', '/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf'],
}
SUBSUP = 0.75           # sub/superscript size relative to the line


class _Approx:
    """Fallback when no Arial-metric font is installed: ~0.55 em per character."""
    def getlength(self, s):
        return 550.0 * len(s)


def _font(bold, italic):
    key = (bool(bold), bool(italic))
    if key not in _FONTS:
        style = {(0, 0): 'Regular', (1, 0): 'Bold', (0, 1): 'Italic', (1, 1): 'BoldItalic'}[(int(key[0]), int(key[1]))]
        found = next((p for p in _FONT_FILES[style] if os.path.exists(p)), None)
        _FONTS[key] = ImageFont.truetype(found, 1000) if found else _Approx()
    return _FONTS[key]


def tw(s, size, bold=False, italic=False):
    return _font(bold, italic).getlength(s) * size / 1000.0


def _runs(runs):
    return [(r, {}) if isinstance(r, str) else r for r in runs]


def rich_w(runs, size, bold=False, italic=False):
    w = 0.0
    for t, o in _runs(runs):
        sz = size * SUBSUP if (o.get('sub') or o.get('sup')) else size
        w += tw(t, sz, o.get('b', bold), o.get('i', italic))
    return w


# ---------------------------------------------------------------- svg primitives
# Kept plain so that PowerPoint's "Convert to Shape" (and Illustrator / Inkscape) reproduce the figure: every shape
# states its fill and stroke, rounded boxes are explicit paths, arrows and small strokes are filled outlines, every
# text run is its own left-anchored text object (sub/superscripts included) and no element carries a transform.
out = []


def add(s):
    out.append(s)


def f(v):
    if isinstance(v, float):
        return f'{v:.1f}'.removesuffix('.0')
    return str(v)


def fs(size):
    return f'{size:.2f}'.rstrip('0').rstrip('.')


def _paint(fill, stroke, sw, cap=None, join=None):
    a = f' fill="{fill}"'
    if not stroke:
        return a + ' stroke="none"'
    a += f' stroke="{stroke}" stroke-width="{f(sw or 1)}"'
    if cap:
        a += f' stroke-linecap="{cap}"'
    if join:
        a += f' stroke-linejoin="{join}"'
    return a


def rrect_d(x, y, w, h, r):
    """Rounded rectangle as path data with cubic corners (an rx attribute gets reinterpreted on conversion)."""
    r = min(r, w / 2, h / 2)
    k = r * (1 - 0.5523)
    x1, y1 = x + w, y + h
    flat_w, flat_h = w - 2 * r > 0.05, h - 2 * r > 0.05
    d = [f'M{f(x + r)} {f(y)}']
    d += [f'L{f(x1 - r)} {f(y)}'] if flat_w else []
    d += [f'C{f(x1 - k)} {f(y)} {f(x1)} {f(y + k)} {f(x1)} {f(y + r)}']
    d += [f'L{f(x1)} {f(y1 - r)}'] if flat_h else []
    d += [f'C{f(x1)} {f(y1 - k)} {f(x1 - k)} {f(y1)} {f(x1 - r)} {f(y1)}']
    d += [f'L{f(x + r)} {f(y1)}'] if flat_w else []
    d += [f'C{f(x + k)} {f(y1)} {f(x)} {f(y1 - k)} {f(x)} {f(y1 - r)}']
    d += [f'L{f(x)} {f(y + r)}'] if flat_h else []
    d += [f'C{f(x)} {f(y + k)} {f(x + k)} {f(y)} {f(x + r)} {f(y)} Z']
    return ' '.join(d)


def rect(x, y, w, h, rx=0, fill='none', stroke=None, sw=None):
    if rx:
        add(f'<path d="{rrect_d(x, y, w, h, rx)}"{_paint(fill, stroke, sw)}/>')
    else:
        add(f'<rect x="{f(x)}" y="{f(y)}" width="{f(w)}" height="{f(h)}"{_paint(fill, stroke, sw)}/>')


def circle(cx, cy, r, fill='none', stroke=None, sw=None):
    add(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(r)}"{_paint(fill, stroke, sw)}/>')


def line(x1, y1, x2, y2, stroke, sw):
    add(f'<line x1="{f(x1)}" y1="{f(y1)}" x2="{f(x2)}" y2="{f(y2)}"{_paint("none", stroke, sw)}/>')


def path(d, fill='none', stroke=None, sw=None, cap='round', join='round'):
    add(f'<path d="{d}"{_paint(fill, stroke, sw, cap if stroke else None, join if stroke else None)}/>')


def _run(x, y, s, size, fill, bold, italic, cls=None):
    """One left-anchored text object; edge spaces become offsets, so no converter can drop or double them."""
    x += tw(s[:len(s) - len(s.lstrip(' '))], size, bold, italic)
    s = s.strip(' ')
    if not s:
        return
    a = f'<text x="{f(x)}" y="{f(y)}" font-family="{FONT}" font-size="{fs(size)}" fill="{fill}"'
    if bold:
        a += ' font-weight="bold"'
    if italic:
        a += ' font-style="italic"'
    if cls:
        a += f' class="{cls}"'
    add(a + f'>{escape(s)}</text>')


def _x0(x, w, anchor):
    return x - w / 2 if anchor == 'middle' else x - w if anchor == 'end' else x


def text(x, y, s, size, fill=INK, bold=False, italic=False, anchor='start'):
    _run(_x0(x, tw(s, size, bold, italic), anchor), y, s, size, fill, bold, italic)


def rich(x, y, runs, size, fill=INK, bold=False, italic=False, anchor='start'):
    """Inline runs: str or (str, {b, i, fill, sub, sup}). Each run is a separate text object at its measured x and
    sub/superscripts sit on their own shifted baseline, so they survive conversion to shapes (no dy, no tspan)."""
    x = _x0(x, rich_w(runs, size, bold, italic), anchor)
    for t, o in _runs(runs):
        sub, sup = o.get('sub'), o.get('sup')
        sz = size * SUBSUP if sub or sup else size
        b, i = o.get('b', bold), o.get('i', italic)
        yy = y + size * 0.24 if sub else y - size * 0.36 if sup else y
        _run(x, yy, t, sz, o.get('fill', fill), b, i, 'sub' if sub else 'sup' if sup else None)
        x += tw(t, sz, b, i)


def group(id_):
    add(f'<g id="{id_}">')


def end():
    add('</g>')


# ---------------------------------------------------------------- strokes drawn as filled outlines
def _unit(p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    n = math.hypot(dx, dy)
    return dx / n, dy / n


def _sides(pts, hw):
    """The two offset lines of a polyline at half-width hw, mitred at the corners."""
    left, right = [], []
    for i, (px, py) in enumerate(pts):
        ua = _unit(pts[i - 1], pts[i]) if i else _unit(pts[0], pts[1])
        ub = _unit(pts[i], pts[i + 1]) if i < len(pts) - 1 else ua
        k = hw / (1 + ua[0] * ub[0] + ua[1] * ub[1])
        ox, oy = (-ua[1] - ub[1]) * k, (ua[0] + ub[0]) * k
        left.append((px + ox, py + oy))
        right.append((px - ox, py - oy))
    return left, right


def _cap(p, u, hw, steps=6):
    """Round end at p facing u: the points between the left (+n) and the right (-n) side."""
    nx, ny = -u[1], u[0]
    angles = [math.pi * j / steps for j in range(1, steps)]
    return [(p[0] + hw * (nx * math.cos(a) + u[0] * math.sin(a)), p[1] + hw * (ny * math.cos(a) + u[1] * math.sin(a)))
            for a in angles]


def stroke_pts(pts, sw, cap='butt'):
    """Outline of an open polyline drawn sw wide (mitred corners, butt or round ends)."""
    left, right = _sides(pts, sw / 2)
    if cap != 'round':
        return left + right[::-1]
    return (left + _cap(pts[-1], _unit(pts[-2], pts[-1]), sw / 2) + right[::-1]
            + _cap(pts[0], _unit(pts[1], pts[0]), sw / 2))


def shape(polys, color):
    """One filled path made of one or more outlines."""
    d = ' '.join('M' + ' L'.join(f'{f(x)} {f(y)}' for x, y in p) + ' Z' for p in polys)
    add(f'<path d="{d}"{_paint(color, None, None)}/>')


def _head(tip, u, hw, kind, length, half):
    """Arrowhead at tip pointing along u, for a shaft of half-width hw: the length of shaft it replaces and its
    outline from the left barb over the tip to the right barb ('tri': flat base; 'stealth': swept back, notched)."""
    nx, ny = -u[1], u[0]
    bx, by = tip[0] - u[0] * length, tip[1] - u[1] * length
    barbs = [(bx + nx * half, by + ny * half), tuple(tip), (bx - nx * half, by - ny * half)]
    return (length * (0.7 + 0.3 * hw / half) if kind == 'stealth' else length), barbs


def arrow(points, color, sw=3.0, head=11.0, both=False, kind='tri', half=None):
    """Arrow along a polyline as ONE filled outline (mitred shaft plus head or heads): no markers, stroke caps or
    loose heads that drift when the SVG is converted to shapes."""
    pts = [tuple(p) for p in points]
    hw, half = sw / 2, half or head * 0.55
    u1 = _unit(pts[-2], pts[-1])
    d1, h1 = _head(pts[-1], u1, hw, kind, head, half)
    pts[-1] = (pts[-1][0] - u1[0] * d1, pts[-1][1] - u1[1] * d1)
    h0 = []
    if both:
        u0 = _unit(pts[1], pts[0])
        d0, h0 = _head(points[0], u0, hw, kind, head, half)
        pts[0] = (pts[0][0] - u0[0] * d0, pts[0][1] - u0[1] * d0)
    left, right = _sides(pts, hw)
    shape([left + h1 + right[::-1] + h0], color)


def bezier_pts(p0, c1, c2, p3, n=32):
    pts = []
    for j in range(n + 1):
        t = j / n
        a, b, c, d = (1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3
        pts.append((a * p0[0] + b * c1[0] + c * c2[0] + d * p3[0], a * p0[1] + b * c1[1] + c * c2[1] + d * p3[1]))
    return pts


def dash_pieces(pts, on, off):
    """Cut a dense polyline into dashes, as stroke-dasharray would."""
    cum = [0.0]
    for p, q in zip(pts, pts[1:]):
        cum.append(cum[-1] + math.dist(p, q))

    def at(s):
        i = min(bisect.bisect_right(cum, s) - 1, len(cum) - 2)
        t = (s - cum[i]) / (cum[i + 1] - cum[i])
        return pts[i][0] + t * (pts[i + 1][0] - pts[i][0]), pts[i][1] + t * (pts[i + 1][1] - pts[i][1])

    pieces, s = [], 0.0
    while s < cum[-1] - 1e-6:
        e = min(s + on, cum[-1])
        piece = [at(s)] + [p for p, c in zip(pts, cum) if s + 1e-6 < c < e - 1e-6] + [at(e)]
        pieces.append(piece)
        s += on + off
    return pieces


def rrect_pts(x, y, w, h, r, n=8):
    """Closed outline of a (rounded) rectangle as points, clockwise from the left end of the top edge."""
    if r <= 0:
        return [(x, y), (x + w, y), (x + w, y + h), (x, y + h), (x, y)]
    pts = [(x + r, y)]
    for cx, cy, a0 in [(x + w - r, y + r, -90), (x + w - r, y + h - r, 0), (x + r, y + h - r, 90), (x + r, y + r, 180)]:
        pts += [(cx + r * math.cos(math.radians(a0 + 90 * j / n)), cy + r * math.sin(math.radians(a0 + 90 * j / n)))
                for j in range(n + 1)]
    return pts


def dashed_rect(x, y, w, h, r, color, sw, on, off):
    """Dashed outline drawn as filled dashes (stroke-dasharray is lost by some shape converters)."""
    shape([stroke_pts(d, sw) for d in dash_pieces(rrect_pts(x, y, w, h, r), on, off)], color)


# ---------------------------------------------------------------- path data with the placement baked in
_TOKENS = re.compile(r'[A-Za-z]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')


def xform_d(d, s, tx, ty):
    """Path data rewritten as absolute M/L/C/Z with (x, y) -> (tx + s x, ty + s y) applied, so no transform is left
    for a converter to reinterpret (H/V become L, quadratic curves become cubic)."""
    tok = _TOKENS.findall(d)
    res, i, cmd = [], 0, ''
    cx = cy = mx = my = 0.0

    def nums(k):
        nonlocal i
        i += k
        return [float(t) for t in tok[i - k:i]]

    def pt(x, y):
        return f'{f(tx + s * x)} {f(ty + s * y)}'

    while i < len(tok):
        if tok[i].isalpha():
            cmd = tok[i]
            i += 1
            if cmd in 'Zz':
                res.append('Z')
                cx, cy = mx, my
                continue
        ox, oy = (cx, cy) if cmd.islower() else (0.0, 0.0)
        c = cmd.upper()
        if c in 'ML':
            x, y = nums(2)
            cx, cy = x + ox, y + oy
            res.append(c + pt(cx, cy))
            if c == 'M':
                mx, my = cx, cy
                cmd = 'l' if cmd == 'm' else 'L'          # further pairs after a moveto are linetos
        elif c == 'H':
            cx = nums(1)[0] + ox
            res.append('L' + pt(cx, cy))
        elif c == 'V':
            cy = nums(1)[0] + oy
            res.append('L' + pt(cx, cy))
        elif c == 'C':
            x1, y1, x2, y2, x, y = nums(6)
            res.append(f'C{pt(x1 + ox, y1 + oy)} {pt(x2 + ox, y2 + oy)} {pt(x + ox, y + oy)}')
            cx, cy = x + ox, y + oy
        elif c == 'Q':
            qx, qy, x, y = (v + o for v, o in zip(nums(4), (ox, oy, ox, oy)))
            res.append(f'C{pt(cx + (qx - cx) * 2 / 3, cy + (qy - cy) * 2 / 3)} '
                       f'{pt(x + (qx - x) * 2 / 3, y + (qy - y) * 2 / 3)} {pt(x, y)}')
            cx, cy = x, y
        else:
            raise ValueError(f'path command {cmd!r} not supported')
    return ' '.join(res)


def mix(c0, c1, t):
    a, b = (int(c0[k:k + 2], 16) for k in (1, 3, 5)), (int(c1[k:k + 2], 16) for k in (1, 3, 5))
    return '#' + ''.join(f'{round(p + (q - p) * t):02X}' for p, q in zip(a, b))


# ---------------------------------------------------------------- icons
def logo(name, x, y, size):
    """Fit a logo (paths taken from the source draft) into a size x size box at (x, y)."""
    items = LOGOS[name]
    x0 = min(i['bb'][0] for i in items)
    y0 = min(i['bb'][1] for i in items)
    x1 = max(i['bb'][2] for i in items)
    y1 = max(i['bb'][3] for i in items)
    s = size / max(x1 - x0, y1 - y0)
    tx = x + (size - (x1 - x0) * s) / 2 - x0 * s
    ty = y + (size - (y1 - y0) * s) / 2 - y0 * s
    for i in items:
        d = xform_d(i['d'], s, tx, ty)
        if i.get('stroke'):
            add(f'<path d="{d}"{_paint("none", i["stroke"], i["sw"] * s, cap="round", join="round")}/>')
        else:
            add(f'<path d="{d}" fill-rule="evenodd"{_paint(i.get("fill", "#000000"), None, None)}/>')


def wordmark(x, baseline, size, anchor='middle'):
    """AutoZT wordmark (outlined Sora) at x (centre or left edge); returns left, right and cap-top coordinates.
    The blue-to-crimson shade of "ZT" is one bounding-box gradient per letter, which converters map one to one."""
    wm = LOGOS['autozt_wordmark']
    s = size / 100
    bx0, by0, bx1, _ = wm['bbox']
    tx = (x - (bx1 - bx0) * s / 2 if anchor == 'middle' else x) - bx0 * s
    add('<g id="autozt-wordmark">')
    add(f'<path d="{xform_d(wm["auto"], s, tx, baseline)}"{_paint(INK, None, None)}/>')
    span = wm['zt_x1'] - wm['zt_x0']
    for k, glyph in enumerate(re.findall(r'M[^M]+', wm['zt'])):
        xs = [float(v) for v in _TOKENS.findall(xform_d(glyph, 1, 0, 0)) if not v.isalpha()][0::2]
        c0, c1 = (mix(BLUE, CRIMSON, (v - wm['zt_x0']) / span) for v in (min(xs), max(xs)))
        add(f'<defs><linearGradient id="wm-zt{k}" x1="0" y1="0" x2="1" y2="0"><stop offset="0" stop-color="{c0}"/>'
            f'<stop offset="1" stop-color="{c1}"/></linearGradient></defs>')
        add(f'<path d="{xform_d(glyph, s, tx, baseline)}"{_paint(f"url(#wm-zt{k})", None, None)}/>')
    add('</g>')
    return tx + bx0 * s, tx + bx1 * s, baseline + by0 * s


def robot_avatar(cx, cy, s=0.58):
    w, h = 60 * s, 46 * s
    rect(cx - w / 2 - 8 * s, cy - 8 * s, 8 * s + 1, 16 * s, rx=3 * s, fill=INK)
    rect(cx + w / 2 - 1, cy - 8 * s, 8 * s + 1, 16 * s, rx=3 * s, fill=INK)
    rect(cx - 3 * s, cy - h / 2 - 13 * s, 6 * s, 14 * s, fill=INK)
    circle(cx, cy - h / 2 - 16 * s, 6.5 * s, fill=INK)
    rect(cx - w / 2, cy - h / 2, w, h, rx=11 * s, fill=INK)
    rect(cx - 15 * s, cy - 9 * s, 11 * s, 18 * s, rx=5 * s, fill='#FFFFFF')
    rect(cx + 4 * s, cy - 9 * s, 11 * s, 18 * s, rx=5 * s, fill='#FFFFFF')


def person(cx, cy, s=0.6):
    circle(cx, cy - 13 * s, 15 * s, fill=INK)
    path(f'M{f(cx - 25 * s)} {f(cy + 26 * s)} C{f(cx - 25 * s)} {f(cy + 4 * s)} {f(cx + 25 * s)} {f(cy + 4 * s)} '
         f'{f(cx + 25 * s)} {f(cy + 26 * s)} Z', fill=INK)


def gate(cx, cy, r=13.0):
    circle(cx, cy, r, fill=GATE, stroke='#FFFFFF', sw=2.6)
    tick = [(cx - r * 0.45, cy + r * 0.02), (cx - r * 0.1, cy + r * 0.38), (cx + r * 0.5, cy - r * 0.36)]
    shape([stroke_pts(tick, r * 0.24, 'round')], '#FFFFFF')


def cursor(x, y, s=0.8):
    pts = [(0, 0), (0, 30), (7, 23.5), (12, 34), (17, 31.8), (12.2, 21.6), (21.5, 21.6)]
    path('M' + ' L'.join(f'{f(x + px * s)} {f(y + py * s)}' for px, py in pts) + ' Z', fill='#FFFFFF', stroke=INK, sw=1.8)


def chip(x, y, w, h, label, bg, fg, size=FS_XS, bold=True):
    rect(x, y, w, h, rx=5, fill=bg)
    text(x + w / 2, y + h / 2 + size * 0.35, label, size, fg, bold=bold, anchor='middle')


def pill_button(x, y, w, label, fill, fg, stroke=None):
    rect(x, y, w, 24, rx=12, fill=fill, stroke=stroke, sw=1.8 if stroke else None)
    text(x + w / 2, y + 12 + FS_XS * 0.36, label, FS_XS, fg, bold=True, anchor='middle')


ZT = [('zT', {'i': True})]
ZTNT = [('zT', {'i': True}), '(', ('n', {'i': True}), ', ', ('T', {'i': True}), ')']
# ================================================================ layout
M, CG = 12, 24                   # outer margin, gap between highlight cards
CW = (W - 2 * M - 2 * CG) / 3    # highlight card width
CX = [M + k * (CW + CG) for k in range(3)]
TOP_Y, TOP_H = 12, 236           # a-c
WEDGE = 70                       # room for the spotlight wedges
BAND_H, MID = 124, 48
PB = (TOP_Y + TOP_H + WEDGE, TOP_Y + TOP_H + WEDGE + BAND_H)   # phonon band
EB = (PB[1] + MID, PB[1] + MID + BAND_H)                        # electron band
BOT_Y, BOT_H = EB[1] + WEDGE, 272                               # d-f
H = BOT_Y + BOT_H + M
BX0, BX1 = 140, 1520             # chain bands
ZX = (1550, W - M)               # zT box
BOX_W, BOX_H = 200, 96
GAPX = (BX1 - BX0 - 28 - 6 * BOX_W) / 5
BXS = [BX0 + 14 + k * (BOX_W + GAPX) for k in range(6)]
Y_PH, Y_EL = (PB[0] + PB[1]) / 2, (EB[0] + EB[1]) / 2
AT_O, AT_M = '#E0A43A', '#3B6FC4'
GPU_BG, GPU_TX = '#DCEFE3', '#1E6B3A'
TAG_BG = '#EEF1F5'


def card(k, y, h, letter, title, col, tx, tag, gpu=False):
    x = CX[k]
    rect(x, y, CW, h, rx=14, fill='#FFFFFF', stroke=col, sw=2)
    text(x + 18, y + 40, letter, FS_LETTER, INK, bold=True)
    rich(x + 48, y + 40, title, FS_H, tx, bold=True)
    if gpu:
        gx = x + 48 + rich_w(title, FS_H, True) + 12
        rect(gx, y + 21, tw('GPU', FS_XS, True) + 16, 22, rx=11, fill=GPU_BG)
        text(gx + 8, y + 32 + FS_XS * 0.36, 'GPU', FS_XS, GPU_TX, bold=True)
    text(x + 48, y + 66, tag, FS_XS, MUTED)
    return x


def wedge(k, y_card, x0, x1, y_box, fill):
    """Spotlight from a step box to its highlight card."""
    cx0, cx1 = CX[k] + 0.2 * CW, CX[k] + 0.8 * CW
    shape([[(cx0, y_card), (cx1, y_card), (x1, y_box), (x0, y_box)]], fill)


def step(k, band, tag, title, sub=None, dashed=False, stack=False):
    x, y = BXS[k], band[0] + 14
    if stack:                    # fan-out step: one job per displaced or strained cell
        rect(x + 16, y - 12, BOX_W - 32, BOX_H, rx=12, fill='#FFFFFF', stroke='#A3ADB8', sw=1.6)
        rect(x + 8, y - 6, BOX_W - 16, BOX_H, rx=12, fill='#FFFFFF', stroke='#7F8A96', sw=1.6)
    rect(x, y, BOX_W, BOX_H, rx=12, fill='#FFFFFF', stroke=None if dashed else INK, sw=2.2)
    if dashed:                   # optional step
        dashed_rect(x, y, BOX_W, BOX_H, 12, INK, 2.2, 10, 7)
    chip(x + 10, y + 10, tw(tag, FS_XS) + 16, 22, tag, TAG_BG, SLATE, bold=False)
    lines = title if isinstance(title[0], list) else [title]
    if sub is None and len(lines) == 2:
        base = [y + 60, y + 84]
    else:
        base = [y + (60 if sub else 66)]
    for runs, by in zip(lines, base):
        rich(x + BOX_W / 2, by, runs, FS_B, INK, bold=True, anchor='middle')
    if sub:
        rich(x + BOX_W / 2, y + 85, sub, FS_S, MUTED, anchor='middle')


def output(k, band, runs, col, bg, tx):
    x, y = BXS[k], band[0] + 14
    rect(x, y, BOX_W, BOX_H, rx=12, fill=bg, stroke=col, sw=2.8)
    rich(x + BOX_W / 2, y + BOX_H / 2 + FS_TITLE * 0.36, runs, FS_TITLE, tx, bold=True, anchor='middle')


def link(k, yc, kind='arrow'):
    x0, x1 = BXS[k] + BOX_W, BXS[k + 1]
    if kind == 'plus':
        text((x0 + x1) / 2, yc + FS_H * 0.36, '+', FS_H, SLATE, bold=True, anchor='middle')
    else:
        arrow([(x0 + 4, yc), (x1 - 4, yc)], SLATE, sw=2.6, head=10)


def skill_tag(x, y, name, col, tx, label):
    w = tw(name, FS_S, True) + 24
    rect(x, y, w, 28, rx=8, fill='#FFFFFF', stroke=col, sw=2.2)
    text(x + w / 2, y + 14 + FS_S * 0.36, name, FS_S, tx, bold=True, anchor='middle')
    text(x + w + 12, y + 14 + FS_B * 0.36, label, FS_B, tx, bold=True)


def dd(k, y, chip_label, runs, x=None):
    """One 3D/2D rule row: dark 3D or green 2D chip, then the rule."""
    x = CX[k] + 18 if x is None else x
    for i, lab in enumerate(chip_label.split()):
        chip(x + i * 46, y - 18, 40, 24, lab, '#4A5563' if lab == '3D' else GREEN, '#FFFFFF')
    rich(x + 46 * len(chip_label.split()) + 6, y, runs, FS_S, INK)


def poscar(cx, cy):
    dw, dh, fold = 84, 104, 22
    dx0, dy0 = cx - dw / 2, cy - dh / 2
    text(cx, dy0 - 12, 'POSCAR', FS_B, INK, bold=True, anchor='middle')
    path(f'M{f(dx0)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw)} {f(dy0 + fold)} '
         f'L{f(dx0 + dw)} {f(dy0 + dh)} L{f(dx0)} {f(dy0 + dh)} Z', fill='#FFFFFF', stroke=INK, sw=2.6)
    path(f'M{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0 + fold)} L{f(dx0 + dw)} {f(dy0 + fold)}',
         stroke=INK, sw=2.2)
    mcx, mcy, mr = cx, cy + 8, 24
    hexp = [(mcx + mr * math.cos(math.radians(a)), mcy - mr * math.sin(math.radians(a)))
            for a in (90, 30, -30, -90, -150, 150)]
    for i in range(6):
        line(*hexp[i], *hexp[(i + 1) % 6], '#8AA0C0', 2.2)
    for i in (0, 2, 4):
        line(mcx, mcy, *hexp[i], '#8AA0C0', 2.2)
    for i, (px, py) in enumerate(hexp):
        circle(px, py, 7, fill='#D9534F' if i % 2 == 0 else '#3B6FC4')
    circle(mcx, mcy, 7, fill='#3B6FC4')


def mini_axes(x, y, w, h, xl, yl):
    arrow([(x, y), (x, y + h), (x + w, y + h)], GATE, sw=1.8, head=7, both=True)
    text(x - 8, y + 14, yl, FS_XS, MUTED, italic=True, anchor='end')
    text(x + w, y + h + 22, xl, FS_XS, MUTED, italic=True, anchor='end')


HM = [['EEB2C6', 'D16388', 'B23C66', '8A1D46', 'B23C66', 'E38AA8'],
      ['EEB2C6', 'E38AA8', 'D16388', 'B23C66', 'D16388', 'E38AA8'],
      ['F6D3DE', 'E38AA8', 'D16388', 'D16388', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'EEB2C6', 'E38AA8', 'E38AA8', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'F6D3DE', 'EEB2C6', 'EEB2C6', 'EEB2C6', 'F6D3DE'],
      ['F6D3DE', 'F6D3DE', 'F6D3DE', 'EEB2C6', 'F6D3DE', 'F6D3DE']]


def zt_map(x0, y0, cell=16.0, cg=1.2, size=FS_B):
    """Schematic zT(n, T) map with its top-left cell at (x0, y0): T up, n to the right (as in Fig. 1)."""
    gw = 6 * cell + 5 * cg
    for r_, row in enumerate(HM):
        for c_, colr in enumerate(row):
            rect(x0 + c_ * (cell + cg), y0 + r_ * (cell + cg), cell, cell, fill='#' + colr)
    k = cell / 16
    axx, axy = x0 - 7 * k, y0 + gw + 7 * k
    arrow([(axx, y0 - 4 * k), (axx, axy), (x0 + gw + 12 * k, axy)], GATE, sw=2.2 * k, head=9 * k, both=True)
    text(axx - 8 * k, y0 + gw / 2 + size * 0.36, 'T', size, MUTED, italic=True, anchor='end')
    text(x0 + gw + 22 * k, axy + size * 0.36, 'n', size, MUTED, italic=True)


def curve(pts, col, sw=2.2):
    path('M' + ' L'.join(f'{f(px)} {f(py)}' for px, py in pts), stroke=col, sw=sw)


# ================================================================ canvas
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{H / 10:.1f}mm" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add('<title>AutoZT zT workflow (Fig. 2)</title>')
rect(0, 0, W, H, fill='#FFFFFF')

# ================================================================ spotlights (drawn first, under everything)
group('spotlights')
for k, bk in enumerate((2, 3, 4)):
    wedge(k, TOP_Y + TOP_H, BXS[bk], BXS[bk] + BOX_W, PB[0], '#FBEEDF')
    wedge(k, BOT_Y, BXS[bk], BXS[bk] + BOX_W, EB[1], '#E6EEFA')
end()

# ================================================================ centre: POSCAR -> two chains -> zT
group('chains')
poscar(60, (PB[1] + EB[0]) / 2 + 10)
PX = 104
arrow([(PX, (PB[1] + EB[0]) / 2 - 6), (122, (PB[1] + EB[0]) / 2 - 6), (122, Y_PH), (BXS[0] - 4, Y_PH)],
      ORANGE, sw=2.8, head=11)
arrow([(PX, (PB[1] + EB[0]) / 2 + 26), (122, (PB[1] + EB[0]) / 2 + 26), (122, Y_EL), (BXS[0] - 4, Y_EL)],
      BLUE, sw=2.8, head=11)
for band, col, bg in ((PB, ORANGE, ORANGE_BG), (EB, BLUE, BLUE_BG)):
    rect(BX0, band[0], BX1 - BX0, BAND_H, rx=16, fill=bg, stroke=col, sw=2)
skill_tag(BX0, PB[0] - 40, 'kl-dft-cpu', ORANGE, ORANGE_TX, 'Phonons')
skill_tag(BX0, EB[1] + 12, 'ke-dft-cpu', BLUE, BLUE_TX, 'Electrons')
KL = [('κ', {'i': True}), ('L', {'sub': True})]
KE = [('κ', {'i': True}), ('e', {'sub': True})]
# phonons: one sequence of jobs
step(0, PB, 'S2_static', ['Static SCF'], ['reference cell'])
step(1, PB, 'S3_nac', ['NAC'], [('ε', {'i': True}), ('∞', {'sub': True}), ', Born charges'], dashed=True)
step(2, PB, 'S4_disp', ['Forces'], ['N displaced cells'], stack=True)
step(3, PB, 'S5_fc', ['Force constants'], ['fc2 + fc3 fit'])
step(4, PB, 'S6_kappa', ['Phonon BTE'], ['3-phonon scattering'])
output(5, PB, KL + ['(', ('T', {'i': True}), ')'], ORANGE, ORANGE_BG, ORANGE_TX)
for k in range(5):
    link(k, Y_PH)
# electrons: four inputs that AMSET combines
step(0, EB, 'S2.x', ['Bands'], ['PBE / HSE'])
step(1, EB, 'S3 + S4', ['Uniform k-grid'], ['wavefunctions'])
step(2, EB, 'S5 + S6', [['Dielectric ', ('ε', {'i': True})], ['Elastic ', ('C', {'i': True})]])
step(3, EB, 'S7_deform', [['Deformation ', ('E', {'i': True}), ('1', {'sub': True})]], ['N strained cells'], stack=True)
step(4, EB, 'S8 · S8.4', ['AMSET'], ['carrier transport'])
output(5, EB, [('σ', {'i': True}), ', ', ('S', {'i': True}), ', '] + KE, BLUE, BLUE_BG, BLUE_TX)
for k in range(5):
    link(k, Y_EL, 'plus' if k < 3 else 'arrow')
# zT
rect(ZX[0], PB[0], ZX[1] - ZX[0], EB[1] - PB[0], rx=14, fill=CRIMSON_BG, stroke=CRIMSON, sw=2.8)
for yc, col in ((Y_PH, ORANGE), (Y_EL, BLUE)):
    arrow([(BXS[5] + BOX_W + 4, yc), (ZX[0] - 4, yc)], col, sw=2.8, head=11)
zcx, eqy, FZ = (ZX[0] + ZX[1]) / 2, PB[0] + 78, 30
lhs = [('zT', {'i': True}), ' =']
num = [('S', {'i': True, 'fill': BLUE}), ('2', {'sup': True, 'fill': BLUE}), ('σ', {'i': True, 'fill': BLUE}),
       ('T', {'i': True})]
den = [('κ', {'i': True, 'fill': BLUE}), ('e', {'sub': True, 'fill': BLUE}), ' + ', ('κ', {'i': True, 'fill': ORANGE}),
       ('L', {'sub': True, 'fill': ORANGE})]
lw_ = rich_w(lhs, FZ)
fw = max(rich_w(num, FZ), rich_w(den, FZ)) + 12
ex0 = zcx - (lw_ + 10 + fw) / 2
rich(ex0, eqy + 10, lhs, FZ, INK)
fx = ex0 + lw_ + 10
rich(fx + fw / 2, eqy - 9, num, FZ, INK, anchor='middle')
line(fx, eqy, fx + fw, eqy, INK, 2.2)
rich(fx + fw / 2, eqy + 31, den, FZ, INK, anchor='middle')
zt_map(zcx - (6 * 17 + 5 * 1.2) / 2 + 12, eqy + 74, cell=17, cg=1.2)
end()

# ================================================================ a-c: phonon highlights
group('phonon-highlights')
# a random displacements: every atom moves, the frame count comes from ALM
x = card(0, TOP_Y, TOP_H, 'a', ['Random displacements'], ORANGE, ORANGE_TX, 'S4_disp · METHOD = alm')
JIT = [(3, -4), (-4, 2), (4, 3), (-2, -4), (-4, -3), (3, 4), (-3, 3), (4, -2),
       (2, 4), (-4, -1), (-3, -4), (4, 1), (-1, 4), (3, -3), (-4, 2), (2, -4)]
for c_, (kind, lab) in enumerate((('finite', 'finite'), ('random', 'random'))):
    gx, gy, gs = x + 22 + c_ * 116, TOP_Y + 88, 96
    rect(gx, gy, gs, gs, rx=6, fill='#FFFFFF', stroke='#A3ADB8', sw=1.6)
    for i in range(4):
        for j in range(4):
            px, py = gx + gs * (i + 0.5) / 4, gy + gs * (j + 0.5) / 4
            dx, dy = (JIT[i * 4 + j] if kind == 'random' else ((8, -5) if (i, j) == (1, 1) else (0, 0)))
            if kind == 'finite' and (i, j) == (1, 1):
                circle(px, py, 7.5, stroke=FAINT, sw=1.4)
            circle(px + dx, py + dy, 7.5, fill=AT_O if (i + j) % 2 else AT_M)
    text(gx + gs / 2, gy + gs + 26, lab, FS_XS, ORANGE_TX if kind == 'random' else MUTED,
         bold=kind == 'random', anchor='middle')
tx_ = x + 254
rich(tx_, TOP_Y + 118, [('Every atom', {'b': True}), ' moves in each frame'], FS_S, INK)
rich(tx_, TOP_Y + 152, [('N frames', {'b': True}), ' = free FCs × oversample'], FS_S, INK)
rich(tx_, TOP_Y + 186, ['free FCs counted by ALM'], FS_S, MUTED)

# b GPU force-constant fit with the imaginary-mode gate
x = card(1, TOP_Y, TOP_H, 'b', ['Force-constant fit'], ORANGE, ORANGE_TX, 'S5_fc · FIT_ENGINE = pheasy', gpu=True)
rich(x + 18, TOP_Y + 118, [('fc2 + fc3', {'b': True}), ' from random frames'], FS_S, INK)
gate(x + 28, TOP_Y + 145, r=10)
rich(x + 46, TOP_Y + 152, [('Gate:', {'b': True}), ' no imaginary modes'], FS_S, INK)
rich(x + 18, TOP_Y + 186, [('GPU = CPU', {'b': True}), ' to ~10', ('−16', {'sup': True})], FS_S, INK)
ax, ay, aw, ah = x + CW - 128, TOP_Y + 92, 96, 100
mini_axes(ax, ay, aw, ah, 'q', 'ω')
for amp, p in ((0.62, 1.0), (0.45, 1.0), (0.32, 2.0)):
    curve([(ax + aw * t, ay + ah - ah * amp * math.sin(math.pi * t / 2) ** p) for t in [i / 16 for i in range(17)]],
          ORANGE, 2)
for lvl in (0.18, 0.28):
    curve([(ax + aw * t, ay + ah * lvl + 6 * t * t) for t in [i / 16 for i in range(17)]], ORANGE, 2)

# c GPU phonon BTE
x = card(2, TOP_Y, TOP_H, 'c', ['Phonon BTE'], ORANGE, ORANGE_TX, 'S6_kappa · SOLVER = shengbte', gpu=True)
rich(x + 18, TOP_Y + 118, [('3-phonon BTE', {'b': True}), ' (RTA), CUDA build'], FS_S, INK)
rich(x + 18, TOP_Y + 152, KL + [' = CPU to 4–6 digits'], FS_S, INK, bold=False)
rich(x + 18, TOP_Y + 186, [('≈ 4.6× faster', {'b': True}), ' on RTX 3090'], FS_S, INK)
ax, ay, aw, ah = x + CW - 128, TOP_Y + 92, 96, 100
mini_axes(ax, ay, aw, ah, 'T', 'κ')
curve([(ax + 8 + (aw - 12) * t, ay + ah - ah * 0.86 / (1 + 4 * t)) for t in [i / 20 for i in range(21)]], ORANGE, 2.4)
end()

# ================================================================ d-f: electron highlights
group('electron-highlights')
# d scattering inputs, 3D and 2D
x = card(0, BOT_Y, BOT_H, 'd', ['Scattering inputs'], BLUE, BLUE_TX, 'S5_dielect · S6_elastic')
rich(x + 18, BOT_Y + 98, ['Dielectric ', ('ε', {'i': True})], FS_S, INK, bold=True)
dd(0, BOT_Y + 126, '3D', ['DFPT ', ('ε', {'i': True}), ('0', {'sub': True}), ', ', ('ε', {'i': True}),
                          ('∞', {'sub': True})])
dd(0, BOT_Y + 154, '2D', ['in-plane slab ', ('ε', {'i': True}), ', 2D screening'])
rich(x + 18, BOT_Y + 190, ['Elastic ', ('C', {'i': True})], FS_S, INK, bold=True)
dd(0, BOT_Y + 218, '3D', [('C', {'i': True}), ('ij', {'sub': True, 'i': True}), ' in GPa (Voigt)'])
dd(0, BOT_Y + 246, '2D', ['in-plane 2 × 2, ', ('C', {'i': True}), ('2D', {'sub': True}), ' = ', ('c', {'i': True}),
                          ' · ', ('C', {'i': True}), ('slab', {'sub': True})])

# e deformation potential: the reference level is the 3D/2D difference
x = card(1, BOT_Y, BOT_H, 'e', ['Deformation potential ', ('E', {'i': True}), ('1', {'sub': True})], BLUE, BLUE_TX,
         'S7_deform → S7.1_read')
dd(1, BOT_Y + 98, '3D', ['reference: core level'])
dd(1, BOT_Y + 126, '2D', ['reference: vacuum level (LOCPOT)'])
vx0, vx1, vy0, vy1 = x + 62, x + CW - 22, BOT_Y + 146, BOT_Y + 220
slab = (vx0 + 0.34 * (vx1 - vx0), vx0 + 0.66 * (vx1 - vx0))
rect(slab[0], vy0, slab[1] - slab[0], vy1 - vy0, fill=GREEN_BG)
text(x + 18, vy0 + 16, 'V(z)', FS_XS, MUTED, italic=True)
planes = [slab[0] + (slab[1] - slab[0]) * q for q in (0.2, 0.5, 0.8)]
vp = [(vx0 + (vx1 - vx0) * t, vy0 + 8 + (vy1 - vy0 - 20) * sum(math.exp(-((vx0 + (vx1 - vx0) * t - z) / 13) ** 2)
                                                           for z in planes)) for t in [i / 160 for i in range(161)]]
curve(vp, INK, 2.2)
core_y = max(py for _, py in vp)
shape([stroke_pts(d_, 1.6) for d_ in dash_pieces([(slab[0] - 10, core_y), (slab[1] + 10, core_y)], 6, 5)], MUTED)
text(slab[1] + 16, core_y + 6, 'core → 3D', FS_XS, MUTED, italic=True)
text(vx1, vy0 - 2, 'vacuum → 2D', FS_XS, GREEN, bold=True, italic=True, anchor='end')
text(x + 18, BOT_Y + 250, 'ions relaxed: all strains (3D), in-plane strains (2D)', FS_XS, MUTED, italic=True)

# f carrier transport: same five mechanisms, 3D or 2D kernels, plus an independent cross-check
x = card(2, BOT_Y, BOT_H, 'f', ['Carrier transport'], BLUE, BLUE_TX, 'S8_kappa · S8.4_amset2d')
mx = x + 18
for name in ('ADP', 'POP', 'IMP', 'PIE', 'ODP'):
    wd = tw(name, FS_XS, True) + 18
    chip(mx, BOT_Y + 84, wd, 26, name, '#E3ECF8', BLUE_TX)
    mx += wd + 8
text(mx + 4, BOT_Y + 97 + FS_XS * 0.36, 'scattering', FS_XS, MUTED, italic=True)
dd(2, BOT_Y + 146, '3D', [('AMSET', {'b': True}), ', 3D kernels'])
dd(2, BOT_Y + 176, '2D', [('AMSET', {'b': True}), ', 2D kernels'])
dd(2, BOT_Y + 214, '3D 2D', [('BoltzTraP2', {'b': True}), ' CRTA × DPT check'])
text(x + 18 + 2 * 46 + 6, BOT_Y + 244, 'S8.1_bt2 + S8.2_dpt, on by default', FS_XS, MUTED)
end()

add('</svg>')
dst = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fig2_zt_workflow_npj.svg')
with open(dst, 'w') as fh:
    fh.write('\n'.join(out) + '\n')
print(f'wrote {dst}: H {H}, box gap {GAPX:.1f}')
