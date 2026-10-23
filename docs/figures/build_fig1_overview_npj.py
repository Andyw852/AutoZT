#!/usr/bin/env python3
"""Fig. 1 of the AutoZT paper (npj Comput. Mater.): overview, 180 mm wide vector SVG.

Layout:
    a (LLM agent session, left column) | b (AutoZT: automated POSCAR -> zT(n, T) workflow,
                                       |    with the skill and engine library it selects from)
    c (one workflow for 3D and 2D materials, full width)
Units: 1 user unit = 0.1 mm, so the 1800-unit width prints at 180 mm (double column).
Type sizes: 1 pt = 3.528 units; nothing but sub/superscripts is set below 5 pt (17.6 units).
PowerPoint: insert the SVG and use Convert to Shape; the primitives below are chosen so that this keeps
sub/superscripts, arrows and rounded boxes as drawn (see 'svg primitives').
Usage: python3 build_fig1_overview_npj.py [out.svg]   (reads fig1_logos.json next to this script; needs Pillow)
Export: open the SVG in Inkscape / Illustrator, or print it to PDF from a browser, at 180 mm width.
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
A_X = (12, 516)                  # a: left column
R_X = (604, 1788)                # b: right column
TITLE_Y, TOP = 42, 60

# panel b, top to bottom: AutoZT mark + pipeline, automation loop, skill/engine library
LOGO_BASE = TOP + 62
Y0 = TOP + 190                   # pipeline centre line
TH = 90                          # transport box height
EL_Y, PH_Y = Y0 - 78, Y0 + 78
ZY = (EL_Y - TH / 2, PH_Y + TH / 2)
AS = (R_X[0] + 20, R_X[1] - 20, ZY[1] + 38, ZY[1] + 38 + 150)
BY0 = AS[3] + 46                 # library boxes
BY1 = BY0 + 176
ROW_BOT = BY1 + 14               # bottom of panels a and b

# panel c (3D and 2D) spans the width below
D_TITLE = ROW_BOT + 40
D_Y = (ROW_BOT + 58, ROW_BOT + 58 + 282)
H = D_Y[1] + 12


def panel(x0, x1, y0, y1):
    rect(x0, y0, x1 - x0, y1 - y0, rx=16, fill=PANEL_BG, stroke=PANEL_LINE, sw=2)


def panel_title(x, y, letter, runs):
    text(x, y, letter, FS_LETTER, INK, bold=True)
    rich(x + 34, y, runs, FS_TITLE, INK, bold=True)


def dag_icon(cx, cy, dx=30.0, dy=14.0, r=8.5):
    """The AutoZT mark: the zT workflow as a five-node DAG, coloured as in panel b."""
    nodes = {'s': (cx - 1.5 * dx, cy, GREEN), 'r': (cx - 0.5 * dx, cy, INK),
             'e': (cx + 0.5 * dx, cy - dy, BLUE), 'p': (cx + 0.5 * dx, cy + dy, ORANGE),
             'z': (cx + 1.5 * dx, cy, CRIMSON)}
    for a_, b_ in [('s', 'r'), ('r', 'e'), ('r', 'p'), ('e', 'z'), ('p', 'z')]:
        line(nodes[a_][0], nodes[a_][1], nodes[b_][0], nodes[b_][1], FAINT, max(2.0, r * 0.26))
    for nx, ny, ncol in nodes.values():
        circle(nx, ny, r, fill=ncol)


def hw_tags(cx, y, kinds):
    widths = [tw(k, FS_XS, True) + 14 for k in kinds]
    x = cx - (sum(widths) + 5 * (len(kinds) - 1)) / 2
    for k, wd in zip(kinds, widths):
        bg, fg = ('#E3E7EC', GATE) if k == 'CPU' else ('#DCEFE3', '#1E6B3A')
        rect(x, y, wd, 19, rx=9.5, fill=bg)
        text(x + wd / 2, y + 9.5 + FS_XS * 0.36, k, FS_XS, fg, bold=True, anchor='middle')
        x += wd + 5


HM = [['EEB2C6', 'D16388', 'B23C66', '8A1D46', 'B23C66', 'E38AA8'],
      ['EEB2C6', 'E38AA8', 'D16388', 'B23C66', 'D16388', 'E38AA8'],
      ['F6D3DE', 'E38AA8', 'D16388', 'D16388', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'EEB2C6', 'E38AA8', 'E38AA8', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'F6D3DE', 'EEB2C6', 'EEB2C6', 'EEB2C6', 'F6D3DE'],
      ['F6D3DE', 'F6D3DE', 'F6D3DE', 'EEB2C6', 'F6D3DE', 'F6D3DE']]


def zt_map(x0, y0, cell=16.0, cg=1.2, size=FS_B):
    """Schematic zT(n, T) map with its top-left cell at (x0, y0): T up, n to the right."""
    gw = 6 * cell + 5 * cg
    for r_, row in enumerate(HM):
        for c_, colr in enumerate(row):
            rect(x0 + c_ * (cell + cg), y0 + r_ * (cell + cg), cell, cell, fill='#' + colr)
    k = cell / 16
    axx, axy = x0 - 7 * k, y0 + gw + 7 * k
    arrow([(axx, y0 - 4 * k), (axx, axy), (x0 + gw + 12 * k, axy)], GATE, sw=2.2 * k, head=9 * k, both=True)
    text(axx - 8 * k, y0 + gw / 2 + size * 0.36, 'T', size, MUTED, italic=True, anchor='end')
    text(x0 + gw + 22 * k, axy + size * 0.36, 'n', size, MUTED, italic=True)


def mini_icon(kind, cx, cy):
    if kind == 'gate':
        gate(cx, cy, r=10)
    elif kind == 'lock':
        rect(cx - 7.5, cy - 2, 15, 11, rx=2.5, fill=SLATE)
        top, ro, ri = cy - 4.5, 5.2, 3.0
        outer = [(cx + ro * math.cos(math.pi * j / 12), top - ro * math.sin(math.pi * j / 12)) for j in range(13)]
        inner = [(cx - ri * math.cos(math.pi * j / 12), top - ri * math.sin(math.pi * j / 12)) for j in range(13)]
        shape([[(cx + ro, cy - 1)] + outer + [(cx - ro, cy - 1), (cx - ri, cy - 1)] + inner + [(cx + ri, cy - 1)]], SLATE)
        circle(cx, cy + 3.2, 1.8, fill='#FFFFFF')
    else:
        path(f'M{f(cx - 6.5)} {f(cy - 8.5)} L{f(cx + 2.5)} {f(cy - 8.5)} L{f(cx + 6.5)} {f(cy - 4.5)} '
             f'L{f(cx + 6.5)} {f(cy + 8.5)} L{f(cx - 6.5)} {f(cy + 8.5)} Z', fill='#FFFFFF', stroke=SLATE, sw=1.6)
        for dy_ in (-1, 3.5):
            line(cx - 3, cy + dy_, cx + 3.5, cy + dy_, SLATE, 1.5)


# ================================================================ canvas
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{H / 10:.1f}mm" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add('<title>AutoZT overview (Fig. 1)</title>')
rect(0, 0, W, H, fill='#FFFFFF')

# ================================================================ a: the LLM agent session
group('panel-a')
panel_title(A_X[0] + 2, TITLE_Y, 'a', ['LLM agent operates AutoZT'])
panel(A_X[0], A_X[1], TOP, ROW_BOT)
CH_X0, CH_X1 = A_X[0] + 14, A_X[1] - 14
AV_X = CH_X0 + 22
CARD_X0, CARD_X1 = CH_X0 + 42, CH_X1 - 8
USER_X1, USER_AV = CH_X1 - 42, CH_X1 - 20
IX0, IX1 = CARD_X0 + 9, CARD_X1 - 9
TB = 35                          # tool row block height
VIA = [('MCP', {'fill': PURPLE}), ('·', {'fill': FAINT}), ('CLI', {'fill': SLATE})]


def user_msg(y, runs):
    w = rich_w(runs, FS_B) + 32
    rect(USER_X1 - w, y, w, 32, rx=11, fill=USER_BG)
    rich(USER_X1 - w / 2, y + 16 + FS_B * 0.35, runs, FS_B, INK, anchor='middle')
    person(USER_AV, y + 15, s=0.56)


def tool_row(y, stages):
    """Stage chips plus a small MCP·CLI tag: each call is both an MCP tool and an agent-CLI command."""
    rect(IX0, y, IX1 - IX0, TB, rx=8, fill=LAV2)
    x = IX0 + 6
    for name in stages:
        wd = tw(name, FS_XS, True) + 14
        chip(x, y + 6, wd, 23, name, PEACH, PEACH_TX)
        x += wd + 6
    rich(IX1 - 8, y + TB / 2 + FS_XS * 0.36, VIA, FS_XS, PURPLE, bold=True, anchor='end')
    return y + TB


def agent_card(y, h):
    robot_avatar(AV_X, y + 22, s=0.54)
    rect(CARD_X0, y, CARD_X1 - CARD_X0, h, rx=12, fill='#FFFFFF', stroke=LAV_LINE, sw=1.6)


def say(y, runs):
    rich(IX0 + 4, y, runs if isinstance(runs, list) else [runs], FS_B, INK, bold=True)


# the chat opens with the agent itself (any coding agent), then request, plan, submit and report
H_MSG = [112, 32, 10 + TB + 76 + 24 + 12, 10 + TB + 36, 10 + TB + 36]
CH_Y0 = TOP + 12
# the "Submitted" card sits level with the automation loop in b, so the arrows between a and b run flat
G1 = ((AS[2] + AS[3]) / 2 - H_MSG[3] / 2 - (CH_Y0 + 12) - sum(H_MSG[:3])) / 3
ys = [CH_Y0 + 12]
for h_ in H_MSG[:-1]:
    ys.append(ys[-1] + h_ + G1)
CH_Y1 = ys[-1] + H_MSG[-1] + 12
rect(CH_X0, CH_Y0, CH_X1 - CH_X0, CH_Y1 - CH_Y0, rx=14, fill=LAV)
# one thread through the robots: every agent message comes from the agent named in the first card
line(AV_X, ys[0] + 22, AV_X, ys[-1] + 22, '#C9C0E4', 2.4)
# 0 the agent: any coding agent
y = ys[0]
agent_card(y, H_MSG[0])
text(IX0 + 4, y + 26, 'Coding agents', FS_S, PURPLE, bold=True)
AGENTS = [[('deepseek', 'DeepSeek'), ('kimi', 'Kimi'), ('glm', 'GLM')],
          [('claude', 'Claude'), ('chatgpt', 'ChatGPT'), (None, '…')]]
LG = 28                          # logo box; the name starts 8 units to its right
col_w = [max(LG + 8 + tw(lab, FS_S, True) if key else LG for key, lab in col) for col in zip(*AGENTS)]
gut = (IX1 - IX0 - 8 - sum(col_w)) / (len(col_w) - 1)
for r_, row in enumerate(AGENTS):
    gy = y + 38 + r_ * 36
    for c_, (key, lab) in enumerate(row):
        gx = IX0 + 4 + sum(col_w[:c_]) + c_ * gut
        if key:
            logo(key, gx, gy, LG)
            text(gx + LG + 8, gy + LG / 2 + FS_S * 0.36, lab, FS_S, INK, bold=True)
        else:
            text(gx + LG / 2, gy + LG / 2 + FS_S * 0.36, '…', FS_H, MUTED, bold=True, anchor='middle')
# 1 user asks in plain language
user_msg(ys[1], ['Compute '] + ZT + [' for this 2D POSCAR.'])
# 2 agent: discover, plan, dry run; the user approves the plan
y = ys[2]
agent_card(y, H_MSG[1])
yy = tool_row(y + 10, ['Discovery', 'Planning', 'Dry run'])
say(yy + 26, 'Plan ready. Proceed?')
x = IX0 + 4
for num, lab, col in [('1', 'structure', GREEN), ('2', 'electrons', BLUE), ('3', 'phonons', ORANGE), ('4', 'zT', CRIMSON)]:
    runs = ZT if lab == 'zT' else [lab]
    wd = rich_w(runs, FS_XS, True) + 30
    rect(x, yy + 41, wd, 21, rx=10.5, fill='#FFFFFF', stroke='#C9C0E4', sw=1.3)
    circle(x + 11, yy + 51.5, 8.5, fill=col)
    text(x + 11, yy + 51.5 + FS_XS * 0.36, num, FS_XS, '#FFFFFF', bold=True, anchor='middle')
    rich(x + 24, yy + 51.5 + FS_XS * 0.36, runs, FS_XS, INK, bold=True)
    x += wd + 6
pill_button(IX1 - 104, yy + 76, 100, 'Proceed', PURPLE, '#FFFFFF')
cursor(IX1 - 104 + 89, yy + 84, s=0.7)  # tip right of the label, so no letter is hidden
# 3 agent: submit; AutoZT now runs the workflow by itself
y = ys[3]
Y_MID = y + H_MSG[3] / 2
agent_card(y, H_MSG[3])
yy = tool_row(y + 10, ['Execution'])
say(yy + 25, 'Submitted; AutoZT runs the jobs.')
# 4 agent: report
y = ys[4]
agent_card(y, H_MSG[4])
yy = tool_row(y + 10, ['Reporting'])
say(yy + 25, ZTNT + [' is ready.'])
# what AutoZT guarantees whichever agent drives it
SG_Y0 = CH_Y1 + 12
rect(CH_X0, SG_Y0, CH_X1 - CH_X0, ROW_BOT - 14 - SG_Y0, rx=14, fill='#FFFFFF', stroke=LAV_LINE, sw=1.6)
text(CH_X0 + 18, SG_Y0 + 28, 'Safeguards', FS_S, PURPLE, bold=True)
for k, (kind, runs) in enumerate([('gate', [('Physical gates', {'b': True}), ' decide, not the model']),
                                  ('lock', [('Destructive steps', {'b': True}), ' need user approval']),
                                  ('doc', [('Every action', {'b': True}), ' is logged for audit'])]):
    cy = SG_Y0 + 56 + k * 28
    mini_icon(kind, CH_X0 + 30, cy)
    rich(CH_X0 + 52, cy + FS_S * 0.36, runs, FS_S, INK)
end()

# ================================================================ b: AutoZT, from POSCAR to zT(n, T)
group('panel-b')
panel_title(R_X[0] + 2, TITLE_Y, 'b', ['From POSCAR to '] + ZTNT)
panel(R_X[0], R_X[1], TOP, ROW_BOT)
WM_SIZE = 58
dag_icon(R_X[0] + 78, LOGO_BASE - 0.365 * WM_SIZE, dx=28, dy=13, r=8)
wordmark(R_X[0] + 136, LOGO_BASE, WM_SIZE, anchor='start')
# pipeline
DCX = R_X[0] + 76
dw, dh = 84, 104
dx0, dy0 = DCX - dw / 2, Y0 - dh / 2
text(DCX, dy0 - 12, 'POSCAR', FS_B, INK, bold=True, anchor='middle')
fold = 22
path(f'M{f(dx0)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw)} {f(dy0 + fold)} L{f(dx0 + dw)} {f(dy0 + dh)} '
     f'L{f(dx0)} {f(dy0 + dh)} Z', fill='#FFFFFF', stroke=INK, sw=2.6)
path(f'M{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0 + fold)} L{f(dx0 + dw)} {f(dy0 + fold)}', stroke=INK, sw=2.2)
mcx, mcy, mr = DCX, Y0 + 8, 24
hexp = [(mcx + mr * math.cos(math.radians(a)), mcy - mr * math.sin(math.radians(a))) for a in (90, 30, -30, -90, -150, 150)]
for i in range(6):
    line(*hexp[i], *hexp[(i + 1) % 6], '#8AA0C0', 2.2)
for i in (0, 2, 4):
    line(mcx, mcy, *hexp[i], '#8AA0C0', 2.2)
for i, (px, py) in enumerate(hexp):
    circle(px, py, 7, fill='#D9534F' if i % 2 == 0 else '#3B6FC4')
circle(mcx, mcy, 7, fill='#3B6FC4')


def step_box(x0, x1, yc, h, title, tcol, stroke, bg, sw=2.6, sub=None, gated=True):
    rect(x0, yc - h / 2, x1 - x0, h, rx=12, fill=bg, stroke=stroke, sw=sw)
    big = h >= 90
    ty = yc - 5 if sub else yc + 8
    text((x0 + x1) / 2, ty, title, FS_H if big else FS_B, tcol, bold=True, anchor='middle')
    if sub:
        rich((x0 + x1) / 2, ty + (30 if big else 26), sub, FS_B, tcol, bold=True, anchor='middle')
    if gated:
        gate(x1 - 4, yc - h / 2 + 2)


DM = (R_X[0] + 156, R_X[0] + 292)
RX = (R_X[0] + 326, R_X[0] + 462)
arrow([(dx0 + dw + 4, Y0), (DM[0] - 3, Y0)], INK, sw=2.8, head=11)
step_box(DM[0], DM[1], Y0, 80, 'Dimension', GREEN, GREEN, GREEN_BG, sub=['2D / 3D'], gated=False)
arrow([(DM[1] + 3, Y0), (RX[0] - 3, Y0)], INK, sw=2.8, head=11)
step_box(RX[0], RX[1], Y0, 80, 'Relaxation', INK, INK, '#FFFFFF')
TR = (R_X[0] + 532, R_X[0] + 788)
brx = RX[1] + 32
arrow([(RX[1] + 2, Y0 - 12), (brx, Y0 - 12), (brx, EL_Y), (TR[0] - 3, EL_Y)], BLUE, sw=2.8, head=11)
arrow([(RX[1] + 2, Y0 + 12), (brx, Y0 + 12), (brx, PH_Y), (TR[0] - 3, PH_Y)], ORANGE, sw=2.8, head=11)
step_box(TR[0], TR[1], EL_Y, TH, 'Electron transport', BLUE_TX, BLUE, BLUE_BG, sw=2.8,
         sub=[('σ', {'i': True}), ', ', ('S', {'i': True}), ', ', ('κ', {'i': True}), ('e', {'sub': True})])
step_box(TR[0], TR[1], PH_Y, TH, 'Phonon transport', ORANGE_TX, ORANGE, ORANGE_BG, sw=2.8,
         sub=[('κ', {'i': True}), ('L', {'sub': True})])
ZX = (TR[1] + 48, R_X[1] - 20)
for yc, col in [(EL_Y, BLUE), (PH_Y, ORANGE)]:
    arrow([(TR[1] + 3, yc), (ZX[0] - 3, yc)], col, sw=2.8, head=11)
rect(ZX[0], ZY[0], ZX[1] - ZX[0], ZY[1] - ZY[0], rx=14, fill=CRIMSON_BG, stroke=CRIMSON, sw=2.8)
zcx = (ZX[0] + ZX[1]) / 2
FZ = 28
eqy = ZY[0] + 50
lhs = [('zT', {'i': True}), ' =']
num = [('S', {'i': True, 'fill': BLUE}), ('2', {'sup': True, 'fill': BLUE}), ('σ', {'i': True, 'fill': BLUE}),
       ('T', {'i': True})]
den = [('κ', {'i': True, 'fill': BLUE}), ('e', {'sub': True, 'fill': BLUE}), ' + ', ('κ', {'i': True, 'fill': ORANGE}),
       ('L', {'sub': True, 'fill': ORANGE})]
lw_ = rich_w(lhs, FZ)
fw = max(rich_w(num, FZ), rich_w(den, FZ)) + 12
ex0 = zcx - (lw_ + 10 + fw) / 2
rich(ex0, eqy + 9, lhs, FZ, INK)
fx = ex0 + lw_ + 10
rich(fx + fw / 2, eqy - 9, num, FZ, INK, anchor='middle')
line(fx, eqy, fx + fw, eqy, INK, 2.2)
rich(fx + fw / 2, eqy + 29, den, FZ, INK, anchor='middle')
zt_map(zcx - (6 * 16 + 5 * 1.2) / 2 + 12, ZY[0] + 96)

# every step runs itself: the automation loop
rect(AS[0], AS[2], AS[1] - AS[0], AS[3] - AS[2], rx=14, fill='#F3F5FA', stroke='#C9D1DA', sw=1.6)
text(AS[0] + 18, AS[2] + 28, 'Every step runs automatically', FS_S, INK, bold=True)
CW_, CG_ = 202, 28
lx0 = (AS[0] + AS[1]) / 2 - (5 * CW_ + 4 * CG_) / 2
ly0, lh = AS[2] + 42, 50


def loop_icon(kind, cx, cy):
    circle(cx, cy, 16, fill='#EEF1F5')
    if kind == 'doc':
        path(f'M{f(cx - 7)} {f(cy - 9)} L{f(cx + 3)} {f(cy - 9)} L{f(cx + 8)} {f(cy - 4)} L{f(cx + 8)} {f(cy + 9)} '
             f'L{f(cx - 7)} {f(cy + 9)} Z', fill='#FFFFFF', stroke=SLATE, sw=1.8)
        for dy_ in (-1, 3.5):
            line(cx - 3.5, cy + dy_, cx + 4.5, cy + dy_, SLATE, 1.6)
    elif kind == 'rack':
        for dy_ in (-8, 1):
            rect(cx - 9, cy + dy_, 18, 7, rx=1.5, fill='#FFFFFF', stroke=SLATE, sw=1.6)
            circle(cx + 5, cy + dy_ + 3.5, 1.4, fill=GREEN)
    elif kind == 'pulse':
        pulse = [(cx - 10, cy), (cx - 4, cy), (cx - 1, cy - 8), (cx + 3, cy + 8), (cx + 6, cy), (cx + 10, cy)]
        shape([stroke_pts(pulse, 2, 'round')], SLATE)
    elif kind == 'gate':
        gate(cx, cy, r=11)
    else:
        for ox in (-6, 2):
            path(f'M{f(cx + ox - 4)} {f(cy - 7)} L{f(cx + ox + 4)} {f(cy)} L{f(cx + ox - 4)} {f(cy + 7)} Z', fill=SLATE)


for k, (t1, kind) in enumerate([('generate inputs', 'doc'), ('submit', 'rack'), ('monitor', 'pulse'),
                                ('check gate', 'gate'), ('advance', 'next')]):
    x = lx0 + k * (CW_ + CG_)
    rect(x, ly0, CW_, lh, rx=12, fill='#FFFFFF', stroke='#C9D1DA', sw=1.4)
    gx = x + (CW_ - 40 - tw(t1, FS_S, True)) / 2       # icon (32) + gap (8) + label, centred in the box
    loop_icon(kind, gx + 16, ly0 + lh / 2)
    text(gx + 40, ly0 + lh / 2 + FS_S * 0.36, t1, FS_S, INK, bold=True)
    if k < 4:
        arrow([(x + CW_ + 3, ly0 + lh / 2), (x + CW_ + CG_ - 3, ly0 + lh / 2)], SLATE, sw=2.4, head=9)
sx = lx0 + 1 * (CW_ + CG_) + CW_ / 2
cx_ = lx0 + 3 * (CW_ + CG_) + CW_ / 2
yb = ly0 + lh + 24
tip = (sx, ly0 + lh + 3)
reach, barbs = _head(tip, (0, -1), 1.2, 'stealth', 12, 5.5)
loop = (bezier_pts((cx_, ly0 + lh + 2), (cx_, yb), (cx_ - 30, yb), ((cx_ + sx) / 2, yb))
        + bezier_pts(((cx_ + sx) / 2, yb), (sx + 30, yb), (sx, yb), (sx, tip[1] + reach))[1:])
shape([stroke_pts(d, 2.4, 'round') for d in dash_pieces(loop, 7, 6)] + [barbs + [(sx, tip[1] + 12 * 0.7)]], RED)
text((cx_ + sx) / 2, yb + 22, 'retry · resume', FS_XS, RED, bold=True, italic=True, anchor='middle')

# library AutoZT picks from: skills (left) and engines (right)
SK = (R_X[0] + 16, R_X[0] + 16 + 566)
EN = (SK[1] + 16, R_X[1] - 16)
rect(SK[0], BY0, SK[1] - SK[0], BY1 - BY0, rx=14, fill=GREEN_BG, stroke=GREEN, sw=2.2)
text(SK[0] + 16, BY0 + 32, 'Workflow skills', FS_H, GREEN, bold=True)
cw = (SK[1] - SK[0] - 28 - 20) / 3
for k, (name, col, tcol) in enumerate([('ke-dft-cpu', BLUE, BLUE_TX), ('kl-dft-cpu', ORANGE, ORANGE_TX),
                                       ('zt-dft-cpu', CRIMSON, CRIMSON), ('band-dft-cpu', GREEN_SOFT, INK),
                                       ('phonon-dft-cpu', GREEN_SOFT, INK), ('fit-fc-thermal', GREEN_SOFT, INK)]):
    x = SK[0] + 14 + (k % 3) * (cw + 10)
    yk = BY0 + 46 + (k // 3) * 62
    feat = k < 3
    rect(x, yk, cw, 54, rx=9, fill='#FFFFFF', stroke=col, sw=2.6 if feat else 1.6)
    text(x + cw / 2, yk + 27 + FS_B * 0.36, name, FS_B, tcol, bold=feat, anchor='middle')
rect(EN[0], BY0, EN[1] - EN[0], BY1 - BY0, rx=14, fill='#EFF1F4', stroke=MUTED, sw=2.2)
text(EN[0] + 16, BY0 + 32, 'Computational engines', FS_H, INK, bold=True)
ce = (EN[1] - EN[0] - 28 - 20) / 3
for k, (name, hw) in enumerate([('VASP', 'CPU GPU'), ('AMSET', 'CPU'), ('ShengBTE', 'CPU GPU'),
                                ('Pheasy', 'CPU GPU'), ('phono3py', 'CPU'), (None, None)]):
    x = EN[0] + 14 + (k % 3) * (ce + 10)
    yk = BY0 + 46 + (k // 3) * 62
    if name is None:
        dashed_rect(x, yk, ce, 54, 9, '#A3ADB8', 1.6, 12, 8)
        text(x + ce / 2, yk + 34, '…', FS_H, MUTED, bold=True, anchor='middle')
        continue
    rect(x, yk, ce, 54, rx=9, fill='#FFFFFF', stroke='#A3ADB8', sw=1.6)
    text(x + ce / 2, yk + 22, name, FS_B, INK, bold=True, anchor='middle')
    hw_tags(x + ce / 2, yk + 29, hw.split())
SEL_X = SK[1] + 8                # between the skill box and the engine box
end()

# library -> workflow: AutoZT picks the skills and engines a material needs
group('b-select-link')
arrow([(SEL_X, BY0 - 3), (SEL_X, AS[3] + 3)], INK, sw=3.2, head=15, kind='stealth', half=6.5)
text(SEL_X + 16, (AS[3] + BY0) / 2 + 7, 'auto-selects', FS_S, INK, bold=True, italic=True)
end()

# a <-> b: tool calls go in, results come back; level with the "Submitted" card and the automation loop
group('a-b-interface')
GX0, GX1 = A_X[1] + 4, R_X[0] - 4
for yv, lab, col, rightward, ly in [(Y_MID - 14, 'tool calls', PURPLE, True, Y_MID - 26),
                                    (Y_MID + 14, 'results', CRIMSON, False, Y_MID + 40)]:
    x0, x1 = (GX0, GX1) if rightward else (GX1, GX0)
    arrow([(x0, yv), (x1, yv)], col, sw=3.2, head=15, kind='stealth', half=6.5)
    text((GX0 + GX1) / 2, ly, lab, FS_XS, col, bold=True, anchor='middle')
end()

# ================================================================ c: one workflow for 3D and 2D materials
group('panel-c')
panel_title(A_X[0] + 2, D_TITLE, 'c', ['One workflow for 3D and 2D materials'])
panel(A_X[0], R_X[1], D_Y[0], D_Y[1])
AT_O, AT_M = '#E0A43A', '#3B6FC4'


def layer(x0, y0, w):
    """X-M-X trilayer: 4 outer atoms per sheet, 3 inner atoms."""
    for k in range(4):
        circle(x0 + w * (k + 0.5) / 4, y0, 9.5, fill=AT_O)
        circle(x0 + w * (k + 0.5) / 4, y0 + 48, 9.5, fill=AT_O)
    for k in range(3):
        circle(x0 + w * (k + 1) / 4, y0 + 24, 11.5, fill=AT_M)


def c_axis(x, y0, y1):
    arrow([(x, y0 + 2), (x, y1 - 2)], GATE, sw=2.4, head=10, both=True)
    text(x + 12, (y0 + y1) / 2 + 8, 'c', FS_B, INK, italic=True)


SKW = 108
bkx, bkh = A_X[0] + 46, 150
slx, sly, slh = bkx + 214, D_Y[0] + 26, 212
bky = sly + slh - bkh            # bottom-aligned with the slab
rect(bkx, bky, SKW, bkh, rx=4, fill='#FFFFFF', stroke=FAINT, sw=2.2)
layer(bkx + 4, bky + 20, SKW - 8)
layer(bkx + 4, bky + 82, SKW - 8)
c_axis(bkx + SKW + 22, bky, bky + bkh)
text(bkx + SKW / 2, bky + bkh + 30, '3D bulk', FS_B, INK, bold=True, anchor='middle')
slab_y0, slab_h = sly + 69, 74
rect(slx, sly, SKW, slh, fill='#FFFFFF')
rect(slx, slab_y0, SKW, slab_h, fill='#E3F2EA')
dashed_rect(slx, sly, SKW, slh, 0, FAINT, 2.2, 12, 9)
text(slx + SKW / 2, sly + 42, 'vacuum', FS_XS, FAINT, italic=True, anchor='middle')
text(slx + SKW / 2, sly + slh - 28, 'vacuum', FS_XS, FAINT, italic=True, anchor='middle')
layer(slx + 4, slab_y0 + 13, SKW - 8)
path(f'M{f(slx - 8)} {f(slab_y0)} L{f(slx - 18)} {f(slab_y0)} L{f(slx - 18)} {f(slab_y0 + slab_h)} '
     f'L{f(slx - 8)} {f(slab_y0 + slab_h)}', stroke=GREEN, sw=2.6, cap='butt', join='miter')
text(slx - 26, slab_y0 + slab_h / 2 + 9, 't', FS_H, GREEN, bold=True, italic=True, anchor='end')
c_axis(slx + SKW + 22, sly, sly + slh)
text(slx + SKW / 2, sly + slh + 30, '2D slab', FS_B, GREEN, bold=True, anchor='middle')

RX0 = A_X[0] + 450
ry = D_Y[0] + 50
chip(RX0, ry - 19, 40, 26, '3D', '#4A5563', '#FFFFFF')
rich(RX0 + 52, ry, [('Bulk', {'b': True}), ': standard path'], FS_B, INK)
x2d = RX0 + 52 + rich_w([('Bulk', {'b': True}), ': standard path'], FS_B) + 70
chip(x2d, ry - 19, 40, 26, '2D', GREEN, '#FFFFFF')
rich(x2d + 52, ry, [('Slab', {'b': True}), ': same path, plus'], FS_B, GREEN)
rich(R_X[1] - 24, ry, [('Dimension', {'b': True, 'fill': GREEN}), ' auto-detected: ≥ 8 Å vacuum along one axis'],
     FS_S, MUTED, anchor='end')
EPS2D = [('ε', {'i': True}), ('2D', {'sub': True})]
GROUPS = [
    (['Structure'], INK, [[('Vacuum', {'b': True}), ' along ', ('c', {'i': True})],
                          [('Relaxation', {'b': True}), ' in-plane'],
                          [('k-mesh', {'b': True}), ' 1 along ', ('c', {'i': True})],
                          [('Thickness', {'b': True}), ' ', ('t', {'i': True}), ' from vdW radii']]),
    (['Electrons'], BLUE_TX, [[('Elastic constant', {'b': True}), ' in N m', ('−1', {'sup': True})],
                              [('Dielectric', {'b': True}), ' without vacuum'],
                              [('Deformation ', {'b': True}), ('E', {'b': True, 'i': True}), ('1', {'sub': True, 'b': True}),
                               ' vs vacuum'],
                              [('Scattering', {'b': True}), ' with 2D kernels']]),
    (['Phonons'], ORANGE_TX, [[('NAC', {'b': True}), ' off'], [('ZA branch', {'b': True}), ' quadratic'],
                              [('Supercell', {'b': True}), ' 1 along ', ('c', {'i': True})],
                              [('Rotational sum rules', {'b': True}), ' on']]),
    (ZT + [' inputs'], CRIMSON, [[('σ', {'i': True, 'b': True}), (', ', {'b': True}), ('κ', {'i': True, 'b': True}),
                                  ('e', {'sub': True, 'b': True}), (', ', {'b': True}), ('κ', {'i': True, 'b': True}),
                                  ('L', {'sub': True, 'b': True}), ' × ', ('c', {'i': True}), ' / ', ('t', {'i': True})]]),
]
GCW = (R_X[1] - 24 - RX0) / 4
gy0 = ry + 54
for gi, (gname, col, items) in enumerate(GROUPS):
    gx = RX0 + gi * GCW
    rich(gx, gy0, gname, FS_S, col, bold=True)
    line(gx + rich_w(gname, FS_S, True) + 10, gy0 - 6, gx + GCW - 26, gy0 - 6, '#C9D1DA', 1.6)
    for ii, runs in enumerate(items):
        rich(gx, gy0 + 38 + ii * 34, runs, FS_B, col)
fgx = RX0 + 3 * GCW
fy0 = gy0 + 38 + 22
rect(fgx, fy0, GCW - 26, 54, rx=10, fill='#F3FAF6', stroke=GREEN, sw=2)
rich(fgx + (GCW - 26) / 2, fy0 + 35, EPS2D + [' = 1 + (', ('c', {'i': True}), ' / ', ('t', {'i': True}), ') (',
                                               ('ε', {'i': True}), ('slab', {'sub': True}), ' − 1)'],
     FS_B, INK, anchor='middle')
rich(fgx, gy0 + 38 + 3 * 34, [('zT', {'i': True, 'b': True}), ' independent of ', ('t', {'i': True})], FS_B, CRIMSON)
end()

add('</svg>')
dst = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fig1_autozt_overview_npj.svg')
with open(dst, 'w') as fh:
    fh.write('\n'.join(out) + '\n')
print(f'wrote {dst}: H {H}, chat gap {G1:.1f}, safeguards card height {ROW_BOT - 14 - SG_Y0:.0f}, '
      f'Submitted card centre {Y_MID:.0f} vs automation loop centre {(AS[2] + AS[3]) / 2:.0f}')
