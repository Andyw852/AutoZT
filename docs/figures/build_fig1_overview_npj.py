#!/usr/bin/env python3
"""Fig. 1 of the AutoZT paper (npj Comput. Mater.): overview, 180 mm wide vector SVG.

Layout:
    a (LLM agent session, left column) | b (AutoZT: automated POSCAR -> zT(n, T) workflow,
                                       |    with the skill and engine library it selects from)
    c (one workflow for 3D and 2D materials, full width)
Units: 1 user unit = 0.1 mm, so the 1800-unit width prints at 180 mm (double column).
Type sizes: 1 pt = 3.528 units; nothing is set below 5 pt (17.6 units).
Usage: python3 build_fig1_overview_npj.py [out.svg]   (reads fig1_logos.json next to this script; needs Pillow)
Export: open the SVG in Inkscape / Illustrator, or print it to PDF from a browser, at 180 mm width.
"""
import json
import math
import os
import sys
from xml.sax.saxutils import escape

from PIL import ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'fig1_logos.json')) as _fh:
    LOGOS = json.load(_fh)

# ---------------------------------------------------------------- tokens
FONT = "Arial, Helvetica, 'Liberation Sans', sans-serif"
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
PINK = '#F8D7E1'
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
out = []


def add(s):
    out.append(s)


def f(v):
    if isinstance(v, float):
        return f'{v:.1f}'.removesuffix('.0')
    return str(v)


def rect(x, y, w, h, rx=0, fill='none', stroke=None, sw=None, dash=None):
    a = f'<rect x="{f(x)}" y="{f(y)}" width="{f(w)}" height="{f(h)}"'
    if rx:
        a += f' rx="{f(rx)}"'
    a += f' fill="{fill}"'
    if stroke:
        a += f' stroke="{stroke}" stroke-width="{f(sw or 1)}"'
    if dash:
        a += f' stroke-dasharray="{dash}"'
    add(a + '/>')


def circle(cx, cy, r, fill='none', stroke=None, sw=None):
    a = f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(r)}" fill="{fill}"'
    if stroke:
        a += f' stroke="{stroke}" stroke-width="{f(sw or 1)}"'
    add(a + '/>')


def line(x1, y1, x2, y2, stroke, sw):
    add(f'<line x1="{f(x1)}" y1="{f(y1)}" x2="{f(x2)}" y2="{f(y2)}" stroke="{stroke}" stroke-width="{f(sw)}"/>')


def path(d, fill='none', stroke=None, sw=None, cap='round', join='round', dash=None):
    a = f'<path d="{d}" fill="{fill}"'
    if stroke:
        a += f' stroke="{stroke}" stroke-width="{f(sw or 1)}" stroke-linecap="{cap}" stroke-linejoin="{join}"'
    if dash:
        a += f' stroke-dasharray="{dash}"'
    add(a + '/>')


def text(x, y, s, size, fill=INK, bold=False, italic=False, anchor='start'):
    a = f'<text x="{f(x)}" y="{f(y)}" font-size="{f(size)}" fill="{fill}"'
    if bold:
        a += ' font-weight="bold"'
    if italic:
        a += ' font-style="italic"'
    if anchor != 'start':
        a += f' text-anchor="{anchor}"'
    add(a + f'>{escape(s)}</text>')


def rich(x, y, runs, size, fill=INK, bold=False, italic=False, anchor='start'):
    """Inline runs: str or (str, {b, i, fill, sub, sup}); sub/sup via dy so editors keep one text flow."""
    parts, cur = [], 0.0
    for t, o in _runs(runs):
        target = size * 0.24 if o.get('sub') else (-size * 0.36 if o.get('sup') else 0.0)
        a = []
        if abs(target - cur) > 1e-6:
            a.append(f'dy="{target - cur:.1f}"')
            cur = target
        if o.get('sub') or o.get('sup'):
            a.append(f'font-size="{size * SUBSUP:.1f}"')
        if 'b' in o and o['b'] != bold:
            a.append(f'font-weight="{"bold" if o["b"] else "normal"}"')
        if 'i' in o and o['i'] != italic:
            a.append(f'font-style="{"italic" if o["i"] else "normal"}"')
        if 'fill' in o and o['fill'] != fill:
            a.append(f'fill="{o["fill"]}"')
        parts.append(f'<tspan {" ".join(a)}>{escape(t)}</tspan>' if a else escape(t))
    a = f'<text x="{f(x)}" y="{f(y)}" font-size="{f(size)}" fill="{fill}"'
    if bold:
        a += ' font-weight="bold"'
    if italic:
        a += ' font-style="italic"'
    if anchor != 'start':
        a += f' text-anchor="{anchor}"'
    add(a + '>' + ''.join(parts) + '</text>')


def group(id_):
    add(f'<g id="{id_}">')


def end():
    add('</g>')


def arrow(points, color, sw=3.0, head=11.0):
    """Polyline with a filled head at the last point (no SVG markers: they import badly into Illustrator/PowerPoint)."""
    (x1, y1), (x2, y2) = points[-2], points[-1]
    dx, dy = x2 - x1, y2 - y1
    n = math.hypot(dx, dy)
    ux, uy = dx / n, dy / n
    bx, by = x2 - ux * head, y2 - uy * head
    pts = points[:-1] + [(bx, by)]
    path('M' + ' L'.join(f'{f(px)} {f(py)}' for px, py in pts), stroke=color, sw=sw, cap='butt')
    hw = head * 0.55
    path(f'M{f(x2)} {f(y2)} L{f(bx - uy * hw)} {f(by + ux * hw)} L{f(bx + uy * hw)} {f(by - ux * hw)} Z', fill=color)


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
    add(f'<g transform="translate({tx:.2f} {ty:.2f}) scale({s:.5f})">')
    for i in items:
        if i.get('stroke'):
            add(f'<path d="{i["d"]}" fill="none" stroke="{i["stroke"]}" stroke-width="{i["sw"]}" '
                f'stroke-linecap="round" stroke-linejoin="round"/>')
        else:
            add(f'<path d="{i["d"]}" fill="{i.get("fill", "#000000")}" fill-rule="evenodd"/>')
    add('</g>')


def stealth(tip, ux, uy, color, length=15.0, half=6.5):
    """Arrowhead with swept-back barbs and a notch, pointing along (ux, uy)."""
    tx, ty = tip
    bx, by = tx - ux * length, ty - uy * length
    nx, ny = -uy, ux
    kx, ky = tx - ux * length * 0.7, ty - uy * length * 0.7
    path(f'M{f(tx)} {f(ty)} L{f(bx + nx * half)} {f(by + ny * half)} L{f(kx)} {f(ky)} '
         f'L{f(bx - nx * half)} {f(by - ny * half)} Z', fill=color)


def curve_arrow(p0, c1, c2, p3, color, sw=3.4):
    """Cubic Bezier ending in a stealth head at p3."""
    dx, dy = p3[0] - c2[0], p3[1] - c2[1]
    n = math.hypot(dx, dy)
    ux, uy = dx / n, dy / n
    end_ = (p3[0] - ux * 10, p3[1] - uy * 10)
    path(f'M{f(p0[0])} {f(p0[1])} C{f(c1[0])} {f(c1[1])} {f(c2[0])} {f(c2[1])} {f(end_[0])} {f(end_[1])}',
         stroke=color, sw=sw)
    stealth(p3, ux, uy, color)


def flow_arrow(x, y0, y1, gid, w0=20.0, w1=11.0, head_w=44.0, head_l=32.0):
    """Tapered downward arrow, shaded like the ZT of the wordmark (blue to crimson)."""
    add(f'<defs><linearGradient id="{gid}" gradientUnits="userSpaceOnUse" x1="0" y1="{f(y0)}" x2="0" y2="{f(y1)}">'
        f'<stop offset="0" stop-color="{BLUE}" stop-opacity="0.35"/><stop offset="1" stop-color="{CRIMSON}"/>'
        '</linearGradient></defs>')
    hy = y1 - head_l
    path(f'M{f(x - w0 / 2)} {f(y0)} L{f(x + w0 / 2)} {f(y0)} L{f(x + w1 / 2)} {f(hy)} L{f(x + head_w / 2)} {f(hy)} '
         f'L{f(x)} {f(y1)} L{f(x - head_w / 2)} {f(hy)} L{f(x - w1 / 2)} {f(hy)} Z', fill=f'url(#{gid})')


def wordmark(x, baseline, size, anchor='middle'):
    """AutoZT wordmark (outlined Sora) at x (centre or left edge); returns left, right and cap-top coordinates."""
    wm = LOGOS['autozt_wordmark']
    s = size / 100
    bx0, by0, bx1, _ = wm['bbox']
    tx = (x - (bx1 - bx0) * s / 2 if anchor == 'middle' else x) - bx0 * s
    add(f'<defs><linearGradient id="wm-zt" gradientUnits="userSpaceOnUse" x1="{wm["zt_x0"]}" y1="0" '
        f'x2="{wm["zt_x1"]}" y2="0"><stop offset="0" stop-color="{BLUE}"/><stop offset="1" stop-color="{CRIMSON}"/>'
        '</linearGradient></defs>')
    add(f'<g id="autozt-wordmark" transform="translate({tx:.2f} {baseline:.2f}) scale({s:.4f})">')
    add(f'<path d="{wm["auto"]}" fill="{INK}"/>')
    add(f'<path d="{wm["zt"]}" fill="url(#wm-zt)"/>')
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
    path(f'M{f(cx - r * 0.45)} {f(cy + r * 0.02)} L{f(cx - r * 0.1)} {f(cy + r * 0.38)} L{f(cx + r * 0.5)} {f(cy - r * 0.36)}',
         stroke='#FFFFFF', sw=r * 0.24)


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
AS = (R_X[0] + 20, R_X[1] - 20, ZY[1] + 26, ZY[1] + 26 + 150)
BY0 = AS[3] + 62                 # library boxes
BY1 = BY0 + 188
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


# ================================================================ canvas
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{H / 10:.1f}mm" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add('<title>AutoZT overview (Fig. 1)</title>')
rect(0, 0, W, H, fill='#FFFFFF')

# ================================================================ a: the LLM agent session
group('panel-a')
panel_title(A_X[0] + 2, TITLE_Y, 'a', ['LLM agent operates AutoZT'])
panel(A_X[0], A_X[1], TOP, ROW_BOT)
text(A_X[0] + 22, TOP + 38, 'Coding agents', FS_S, PURPLE, bold=True)
cwl = (A_X[1] - A_X[0] - 40) / 3
for k, (key, lab) in enumerate([('deepseek', 'DeepSeek'), ('kimi', 'Kimi'), ('glm', 'GLM'),
                                ('claude', 'Claude'), ('chatgpt', 'ChatGPT'), (None, '…')]):
    gx = A_X[0] + 20 + (k % 3) * cwl
    gy = TOP + 58 + (k // 3) * 52
    if key:
        lw = tw(lab, FS_S, True)
        x0 = gx + (cwl - (34 + 8 + lw)) / 2
        logo(key, x0, gy, 34)
        text(x0 + 42, gy + 17 + FS_S * 0.36, lab, FS_S, INK, bold=True)
    else:
        text(gx + cwl / 2, gy + 23, '…', FS_H, MUTED, bold=True, anchor='middle')

CH_X0, CH_X1 = A_X[0] + 14, A_X[1] - 14
CH_Y0, CH_Y1 = TOP + 176, ROW_BOT - 14
rect(CH_X0, CH_Y0, CH_X1 - CH_X0, CH_Y1 - CH_Y0, rx=14, fill=LAV)
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
        bg, fg = (PINK, RED) if name == 'Approval' else (PEACH, PEACH_TX)
        wd = tw(name, FS_XS, True) + 14
        chip(x, y + 6, wd, 23, name, bg, fg)
        x += wd + 6
    rich(IX1 - 8, y + TB / 2 + FS_XS * 0.36, VIA, FS_XS, PURPLE, bold=True, anchor='end')
    return y + TB


def agent_card(y, h):
    robot_avatar(AV_X, y + 22, s=0.54)
    rect(CARD_X0, y, CARD_X1 - CARD_X0, h, rx=12, fill='#FFFFFF', stroke=LAV_LINE, sw=1.6)


def say(y, runs):
    rich(IX0 + 4, y, runs if isinstance(runs, list) else [runs], FS_B, INK, bold=True)


H_MSG = [32, 10 + TB + 34 + 30 + 10, 10 + TB + 36, 32, 10 + TB + 72, 10 + TB + 36]
# uniform gaps; the arrows to panel b are level with the "Submitted" card (message 3)
GAP = (CH_Y1 - CH_Y0 - 24 - sum(H_MSG)) / (len(H_MSG) - 1)
ys = [CH_Y0 + 12]
for h_ in H_MSG[:-1]:
    ys.append(ys[-1] + h_ + GAP)
# 1 user asks in plain language
user_msg(ys[0], ['Compute '] + ZT + [' for this 2D POSCAR.'])
# 2 agent: discover, plan, dry run; the user approves the plan
y = ys[1]
agent_card(y, H_MSG[1])
yy = tool_row(y + 10, ['Discovery', 'Planning', 'Dry run'])
say(yy + 25, 'Plan ready. Proceed?')
pbx = IX1 - 8 - 100
pill_button(pbx, yy + 7, 100, 'Proceed', PURPLE, '#FFFFFF')
cursor(pbx + 89, yy + 15, s=0.7)  # tip right of the label, so no letter is hidden
x = IX0 + 4
for num, lab, col in [('1', 'structure', GREEN), ('2', 'electrons', BLUE), ('3', 'phonons', ORANGE), ('4', 'zT', CRIMSON)]:
    runs = ZT if lab == 'zT' else [lab]
    wd = rich_w(runs, FS_XS, True) + 30
    rect(x, yy + 40, wd, 21, rx=10.5, fill='#FFFFFF', stroke='#C9C0E4', sw=1.3)
    circle(x + 11, yy + 50.5, 8.5, fill=col)
    text(x + 11, yy + 50.5 + FS_XS * 0.36, num, FS_XS, '#FFFFFF', bold=True, anchor='middle')
    rich(x + 24, yy + 50.5 + FS_XS * 0.36, runs, FS_XS, INK, bold=True)
    x += wd + 6
# 3 agent: submit; AutoZT now runs the workflow by itself
y = ys[2]
Y_MID = y + H_MSG[2] / 2
agent_card(y, H_MSG[2])
yy = tool_row(y + 10, ['Execution'])
say(yy + 25, 'Submitted; AutoZT runs the jobs.')
# 4 user asks for a destructive action
user_msg(ys[3], ['Rerun the deformation step.'])
# 5 agent: approval in chat
y = ys[4]
agent_card(y, H_MSG[4])
yy = tool_row(y + 10, ['Approval'])
say(yy + 25, 'Deletes its folder. Approve?')
pill_button(IX0 + 4, yy + 36, 88, 'Approve', RED, '#FFFFFF')
pill_button(IX0 + 100, yy + 36, 88, 'Deny', '#FFFFFF', RED, stroke=RED)
# 6 agent: report
y = ys[5]
agent_card(y, H_MSG[5])
yy = tool_row(y + 10, ['Reporting'])
say(yy + 25, ZTNT + [' is ready.'])
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
HM = [['EEB2C6', 'D16388', 'B23C66', '8A1D46', 'B23C66', 'E38AA8'],
      ['EEB2C6', 'E38AA8', 'D16388', 'B23C66', 'D16388', 'E38AA8'],
      ['F6D3DE', 'E38AA8', 'D16388', 'D16388', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'EEB2C6', 'E38AA8', 'E38AA8', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'F6D3DE', 'EEB2C6', 'EEB2C6', 'EEB2C6', 'F6D3DE'],
      ['F6D3DE', 'F6D3DE', 'F6D3DE', 'EEB2C6', 'F6D3DE', 'F6D3DE']]
cell, cg = 16, 1.2
gw = 6 * cell + 5 * cg
hmx0, hmy0 = zcx - gw / 2 + 12, ZY[0] + 96
for r_, row in enumerate(HM):
    for c_, colr in enumerate(row):
        rect(hmx0 + c_ * (cell + cg), hmy0 + r_ * (cell + cg), cell, cell, fill='#' + colr)
axx, axy = hmx0 - 7, hmy0 + gw + 7
arrow([(axx, axy), (axx, hmy0 - 4)], GATE, sw=2.2, head=9)
arrow([(axx, axy), (hmx0 + gw + 12, axy)], GATE, sw=2.2, head=9)
text(axx - 8, hmy0 + gw / 2 + 7, 'T', FS_B, MUTED, italic=True, anchor='end')
text(hmx0 + gw + 22, axy + 7, 'n', FS_B, MUTED, italic=True)

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
        path(f'M{f(cx - 10)} {f(cy)} L{f(cx - 4)} {f(cy)} L{f(cx - 1)} {f(cy - 8)} L{f(cx + 3)} {f(cy + 8)} '
             f'L{f(cx + 6)} {f(cy)} L{f(cx + 10)} {f(cy)}', stroke=SLATE, sw=2)
    elif kind == 'gate':
        gate(cx, cy, r=11)
    else:
        for ox in (-6, 2):
            path(f'M{f(cx + ox - 4)} {f(cy - 7)} L{f(cx + ox + 4)} {f(cy)} L{f(cx + ox - 4)} {f(cy + 7)} Z', fill=SLATE)


for k, (t1, kind) in enumerate([('generate inputs', 'doc'), ('submit', 'rack'), ('monitor', 'pulse'),
                                ('check gate', 'gate'), ('advance', 'next')]):
    x = lx0 + k * (CW_ + CG_)
    rect(x, ly0, CW_, lh, rx=12, fill='#FFFFFF', stroke='#C9D1DA', sw=1.4)
    loop_icon(kind, x + 28, ly0 + lh / 2)
    text(x + 52, ly0 + lh / 2 + FS_S * 0.36, t1, FS_S, INK, bold=True)
    if k < 4:
        arrow([(x + CW_ + 3, ly0 + lh / 2), (x + CW_ + CG_ - 3, ly0 + lh / 2)], SLATE, sw=2.4, head=9)
sx = lx0 + 1 * (CW_ + CG_) + CW_ / 2
cx_ = lx0 + 3 * (CW_ + CG_) + CW_ / 2
yb = ly0 + lh + 24
path(f'M{f(cx_)} {f(ly0 + lh + 2)} C{f(cx_)} {f(yb)} {f(cx_ - 30)} {f(yb)} {f((cx_ + sx) / 2)} {f(yb)} '
     f'C{f(sx + 30)} {f(yb)} {f(sx)} {f(yb)} {f(sx)} {f(ly0 + lh + 14)}', stroke=RED, sw=2.4, dash='7 6')
stealth((sx, ly0 + lh + 3), 0, -1, RED, length=12, half=5.5)
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
    yk = BY0 + 48 + (k // 3) * 64
    feat = k < 3
    rect(x, yk, cw, 54, rx=9, fill='#FFFFFF', stroke=col, sw=2.6 if feat else 1.6)
    text(x + cw / 2, yk + 27 + FS_B * 0.36, name, FS_B, tcol, bold=feat, anchor='middle')
rect(EN[0], BY0, EN[1] - EN[0], BY1 - BY0, rx=14, fill='#EFF1F4', stroke=MUTED, sw=2.2)
text(EN[0] + 16, BY0 + 32, 'Computational engines', FS_H, INK, bold=True)
ce = (EN[1] - EN[0] - 28 - 20) / 3
for k, (name, hw) in enumerate([('VASP', 'CPU GPU'), ('AMSET', 'CPU'), ('ShengBTE', 'CPU GPU'),
                                ('Pheasy', 'CPU GPU'), ('phono3py', 'CPU'), (None, None)]):
    x = EN[0] + 14 + (k % 3) * (ce + 10)
    yk = BY0 + 48 + (k // 3) * 64
    if name is None:
        rect(x, yk, ce, 54, rx=9, stroke='#A3ADB8', sw=1.6, dash='12 8')
        text(x + ce / 2, yk + 34, '…', FS_H, MUTED, bold=True, anchor='middle')
        continue
    rect(x, yk, ce, 54, rx=9, fill='#FFFFFF', stroke='#A3ADB8', sw=1.6)
    text(x + ce / 2, yk + 22, name, FS_B, INK, bold=True, anchor='middle')
    hw_tags(x + ce / 2, yk + 29, hw.split())
SEL_X = SK[1] + 8                # between the skill box and the engine box
end()

# library -> workflow: AutoZT picks the skills and engines a material needs
group('b-select-link')
curve_arrow((SEL_X, BY0 - 3), (SEL_X, BY0 - 24), (SEL_X, AS[3] + 24), (SEL_X, AS[3] + 3), INK, sw=3.2)
text(SEL_X + 16, (AS[3] + BY0) / 2 + 7, 'auto-selects', FS_S, INK, bold=True, italic=True)
end()

# a <-> b: tool calls go in, results come back; level with the "Submitted" card and the automation loop
group('a-b-interface')
GX0, GX1 = A_X[1] + 4, R_X[0] - 4
for yv, lab, col, rightward, ly in [(Y_MID - 14, 'tool calls', PURPLE, True, Y_MID - 26),
                                    (Y_MID + 14, 'results', CRIMSON, False, Y_MID + 40)]:
    x0, x1 = (GX0, GX1) if rightward else (GX1, GX0)
    curve_arrow((x0, yv), (x0 + (x1 - x0) * 0.3, yv), (x0 + (x1 - x0) * 0.7, yv), (x1, yv), col, sw=3.2)
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
    ym = (y0 + y1) / 2
    arrow([(x, ym), (x, y0 + 2)], GATE, sw=2.4, head=10)
    arrow([(x, ym), (x, y1 - 2)], GATE, sw=2.4, head=10)
    text(x + 12, ym + 8, 'c', FS_B, INK, italic=True)


SKW = 108
bkx, bky, bkh = A_X[0] + 46, D_Y[0] + 26, 150
rect(bkx, bky, SKW, bkh, rx=4, fill='#FFFFFF', stroke=FAINT, sw=2.2)
layer(bkx + 4, bky + 20, SKW - 8)
layer(bkx + 4, bky + 82, SKW - 8)
c_axis(bkx + SKW + 22, bky, bky + bkh)
text(bkx + SKW / 2, bky + bkh + 30, '3D bulk', FS_B, INK, bold=True, anchor='middle')
slx, sly, slh = bkx + 214, bky, 212
slab_y0, slab_h = sly + 69, 74
rect(slx, sly, SKW, slh, fill='#FFFFFF')
rect(slx, slab_y0, SKW, slab_h, fill='#E3F2EA')
rect(slx, sly, SKW, slh, stroke=FAINT, sw=2.2, dash='12 9')
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
EPS2D = [('ε', {'i': True}), ('2D', {'sub': True})]
GROUPS = [
    (['Structure'], INK, [[('Vacuum', {'b': True}), ' along ', ('c', {'i': True})],
                          [('Relaxation', {'b': True}), ' in-plane']]),
    (['Electrons'], BLUE_TX, [[('Elastic constant', {'b': True}), ' in N m', ('−1', {'sup': True})],
                              [('Dielectric', {'b': True}), ' without vacuum'],
                              [('Deformation ', {'b': True}), ('E', {'b': True, 'i': True}), ('1', {'sub': True, 'b': True}),
                               ' vs vacuum']]),
    (['Phonons'], ORANGE_TX, [[('NAC', {'b': True}), ' off'], [('ZA branch', {'b': True}), ' quadratic']]),
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
rich(fgx, fy0 + 82, [('t', {'i': True}), ': van der Waals thickness'], FS_XS, MUTED)
end()

add('</svg>')
dst = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fig1_autozt_overview_npj.svg')
with open(dst, 'w') as fh:
    fh.write('\n'.join(out) + '\n')
print(f'wrote {dst}: H {H}, chat gap {GAP:.1f}, Submitted card centre {Y_MID:.0f} vs automation loop centre '
      f'{(AS[2] + AS[3]) / 2:.0f}, panel c bottom {D_Y[1]}')
