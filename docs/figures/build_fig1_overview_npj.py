#!/usr/bin/env python3
"""Fig. 1 of the AutoZT paper (npj Comput. Mater.): overview, 180 mm wide vector SVG.

Layout follows the reference draft:
    a (agent session) | centre (coding agents -> AutoZT; loads b, runs c)           | b (skills, engines)
    c (POSCAR -> zT(n, T) pipeline)                                               | d (3D vs 2D rules)
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
H = 1480

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

# ================================================================ canvas
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{H / 10:.1f}mm" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add('<title>AutoZT overview (Fig. 1)</title>')
rect(0, 0, W, H, fill='#FFFFFF')

COL_A = (12, 614)
COL_M = (630, 1080)
COL_B = (1098, 1788)
R1_TOP, R1_BOT = 60, 784
R2_TOP, R2_BOT = 858, H - 12
TITLE_Y1, TITLE_Y2 = 42, 840


# panel b: skill box on top, engine box below; the hub in the centre column is level with the gap
SY0 = R1_TOP + 16
FEAT_H, SUP_H, MORE_H = 132, 88, 44
SH = 56 + FEAT_H + 14 + SUP_H + 14 + MORE_H + 32
EY0 = SY0 + SH + 16


def panel(x0, x1, y0, y1):
    rect(x0, y0, x1 - x0, y1 - y0, rx=16, fill=PANEL_BG, stroke=PANEL_LINE, sw=2)


def panel_title(x, y, letter, runs):
    text(x, y, letter, FS_LETTER, INK, bold=True)
    rich(x + 34, y, runs, FS_TITLE, INK, bold=True)


# ================================================================ a: agent session
group('panel-a-session')
panel_title(COL_A[0] + 2, TITLE_Y1, 'a', ['LLM operates AutoZT through MCP / CLI'])
panel(COL_A[0], COL_A[1], R1_TOP, R1_BOT)
text(COL_A[0] + 22, R1_TOP + 36, 'Coding agent session', FS_H, PURPLE, bold=True)
CH_X0, CH_X1 = COL_A[0] + 14, COL_A[1] - 14
CH_Y0, CH_Y1 = R1_TOP + 52, R1_BOT - 14
rect(CH_X0, CH_Y0, CH_X1 - CH_X0, CH_Y1 - CH_Y0, rx=14, fill=LAV)

AV_X = CH_X0 + 24                      # agent avatar centre
CARD_X0, CARD_X1 = CH_X0 + 48, CH_X1 - 10
USER_X1 = CH_X1 - 46                   # right edge of user bubbles
USER_AV = CH_X1 - 22
IX0, IX1 = CARD_X0 + 10, CARD_X1 - 10  # inner content edges of agent cards
GAP = 10
T_RH, T_PAD, T_GAP = 23, 6, 4


def tool_h(n):
    return 2 * T_PAD + n * T_RH + (n - 1) * T_GAP


def user_msg(y, runs):
    w = rich_w(runs, FS_B) + 36
    h = 32
    rect(USER_X1 - w, y, w, h, rx=11, fill=USER_BG)
    rich(USER_X1 - w / 2, y + h / 2 + FS_B * 0.35, runs, FS_B, INK, anchor='middle')
    person(USER_AV, y + h / 2 - 1, s=0.58)
    return y + h + GAP


VIA = [('MCP', {'fill': PURPLE}), (' · ', {'fill': FAINT, 'b': False}), ('CLI', {'fill': SLATE})]


def via(y):
    """Small right-aligned tag: the same call is available as an MCP tool and an agent-CLI command."""
    rich(IX1 - 10, y, VIA, FS_XS, PURPLE, bold=True, anchor='end')


def tool_block(y, rows):
    rect(IX0, y, IX1 - IX0, tool_h(len(rows)), rx=8, fill=LAV2)
    yy = y + T_PAD
    for stage, bg, fg, tools in rows:
        chip(IX0 + 7, yy, 96, T_RH, stage, bg, fg)
        rich(IX0 + 114, yy + T_RH / 2 + FS_XS * 0.36, tools, FS_XS, CODE)
        via(yy + T_RH / 2 + FS_XS * 0.36)
        yy += T_RH + T_GAP
    return y + tool_h(len(rows))


def agent_card(y, h):
    robot_avatar(AV_X, y + 24)
    rect(CARD_X0, y, CARD_X1 - CARD_X0, h, rx=12, fill='#FFFFFF', stroke=LAV_LINE, sw=1.6)


def stage(name):
    return (name, PINK, RED) if name == 'Approval' else (name, PEACH, PEACH_TX)


y = CH_Y0 + 12
# 1 user asks
y = user_msg(y, ['Compute '] + ZT + [' for this 2D POSCAR.'])
# 2 agent: discover, plan, dry run -> plan card
PLAN_H = 90
h2 = 10 + tool_h(3) + 34 + PLAN_H + 10
agent_card(y, h2)
yy = tool_block(y + 10, [stage('Discovery') + (['list_skills · check_env'],),
                         stage('Planning') + (['research_plan · preflight'],),
                         stage('Execution') + (['cycle (dry run)'],)])
text(IX0 + 4, yy + 25, 'Plan ready; the dry run is clean. Proceed?', FS_B, INK, bold=True)
py0 = yy + 34
rect(IX0, py0, IX1 - IX0, PLAN_H, rx=10, fill='#FBFAFE', stroke='#C9C0E4', sw=1.6)
rich(IX0 + 14, py0 + 23, ['2D thermoelectric '] + ZT + [' workflow'], FS_S, INK, bold=True)
pbx = IX1 - 12 - 92
pill_button(pbx, py0 + 7, 92, 'Proceed', PURPLE, '#FFFFFF')
pw = (IX1 - IX0 - 14 * 2 - 10) / 2
for k, (num, lab, col) in enumerate([('1', 'structure', GREEN), ('2', 'electron transport', BLUE),
                                     ('3', 'phonon transport', ORANGE), ('4', 'combine', CRIMSON)]):
    px = IX0 + 14 + (k % 2) * (pw + 10)
    pyy = py0 + 38 + (k // 2) * 26
    rect(px, pyy, pw, 22, rx=11, fill='#FFFFFF', stroke='#C9C0E4', sw=1.4)
    circle(px + 12, pyy + 11, 10, fill=col)
    text(px + 12, pyy + 11 + FS_XS * 0.36, num, FS_XS, '#FFFFFF', bold=True, anchor='middle')
    text(px + 29, pyy + 11 + FS_XS * 0.36, lab, FS_XS, INK, bold=True)
cursor(pbx + 92 - 8, py0 + 24)
y += h2 + GAP
# 3 agent: submit; AutoZT owns the job loop
h3 = 10 + tool_h(1) + 40
agent_card(y, h3)
yy = tool_block(y + 10, [stage('Execution') + (['register_material · cycle'],)])
text(IX0 + 4, yy + 25, 'Submitted; AutoZT advances and resumes jobs.', FS_B, INK, bold=True)
y += h3 + GAP
# 4 user asks for a destructive action
y = user_msg(y, ['Rerun the deformation step.'])
# 5 agent: approval gate in chat
h5 = 10 + tool_h(1) + 72
agent_card(y, h5)
yy = tool_block(y + 10, [stage('Approval') + (['request_destructive_action'],)])
text(IX0 + 4, yy + 25, 'This deletes its folder. Approve?', FS_B, INK, bold=True)
pill_button(IX0 + 4, yy + 36, 92, 'Approve', RED, '#FFFFFF')
pill_button(IX0 + 106, yy + 36, 92, 'Deny', '#FFFFFF', RED, stroke=RED)
rich(IX0 + 218, yy + 48 + FS_XS * 0.36, [('→ ', {'fill': MUTED}), 'approve_request'], FS_XS, CODE)
via(yy + 48 + FS_XS * 0.36)
y += h5 + GAP
# 6 agent: report
h6 = 10 + tool_h(1) + 40
agent_card(y, h6)
yy = tool_block(y + 10, [stage('Reporting') + (['results · task_graph'],)])
rich(IX0 + 4, yy + 25, ZTNT + [' is ready, with units and provenance.'], FS_B, INK, bold=True)
SESSION_END = y + h6
end()

# ================================================================ centre: coding agents -> AutoZT; b -> AutoZT -> c
group('panel-a-hub')
MX = (COL_M[0] + COL_M[1]) / 2
HUB_W, HUB_H = 300, 220
HUB_CY = EY0 - 8                 # level with the gap between the two boxes of panel b
HUB_CX = 818                     # over the gap between the transport boxes and the zT box of panel c
HUB_X0, HUB_Y0 = HUB_CX - HUB_W / 2, HUB_CY - HUB_H / 2
# coding agents: speech bubble whose tail points at AutoZT
BX0, BX1 = COL_M[0] + 16, COL_M[1] - 16
BY0 = R1_TOP
BH = HUB_Y0 - 22 - BY0
rect(BX0, BY0, BX1 - BX0, BH, rx=18, fill=LAV)
path(f'M{f(HUB_CX - 20)} {f(BY0 + BH - 1)} L{f(HUB_CX + 20)} {f(BY0 + BH - 1)} L{f(HUB_CX)} {f(BY0 + BH + 16)} Z',
     fill=LAV)
text(MX, BY0 + 42, 'Coding agents', FS_H, PURPLE, bold=True, anchor='middle')
cw = (BX1 - BX0 - 24) / 3
LOGO = 48
for k, (key, lab) in enumerate([('deepseek', 'DeepSeek'), ('kimi', 'Kimi'), ('glm', 'GLM'),
                                ('claude', 'Claude'), ('chatgpt', 'ChatGPT'), (None, '…')]):
    ccx = BX0 + 12 + (k % 3 + 0.5) * cw
    gy = BY0 + 70 + (k // 3) * 100
    if key:
        logo(key, ccx - LOGO / 2, gy, LOGO)
        text(ccx, gy + LOGO + 22, lab, FS_B, INK, bold=True, anchor='middle')
    else:
        text(ccx, gy + LOGO / 2 + 8, '…', FS_H, MUTED, bold=True, anchor='middle')
# AutoZT hub: a small DAG in the colours of panel c, then the name
rect(HUB_X0, HUB_Y0, HUB_W, HUB_H, rx=20, fill='#FFFFFF', stroke=INK, sw=5)
dcx, dcy, ddx, ddy = HUB_CX, HUB_Y0 + 62, 52, 22
nodes = {'s': (dcx - 1.5 * ddx, dcy, GREEN), 'r': (dcx - 0.5 * ddx, dcy, INK),
         'e': (dcx + 0.5 * ddx, dcy - ddy, BLUE), 'p': (dcx + 0.5 * ddx, dcy + ddy, ORANGE),
         'z': (dcx + 1.5 * ddx, dcy, CRIMSON)}
for a_, b_ in [('s', 'r'), ('r', 'e'), ('r', 'p'), ('e', 'z'), ('p', 'z')]:
    line(nodes[a_][0], nodes[a_][1], nodes[b_][0], nodes[b_][1], FAINT, 2.6)
for nx, ny, ncol in nodes.values():
    circle(nx, ny, 11, fill=ncol)
text(HUB_CX, HUB_Y0 + 150, 'AutoZT', 44, INK, bold=True, anchor='middle')
text(HUB_CX, HUB_Y0 + 190, 'workflow engine', FS_S, MUTED, anchor='middle')
# AutoZT -> c
arrow([(HUB_CX, HUB_Y0 + HUB_H + 4), (HUB_CX, R2_TOP - 2)], INK, sw=3.2, head=13)
text(HUB_CX + 14, (HUB_Y0 + HUB_H + R2_TOP) / 2 + 6, 'runs', FS_S, INK, bold=True, italic=True)
end()

# ================================================================ b: skills and engines
group('panel-b')
panel_title(COL_B[0] + 2, TITLE_Y1, 'b', ['Workflow skills and engines'])
panel(COL_B[0], COL_B[1], R1_TOP, R1_BOT)
SX0, SX1 = COL_B[0] + 16, COL_B[1] - 16
GX0, GX1 = SX0 + 18, SX1 - 18
CW3 = (GX1 - GX0 - 2 * 12) / 3
rect(SX0, SY0, SX1 - SX0, SH, rx=14, fill=GREEN_BG, stroke=GREEN, sw=2.2)
text(SX0 + 20, SY0 + 35, 'Workflow skills', FS_H, GREEN, bold=True)
text(SX1 - 20, SY0 + 35, 'skill.yaml: steps + acceptance checks', FS_XS, MUTED, italic=True, anchor='end')
fy = SY0 + 56
for k, (name, col, tcol, role, q) in enumerate([
        ('ke-dft-cpu', BLUE, BLUE_TX, 'electron transport',
         [('σ', {'i': True}), ', ', ('S', {'i': True}), ', ', ('κ', {'i': True}), ('e', {'sub': True})]),
        ('kl-dft-cpu', ORANGE, ORANGE_TX, 'phonon transport', [('κ', {'i': True}), ('L', {'sub': True})]),
        ('zt-dft-cpu', CRIMSON, CRIMSON, 'combine', ZTNT)]):
    x = GX0 + k * (CW3 + 12)
    rect(x, fy, CW3, FEAT_H, rx=10, fill='#FFFFFF', stroke=col, sw=3)
    text(x + CW3 / 2, fy + 42, name, FS_H, tcol, bold=True, anchor='middle')
    text(x + CW3 / 2, fy + 72, role, FS_XS, MUTED, anchor='middle')
    rich(x + CW3 / 2, fy + 106, q, FS_B, tcol, bold=True, anchor='middle')
sy = fy + FEAT_H + 14
for k, (name, role) in enumerate([('band-dft-cpu', 'band structure'), ('phonon-dft-cpu', 'phonon dispersion'),
                                  ('fit-fc-thermal', 'force-constant fitting')]):
    x = GX0 + k * (CW3 + 12)
    rect(x, sy, CW3, SUP_H, rx=10, fill='#FFFFFF', stroke=GREEN_SOFT, sw=1.8)
    text(x + CW3 / 2, sy + 38, name, FS_B, INK, anchor='middle')
    text(x + CW3 / 2, sy + 65, role, FS_XS, MUTED, anchor='middle')
sy += SUP_H + 14
rect(GX0, sy, GX1 - GX0, MORE_H, rx=10, stroke=GREEN_SOFT, sw=1.8, dash='14 9')
text((GX0 + GX1) / 2, sy + MORE_H / 2 + FS_XS * 0.36, '…  plus MLFF, defect, elastic and screening skills',
     FS_XS, MUTED, anchor='middle')

EH = R1_BOT - 16 - EY0
rect(SX0, EY0, SX1 - SX0, EH, rx=14, fill='#EFF1F4', stroke=MUTED, sw=2.2)
text(SX0 + 20, EY0 + 35, 'Computational engines', FS_H, INK, bold=True)
text(SX1 - 20, EY0 + 35, 'on CPU / GPU clusters', FS_XS, MUTED, italic=True, anchor='end')
ey = EY0 + 54
eh = (EY0 + EH - 16 - ey - 12) / 2
for k, (name, role, rcol, hw) in enumerate([('VASP', 'DFT', MUTED, 'CPU GPU'),
                                            ('AMSET', 'carrier transport', BLUE_TX, 'CPU'),
                                            ('ShengBTE', 'phonon BTE', ORANGE_TX, 'CPU GPU'),
                                            ('Pheasy', 'force constants', ORANGE_TX, 'CPU GPU'),
                                            ('phono3py', 'anharmonic phonons', ORANGE_TX, 'CPU'),
                                            (None, None, None, None)]):
    x = GX0 + (k % 3) * (CW3 + 12)
    yk = ey + (k // 3) * (eh + 12)
    if name is None:
        rect(x, yk, CW3, eh, rx=10, stroke='#A3ADB8', sw=1.8, dash='14 9')
        text(x + CW3 / 2, yk + eh / 2 + 8, '…', FS_H, MUTED, bold=True, anchor='middle')
        continue
    rect(x, yk, CW3, eh, rx=10, fill='#FFFFFF', stroke='#A3ADB8', sw=1.8)
    mid = yk + eh / 2
    text(x + CW3 / 2, mid - 18, name, FS_H, INK, bold=True, anchor='middle')
    text(x + CW3 / 2, mid + 6, role, FS_XS, rcol, anchor='middle')
    tags = hw.split()
    widths = [tw(t, FS_XS, True) + 16 for t in tags]
    tx0 = x + CW3 / 2 - (sum(widths) + 6 * (len(tags) - 1)) / 2
    for t, wd in zip(tags, widths):
        bg, fg = ('#E3E7EC', GATE) if t == 'CPU' else ('#DCEFE3', '#1E6B3A')
        rect(tx0, mid + 16, wd, 22, rx=11, fill=bg)
        text(tx0 + wd / 2, mid + 27 + FS_XS * 0.36, t, FS_XS, fg, bold=True, anchor='middle')
        tx0 += wd + 6
end()

group('panel-a-hub-link')
arrow([(COL_B[0] + 26, HUB_CY), (HUB_X0 + HUB_W + 4, HUB_CY)], INK, sw=3.2, head=13)
text((HUB_X0 + HUB_W + COL_B[0]) / 2 + 2, HUB_CY - 12, 'loads', FS_S, INK, bold=True, italic=True, anchor='middle')
end()

# ================================================================ c: POSCAR -> zT(n, T)
group('panel-c')
CX0, CX1 = COL_A[0], COL_M[1]
panel_title(CX0 + 2, TITLE_Y2, 'c', ['From POSCAR to '] + ZTNT)
panel(CX0, CX1, R2_TOP, R2_BOT)
# POSCAR -> Dimension (left column)
RX_Y = R2_TOP + 300                    # centre line of the relaxation step
TH_ = 104
EL_Y, PH_Y = RX_Y - 160, RX_Y + 160
DCX = CX0 + 124
dw, dh = 116, 140
dx0, dy0 = DCX - dw / 2, R2_TOP + 64
text(DCX, R2_TOP + 46, 'POSCAR', FS_H, INK, bold=True, anchor='middle')
fold = 28
path(f'M{f(dx0)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw)} {f(dy0 + fold)} L{f(dx0 + dw)} {f(dy0 + dh)} '
     f'L{f(dx0)} {f(dy0 + dh)} Z', fill='#FFFFFF', stroke=INK, sw=2.8)
path(f'M{f(dx0 + dw - fold)} {f(dy0)} L{f(dx0 + dw - fold)} {f(dy0 + fold)} L{f(dx0 + dw)} {f(dy0 + fold)}', stroke=INK, sw=2.4)
mcx, mcy, mr = DCX, dy0 + dh / 2 + 10, 32
hexp = [(mcx + mr * math.cos(math.radians(a)), mcy - mr * math.sin(math.radians(a))) for a in (90, 30, -30, -90, -150, 150)]
for i in range(6):
    line(*hexp[i], *hexp[(i + 1) % 6], '#8AA0C0', 2.6)
for i in (0, 2, 4):
    line(mcx, mcy, *hexp[i], '#8AA0C0', 2.6)
for i, (px, py) in enumerate(hexp):
    circle(px, py, 9, fill='#D9534F' if i % 2 == 0 else '#3B6FC4')
circle(mcx, mcy, 9, fill='#3B6FC4')
DM_Y = RX_Y                            # Dimension -> Relaxation on one line
dmw, dmh = 188, 100
arrow([(DCX, dy0 + dh + 6), (DCX, DM_Y - dmh / 2 - 2)], GREEN, sw=3, head=12)
rect(DCX - dmw / 2, DM_Y - dmh / 2, dmw, dmh, rx=12, fill=GREEN_BG, stroke=GREEN, sw=2.6)
text(DCX, DM_Y - 10, 'Dimension', FS_H, GREEN, bold=True, anchor='middle')
text(DCX, DM_Y + 28, '2D / 3D', FS_H, INK, bold=True, anchor='middle')
# Relaxation
RX_X, rxw, rxh = CX0 + 372, 180, 92
arrow([(DCX + dmw / 2 + 2, DM_Y), (RX_X - rxw / 2 - 3, RX_Y)], INK, sw=3, head=12)
rect(RX_X - rxw / 2, RX_Y - rxh / 2, rxw, rxh, rx=12, fill='#FFFFFF', stroke=INK, sw=2.8)
text(RX_X, RX_Y + 9, 'Relaxation', FS_H, INK, bold=True, anchor='middle')
gate(RX_X + rxw / 2 - 4, RX_Y - rxh / 2 + 2)
# transport branches; the label above each box names the skill from panel b
TX0, TW_ = CX0 + 528, 256
brx = RX_X + rxw / 2 + 36
arrow([(RX_X + rxw / 2 + 2, RX_Y - 12), (brx, RX_Y - 12), (brx, EL_Y), (TX0 - 3, EL_Y)], BLUE, sw=3, head=12)
arrow([(RX_X + rxw / 2 + 2, RX_Y + 12), (brx, RX_Y + 12), (brx, PH_Y), (TX0 - 3, PH_Y)], ORANGE, sw=3, head=12)
for yc, col, bg, tcol, title, skill, q in [
        (PH_Y, ORANGE, ORANGE_BG, ORANGE_TX, 'Phonon transport', 'kl-dft-cpu', [('κ', {'i': True}), ('L', {'sub': True})]),
        (EL_Y, BLUE, BLUE_BG, BLUE_TX, 'Electron transport', 'ke-dft-cpu',
         [('σ', {'i': True}), ', ', ('S', {'i': True}), ', ', ('κ', {'i': True}), ('e', {'sub': True})])]:
    rect(TX0, yc - TH_ / 2, TW_, TH_, rx=12, fill=bg, stroke=col, sw=2.8)
    text(TX0 + 2, yc - TH_ / 2 - 9, skill, FS_XS, tcol, bold=True)
    text(TX0 + TW_ / 2, yc - 8, title, FS_H, tcol, bold=True, anchor='middle')
    rich(TX0 + TW_ / 2, yc + 28, q, FS_B, tcol, bold=True, anchor='middle')
    gate(TX0 + TW_ - 4, yc - TH_ / 2 + 2)
# zT box (spans both branches): formula + zT(n, T) map
ZX0, ZX1 = TX0 + TW_ + 44, CX1 - 22
ZY0, ZY1 = EL_Y - TH_ / 2, PH_Y + TH_ / 2
for yc, col in [(PH_Y, ORANGE), (EL_Y, BLUE)]:
    arrow([(TX0 + TW_ + 3, yc), (ZX0 - 3, yc)], col, sw=3, head=12)
rect(ZX0, ZY0, ZX1 - ZX0, ZY1 - ZY0, rx=14, fill=CRIMSON_BG, stroke=CRIMSON, sw=2.8)
text(ZX0 + 2, ZY0 - 9, 'zt-dft-cpu', FS_XS, CRIMSON, bold=True)
zcx = (ZX0 + ZX1) / 2
FZ = 28
eqy = ZY0 + 74
lhs = [('zT', {'i': True}), ' =']
num = [('S', {'i': True, 'fill': BLUE}), ('2', {'sup': True, 'fill': BLUE}), ('σ', {'i': True, 'fill': BLUE}),
       ('T', {'i': True})]
den = [('κ', {'i': True, 'fill': BLUE}), ('e', {'sub': True, 'fill': BLUE}), ' + ', ('κ', {'i': True, 'fill': ORANGE}),
       ('L', {'sub': True, 'fill': ORANGE})]
lw_ = rich_w(lhs, FZ)
fw = max(rich_w(num, FZ), rich_w(den, FZ)) + 14
ex0 = zcx - (lw_ + 10 + fw) / 2
rich(ex0, eqy + 10, lhs, FZ, INK)
fx = ex0 + lw_ + 10
rich(fx + fw / 2, eqy - 10, num, FZ, INK, anchor='middle')
line(fx, eqy, fx + fw, eqy, INK, 2.4)
rich(fx + fw / 2, eqy + 33, den, FZ, INK, anchor='middle')
HM = [['EEB2C6', 'D16388', 'B23C66', '8A1D46', 'B23C66', 'E38AA8'],
      ['EEB2C6', 'E38AA8', 'D16388', 'B23C66', 'D16388', 'E38AA8'],
      ['F6D3DE', 'E38AA8', 'D16388', 'D16388', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'EEB2C6', 'E38AA8', 'E38AA8', 'E38AA8', 'EEB2C6'],
      ['F6D3DE', 'F6D3DE', 'EEB2C6', 'EEB2C6', 'EEB2C6', 'F6D3DE'],
      ['F6D3DE', 'F6D3DE', 'F6D3DE', 'EEB2C6', 'F6D3DE', 'F6D3DE']]
cell, cg = 25, 1.5
gw = 6 * cell + 5 * cg
hmx0 = zcx - gw / 2 + 10
hmy0 = ZY0 + 190
rich(zcx, hmy0 - 22, ZTNT, FS_H, CRIMSON, bold=True, anchor='middle')
for r_, row in enumerate(HM):
    for c_, col in enumerate(row):
        rect(hmx0 + c_ * (cell + cg), hmy0 + r_ * (cell + cg), cell, cell, fill='#' + col)
axx, axy = hmx0 - 9, hmy0 + gw + 9
arrow([(axx, axy), (axx, hmy0 - 4)], GATE, sw=2.4, head=10)
arrow([(axx, axy), (hmx0 + gw + 14, axy)], GATE, sw=2.4, head=10)
text(axx - 10, hmy0 + gw / 2 + 8, 'T', FS_B, MUTED, italic=True, anchor='end')
text(hmx0 + gw / 2, axy + 30, 'n', FS_B, MUTED, italic=True, anchor='middle')
# legend
gate(CX0 + 36, R2_BOT - 30, r=12)
text(CX0 + 56, R2_BOT - 30 + FS_XS * 0.36, 'acceptance gate: checks must pass before the next stage runs', FS_XS, MUTED)
end()

# ================================================================ d: 3D and 2D
group('panel-d')
DX0, DX1 = COL_B
panel_title(DX0 + 2, TITLE_Y2, 'd', ['One workflow for 3D and 2D materials'])
panel(DX0, DX1, R2_TOP, R2_BOT)
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


SKX, SKW = DX0 + 60, 108
bk_y0, bk_h = R2_TOP + 30, 160
rect(SKX, bk_y0, SKW, bk_h, rx=4, fill='#FFFFFF', stroke=FAINT, sw=2.2)
layer(SKX + 4, bk_y0 + 22, SKW - 8)
layer(SKX + 4, bk_y0 + 92, SKW - 8)
c_axis(SKX + SKW + 22, bk_y0, bk_y0 + bk_h)
text(SKX + SKW / 2, bk_y0 + bk_h + 30, '3D bulk', FS_B, INK, bold=True, anchor='middle')
sl_y0, sl_h = bk_y0 + bk_h + 58, 250
slab_y0, slab_h = sl_y0 + 88, 74
rect(SKX, sl_y0, SKW, sl_h, fill='#FFFFFF')
rect(SKX, slab_y0, SKW, slab_h, fill='#E3F2EA')
rect(SKX, sl_y0, SKW, sl_h, stroke=FAINT, sw=2.2, dash='12 9')
text(SKX + SKW / 2, sl_y0 + 50, 'vacuum', FS_XS, FAINT, italic=True, anchor='middle')
text(SKX + SKW / 2, sl_y0 + sl_h - 38, 'vacuum', FS_XS, FAINT, italic=True, anchor='middle')
layer(SKX + 4, slab_y0 + 13, SKW - 8)
path(f'M{f(SKX - 8)} {f(slab_y0)} L{f(SKX - 18)} {f(slab_y0)} L{f(SKX - 18)} {f(slab_y0 + slab_h)} '
     f'L{f(SKX - 8)} {f(slab_y0 + slab_h)}', stroke=GREEN, sw=2.6, cap='butt', join='miter')
text(SKX - 26, slab_y0 + slab_h / 2 + 9, 't', FS_H, GREEN, bold=True, italic=True, anchor='end')
c_axis(SKX + SKW + 22, sl_y0, sl_y0 + sl_h)
text(SKX + SKW / 2, sl_y0 + sl_h + 30, '2D slab', FS_B, GREEN, bold=True, anchor='middle')
rich(DX0 + 22, R2_BOT - 20, [('t', {'i': True}), ': van der Waals thickness'], FS_XS, MUTED)

TX = DX0 + 236           # item column
TX2 = TX + 196           # treatment column
TXE = DX1 - 20
yy = R2_TOP + 52
chip(TX, yy - 19, 40, 26, '3D', '#4A5563', '#FFFFFF')
text(TX + 52, yy, 'Bulk crystals', FS_B, INK, bold=True)
text(TX2, yy, 'standard path', FS_B, INK)
yy += 40
chip(TX, yy - 19, 40, 26, '2D', GREEN, '#FFFFFF')
text(TX + 52, yy, 'Slabs: same path, plus:', FS_B, GREEN, bold=True)
GROUPS = [
    (['Structure'], INK, [(['Vacuum axis'], ['normalised to ', ('c', {'i': True})]),
                          (['Relaxation'], ['in-plane only'])]),
    (['Electrons'], BLUE_TX, [(['Elastic constant'], ['in N m', ('−1', {'sup': True})]),
                              (['Dielectric ', ('ε', {'i': True})], ['vacuum removed']),
                              (['Deformation ', ('E', {'i': True}), ('1', {'sub': True})], ['vacuum-level reference'])]),
    (['Phonons'], ORANGE_TX, [(['NAC'], ['switched off']),
                              (['ZA branch'], ['quadratic, gated'])]),
    (ZT + [' inputs'], CRIMSON, [([('σ', {'i': True}), ', ', ('κ', {'i': True}), ('e', {'sub': True}), ', ',
                                   ('κ', {'i': True}), ('L', {'sub': True})],
                                  ['in-plane, × ', ('c', {'i': True}), ' / ', ('t', {'i': True})])]),
]
yy += 14
for gname, col, items in GROUPS:
    yy += 34
    rich(TX, yy, gname, FS_XS, col, bold=True)
    line(TX + rich_w(gname, FS_XS, True) + 10, yy - 6, TXE, yy - 6, '#C9D1DA', 1.6)
    for left, right in items:
        yy += 31
        rich(TX + 12, yy, left, FS_B, col, bold=True)
        rich(TX2, yy, right, FS_B, INK)
fy0 = yy + 24
rect(TX, fy0, TXE - TX, 56, rx=10, fill='#F3FAF6', stroke=GREEN, sw=2.2)
rich((TX + TXE) / 2, fy0 + 37, [('ε', {'i': True}), ('2D', {'sub': True}), ' = 1 + (', ('c', {'i': True}), ' / ',
                                 ('t', {'i': True}), ') (', ('ε', {'i': True}), ('slab', {'sub': True}), ' − 1)'],
     FS_H, INK, anchor='middle')
D_END = fy0 + 56
end()

add('</svg>')
dst = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fig1_autozt_overview_npj.svg')
with open(dst, 'w') as fh:
    fh.write('\n'.join(out) + '\n')
print(f'wrote {dst}: session ends {SESSION_END:.0f} (row bottom {R1_BOT}), hub {HUB_Y0:.0f}-{HUB_Y0 + HUB_H:.0f}, '
      f'skills box {SY0}-{SY0 + SH}, d table ends {D_END:.0f} (row bottom {R2_BOT})')
