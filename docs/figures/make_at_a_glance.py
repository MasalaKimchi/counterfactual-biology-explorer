"""
Build the 'at-a-glance' verdict scoreboard for Cell-State Reachability.

A complementary hero to fig_central_illustration.png (which is the method schematic).
This one is the 5-second scoreboard a reviewer sees first: the reframe, the flagship
verdict + headline numbers, the signed decomposition, the four outputs, and the honest
scope line. Every number is taken verbatim from manuscript_facts.json.

Output: fig_at_a_glance.png (+ .pdf), matplotlib to match the repo's figure toolchain.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

# ---- palette (matched to fig_central_illustration.png) ----------------------
NAVY   = "#1F3A4D"   # primary text / rules
INK    = "#2C4A5E"   # secondary text
MUTE   = "#6B7C86"   # muted captions
TEAL   = "#2E7D6B"   # LOF / reachable
GOLD   = "#C0902F"   # GOF / activation-required
GRAY   = "#B9B3A6"   # neither / residual
PANEL  = "#EAF1F4"   # light panel fill
PANELB = "#D3E1E8"   # panel border
CREAM  = "#F6F3EC"   # warm band
WHITE  = "#FFFFFF"
GREENBG= "#E7F1EE"   # transfer badge fill
GREENBD= "#8FBFB1"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

W, H = 14.0, 8.0
fig = plt.figure(figsize=(W, H), dpi=200)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.axis("off")
fig.patch.set_facecolor(WHITE)


def rbox(x, y, w, h, *, fc=WHITE, ec=PANELB, lw=1.4, r=0.9, z=1):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        fc=fc, ec=ec, lw=lw, zorder=z, mutation_aspect=H / W))


def txt(x, y, s, *, size=11, color=NAVY, weight="normal", ha="left", va="center",
        style="normal", z=5, spacing=None):
    t = ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
                style=style, zorder=z)
    if spacing:
        t.set_linespacing(spacing)
    return t


# ============================ TITLE BAND =====================================
txt(4.0, 94.4, "Cell-State Reachability", size=27, color=NAVY, weight="bold")
txt(4.2, 89.0, "Can knockdown point a cell where you want it to go?",
    size=15.5, color=TEAL, weight="bold", style="italic")
txt(96.0, 94.2, "Built with Claude", size=10.5, color=MUTE, ha="right", weight="bold")
txt(96.0, 90.6, "Life-Sciences Hackathon · Research / Lab track",
    size=9.5, color=MUTE, ha="right")
ax.add_line(Line2D([4.0, 96.0], [86.4, 86.4], color=PANELB, lw=1.3, zorder=1))

# ============================ REFRAME LINE ===================================
rbox(4.0, 77.6, 92.0, 7.6, fc=CREAM, ec="#E4DCC8", lw=1.2, r=0.8)
txt(6.0, 83.2, "THE REFRAME", size=9.5, color=GOLD, weight="bold")
txt(6.0, 80.55,
    "Most tools rank what a perturbation might do.  This asks an earlier question — can non-negative combinations of a screen's",
    size=11.3, color=INK)
txt(6.0, 78.35,
    "measured effects point the transcriptome toward the target — and it certifies when they provably cannot.",
    size=11.3, color=INK)

# ============================ FLAGSHIP (left) ================================
FX, FW = 4.0, 55.0
rbox(FX, 40.0, FW, 35.2, fc=WHITE, ec=PANELB, lw=1.5, r=0.9)
txt(FX + 2.2, 71.9, "FLAGSHIP  ·  Th2 → Th1,  resting primary human CD4⁺ T cells",
    size=12.0, color=NAVY, weight="bold")
ax.add_line(Line2D([FX + 2.2, FX + FW - 2.2], [69.4, 69.4], color=PANELB, lw=1.1))

# big verdict stat
txt(FX + 2.4, 62.6, "0.448", size=40, color=TEAL, weight="bold", va="center")
txt(FX + 20.5, 65.4, "held-out cosine", size=12.5, color=NAVY, weight="bold")
txt(FX + 20.5, 61.7, "in-sample 0.627 · the honest,", size=10.3, color=MUTE)
txt(FX + 20.5, 58.9, "cross-validated number", size=10.3, color=MUTE)

# verdict chip
rbox(FX + 2.4, 49.6, 34.0, 5.6, fc=GREENBG, ec=GREENBD, lw=1.3, r=0.7)
txt(FX + 4.0, 52.4, "VERDICT:  PARTIALLY REACHABLE", size=11.5, color="#1E6B58", weight="bold")

# two small stat tiles
rbox(FX + 2.4, 41.2, 25.2, 6.9, fc=PANEL, ec=PANELB, lw=1.1, r=0.6)
txt(FX + 3.6, 45.7, "> all 60 shuffled targets", size=10.2, color=NAVY, weight="bold")
txt(FX + 3.6, 43.0, "empirical p = 1/61  ·  z ≈ 24", size=9.6, color=MUTE)

rbox(FX + 29.2, 41.2, 23.4, 6.9, fc=PANEL, ec=PANELB, lw=1.1, r=0.6)
txt(FX + 30.4, 45.7, "KKT / Farkas residual", size=10.2, color=NAVY, weight="bold")
txt(FX + 30.4, 43.0, "1.1 × 10⁻¹¹  (cone optimality)", size=9.6, color=MUTE)

# ============================ DECOMPOSITION (right) ==========================
DX, DW = 61.5, 34.5
rbox(DX, 40.0, DW, 35.2, fc=WHITE, ec=PANELB, lw=1.5, r=0.9)
txt(DX + 2.2, 71.9, "SIGNED DECOMPOSITION", size=12.0, color=NAVY, weight="bold")
txt(DX + 2.2, 68.6, "of the target direction", size=10.2, color=MUTE)

# stacked horizontal bar  (exact fractions: 0.393 + 0.253 + 0.354 = 1.000)
bx, by, bw, bh = DX + 2.4, 59.4, DW - 4.8, 6.4
segs = [(0.393, TEAL, "39%"), (0.253, GOLD, "25%"), (0.354, GRAY, "35%")]
cur = bx
for frac, col, lab in segs:
    ax.add_patch(Rectangle((cur, by), bw * frac, bh, fc=col, ec=WHITE, lw=1.6, zorder=3))
    tc = WHITE if col != GRAY else NAVY
    txt(cur + bw * frac / 2, by + bh / 2, lab, size=12.5, color=tc, weight="bold", ha="center")
    cur += bw * frac

# legend rows (single line each)
def legrow(y, col, name):
    ax.add_patch(Rectangle((DX + 2.4, y - 0.95), 2.0, 2.0, fc=col, ec=WHITE, lw=1, zorder=3))
    txt(DX + 5.6, y, name, size=10.4, color=NAVY, va="center")

legrow(54.2, TEAL, "LOF 39% — reachable by knockdown")
legrow(49.6, GOLD, "GOF 25% — needs activation (CRISPRa)")
legrow(45.0, GRAY, "neither 35% — irreducible residual")

txt(DX + 2.2, 41.4, "knockdown is never the majority modality  (atlas mean LOF 0.34)",
    size=9.3, color=MUTE, style="italic")

# ============================ FOUR OUTPUTS ===================================
txt(4.0, 35.7, "FOUR OUTPUTS, EVERY RUN", size=11.5, color=NAVY, weight="bold")
cards = [
    ("①  Directional verdict", "held-out + null-calibrated,\nnot a similarity ranking"),
    ("②  Ranked sparse panel", "greedy mixture → compact\nknockdown candidate set"),
    ("③  Dual certificate", "Farkas/KKT proof when the\ntarget is outside the cone"),
    ("④  Unmet-readout list", "ranked CRISPRa / de-repression\nhypotheses — to be tested"),
]
cw, gap = 21.9, 1.45
x0 = 4.0
for i, (head, body) in enumerate(cards):
    x = x0 + i * (cw + gap)
    rbox(x, 24.0, cw, 9.7, fc=PANEL, ec=PANELB, lw=1.2, r=0.7)
    txt(x + 1.5, 31.0, head, size=11.0, color=TEAL, weight="bold")
    txt(x + 1.5, 27.0, body, size=9.4, color=INK, spacing=1.28, va="center")

# ============================ TRANSFER BADGE =================================
rbox(4.0, 14.6, 92.0, 6.9, fc=GREENBG, ec=GREENBD, lw=1.3, r=0.8)
txt(6.0, 19.1, "SAME OPERATOR, NO RETUNING", size=9.8, color="#1E6B58", weight="bold")
txt(6.0, 16.3,
    "Transfers unchanged to a Norman K562 CRISPRa screen (held-out CEBPA state, cosine 0.878).   "
    "Additivity guardrail: 126 doubles, median cos 0.71 (bounded, not exact).",
    size=10.4, color=INK)

# ============================ SCOPE FOOTER ===================================
ax.add_line(Line2D([4.0, 96.0], [11.2, 11.2], color=PANELB, lw=1.1, zorder=1))
txt(4.0, 8.2, "SCOPE  ", size=9.5, color=GOLD, weight="bold", va="center")
txt(11.2, 8.3,
    "Directional feasibility relative to one measured screen — not functional rescue.  Combinations assume bounded additivity.",
    size=10.2, color=MUTE, va="center")
txt(11.2, 5.4,
    "Recovered regulators (GATA3↓, TBX21) are positive controls, not discoveries.  Every nomination is a wet-lab hypothesis, never a validated target.",
    size=10.2, color=MUTE, va="center")

# source data footnote
txt(96.0, 2.2, "Source: Zhu et al. 2025, genome-scale CRISPRi Perturb-seq in CD4⁺ T cells (CZI VCP)  ·  numbers verbatim from manuscript_facts.json",
    size=8.2, color="#9AA6AD", ha="right", va="center")

import os
out = os.path.dirname(os.path.abspath(__file__))
fig.savefig(f"{out}/fig_at_a_glance.png", dpi=200, facecolor=WHITE,
            bbox_inches="tight", pad_inches=0.12)
fig.savefig(f"{out}/fig_at_a_glance.pdf", facecolor=WHITE,
            bbox_inches="tight", pad_inches=0.12)
print("wrote fig_at_a_glance.png / .pdf")
