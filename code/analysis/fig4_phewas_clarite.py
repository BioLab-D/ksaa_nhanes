# Figure 4: draws panel A with clarite.plot.manhattan and combines it with panel B (drawn by fig4_phewas.R).
#   CLARITE defaults are kept (alphabetical categories, alternating gray/white bands, FDR 5% line, labels right of points);
#   only the layout is adjusted (font size, axis ranges, label overlap, category colors matching panel B).
# Input: results/figures/Figure4A.csv and Figure4B.{png,pdf} from fig4_phewas.R
# Output: results/figures/Figure4A.{png,pdf}, Figure4.{png,pdf}
# Run: $CLARITE_PY code/analysis/fig4_phewas_clarite.py . (CLARITE 2.3.6 needs Python <= 3.10, plus Pillow and PyMuPDF)
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont
import clarite
from clarite.modules.analyze import result_columns, corrected_pvalue_columns

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
FIG = f"{ROOT}/results/figures"
for f in ["/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Supplemental/Arial Bold.ttf"]:
    font_manager.fontManager.addfont(f)
plt.rcParams.update({"pdf.fonttype": 42, "font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.4, "ytick.major.width": 0.4})

A = pd.read_csv(f"{FIG}/Figure4A.csv")
n = len(A)
# CLARITE sorts variables by name within a category, so it gets the NHANES codes; display labels are set after plotting.
df = pd.DataFrame({c: np.nan for c in result_columns + corrected_pvalue_columns}, index=pd.MultiIndex.from_arrays([A["variable"], ["KSAA"] * n], names=["Variable", "Outcome"]))
df["pvalue"] = A["p"].values; df["pvalue_fdr"] = A["fdr"].values; df["pvalue_bonferroni"] = np.minimum(1, A["p"].values * n)
df["N"] = n; df["Converged"] = True; df["Variable_type"] = "continuous"
cats = dict(zip(A["variable"], A["category"])); SHOW = dict(zip(A["variable"], A["label"]))
cat_order = sorted(A["category"].unique())                       # CLARITE sorts by category name
colors = [A.loc[A["category"] == c, "colour"].iloc[0] for c in cat_order]
bg = ["#EBEBEB" if i % 2 == 0 else "#FFFFFF" for i in range(len(cat_order))]
LABEL_CRP = False   # True = annotate the reference item C-reactive protein (label, leader line, black ring)
labels = list(A.loc[A["is_top"], "variable"]) + (["CRP"] if LABEL_CRP else [])

# Layout after plotting: labels are moved vertically to the nearest spot where they cover no point or label; in the three
#   rightmost categories they form a right-aligned column with leader lines. Category legend at the bottom instead of
#   x-axis names; y ticks every 5; font size 5.2 pt as in Figures 2 and 3.
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
FS = 5.2
# figsize: CLARITE adds its own margins; about 3.85 x 4.2 in after saving with a tight bounding box
# FDR 5% line only (red dotted); no Bonferroni line
fig = clarite.plot.manhattan({"": df}, categories=cats, bonferroni=None, fdr=0.05, num_labeled=0, label_vars=labels,
                             figsize=(4.15, 4.95), dpi=300, colors=colors, background_colors=bg, return_figure=True)   # so the full figure fits the 183 mm full width
ax = fig.axes[0]
ymax = -np.log10(A["p"].min())
ax.set_ylim(-0.4, ymax * 1.10); ax.set_title(""); fig.suptitle("")
ax.set_ylabel("−log$_{10}$ P", fontsize=FS, labelpad=2); ax.yaxis.set_major_locator(MultipleLocator(5))
ax.tick_params(axis="y", labelsize=FS, labelrotation=0, length=2, width=0.4, pad=1.5)
ax.set_xticks([]); ax.set_xlabel("Phenotype", fontsize=FS, labelpad=3)
for sp in ax.spines.values(): sp.set_linewidth(0.5)
for coll in ax.collections:
    try: coll.set_sizes([2])
    except Exception: pass
leg = ax.get_legend()
if leg is not None:
    leg.set_bbox_to_anchor((1.0, 1.0)); leg.get_frame().set_linewidth(0.4); leg.shadow = False
    for t in leg.get_texts(): t.set_fontsize(FS)
    for h in leg.get_lines(): h.set_linewidth(1.0)
    leg._fontsize = FS; ax.add_artist(leg)
ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=c, markersize=2.6) for c in colors], labels=cat_order, loc="upper center", bbox_to_anchor=(0.5, -0.075),
          ncol=3, frameon=False, fontsize=FS, handletextpad=0.1, columnspacing=1.2, labelspacing=0.5, borderpad=0)

# ---- label placement ----
fig.canvas.draw(); RD = fig.canvas.get_renderer()
PTS = ax.transData.transform(np.vstack([np.asarray(c.get_offsets()) for c in ax.collections if len(c.get_offsets())]))   # all data points (display coordinates, px)
placed = []                                                     # bounding boxes of labels already placed
def clear(t, pad=4):
    bb = t.get_window_extent(RD); x0, y0, x1, y1 = bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad
    if ((PTS[:, 0] > x0) & (PTS[:, 0] < x1) & (PTS[:, 1] > y0) & (PTS[:, 1] < y1)).any(): return False
    return not any(x0 - 14 < b.x1 and x1 + 14 > b.x0 and y0 < b.y1 and y1 > b.y0 for b in placed)   # keep about 1 mm between labels horizontally
def settle(t, x, y, step, ylo, yhi, up_first=True):             # search vertically from y for the nearest free spot
    for k in range(0, 80):
        for sgn in ((1, -1) if up_first else (-1, 1)):
            yy = y + sgn * k * step
            if yy < ylo or yy > yhi: continue
            t.set_position((x, yy))
            if clear(t): placed.append(t.get_window_extent(RD)); return yy
            if k == 0: break
    t.set_position((x, y)); placed.append(t.get_window_extent(RD)); return y
texts = list(ax.texts); anchor = {id(t): t.get_position() for t in texts}; code = {id(t): t.get_text().strip() for t in texts}
for t in texts: t.set_text(SHOW[code[id(t)]]); t.set_fontsize(FS); t.set_va("center")   # variable code -> display label
gap = ymax * 1.10 / 30; top = ymax * 1.07; lines = []
COLCATS = {"Physical function", "Self-/physician-rated health", "Symptoms"}   # the three rightmost categories: no room right of the points, so labels form a right-aligned column
crp = [t for t in texts if code[id(t)] == "CRP"]
column = sorted([t for t in texts if cats[code[id(t)]] in COLCATS and t not in crp], key=lambda t: anchor[id(t)][1])
single = sorted([t for t in texts if t not in column and t not in crp], key=lambda t: -anchor[id(t)][1])
if column:                                                      # (1) right-aligned column first, bottom-up with spacing >= gap
    xr = min(anchor[id(t)][0] for t in column) - 14; prev = -1e9
    y_start = min(anchor[id(column[0])][1], top - gap * (len(column) - 1))
    for i, t in enumerate(column):
        t.set_ha("right"); y = max(y_start + i * gap, prev + gap)
        y = settle(t, xr, y, gap / 6, prev + gap, top + gap, up_first=True); prev = y
        lines.append(((xr + 1.5, y), anchor[id(t)]))
def best_side(t, x0, y0, ylo, yhi):                             # (2) other labels: right or left of the point, whichever needs less vertical shift
    cand = []
    for ha, dx in (("left", 9), ("right", -9)):
        t.set_ha(ha)
        for k in range(0, 60):
            hit = None
            for sgn in (1, -1):
                yy = y0 + sgn * k * gap / 6
                if ylo <= yy <= yhi:
                    t.set_position((x0 + dx, yy))
                    if clear(t): hit = yy; break
                if k == 0: break
            if hit is not None: cand.append((abs(hit - y0), 0 if ha == "left" else 1, ha, dx, hit)); break
    if not cand: t.set_ha("left"); t.set_position((x0 + 9, y0)); placed.append(t.get_window_extent(RD)); return 9, y0
    _, _, ha, dx, yy = min(cand); t.set_ha(ha); t.set_position((x0 + dx, yy)); placed.append(t.get_window_extent(RD)); return dx, yy
for t in single:
    x0, y0 = anchor[id(t)]; dx, y = best_side(t, x0, y0, 0.8, top); lines.append(((x0 + dx - (1 if dx > 0 else -1), y), (x0, y0)))
for t in crp:                                                   # (3) reference item: move up into the empty area
    x0, y0 = anchor[id(t)]; t.set_ha("left")
    y = settle(t, x0 + 10, 5.2, gap / 6, 4.6, 9.0, up_first=True); lines.append(((x0 + 9, y - 0.25), (x0, y0)))
# CLARITE puts labels at xpos + 1, so leader lines end at xb - 1 (the point center)
for (xa, ya), (xb, yb) in lines: ax.plot([xa, xb - 1], [ya, yb], color="#898781", linewidth=0.3, zorder=3, solid_capstyle="round")
# the reference CRP point sits in the bottom cloud; add a thin black ring to identify it
for t in crp: ax.scatter([anchor[id(t)][0] - 1], [anchor[id(t)][1]], s=11, facecolors="none", edgecolors="#0b0b0b", linewidths=0.5, zorder=6)
fig.savefig(f"{FIG}/Figure4A.png", dpi=300, bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{FIG}/Figure4A.pdf", bbox_inches="tight", pad_inches=0.03)   # used for Figure4.pdf below

# ---- assembly: A left, B right ----
a = Image.open(f"{FIG}/Figure4A.png").convert("RGB"); b = Image.open(f"{FIG}/Figure4B.png").convert("RGB")
PADX, PADY, LET = 12, 40, 33
H = max(a.height, b.height)
canvas = Image.new("RGB", (a.width + b.width + 3 * PADX, H + PADY), "white")
canvas.paste(a, (PADX, PADY)); canvas.paste(b, (a.width + 2 * PADX, PADY))   # panels top-aligned
d = ImageDraw.Draw(canvas); fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", LET)
d.text((6, 2), "A", fill="#0b0b0b", font=fnt); d.text((a.width + 2 * PADX - 6, 2), "B", fill="#0b0b0b", font=fnt)
canvas.save(f"{FIG}/Figure4.png", dpi=(300, 300))
print(f"-> Figure4.png {canvas.width}x{canvas.height}px ({canvas.width/300*25.4:.0f} x {canvas.height/300*25.4:.0f} mm)")

# ---- vector PDF: same layout as the PNG ----
import sys, fitz
from PIL import Image
ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
FIG = f"{ROOT}/results/figures"
k = 72 / 300                                            # pixels -> pt
a_png = Image.open(f"{FIG}/Figure4A.png"); b_png = Image.open(f"{FIG}/Figure4B.png")
aw, ah = a_png.size; bw, bh = b_png.size
PADX, PADY, LET = 12, 40, 33   # same as the PNG assembly above
W, H = (aw + bw + 3 * PADX) * k, (max(ah, bh) + PADY) * k
out = fitz.open(); page = out.new_page(width=W, height=H)
for f, x0, y0, w, h in [("Figure4A.pdf", PADX, PADY, aw, ah), ("Figure4B.pdf", aw + 2 * PADX, PADY, bw, bh)]:
    src = fitz.open(f"{FIG}/{f}")
    page.show_pdf_page(fitz.Rect(x0 * k, y0 * k, (x0 + w) * k, (y0 + h) * k), src, 0, keep_proportion=True)
font = fitz.Font(fontfile="/System/Library/Fonts/Supplemental/Arial Bold.ttf"); page.insert_font(fontname="arialbd", fontbuffer=font.buffer)
fs = LET * k; base = (2 + font.ascender * LET) * k     # PIL positions by glyph top, PDF by baseline
page.insert_text((6 * k, base), "A", fontname="arialbd", fontsize=fs, color=(11 / 255,) * 3)
page.insert_text(((aw + 2 * PADX - 6) * k, base), "B", fontname="arialbd", fontsize=fs, color=(11 / 255,) * 3)
out.save(f"{FIG}/Figure4.pdf", garbage=4, deflate=True)
print(f"-> Figure4.pdf {W:.1f} x {H:.1f} pt ({W / 72:.2f} x {H / 72:.2f} in)")
