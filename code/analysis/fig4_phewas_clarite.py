# Figure 4A:用 CLARITE 的 clarite.plot.manhattan(matplotlib)作图,再与 R 画的 B 拼成整图。
# 输入 results/figures/Figure4A.csv(fig4_phewas.R 导出)、Figure4B.png;输出 Figure4A.{png,pdf}、Figure4.png
# CLARITE 的约定原样保留(类别按名称字母序、灰白交替底色、FDR 5% 线、标签在点右侧);出图后只做排版微调(字号、坐标范围、标签防重叠、类别色 = B 的健康类别色)。
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
# 传给 CLARITE 的变量名用 NHANES 变量代码(CLARITE 在类别内按变量名排序,带排名前缀的显示标签会排到最前);出图后再换回显示标签。
df = pd.DataFrame({c: np.nan for c in result_columns + corrected_pvalue_columns}, index=pd.MultiIndex.from_arrays([A["variable"], ["KSAA"] * n], names=["Variable", "Outcome"]))
df["pvalue"] = A["p"].values; df["pvalue_fdr"] = A["fdr"].values; df["pvalue_bonferroni"] = np.minimum(1, A["p"].values * n)
df["N"] = n; df["Converged"] = True; df["Variable_type"] = "continuous"
cats = dict(zip(A["variable"], A["category"])); SHOW = dict(zip(A["variable"], A["label"]))
cat_order = sorted(A["category"].unique())                       # CLARITE 按类别名排序
colors = [A.loc[A["category"] == c, "colour"].iloc[0] for c in cat_order]
bg = ["#EBEBEB" if i % 2 == 0 else "#FFFFFF" for i in range(len(cat_order))]
LABEL_CRP = False   # True = 标注参照项 C-reactive protein(标签、引线、黑圈)
labels = list(A.loc[A["is_top"], "variable"]) + (["CRP"] if LABEL_CRP else [])

# 出图后的排版:标签不压点(与数据点及已放好的标签比对,沿纵向挪到最近空位);图右端三类的标签排成右对齐一列,引线呈扇形连到各点;
#   横轴不写类别名,用底部色块图例 + 标题 "Phenotype";纵轴每 5 一格;字号 5.2 pt(同图 2、3)。
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
FS = 5.2
# figsize:CLARITE 自带边距,按紧边界存出后约 3.85 × 4.2 in
# 不画 Bonferroni 线,只留 FDR 5%(红点线)
fig = clarite.plot.manhattan({"": df}, categories=cats, bonferroni=None, fdr=0.05, num_labeled=0, label_vars=labels,
                             figsize=(4.15, 4.95), dpi=300, colors=colors, background_colors=bg, return_figure=True)   # 使整图收在通栏 183 mm
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

# ---- 标签排布 ----
fig.canvas.draw(); RD = fig.canvas.get_renderer()
PTS = ax.transData.transform(np.vstack([np.asarray(c.get_offsets()) for c in ax.collections if len(c.get_offsets())]))   # 全部数据点(显示坐标,px)
placed = []                                                     # 已放好的标签文字框
def clear(t, pad=4):
    bb = t.get_window_extent(RD); x0, y0, x1, y1 = bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad
    if ((PTS[:, 0] > x0) & (PTS[:, 0] < x1) & (PTS[:, 1] > y0) & (PTS[:, 1] < y1)).any(): return False
    return not any(x0 - 14 < b.x1 and x1 + 14 > b.x0 and y0 < b.y1 and y1 > b.y0 for b in placed)   # 标签之间左右至少留约 1 mm,免得两条标签首尾相接
def settle(t, x, y, step, ylo, yhi, up_first=True):             # 从 y 起沿纵向找最近的空位
    for k in range(0, 80):
        for sgn in ((1, -1) if up_first else (-1, 1)):
            yy = y + sgn * k * step
            if yy < ylo or yy > yhi: continue
            t.set_position((x, yy))
            if clear(t): placed.append(t.get_window_extent(RD)); return yy
            if k == 0: break
    t.set_position((x, y)); placed.append(t.get_window_extent(RD)); return y
texts = list(ax.texts); anchor = {id(t): t.get_position() for t in texts}; code = {id(t): t.get_text().strip() for t in texts}
for t in texts: t.set_text(SHOW[code[id(t)]]); t.set_fontsize(FS); t.set_va("center")   # 变量代码 → 显示标签
gap = ymax * 1.10 / 30; top = ymax * 1.07; lines = []
COLCATS = {"Physical function", "Self-/physician-rated health", "Symptoms"}   # 图右端的三类:标签放不到点的右侧,排成右对齐的一列
crp = [t for t in texts if code[id(t)] == "CRP"]
column = sorted([t for t in texts if cats[code[id(t)]] in COLCATS and t not in crp], key=lambda t: anchor[id(t)][1])
single = sorted([t for t in texts if t not in column and t not in crp], key=lambda t: -anchor[id(t)][1])
if column:                                                      # ① 先排这一列:自下而上、间距不小于 gap、不压点;引线从标签右端呈扇形连到各点
    xr = min(anchor[id(t)][0] for t in column) - 14; prev = -1e9
    y_start = min(anchor[id(column[0])][1], top - gap * (len(column) - 1))
    for i, t in enumerate(column):
        t.set_ha("right"); y = max(y_start + i * gap, prev + gap)
        y = settle(t, xr, y, gap / 6, prev + gap, top + gap, up_first=True); prev = y
        lines.append(((xr + 1.5, y), anchor[id(t)]))
def best_side(t, x0, y0, ylo, yhi):                             # ② 其余标签:点的右侧与左侧各沿纵向找空位,取纵向位移最小的一个(相同则右侧)
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
for t in crp:                                                   # ③ 参照项:提到上方空白处
    x0, y0 = anchor[id(t)]; t.set_ha("left")
    y = settle(t, x0 + 10, 5.2, gap / 6, 4.6, 9.0, up_first=True); lines.append(((x0 + 9, y - 0.25), (x0, y0)))
# CLARITE 把标签放在 xpos + 1 处,数据点本身在 xpos:引线终点取 xb - 1,正对点心
for (xa, ya), (xb, yb) in lines: ax.plot([xa, xb - 1], [ya, yb], color="#898781", linewidth=0.3, zorder=3, solid_capstyle="round")
# 参照项 CRP 的点在底部点云里,加黑色细圈以便辨认
for t in crp: ax.scatter([anchor[id(t)][0] - 1], [anchor[id(t)][1]], s=11, facecolors="none", edgecolors="#0b0b0b", linewidths=0.5, zorder=6)
fig.savefig(f"{FIG}/Figure4A.png", dpi=300, bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{FIG}/Figure4A.pdf", bbox_inches="tight", pad_inches=0.03)   # 矢量版,供 fig4_phewas_clarite.py 拼整图 PDF

# ---- 拼图:A 左、B 右 ----
a = Image.open(f"{FIG}/Figure4A.png").convert("RGB"); b = Image.open(f"{FIG}/Figure4B.png").convert("RGB")
PADX, PADY, LET = 12, 40, 33
H = max(a.height, b.height)
canvas = Image.new("RGB", (a.width + b.width + 3 * PADX, H + PADY), "white")
canvas.paste(a, (PADX, PADY)); canvas.paste(b, (a.width + 2 * PADX, PADY))   # 两面板顶对齐
d = ImageDraw.Draw(canvas); fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", LET)
d.text((6, 2), "A", fill="#0b0b0b", font=fnt); d.text((a.width + 2 * PADX - 6, 2), "B", fill="#0b0b0b", font=fnt)
canvas.save(f"{FIG}/Figure4.png", dpi=(300, 300))
print(f"-> Figure4.png {canvas.width}x{canvas.height}px ({canvas.width/300*25.4:.0f} x {canvas.height/300*25.4:.0f} mm)")

# ---- 矢量 PDF:按上面拼 PNG 的同一版式把 A(CLARITE)与 B(R)拼成 Figure4.pdf ----
import sys, fitz
from PIL import Image
ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
FIG = f"{ROOT}/results/figures"
k = 72 / 300                                            # 像素 → pt
a_png = Image.open(f"{FIG}/Figure4A.png"); b_png = Image.open(f"{FIG}/Figure4B.png")
aw, ah = a_png.size; bw, bh = b_png.size
PADX, PADY, LET = 12, 40, 33   # 与 fig4_phewas_clarite.py 的拼图段同步
W, H = (aw + bw + 3 * PADX) * k, (max(ah, bh) + PADY) * k
out = fitz.open(); page = out.new_page(width=W, height=H)
for f, x0, y0, w, h in [("Figure4A.pdf", PADX, PADY, aw, ah), ("Figure4B.pdf", aw + 2 * PADX, PADY, bw, bh)]:
    src = fitz.open(f"{FIG}/{f}")
    page.show_pdf_page(fitz.Rect(x0 * k, y0 * k, (x0 + w) * k, (y0 + h) * k), src, 0, keep_proportion=True)
font = fitz.Font(fontfile="/System/Library/Fonts/Supplemental/Arial Bold.ttf"); page.insert_font(fontname="arialbd", fontbuffer=font.buffer)
fs = LET * k; base = (2 + font.ascender * LET) * k     # PIL 以字形上沿定位,PDF 以基线定位
page.insert_text((6 * k, base), "A", fontname="arialbd", fontsize=fs, color=(11 / 255,) * 3)
page.insert_text(((aw + 2 * PADX - 6) * k, base), "B", fontname="arialbd", fontsize=fs, color=(11 / 255,) * 3)
out.save(f"{FIG}/Figure4.pdf", garbage=4, deflate=True)
print(f"-> Figure4.pdf {W:.1f} x {H:.1f} pt ({W / 72:.2f} x {H / 72:.2f} in)")
# 跑法:$CLARITE_PY code/analysis/fig4_phewas_clarite.py .(CLARITE 2.3.6 需 Python ≤ 3.10;该环境另装 Pillow 与 PyMuPDF)
