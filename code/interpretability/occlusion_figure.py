# Fig. S2: occlusion sensitivity maps of the six example right knees (occlusion_map.py), in the layout of the Grad-CAM examples.
#   Maps are lightly smoothed for display (sigma 2 grid cells = 8 px at 448) and share one color scale (99.5th percentile).
# Input: results/analysis/saliency/{occlusion_map_examples.npz, ksaa_for_saliency.csv}, ROI images
# Output: results/figures/FigureS2.{png,pdf}
# Run: python code/interpretability/occlusion_figure.py
import os, sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import KB, OUT, roi_file, pre
import torch, torch.nn.functional as F
import torchvision.transforms.functional as TF
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
FIGD = os.path.join(KB, "results", "figures"); os.makedirs(FIGD, exist_ok=True)
INK, Q1C, Q4C = "#1d1d1d", "#9fb3c8", "#c8382c"
CM = LinearSegmentedColormap.from_list("sal", ["#000000", "#5b1a8a", "#d1335b", "#f8961e", "#fff3a0"])
mi = lambda s: s.replace("-", "−")

Z = np.load(os.path.join(OUT, "occlusion_map_examples.npz")); SEQ = list(Z["seqn"]); assert str(Z["side"]) == "R"
M = TF.gaussian_blur(torch.from_numpy(Z["abs"]).float()[:, None], kernel_size=13, sigma=2.0)[:, 0].numpy()   # smooth out the 16-pixel stride
info = {r["SEQN"]: r for r in csv.DictReader(open(os.path.join(OUT, "ksaa_for_saliency.csv")))}
pairs = [(0, 1), (2, 3), (4, 5)]                                     # (younger-appearing, older-appearing), in the order of occlusion_map.py --seqn
VMAX = float(np.percentile(M, 99.5))                                 # shared color scale

WF, SF, GX, GY, XL, YT = 160.0, 47.0, 3.0, 9.0, 9.0, 8.0; CBW = 20.0; HF = YT + 2 * SF + GY + 2.0
fig = plt.figure(figsize=((WF + CBW) / 25.4, HF / 25.4)); WT = WF + CBW
for k, (j, i) in enumerate(pairs):
    for m_, (p, col) in enumerate([(j, Q1C), (i, Q4C)]):
        s = SEQ[p]; x_, y_ = XL + k * (SF + GX), YT + m_ * (SF + GY)
        ax = fig.add_axes([x_ / WT, 1 - (y_ + SF) / HF, SF / WT, SF / HF]); u = pre(roi_file(s, "R"), "R").numpy()
        mp = F.interpolate(torch.from_numpy(M[p])[None, None], size=u.shape, mode="bilinear", align_corners=False)[0, 0].numpy()
        ax.imshow(u, cmap="gray", vmin=0, vmax=1, interpolation="bilinear")
        a = np.clip(mp / VMAX, 0, 1); rgba = CM(a); rgba[..., 3] = 0.75 * a ** 0.8; ax.imshow(rgba, interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([])
        for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.3)
        r = info[s]; ax.set_title(mi(f"{'Woman' if r['sex'] == '2' else 'Man'}, {float(r['age']):.0f} years, KSAA {float(r['gap']):+.1f} years"), fontsize=7, color=INK, pad=2.5)
        if k == 0:
            for x0, t0, h0 in [(0.04, "L", "left"), (0.96, "M", "right")]: ax.text(x0, 0.96, t0, transform=ax.transAxes, color="white", fontsize=6.5, va="top", ha=h0)
for m_, (t_, col) in enumerate([("Younger-appearing", Q1C), ("Older-appearing", Q4C)]):
    fig.text((XL - 3.5) / WT, 1 - (YT + m_ * (SF + GY) + SF / 2) / HF, t_, rotation=90, fontsize=7.5, fontweight="bold", color=INK, ha="center", va="center")
cax = fig.add_axes([(WF + 3) / WT, 1 - (YT + 2 * SF + GY - 10) / HF, 2.5 / WT, (2 * SF + GY - 20) / HF])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, VMAX), cmap=CM), cax=cax); cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.set_label("Change in predicted age\nwhen the region is blurred (years)", fontsize=6.5, color=INK, labelpad=3)
for ext in ("png", "pdf"): fig.savefig(os.path.join(FIGD, f"FigureS2.{ext}"), dpi=300)
print(f"common scale: 0 to {VMAX:.2f} years (99.5th percentile over the six maps); per-map max |Δ| " + ", ".join(f"{M[p].max():.2f}" for p in range(len(SEQ))))
print("-> results/figures/FigureS2.{png,pdf}")
