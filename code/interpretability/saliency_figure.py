# Interpretability figures and quantification (Fig. S2); inputs come from saliency_images.py and saliency.py --heldout [--occlude].
#   Default layout: 3 same-sex example pairs from the internal validation set, age difference <= 1 year (KSAA near the cohort's 90th / 10th percentile), right-knee image + Grad-CAM.
#   Quantitative output (saliency_summary.txt / saliency_regions.csv): correlation of the age and KSAA difference maps (null = KSAA permuted within age strata),
#   enrichment ratios of the 6 anatomical regions and MAE increase after occlusion (bootstrap CI).
# Run (from the repository root): python code/interpretability/saliency_figure.py
import os, sys, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec, torch, torch.nn.functional as F
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import KB, OUT, R, REG, roi_file, pre, regions, otsu
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.4, "pdf.fonttype": 42, "svg.fonttype": "none"})
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"; Q1C, Q4C = "#9DB4C8", "#C8372D"      # quartile colors as in the main figures
import argparse; ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="_heldout"); ap.add_argument("--out", default="")
ap.add_argument("--layout", default="examples", choices=["examples", "full"]); A = ap.parse_args()
# --layout full: six-panel version (not in the manuscript); --tag _n40 --out <dir>: layout debugging only
TEST = A.tag != "_heldout"; FIGD = A.out or os.path.join(KB, "results", "figures"); SUM = A.out or OUT
rng = np.random.default_rng(20260922); NB = 100 if TEST else 500
mi = lambda s: s.replace("-", "−")
f1 = lambda x: mi(f"{x:.1f}"); f2 = lambda x: mi(f"{x:.2f}")

# ---------- whole cohort: mean and difference maps ----------
rows = list(csv.DictReader(open(os.path.join(OUT, "ksaa_for_saliency.csv"))))
assert [r["SEQN"] for r in rows] == open(os.path.join(OUT, "img_all_seqn.txt")).read().split()
I = np.load(os.path.join(OUT, "img_all.npy")).astype(np.float32).mean(1)                 # mean of both knees [N,R,R]
age = np.array([float(r["age"]) for r in rows]); gap = np.array([float(r["gap"]) for r in rows]); sex = np.array([int(r["sex"]) for r in rows])
qc = np.quantile(gap, [.25, .5, .75]); Q = np.digitize(gap, qc)                           # same as R quantile (type 7)
band = np.minimum(((age - 60) // 5).astype(int), 5)                                     # 60–64 … 85–90
yng, old = age < 65, age >= 80
def std_means(a, b, idx=None):
    """Means of groups a and b within 5-year age strata, weighted by the combined stratum size (age-standardized). idx: bootstrap sample indices."""
    ii = np.arange(len(age)) if idx is None else idx; ma = mb = 0.0; wt = 0
    for s in range(6):
        ia = ii[(band[ii] == s) & a[ii]]; ib = ii[(band[ii] == s) & b[ii]]
        if len(ia) == 0 or len(ib) == 0: continue
        w = len(ia) + len(ib); ma = ma + w * I[ia].mean(0); mb = mb + w * I[ib].mean(0); wt += w
    return ma / wt, mb / wt
MY, MO = I[yng].mean(0), I[old].mean(0); DA = MO - MY
K1, K4 = std_means(Q == 0, Q == 3); DK = K4 - K1
C = slice(8, R - 8); pr = lambda a, b: np.corrcoef(a[C, C].ravel(), b[C, C].ravel())[0, 1]
r_map = pr(DA, DK); perm = []                                                            # null: KSAA permuted within 5-year age strata (the bootstrap would systematically underestimate r)
for _ in range(NB):
    gp = gap.copy()
    for s_ in range(6): k_ = np.where(band == s_)[0]; gp[k_] = gp[rng.permutation(k_)]
    Qp = np.digitize(gp, qc); a_, b_ = std_means(Qp == 3, Qp == 0); perm.append(pr(DA, a_ - b_))
r_null = np.percentile(np.abs(perm), 95)
slope = np.polyfit(DA[C, C].ravel(), DK[C, C].ravel(), 1)[0]                        # magnitude: slope of the KSAA difference map on the age difference map (1 = equally large)
# regional summary of the difference maps (same anatomical regions, joint line of the cohort mean image)
REGM, jm = regions(I.mean(0))
def reg_mean(D): return [float(D[m].mean()) for m in REGM]

# ---------- internal validation set: Grad-CAM ----------
ho = list(csv.DictReader(open(os.path.join(OUT, f"saliency_pred{A.tag}.csv"))))
CAM = np.load(os.path.join(OUT, f"cam{A.tag}.npy")).astype(np.float32); IH = np.load(os.path.join(OUT, f"img{A.tag}.npy")).astype(np.float32)
n_ho = len(ho); assert CAM.shape[0] == n_ho and (TEST or n_ho == 487)
aH = np.array([float(r["age"]) for r in ho]); gH = np.array([float(r["gap"]) for r in ho]); pH = np.array([float(r["pred_full_ens"]) for r in ho])
oH = np.array([float(r["pred_oof"]) for r in ho]); sH = {r["SEQN"]: int(x["sex"]) for x in rows for r in [x]}
share = np.zeros((n_ho, 2, 6)); area = np.zeros((n_ho, 2, 6)); bone_sh = np.zeros((n_ho, 2)); bone_ar = np.zeros((n_ho, 2)); JL = np.zeros((n_ho, 2), int)
for i in range(n_ho):
    for v in range(2):
        M, JL[i, v] = regions(IH[i, v]); c = CAM[i, v] / CAM[i, v].sum()
        share[i, v] = [c[m].sum() for m in M]; area[i, v] = [m.mean() for m in M]
        b = IH[i, v] > otsu(IH[i, v]); b &= ~(M[2] | M[3]); bone_sh[i, v] = c[b].sum(); bone_ar[i, v] = b.mean()   # bone = above Otsu and outside the joint band
sh, ar = share.mean(1), area.mean(1)                                                     # per participant, mean of both knees
def boot_ratio(num, den):
    est = num.mean(0) / den.mean(0); b = []
    for _ in range(NB):
        ii = rng.integers(0, len(num), len(num)); b.append(num[ii].mean(0) / den[ii].mean(0))
    return est, np.percentile(b, [2.5, 97.5], axis=0)
enr, enr_ci = boot_ratio(sh, ar)
bone_enr, bone_ci = boot_ratio(bone_sh.mean(1)[:, None], bone_ar.mean(1)[:, None])
MCAM = CAM.mean((0, 1))

# ---------- occlusion test ----------
labs = ["None"] + REG
if TEST:                                                                                 # debugging: unoccluded prediction + noise instead
    E = np.repeat(pH[:, None], 7, 1) + rng.normal(0, 1, (n_ho, 7)) * (np.arange(7) > 0); occ = [{"region": labs[g], "SEQN": ho[i]["SEQN"]} for i in range(n_ho) for g in range(7)]
else:
    occ = list(csv.DictReader(open(os.path.join(OUT, f"occlusion{A.tag}.csv")))); E = np.array([float(r["pred_full_ens"]) for r in occ]).reshape(n_ho, 7)
assert [r["region"] for r in occ[:7]] == labs and [occ[7 * i]["SEQN"] for i in range(n_ho)] == [r["SEQN"] for r in ho]
assert np.allclose(E[:, 0], pH, atol=0.05)                                                # unoccluded = forward pass of the Grad-CAM run
AE = np.abs(E - aH[:, None]); dMAE = AE[:, 1:] - AE[:, [0]]; shift = E[:, 1:] - E[:, [0]]
def boot_mean(x):
    b = [x[rng.integers(0, len(x), len(x))].mean(0) for _ in range(NB)]; return x.mean(0), np.percentile(b, [2.5, 97.5], axis=0)
dmae, dmae_ci = boot_mean(dMAE); sft, sft_ci = boot_mean(shift)
r_occ = [np.corrcoef(E[:, g], aH)[0, 1] for g in range(7)]

# ---------- example pairs ----------
p10, p90 = np.percentile(gap, [10, 90]); seqH = [r["SEQN"] for r in ho]; sexH = np.array([sH[s] for s in seqH])
pairs = []
for lo, hi in [(60, 70), (70, 80), (80, 91)]:
    cand = [i for i in np.argsort(np.abs(gH - p90)) if lo <= aH[i] < hi]
    for i in cand:
        js = [j for j in range(n_ho) if sexH[j] == sexH[i] and abs(aH[j] - aH[i]) <= 1 and j != i]
        if js:
            j = min(js, key=lambda j: abs(gH[j] - p10))
            if gH[j] < 0: pairs.append((j, i)); break

# ---------- numeric output ----------
L = []
L.append(f"Mean images: age 60–64 n={yng.sum()}, age >=80 n={old.sum()}; KSAA Q1 n={(Q == 0).sum()} (mean KSAA {gap[Q == 0].mean():.2f}, mean age {age[Q == 0].mean():.1f}), Q4 n={(Q == 3).sum()} (mean KSAA {gap[Q == 3].mean():.2f}, mean age {age[Q == 3].mean():.1f})")
L.append(f"Pixelwise r between age-difference and age-standardised KSAA-difference maps: {r_map:.3f} (95th percentile of |r| when KSAA is permuted within 5-year age bands: {r_null:.3f}; {NB} permutations)")
L.append(f"Magnitude: slope of the KSAA-difference map on the age-difference map {slope:.3f}; mean age {age[old].mean():.1f} (>=80) vs {age[yng].mean():.1f} (60-64)")
L.append("Regional mean of difference maps (normalised intensity; + = brighter in older / higher KSAA):")
for g, (a_, k_) in enumerate(zip(reg_mean(DA), reg_mean(DK))): L.append(f"  {REG[g]:<22} age {a_:+.4f}  KSAA {k_:+.4f}")
L.append(f"Joint line (mean image) row {jm} of {R}; per-knee joint line median {np.median(JL):.0f} (IQR {np.percentile(JL, 25):.0f}-{np.percentile(JL, 75):.0f})")
L.append(f"Held-out n={n_ho}: mean age {aH.mean():.1f}; full-model ensemble vs age MAE {np.mean(np.abs(pH - aH)):.2f}, r {np.corrcoef(pH, aH)[0, 1]:.3f}; vs out-of-fold prediction r {np.corrcoef(pH, oH)[0, 1]:.3f}; out-of-fold MAE in same people {np.mean(np.abs(oH - aH)):.2f}")
L.append("Grad-CAM saliency by region (share of saliency / share of area = enrichment, 95% bootstrap CI):")
for g in range(6): L.append(f"  {REG[g]:<22} saliency {sh[:, g].mean():.3f}  area {ar[:, g].mean():.3f}  enrichment {enr[g]:.2f} ({enr_ci[0, g]:.2f}-{enr_ci[1, g]:.2f})")
L.append(f"  medial total saliency {sh[:, 1::2].sum(1).mean():.3f} vs area {ar[:, 1::2].sum(1).mean():.3f}; joint band total saliency {sh[:, 2:4].sum(1).mean():.3f} vs area {ar[:, 2:4].sum(1).mean():.3f}")
L.append(f"  bone (Otsu, outside joint band) saliency {bone_sh.mean():.3f} area {bone_ar.mean():.3f} enrichment {bone_enr[0]:.2f} ({bone_ci[0, 0]:.2f}-{bone_ci[1, 0]:.2f})")
L.append("Occlusion (Gaussian blur of one region in both knees): increase in MAE, mean shift in predicted age, r with age:")
L.append(f"  None                   MAE {AE[:, 0].mean():.2f}  r {r_occ[0]:.3f}")
for g in range(6): L.append(f"  {REG[g]:<22} dMAE {dmae[g]:+.2f} ({dmae_ci[0, g]:+.2f} to {dmae_ci[1, g]:+.2f})  shift {sft[g]:+.2f} ({sft_ci[0, g]:+.2f} to {sft_ci[1, g]:+.2f})  r {r_occ[g + 1]:.3f}  area {ar[:, g].mean():.3f}")
L.append(f"Exemplars (KSAA P10 {p10:.2f}, P90 {p90:.2f}):")
for j, i in pairs: L.append(f"  sex {sexH[i]}: low SEQN {seqH[j]} age {aH[j]:.0f} KSAA {gH[j]:+.2f} | high SEQN {seqH[i]} age {aH[i]:.0f} KSAA {gH[i]:+.2f}")
open(os.path.join(SUM, "saliency_summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
with open(os.path.join(SUM, "saliency_regions.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["region", "area_share", "saliency_share", "enrichment", "enr_lo", "enr_hi", "occl_dMAE", "dMAE_lo", "dMAE_hi", "occl_shift", "shift_lo", "shift_hi", "diff_age", "diff_ksaa"])
    for g in range(6): w.writerow([REG[g], f"{ar[:, g].mean():.4f}", f"{sh[:, g].mean():.4f}", f"{enr[g]:.3f}", f"{enr_ci[0, g]:.3f}", f"{enr_ci[1, g]:.3f}", f"{dmae[g]:.3f}", f"{dmae_ci[0, g]:.3f}", f"{dmae_ci[1, g]:.3f}", f"{sft[g]:.3f}", f"{sft_ci[0, g]:.3f}", f"{sft_ci[1, g]:.3f}", f"{reg_mean(DA)[g]:.5f}", f"{reg_mean(DK)[g]:.5f}"])

# ---------- default layout: paired examples ----------
CM = LinearSegmentedColormap.from_list("sal", ["#000000", "#5b1a8a", "#d1335b", "#f8961e", "#fff3a0"])
if A.layout == "examples":   # 2 rows x 3 columns, one pair per column (top = younger-appearing, bottom = older-appearing), 160 mm wide
    WF, SF, GX, GY, XL, YT = 160.0, 47.0, 3.0, 9.0, 9.0, 8.0; HF = YT + 2 * SF + GY + 2.0
    fig = plt.figure(figsize=(WF / 25.4, HF / 25.4))
    for k, (j, i) in enumerate(pairs):
        for m_, (p, col) in enumerate([(j, Q1C), (i, Q4C)]):
            x_, y_ = XL + k * (SF + GX), YT + m_ * (SF + GY)
            ax = fig.add_axes([x_ / WF, 1 - (y_ + SF) / HF, SF / WF, SF / HF]); u = pre(roi_file(seqH[p], "R"), "R").numpy()
            cm_ = F.interpolate(torch.from_numpy(CAM[p, 1])[None, None], size=u.shape, mode="bilinear", align_corners=False)[0, 0].numpy()
            ax.imshow(u, cmap="gray", vmin=0, vmax=1, interpolation="bilinear")
            a = np.clip(cm_ / np.percentile(cm_, 99.5), 0, 1); rgba = CM(a); rgba[..., 3] = 0.75 * a ** 0.8; ax.imshow(rgba, interpolation="bilinear")
            ax.set_xticks([]); ax.set_yticks([])
            for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.3)
            ax.set_title(mi(f"{'Woman' if sexH[p] == 2 else 'Man'}, {aH[p]:.0f} years, KSAA {gH[p]:+.1f} years"), fontsize=7, color=INK, pad=2.5)
            if k == 0:
                for x0, t0, h0 in [(0.04, "L", "left"), (0.96, "M", "right")]: ax.text(x0, 0.96, t0, transform=ax.transAxes, color="white", fontsize=6.5, va="top", ha=h0)
    for m_, (t_, col) in enumerate([("Younger-appearing", Q1C), ("Older-appearing", Q4C)]):
        fig.text((XL - 3.5) / WF, 1 - (YT + m_ * (SF + GY) + SF / 2) / HF, t_, rotation=90, fontsize=7.5, fontweight="bold", color=INK, ha="center", va="center")
    for ext in ("png", "pdf"): fig.savefig(os.path.join(FIGD, f"gradcam_examples.{ext}"), dpi=300)
    print("figure written (examples only)"); sys.exit(0)

# ---------- --layout full six-panel version (180 x 146 mm, positioned in mm) ----------
CM = LinearSegmentedColormap.from_list("sal", ["#000000", "#5b1a8a", "#d1335b", "#f8961e", "#fff3a0"])
WF, HF = 180.0, 146.0
fig = plt.figure(figsize=(WF / 25.4, HF / 25.4))
def axmm(x, y, w, h): return fig.add_axes([x / WF, 1 - (y + h) / HF, w / WF, h / HF])      # (x, y) = top-left corner, measured from the page's top-left
def letter(x, y, t): fig.text(x / WF, 1 - y / HF, t, fontsize=9, fontweight="bold", va="top", ha="left", color=INK)
def title(ax, t): ax.set_title(t, fontsize=6.5, color=INK, pad=2, linespacing=1.1)
def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_linewidth(0.3); s_.set_edgecolor(MUTED)
def gray(ax, img):
    ax.imshow(img, cmap="gray", vmin=0, vmax=1, interpolation="bilinear"); bare(ax)
def overlay(ax, img, cam, vmax):
    ax.imshow(img, cmap="gray", vmin=0, vmax=1, interpolation="bilinear")
    a = np.clip(cam / vmax, 0, 1); rgba = CM(a); rgba[..., 3] = 0.75 * a ** 0.8
    ax.imshow(rgba, interpolation="bilinear"); bare(ax)
def sides(ax):
    for x_, t_, h_ in [(0.04, "L", "left"), (0.96, "M", "right")]: ax.text(x_, 0.96, t_, transform=ax.transAxes, color="white", fontsize=5.5, va="top", ha=h_)
SZ, G = 32.0, 1.5; X0 = 5.0; XS = [X0, X0 + SZ + G, X0 + 2 * (SZ + G)]; YA, YB = 9.0, 53.0
vd = float(max(np.abs(DA[C, C]).max(), np.abs(DK[C, C]).max()))
letter(1, 1.5, "A")
for k_, (img, t_) in enumerate([(MY, f"Age 60\u201364 years\n(n = {yng.sum():,})"), (MO, f"Age \u226580 years\n(n = {old.sum():,})")]):
    ax = axmm(XS[k_], YA, SZ, SZ); gray(ax, img); title(ax, t_)
    if k_ == 0: sides(ax)
ax = axmm(XS[2], YA, SZ, SZ); im = ax.imshow(DA, cmap="RdBu_r", vmin=-vd, vmax=vd, interpolation="bilinear"); bare(ax); title(ax, "Difference\n(\u226580 minus 60\u201364 years)")
letter(1, 45.5, "B")
for k_, (img, t_) in enumerate([(K1, f"KSAA Q1, youngest-appearing\n(n = {(Q == 0).sum():,})"), (K4, f"KSAA Q4, oldest-appearing\n(n = {(Q == 3).sum():,})")]):
    ax = axmm(XS[k_], YB, SZ, SZ); gray(ax, img); title(ax, t_)
ax = axmm(XS[2], YB, SZ, SZ); ax.imshow(DK, cmap="RdBu_r", vmin=-vd, vmax=vd, interpolation="bilinear"); bare(ax); title(ax, "Difference (Q4 minus Q1),\nage-standardised")
cax = axmm(XS[2] + 3, YB + SZ + 3.2, SZ - 6, 1.6); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.outline.set_linewidth(0.3)
cb.ax.tick_params(labelsize=5.5, width=0.3, length=1.5, pad=1); cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: mi(f"{x:.2f}")))
cb.set_label("Difference in normalised intensity\n(red, brighter in older knees)", fontsize=5.5, color=INK, labelpad=1.5, linespacing=1.1)
fig.text((XS[0] + XS[1] + SZ) / 2 / WF, 1 - (YB + SZ + 3.5) / HF, mi(f"Pixelwise correlation of the two difference maps: r = {r_map:.2f}"), fontsize=6, color=INK, ha="center", va="top")
# C mean Grad-CAM
XR = 116.0; letter(XR - 4, 1.5, "C")
ax = axmm(XR + 12, YA, 40, 40); overlay(ax, I.mean(0), MCAM, np.percentile(MCAM, 99.5)); sides(ax)
for jl in (jm - 8 - 0.5, jm + 8 + 0.5): ax.axhline(jl, color="white", lw=0.4, ls=(0, (2, 2)))
ax.axvline(R / 2 - 0.5, color="white", lw=0.4, ls=(0, (2, 2))); title(ax, f"Mean Grad-CAM, internal validation set\n(n = {n_ho})")
for t_, y_ in [("Femur", (jm - 8) / 2), ("Joint band", jm), ("Tibia", (jm + 8 + R) / 2)]: ax.text(-2, y_, t_, fontsize=5.5, color=INK, ha="right", va="center")
# D enrichment ratios, E occlusion
lab = ["Femur, lateral", "Femur, medial", "Joint band, lateral", "Joint band, medial", "Tibia, lateral", "Tibia, medial"]; yy = np.arange(6)[::-1]
YD, HD, WD = 66.0, 27.0, 18.0; XD = XR + 17
def dots(ax, est, lo, hi, ref, col, xl, t_):
    ax.axvline(ref, color=MUTED, lw=0.4, ls=(0, (2, 2))); ax.hlines(yy, lo, hi, color=col, lw=0.7); ax.plot(est, yy, "o", ms=2.6, color=col, mec="none")
    ax.set_ylim(-0.6, 5.6); ax.set_yticks(yy); ax.set_yticklabels([]); ax.tick_params(labelsize=5.5, width=0.3, length=1.5, pad=1)
    for k_ in ("top", "right"): ax.spines[k_].set_visible(False)
    for k_ in ("left", "bottom"): ax.spines[k_].set_linewidth(0.3)
    ax.set_xlabel(xl, fontsize=5.5, color=INK, labelpad=1.5, linespacing=1.1); title(ax, t_)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: mi(f"{x:g}"))); ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
axd = axmm(XD, YD, WD, HD); dots(axd, enr, enr_ci[0], enr_ci[1], 1, "#2a78d6", "Share of saliency /\nshare of area", "Saliency"); axd.set_xlim(0.78, 1.28); axd.set_xticks([0.8, 1.0, 1.2])
axd.set_yticklabels(lab, fontsize=5.5); letter(XR - 4, YD - 7, "D")
axe = axmm(XD + WD + 6, YD, WD, HD); dots(axe, dmae, dmae_ci[0], dmae_ci[1], 0, "#eb6834", "Increase in MAE\nwhen blurred (years)", "Occlusion")
letter(XD + WD + 2, YD - 7, "E")
# F example pairs
YF = 112.0; SF, GP, GQ = 26.5, 1.2, 4.8; letter(1, YF - 10, "F")
vx = [X0 + k * (2 * SF + GP + GQ) + m * (SF + GP) for k in range(3) for m in range(2)]
for k, (j, i) in enumerate(pairs):
    for m_, (p, col) in enumerate([(j, Q1C), (i, Q4C)]):
        ax = axmm(vx[2 * k + m_], YF, SF, SF); u = pre(roi_file(seqH[p], "R"), "R").numpy()
        cm_ = F.interpolate(torch.from_numpy(CAM[p, 1])[None, None], size=u.shape, mode="bilinear", align_corners=False)[0, 0].numpy()
        overlay(ax, u, cm_, np.percentile(cm_, 99.5))
        for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.3)
        title(ax, mi(f"{'Woman' if sexH[p] == 2 else 'Man'}, {aH[p]:.0f} years\nKSAA {gH[p]:+.1f} years"))
for ext in ("png", "pdf"): fig.savefig(os.path.join(FIGD, f"saliency_full.{ext}"), dpi=300)
print("figure written")
