# Corner occlusion test: Grad-CAM weight was concentrated in the two upper image corners, so check whether predictions depend on them.
#   Same models, preprocessing, blur (sigma 16 px at 448) and ensemble weights as saliency.py --occlude; internal validation set (487).
#   Conditions (both knees, 64 x 64 px): none; each corner; both upper corners; all four; control = left and right mid-edge patches.
# Output: results/analysis/saliency/occlusion_corners_heldout.csv, corner_occlusion_summary.txt (change in MAE and mean shift,
#   participant-level bootstrap 95% CI)
# Run: python code/interpretability/corner_occlusion.py --bs 1
import os, sys, csv, time, argparse, numpy as np
os.environ["N_VIEWS"] = "2"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import OUT, MODELS, SNAP, S, roi_file, pre
sys.path.insert(0, SNAP)
import torch, torch.nn.functional as F, torchvision.models as tvm
import torchvision.transforms.functional as TF
for _a in ["resnet18", "resnet50", "efficientnet_b0", "convnext_tiny", "densenet121"]:   # weights come from the checkpoints; do not download ImageNet weights
    _f = getattr(tvm, _a); setattr(tvm, _a, (lambda f: (lambda *a, **k: f(weights=None)))(_f)); tvm.__dict__[_a] = getattr(tvm, _a)
from mammo_age_model import Mammo_AGE

ap = argparse.ArgumentParser(); ap.add_argument("--bs", type=int, default=8); ap.add_argument("--n", type=int, default=0); A = ap.parse_args()
DEV = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
ARCH = {"resnet18": "resnet18", "resnet50": "resnet50", "efficientnetb0": "efficientnet_b0", "convnexttiny": "convnext_tiny", "densenet121": "densenet121"}
rows = list(csv.DictReader(open(os.path.join(OUT, "ksaa_for_saliency.csv"))))
HO = set(open(os.path.join(OUT, "heldout_seqn.txt")).read().split()); rows = [r for r in rows if r["SEQN"] in HO]
if A.n: rows = rows[:A.n]

K = S // 7                                                           # 64 pixels
def box(r0, c0): m = torch.zeros(S, S, dtype=torch.bool); m[r0:r0 + K, c0:c0 + K] = True; return m
TL, TR, BL, BR = box(0, 0), box(0, S - K), box(S - K, 0), box(S - K, S - K)
ML, MR = box(S // 2 - K // 2, 0), box(S // 2 - K // 2, S - K)
COND = {"None": None, "Top-left corner": TL, "Top-right corner": TR, "Bottom-left corner": BL, "Bottom-right corner": BR,
        "Both top corners": TL | TR, "All four corners": TL | TR | BL | BR, "Control: mid-edge left + right": ML | MR}

models, wts = {}, {}
for k, arch in ARCH.items():
    ck = torch.load(os.path.join(MODELS, f"m5_full_{k}_full.pt"), map_location="cpu", weights_only=False)
    a = ck["args"]; assert a["norm"] == "pct" and a["fusion"] == "tf" and a["views"] == "both" and not a["no_mirror"] and ck["views"] == ["L", "R"] and ck["min_age"] == 0 and ck["OFF"] == 0.0
    m = Mammo_AGE(arch=arch, num_output=a["num_output"], nblock=a["nblock"], hidden_size=a["hidden"]); m.load_state_dict(ck["state_dict"], strict=True)
    m.eval().to(DEV); models[k] = (m, ck["par"]); wts[k] = 1.0 / ck["inner_val_mae"]
W = {k: v / sum(wts.values()) for k, v in wts.items()}
print(f"{len(rows)} participants, device {DEV}, ensemble weights", {k: round(v, 3) for k, v in W.items()}, flush=True)

def variants(r):
    u = torch.stack([pre(roi_file(r["SEQN"], sd), sd) for sd in "LR"]).float().to(DEV)     # [2,S,S], L already mirrored; float32 (a 97 x 97 blur in float64 is very slow on CPU)
    bl = TF.gaussian_blur(u[:, None], kernel_size=97, sigma=16.0)[:, 0]
    return torch.stack([u if M is None else torch.where(M.to(DEV), bl, u) for M in COND.values()])  # [NV,2,S,S]

NV = len(COND); N = len(rows); P = {k: np.zeros((N, NV)) for k in ARCH}; t0 = time.time()
AR = torch.arange(100, dtype=torch.float32, device=DEV)
for b0 in range(0, N, A.bs):
    batch = rows[b0:b0 + A.bs]; xu = torch.cat([variants(r) for r in batch])             # [B*NV,2,S,S]
    for k in ARCH:
        m, par = models[k]; x = ((xu - par["mu"]) / par["sd"])[:, :, None].float().to(DEV)
        with torch.no_grad():
            _, loc, *_ = m(x, use_sto=False)
            p = torch.stack([(F.softmax(l, -1) * AR).sum(-1) for l in loc], 1).mean(1).cpu().numpy()
        P[k][b0:b0 + len(batch)] = p.reshape(len(batch), NV)
    if (b0 // A.bs) % 10 == 0 or b0 < 3: print(f"{b0 + len(batch)}/{N}  {time.time() - t0:.0f}s", flush=True)
E = sum(P[k] * W[k] for k in ARCH)                                                         # [N,NV] ensemble predicted age
age = np.array([float(r["age"]) for r in rows])

with open(os.path.join(OUT, "occlusion_corners_heldout.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["SEQN", "age", "gap", "condition"] + [f"pred_{k}" for k in ARCH] + ["pred_full_ens"])
    for i, r in enumerate(rows):
        for g, c in enumerate(COND): w.writerow([r["SEQN"], r["age"], r["gap"], c] + [f"{P[k][i, g]:.4f}" for k in ARCH] + [f"{E[i, g]:.4f}"])

rng = np.random.default_rng(20261016); B = 2000; idx = rng.integers(0, N, (B, N))
ae = np.abs(E - age[:, None]); lines = []
for g, c in enumerate(COND):
    dm = ae[:, g] - ae[:, 0]; sh = E[:, g] - E[:, 0]
    dm_ci = np.percentile(dm[idx].mean(1), [2.5, 97.5]); sh_ci = np.percentile(sh[idx].mean(1), [2.5, 97.5])
    lines.append(f"{c:<32} MAE {ae[:, g].mean():.2f}  dMAE {dm.mean():+.3f} ({dm_ci[0]:+.3f} to {dm_ci[1]:+.3f})  "
                 f"shift {sh.mean():+.3f} ({sh_ci[0]:+.3f} to {sh_ci[1]:+.3f})  mean |shift| {np.abs(sh).mean():.3f}  r {np.corrcoef(E[:, g], age)[0, 1]:.3f}")
open(os.path.join(OUT, "corner_occlusion_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines)); print(f"done {time.time() - t0:.0f}s", flush=True)
