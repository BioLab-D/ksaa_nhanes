import os, sys, csv, time, argparse, numpy as np
os.environ["N_VIEWS"] = "2"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import OUT, MODELS, SNAP, S, R, roi_file, pre
sys.path.insert(0, SNAP)
import torch, torch.nn.functional as F, torchvision.models as tvm
import torchvision.transforms.functional as TF
# weights come from the checkpoints; do not download ImageNet weights
for _a in ["resnet18", "resnet50", "efficientnet_b0", "convnext_tiny", "densenet121"]:
    _f = getattr(tvm, _a)
    setattr(tvm, _a, (lambda f: (lambda *a, **k: f(weights=None)))(_f))
    tvm.__dict__[_a] = getattr(tvm, _a)
from mammo_age_model import Mammo_AGE

ap = argparse.ArgumentParser()
ap.add_argument("--seqn", required=True, help="comma-separated; 'heldout' = the 487 internal-validation participants")
ap.add_argument("--side", default="R", choices=["L", "R"])
ap.add_argument("--patch", type=int, default=64)
ap.add_argument("--stride", type=int, default=16)
ap.add_argument("--bs", type=int, default=64)
ap.add_argument("--tag", default="")
A = ap.parse_args()
DEV = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
ARCH = {"resnet18": "resnet18", "resnet50": "resnet50", "efficientnetb0": "efficientnet_b0",
        "convnexttiny": "convnext_tiny", "densenet121": "densenet121"}
SEQ = open(os.path.join(OUT, "heldout_seqn.txt")).read().split() if A.seqn == "heldout" else A.seqn.split(",")

models, wts = {}, {}
for k, arch in ARCH.items():
    ck = torch.load(os.path.join(MODELS, f"full_{k}.pt"), map_location="cpu", weights_only=False)
    a = ck["args"]
    assert a["norm"] == "pct" and a["fusion"] == "tf" and a["views"] == "both" and not a["no_mirror"] and ck["views"] == ["L", "R"]
    m = Mammo_AGE(arch=arch, num_output=a["num_output"], nblock=a["nblock"], hidden_size=a["hidden"])
    m.load_state_dict(ck["state_dict"], strict=True)
    m.eval().to(DEV)
    models[k] = (m, ck["par"])
    wts[k] = 1.0 / ck["inner_val_mae"]
W = {k: v / sum(wts.values()) for k, v in wts.items()}
V = ["L", "R"].index(A.side); P0 = list(range(0, S - A.patch + 1, A.stride))
POS = [(r0, c0) for r0 in P0 for c0 in P0]
print(f"{len(SEQ)} participants, side {A.side}, {len(POS)} patches of {A.patch}px (stride {A.stride}), device {DEV}", flush=True)
AR = torch.arange(100, dtype=torch.float32, device=DEV)

def predict(xu):                                                     # xu: [B,2,S,S] unit gray -> ensemble predicted age [B]
    out = 0
    for k, (m, par) in models.items():
        x = ((xu - par["mu"]) / par["sd"])[:, :, None]
        with torch.no_grad():
            _, loc, *_ = m(x, use_sto=False)
            out = out + W[k] * torch.stack([(F.softmax(l, -1) * AR).sum(-1) for l in loc], 1).mean(1)
    return out

ABS = np.zeros((len(SEQ), R, R), np.float32)
SGN = np.zeros((len(SEQ), R, R), np.float32)
BASE = np.zeros(len(SEQ))
t0 = time.time()
for n, s in enumerate(SEQ):
    u = torch.stack([pre(roi_file(s, sd), sd) for sd in "LR"]).float().to(DEV)          # [2,S,S], L already mirrored
    bl = TF.gaussian_blur(u[V][None, None], kernel_size=97, sigma=16.0)[0, 0]
    base = float(predict(u[None])[0]); BASE[n] = base
    acc_a = torch.zeros(S, S, device=DEV)
    acc_s = torch.zeros(S, S, device=DEV)
    cnt = torch.zeros(S, S, device=DEV)
    for b0 in range(0, len(POS), A.bs):
        pos = POS[b0:b0 + A.bs]; xu = u[None].repeat(len(pos), 1, 1, 1)
        for i, (r0, c0) in enumerate(pos):
            xu[i, V, r0:r0 + A.patch, c0:c0 + A.patch] = bl[r0:r0 + A.patch, c0:c0 + A.patch]
        d = predict(xu) - base
        for i, (r0, c0) in enumerate(pos):
            acc_a[r0:r0 + A.patch, c0:c0 + A.patch] += d[i].abs()
            acc_s[r0:r0 + A.patch, c0:c0 + A.patch] += d[i]
            cnt[r0:r0 + A.patch, c0:c0 + A.patch] += 1
    pool = lambda t: F.adaptive_avg_pool2d((t / cnt)[None, None], R)[0, 0].cpu().numpy()
    ABS[n], SGN[n] = pool(acc_a), pool(acc_s)
    print(f"{n + 1}/{len(SEQ)} SEQN {s}: baseline {base:.2f}, max |change| {ABS[n].max():.3f} years, {time.time() - t0:.0f}s", flush=True)
np.savez(os.path.join(OUT, f"occlusion_map{A.tag}.npz"), seqn=np.array(SEQ), side=A.side, patch=A.patch,
         stride=A.stride, abs=ABS, signed=SGN, base=BASE)
print(f"-> occlusion_map{A.tag}.npz  done {time.time() - t0:.0f}s", flush=True)
