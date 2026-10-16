# 解释性分析 · 第 2 步:用 5 个全数据最终模型对预测年龄做 Grad-CAM(只作可视化;KSAA 仍是交叉拟合 OOF 值)。
#   预处理同 model/knee_v2.py(pct 归一化、左膝镜像、双线性缩放到 448、存档 mu/sd 标准化);读出 = 两膝局部头期望年龄均值,确定性前向。
#   Grad-CAM 取各骨干最后一个特征图,每膝归一化为和为 1,5 个骨干按内部验证 MAE 倒数加权平均,存 112 × 112。
#   --occlude:6 个解剖区逐一高斯模糊(σ = 16 像素 @448)后重新预测;模糊比常数填充少引入人为边缘。
# 全数据模型在训练过的人身上几乎记住了年龄,故正式结果只用 --heldout(487 人,未参与训练)。
# 跑法(仓库根目录下):python code/interpretability/saliency.py --heldout --bs 8 [--occlude];结果写到 results/analysis/saliency/。
import os, sys, csv, time, argparse, numpy as np
os.environ["N_VIEWS"] = "2"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import OUT, MODELS, SNAP, S, R, REG, roi_file, pre, pool, regions
sys.path.insert(0, SNAP)
import torch, torch.nn.functional as F, torchvision.models as tvm
import torchvision.transforms.functional as TF
for _a in ["resnet18", "resnet50", "efficientnet_b0", "convnext_tiny", "densenet121"]:   # 权重来自存档,不下载 ImageNet 权重
    _f = getattr(tvm, _a); setattr(tvm, _a, (lambda f: (lambda *a, **k: f(weights=None)))(_f)); tvm.__dict__[_a] = getattr(tvm, _a)
from mammo_age_model import Mammo_AGE

ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=0); ap.add_argument("--bs", type=int, default=4)
ap.add_argument("--heldout", action="store_true", help="只算全数据模型的内部验证集(487 人;5 个模型相同,未参与训练)")
ap.add_argument("--occlude", action="store_true", help="区域遮挡检验(只前向)"); A = ap.parse_args()
DEV = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
ARCH = {"resnet18": "resnet18", "resnet50": "resnet50", "efficientnetb0": "efficientnet_b0", "convnexttiny": "convnext_tiny", "densenet121": "densenet121"}
rows = list(csv.DictReader(open(os.path.join(OUT, "ksaa_for_saliency.csv"))))
if A.heldout:
    HO = set(open(os.path.join(OUT, "heldout_seqn.txt")).read().split()); rows = [r for r in rows if r["SEQN"] in HO]
if A.n: rows = rows[:A.n]
print(f"{len(rows)} participants, device {DEV}", flush=True)
UNIT = {}                                                          # SEQN → [2,S,S] 单位灰度(L 已镜像)
def unit(r):
    if r["SEQN"] not in UNIT: UNIT[r["SEQN"]] = torch.stack([pre(roi_file(r["SEQN"], sd), sd) for sd in "LR"])
    return UNIT[r["SEQN"]]

models, wts = {}, {}
for k, arch in ARCH.items():
    ck = torch.load(os.path.join(MODELS, f"m5_full_{k}_full.pt"), map_location="cpu", weights_only=False)
    a = ck["args"]; assert a["norm"] == "pct" and a["fusion"] == "tf" and a["views"] == "both" and not a["no_mirror"] and ck["views"] == ["L", "R"] and ck["min_age"] == 0 and ck["OFF"] == 0.0
    m = Mammo_AGE(arch=arch, num_output=a["num_output"], nblock=a["nblock"], hidden_size=a["hidden"]); m.load_state_dict(ck["state_dict"], strict=True)
    dev = "cpu" if k == "convnexttiny" and DEV == "mps" and not A.occlude else DEV   # ConvNeXt 反向传播在 MPS 上报错,改用 CPU
    m.eval().to(dev); models[k] = (m, ck["par"], dev); wts[k] = 1.0 / ck["inner_val_mae"]
W = {k: v / sum(wts.values()) for k, v in wts.items()}; print("ensemble weights", {k: round(v, 3) for k, v in W.items()}, flush=True)
def run(m, par, dev, batch):
    AR = torch.arange(100, dtype=torch.float32, device=dev)
    feats = []
    x = ((torch.stack([unit(r) for r in batch]) - par["mu"]) / par["sd"])[:, :, None].float().to(dev)   # [B,2,1,S,S]
    def _hook(mod, i, o):                                              # 不返回值(返回值会替换该层输出)
        o.retain_grad(); feats.append(o)
    h = m.backbone_model.model.register_forward_hook(_hook)
    g, loc, *_ = m(x, use_sto=False); h.remove()
    ages = torch.stack([(F.softmax(l, -1) * AR).sum(-1) for l in loc], 1)   # [B,2] 两膝期望年龄
    pred = ages.mean(1); m.zero_grad(); pred.sum().backward()
    cams = []
    for A_ in feats:                                                   # 两个视图各一次骨干前向:L、R
        w_ = A_.grad.mean((2, 3), keepdim=True); c = F.relu((w_ * A_).sum(1, keepdim=True))
        c = F.interpolate(c, size=(S, S), mode="bilinear", align_corners=False)[:, 0]
        c = c / c.sum((1, 2), keepdim=True).clamp_min(1e-12); cams.append(c.detach())
    return pred.detach().cpu().numpy(), torch.stack(cams, 1).cpu()     # [B], [B,2,S,S]

def occ_inputs(batch):
    """每人 7 个变体(0 = 不遮挡,1–6 = REG 各区模糊),[B*7,2,S,S] 单位灰度;区域按每膝自己的关节线划分。"""
    out = []
    for r in batch:
        u = unit(r); M = [torch.from_numpy(regions(pool(u[v]).numpy())[0]).repeat_interleave(S // R, 1).repeat_interleave(S // R, 2).to(DEV) for v in range(2)]   # [6,S,S]
        u = u.to(DEV); bl = TF.gaussian_blur(u[:, None], kernel_size=97, sigma=16.0)[:, 0]
        out.append(u)
        for g in range(len(REG)): out.append(torch.stack([torch.where(M[v][g], bl[v], u[v]) for v in range(2)]))
    return torch.stack(out)
def run_occ(m, par, dev, xu):
    AR = torch.arange(100, dtype=torch.float32, device=dev); x = ((xu - par["mu"]) / par["sd"])[:, :, None].float().to(dev)
    with torch.no_grad():
        _, loc, *_ = m(x, use_sto=False)
        return torch.stack([(F.softmax(l, -1) * AR).sum(-1) for l in loc], 1).mean(1).cpu().numpy()
tag = ("_heldout" if A.heldout else "") + (f"_n{len(rows)}" if A.n else "")
if A.occlude:
    NV = 1 + len(REG); N = len(rows); P = {k: np.zeros(N * NV) for k in ARCH}; t0 = time.time()
    for b0 in range(0, N, A.bs):
        batch = rows[b0:b0 + A.bs]; xu = occ_inputs(batch)
        for k in ARCH:
            m, par, dev = models[k]; P[k][b0 * NV:(b0 + len(batch)) * NV] = run_occ(m, par, dev, xu)
        for r in batch: UNIT.pop(r["SEQN"], None)
        if (b0 // A.bs) % 25 == 0: print(f"{b0 + len(batch)}/{N}  {time.time() - t0:.0f}s", flush=True)
    ens = sum(P[k] * W[k] for k in ARCH)
    with open(os.path.join(OUT, f"occlusion{tag}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["SEQN", "age", "gap", "region"] + [f"pred_{k}" for k in ARCH] + ["pred_full_ens"])
        for i, r in enumerate(rows):
            for g in range(NV): w.writerow([r["SEQN"], r["age"], r["gap"], "None" if g == 0 else REG[g - 1]] + [f"{P[k][i * NV + g]:.4f}" for k in ARCH] + [f"{ens[i * NV + g]:.4f}"])
    age = np.array([float(r["age"]) for r in rows]); E = ens.reshape(N, NV)
    for g in range(NV): print(f"{'None' if g == 0 else REG[g - 1]:<22} MAE {np.mean(np.abs(E[:, g] - age)):.2f}  r {np.corrcoef(E[:, g], age)[0, 1]:.3f}  mean shift {np.mean(E[:, g] - E[:, 0]):+.2f}")
    print(f"done {time.time() - t0:.0f}s", flush=True); sys.exit()

N = len(rows); P = {k: np.zeros(N) for k in ARCH}; CAM = np.zeros((N, 2, R, R), np.float16); IMG = np.zeros((N, 2, R, R), np.float16)
t0 = time.time()
for b0 in range(0, N, A.bs):
    batch = rows[b0:b0 + A.bs]; acc = None
    for k in ARCH:
        m, par, dev = models[k]; p, c = run(m, par, dev, batch); P[k][b0:b0 + len(batch)] = p
        acc = c * W[k] if acc is None else acc + c * W[k]
    CAM[b0:b0 + len(batch)] = F.adaptive_avg_pool2d(acc.flatten(0, 1)[:, None], R)[:, 0].view(len(batch), 2, R, R).numpy() * (S * S / (R * R))
    for j, r in enumerate(batch):
        for v in range(2): IMG[b0 + j, v] = pool(unit(r)[v]).numpy()
        UNIT.pop(r["SEQN"], None)
    if (b0 // A.bs) % 25 == 0: print(f"{b0 + len(batch)}/{N}  {time.time() - t0:.0f}s", flush=True)
ens = sum(P[k] * W[k] for k in ARCH)
age = np.array([float(r["age"]) for r in rows]); oof = np.array([float(r["pred_raw"]) for r in rows])
print(f"full-model ensemble vs chronological age: MAE {np.mean(np.abs(ens - age)):.2f}, r {np.corrcoef(ens, age)[0, 1]:.3f}; vs out-of-fold prediction r {np.corrcoef(ens, oof)[0, 1]:.3f}", flush=True)
with open(os.path.join(OUT, f"saliency_pred{tag}.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["SEQN", "age", "gap", "pred_oof"] + [f"pred_{k}" for k in ARCH] + ["pred_full_ens"])
    for i, r in enumerate(rows): w.writerow([r["SEQN"], r["age"], r["gap"], r["pred_raw"]] + [f"{P[k][i]:.4f}" for k in ARCH] + [f"{ens[i]:.4f}"])
np.save(os.path.join(OUT, f"cam{tag}.npy"), CAM); np.save(os.path.join(OUT, f"img{tag}.npy"), IMG)
print(f"done {time.time() - t0:.0f}s", flush=True)
