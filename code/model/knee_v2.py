"""膝适配版训练脚本 v2(2026-09-10)。

以官方移植版 mammo_official.py 为起点,把核查中提出的问题逐项做成开关;**所有开关的默认值 = v1(官方移植版)**,
不带开关运行即复现 v1,筛选时一次只改一项。
  --aug official|knee      knee:删弹性/网格/光学畸变与中值模糊;几何(旋转 ±5°、平移 ≤5%、缩放 0.95–1.05,p 0.5)左右膝各自独立;
                           光度(亮度/对比度 ±10% p 0.5、gamma 0.9–1.1 p 0.3、3×3 高斯模糊 p 0.15)同一张片左右共用一组参数;
                           高斯噪声(SD 0.007–0.012,即 v1 的 var 3–9 / 255²,p 0.2)强度共用、各自抽样;
                           CoarseDropout(8 块、边长 10%,p 0.2,同 v1 的 albumentations 默认行为)各自独立
  --norm minmax8|pct|global  minmax8:逐图 min–max → uint8(v1);pct:逐图 [P0.5, P99.5] 线性映射到 float32 [0,1];
                           global:训练折像素池化的 [P0.5, P99.5] 映射到 [0,1],保留图间灰度差异;pct/global 再用训练折像素均值/SD 标准化
  --sched step|cosine      step:Adam、wd 0、第 20/40/60/80 轮 ×0.1(v1);cosine:AdamW(--wd)、--warmup 轮线性升温后逐步余弦
  --fusion tf|concat|mean|gate  tf:官方 instance-bag transformer(--nblock,v1 = 10);concat:官方 POE 基线(各视图池化特征拼接,无局部头);
                           mean:双膝池化特征平均;gate:门控注意力汇聚(Ilse 2018)。mean/gate 另有两膝共享的单膝局部头
  --lam_local λ            局部头损失权重(v1 = 1,即 L_global + L_R + L_L)
  --primary local|global   早停与 pred_age 用哪个读出(v1 = local:两膝局部头期望年龄的均值)
  --min_age/--num_output   POE 类别范围(v1 = 0/100;POE-31 = 60/31)
  --reg [--no_poe]         标量回归(L1,目标以训练折平均年龄为中心);加 --no_poe 即无 POE 采样的直接 L1 回归
  --no_mirror              左膝不镜像成右膝方向
  --views both|R|L         单膝模型(只用右膝或左膝训练和预测,须配 mean/gate);默认双膝
  --no_cross               transformer 仅自注意力(交叉注意力输出置零;配 --fusion tf 且 --nblock ≥ 1)
  --no_density             官方密度分支输出置零。该分支由各视图拼接特征算出、拼进每个局部头,不置零时"局部头"并非单膝;
                           本任务无密度标签(λ = 0),它只是个无监督的双膝瓶颈。要让局部通路真正只看单膝(N = 0 或仅自注意力),必须加此开关
  --roi_dir                ROI 目录(默认 1.15 倍骨宽的 knee_roi)
  --folds                  只跑部分外层折(筛选用 0,1)
  --smoke N                冒烟测试:训练 N 人、1 轮,不写任何输出
折划分固定用 --split_seed(= v1 的 20260831),与训练随机种子 --seed 解耦(v1 两者是同一个,换种子会连折划分一起换)。
只纳入左右膝都有 ROI 的人(v1 对单侧者复制成两份;实际 0 人)。v1 的文件名过滤只认 K1 片,漏掉了只有 K2 片的 3 人(有年龄和死亡随访);
v2 纳入他们,并保证原 4,872 人的外层折与内部验证集不变(见"折"一段)。
输出(逐折追加):/home/data/knee_yolo/v2/{tag}_folds.csv(fold、inner_val_mae、best_ep、stop_ep、minutes)——筛选只看这个;
留出折盲化:blind/{tag}_oof.csv(SEQN、fold、age、pred_age = 主读出、pred_global、pred_local、pred_L、pred_R)与 blind/{tag}_holdout.csv(各读出的留出 MAE),
日志不打印留出 MAE,模型定稿前不打开 blind/。冒烟测试(--smoke)写到 v2/smoke/。
"""
import os, re, csv, glob, math, random, argparse, time, numpy as np, cv2
os.environ.setdefault("N_VIEWS", "2")
import torch, torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mammo_age_model import Mammo_AGE, POE_Baselie_Model
from mean_variance_loss import MeanVarianceLoss
from probordiloss import ProbOrdiLoss
from metrics import cal_mae_acc_cls, cal_mae_acc_reg
cv2.setNumThreads(1)

ap = argparse.ArgumentParser()
ap.add_argument("--arch", required=True, choices=["resnet18", "resnet50", "efficientnet_b0", "convnext_tiny", "densenet121"])
ap.add_argument("--tag", required=True); ap.add_argument("--folds", default="0,1,2,3,4"); ap.add_argument("--smoke", type=int, default=0)
ap.add_argument("--seed", type=int, default=20260831); ap.add_argument("--split_seed", type=int, default=20260831)
ap.add_argument("--img", type=int, default=448); ap.add_argument("--bs", type=int, default=16); ap.add_argument("--lr", type=float, default=1e-4)
ap.add_argument("--epochs", type=int, default=100); ap.add_argument("--min_ep", type=int, default=30); ap.add_argument("--patience", type=int, default=10)
ap.add_argument("--sched", default="step", choices=["step", "cosine"]); ap.add_argument("--wd", type=float, default=0.); ap.add_argument("--warmup", type=int, default=5)
ap.add_argument("--aug", default="official", choices=["official", "knee"]); ap.add_argument("--norm", default="minmax8", choices=["minmax8", "pct", "global"])
ap.add_argument("--fusion", default="tf", choices=["tf", "concat", "mean", "gate"]); ap.add_argument("--nblock", type=int, default=10); ap.add_argument("--hidden", type=int, default=128)
ap.add_argument("--lam_local", type=float, default=1.0); ap.add_argument("--primary", default="local", choices=["local", "global"])
ap.add_argument("--min_age", type=int, default=0); ap.add_argument("--num_output", type=int, default=100); ap.add_argument("--max_t", type=int, default=50)
ap.add_argument("--reg", action="store_true"); ap.add_argument("--no_poe", action="store_true"); ap.add_argument("--no_mvl", action="store_true"); ap.add_argument("--no_sto", action="store_true")
ap.add_argument("--no_mirror", action="store_true"); ap.add_argument("--roi_dir", default="/home/data/knee_roi"); ap.add_argument("--save_dir", default="")
ap.add_argument("--views", default="both", choices=["both", "R", "L"]); ap.add_argument("--no_cross", action="store_true"); ap.add_argument("--no_density", action="store_true")
# --full:全部人训练最终模型并存档(内部验证 10% 仅作早停,无留出)
ap.add_argument("--full", action="store_true")
# --target phenoage:监督目标换成 --label_csv 中的该列;折与膝龄模型相同,无目标值者剔除;目标取整并截到 POE 类别范围
ap.add_argument("--target", default="age", choices=["age", "phenoage"]); ap.add_argument("--label_csv", default="/home/data/phenoage_labels.csv")
Aa = ap.parse_args(); DEV = "cuda"; OUT = "/home/data/knee_yolo/v2"; os.makedirs(OUT, exist_ok=True)
assert not Aa.full or Aa.save_dir, "--full 必须配 --save_dir(它的唯一产物就是存档)"
if Aa.aug == "official": assert Aa.norm == "minmax8", "官方增强(GaussNoise 等)按 uint8 灰阶设参,只能配 minmax8"
if Aa.fusion == "concat": assert Aa.primary == "global", "concat 没有局部头,--primary 须为 global"
if Aa.views != "both": assert Aa.fusion in ("mean", "gate"), "单膝模型只能配 mean/gate(官方模块按 N_VIEWS=2 建层)"
VS = ["L", "R"] if Aa.views == "both" else [Aa.views]
if Aa.no_cross:
    assert Aa.fusion == "tf" and Aa.nblock >= 1, "--no_cross 只用于官方 transformer 且 N ≥ 1"
    import transformer as _T
    def _self_only(self, local_feats, global_feat):
        """仅自注意力:每块的交叉注意力输出置零(FFN 输入维度不变),两膝在局部通路里不交互;其余同官方 forward。"""
        out = []
        for feat in local_feats:
            x = self.local_proj(feat); shp = x.shape
            for bi in range(self.num_blocks):
                sa = self._to_spatial(self.self_attn_layers[bi](self._to_sequence(x)), shp)
                x = x + self.ffn_layers[bi](torch.cat([x, sa, torch.zeros_like(sa)], dim=1))
            out.append(torch.cat([x, feat], dim=1))
        return out
    _T.Instance_bag_transformer.forward = _self_only
torch.manual_seed(Aa.seed); np.random.seed(Aa.seed)
S = Aa.img; SPLIT = Aa.split_seed; use_sto = not (Aa.no_sto or Aa.no_poe)

# ---------------- 数据与折 ----------------
lab = {r["SEQN"]: float(r["age"]) for r in csv.DictReader(open("/tmp/surv_labels.csv"))}
PAT = re.compile(r"^(\d+)K(\d)_([RL])$"); films = {}
for f in sorted(glob.glob(f"{Aa.roi_dir}/*.png")):
    m_ = PAT.match(os.path.basename(f)[:-4])
    if m_ and m_.group(1) in lab: films.setdefault(m_.group(1), {}).setdefault(int(m_.group(2)), {})[m_.group(3)] = f
bags, single = {}, []
for s, fs in films.items():
    ok = [n for n in sorted(fs) if len(fs[n]) == 2]
    if ok: bags[s] = fs[ok[0]]                       # 取编号最小、左右齐全的一张片
    else: single.append(s)
# 折:有 K1 片者按 split_seed 分五折;只有 K2 片者另行分折,只并入训练部分,不进内部验证集
base = sorted(s for s in bags if len(films[s].get(1, {})) == 2); bset = set(base); extra = sorted(s for s in bags if s not in bset)
random.seed(SPLIT); random.shuffle(base); random.Random(SPLIT).shuffle(extra)
fold_of = {s: i % 5 for i, s in enumerate(base)}; fold_of.update({s: i % 5 for i, s in enumerate(extra)})
random.seed(Aa.seed); OFF = 0.0
print(f"{len(bags)} 人(v1 队列 {len(base)} + 仅 K2 片 {len(extra)};无双侧 ROI 排除 {len(single)})  " + "  ".join(f"{k}={v}" for k, v in vars(Aa).items()), flush=True)
AGE, raw_t = lab, {}
if Aa.target != "age":   # 折已按年龄队列定好,此后 lab = 训练目标(只含有值者)
    raw_t = {r["SEQN"]: float(r[Aa.target]) for r in csv.DictReader(open(Aa.label_csv)) if r[Aa.target] not in ("", "NA")}
    lo_c, hi_c = Aa.min_age, Aa.min_age + Aa.num_output - 1
    lab = {s: float(min(max(round(v), lo_c), hi_c)) for s, v in raw_t.items() if s in bags}
    print(f"目标 {Aa.target}:队列内有值 {len(lab)} / {len(bags)} 人;截到 [{lo_c}, {hi_c}] 的 {sum(1 for s in lab if round(raw_t[s]) != lab[s])} 人", flush=True)

def to_unit(im, par):
    """16-bit ROI → 模型灰度。minmax8 返回 uint8(v1);pct/global 返回 float32 [0,1]。"""
    if Aa.norm == "minmax8":
        mn, mx = float(im.min()), float(im.max()); return ((im - mn) / max(mx - mn, 1.0) * 255.0).astype(np.uint8)
    lo, hi = np.percentile(im[::2, ::2], [0.5, 99.5]) if Aa.norm == "pct" else (par["lo"], par["hi"])
    return np.clip((im - lo) / max(hi - lo, 1.0), 0.0, 1.0).astype(np.float32)

def fit_norm(tr, k):
    """只在训练折内估计归一化参数:固定种子抽 ≤400 人(两膝),每图取 1/64 像素。"""
    if Aa.norm == "minmax8": return {"mu": 0.5, "sd": 0.5}
    pick = np.random.default_rng(SPLIT + k).choice(len(tr), min(400, len(tr)), replace=False); raw, unit = [], []
    for i in pick:
        for f in bags[tr[i]].values():
            im = cv2.imread(f, cv2.IMREAD_ANYDEPTH).astype(np.float32)
            (unit if Aa.norm == "pct" else raw).append(to_unit(im, None)[::8, ::8].ravel() if Aa.norm == "pct" else im[::8, ::8].ravel())
    par = {}
    if Aa.norm == "global":
        raw = np.concatenate(raw); lo, hi = np.percentile(raw, [0.5, 99.5]); par.update(lo=float(lo), hi=float(hi))
        px = np.clip((raw - lo) / max(hi - lo, 1.0), 0, 1)
    else: px = np.concatenate(unit)
    par.update(mu=float(px.mean()), sd=float(px.std())); return par

if Aa.aug == "official":
    import albumentations as A
    from albumentations.pytorch import ToTensorV2
    TF_TR = A.Compose([A.Resize(S, S, interpolation=cv2.INTER_LINEAR),
        A.RandomBrightnessContrast(brightness_limit=(-0.1, 0.1), contrast_limit=(-0.1, 0.1), p=0.5), A.ShiftScaleRotate(rotate_limit=(-5, 5), p=0.3),
        A.OneOf([A.MedianBlur(blur_limit=3, p=0.3), A.GaussianBlur(blur_limit=3, p=0.3), A.GaussNoise(var_limit=(3.0, 9.0), p=0.3)], p=0.3),
        A.OneOf([A.ElasticTransform(p=0.3), A.GridDistortion(num_steps=5, distort_limit=1., p=0.3), A.OpticalDistortion(distort_limit=1., p=0.3)], p=0.3),
        A.CoarseDropout(max_height=int(S * 0.1), max_width=int(S * 0.1), max_holes=8, p=0.2), A.Normalize(mean=0.5, std=0.5), ToTensorV2()])
    TF_TE = A.Compose([A.Resize(S, S, interpolation=cv2.INTER_LINEAR), A.Normalize(mean=0.5, std=0.5), ToTensorV2()])

def knee_aug(xs, rng):
    """xs:同一人两膝的 float32 [0,1] 图(已缩放到 S)。光度参数左右共用;几何、噪声实现、遮挡各自独立。"""
    bc = (rng.uniform(0.9, 1.1), rng.uniform(-0.1, 0.1)) if rng.random() < 0.5 else None
    gm = rng.uniform(0.9, 1.1) if rng.random() < 0.3 else None
    bl = rng.random() < 0.15
    ns = rng.uniform(0.007, 0.012) if rng.random() < 0.2 else None
    out = []
    for x in xs:
        if bc: x = np.clip(x * bc[0] + bc[1], 0, 1)
        if gm: x = np.power(x, gm)
        if rng.random() < 0.5:
            M = cv2.getRotationMatrix2D((S / 2, S / 2), rng.uniform(-5, 5), rng.uniform(0.95, 1.05)); M[:, 2] += rng.uniform(-0.05, 0.05, 2) * S
            x = cv2.warpAffine(x.astype(np.float32), M, (S, S), flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT_101)
        if bl: x = cv2.GaussianBlur(x.astype(np.float32), (3, 3), 0)
        if ns: x = np.clip(x + rng.normal(0, ns, x.shape), 0, 1)
        if rng.random() < 0.2:
            x = x.copy(); h = int(S * 0.1)
            for _ in range(8):
                y0, x0 = rng.integers(0, S - h + 1, 2); x[y0:y0 + h, x0:x0 + h] = 0
        out.append(x.astype(np.float32))
    return out

class BagDS(Dataset):
    def __init__(s, sq, train, par): s.sq, s.train, s.par = sq, train, par
    def __len__(s): return len(s.sq)
    def __getitem__(s, i):
        sq = s.sq[i]; sides = list(VS)
        if s.train: random.shuffle(sides)          # 训练时视图顺序随机,预测时固定
        ims = []
        for sd in sides:
            x = to_unit(cv2.imread(bags[sq][sd], cv2.IMREAD_ANYDEPTH).astype(np.float32), s.par)
            if sd == "L" and not Aa.no_mirror: x = x[:, ::-1]
            ims.append(np.ascontiguousarray(x))
        if Aa.aug == "official":
            tf = TF_TR if s.train else TF_TE; t = torch.stack([tf(image=x)["image"] for x in ims])
        else:
            ims = [cv2.resize(x, (S, S), interpolation=cv2.INTER_LINEAR) if x.shape != (S, S) else x for x in ims]
            ims = [x.astype(np.float32) / 255.0 if x.dtype == np.uint8 else x for x in ims]
            if s.train: ims = knee_aug(ims, np.random.default_rng(int(torch.randint(0, 2**31 - 1, (1,)).item())))
            t = torch.stack([torch.from_numpy((x - s.par["mu"]) / s.par["sd"]).float()[None] for x in ims])
        return {"img": t, "age": torch.tensor(lab[sq], dtype=torch.float32)}

# ---------------- 模型 ----------------
class Head(nn.Module):
    """POE 头:同官方 POE_Baselie_Model 的 emd/var/drop/final(去掉无标签的密度分支)。"""
    def __init__(s, din, nout):
        super().__init__()
        s.emd = nn.Sequential(nn.Linear(din, din), nn.ReLU(True)); s.var = nn.Sequential(nn.Linear(din, din), nn.BatchNorm1d(din, eps=0.001, affine=False))
        s.drop = nn.Dropout(0.1); s.final = nn.Linear(din, nout)
        for l in (s.emd[0], s.var[0], s.final): nn.init.kaiming_normal_(l.weight, mode="fan_out", nonlinearity="relu"); nn.init.zeros_(l.bias)
    def forward(s, x, max_t, use_sto):
        emb = s.emd(x); lv = s.var(x)
        if use_sto:
            sto = emb[None].expand(max_t, *emb.shape) + torch.exp(lv * 0.5)[None].expand(max_t, *lv.shape) * torch.randn(max_t, *emb.shape, device=emb.device)
            return s.final(s.drop(sto)), emb, lv
        return s.final(s.drop(emb)), emb, lv

class NoTF(nn.Module):
    def __init__(s, nout): super().__init__(); s.m = POE_Baselie_Model(Aa.arch, nout, backbone_model=False)
    def forward(s, x, max_t=50, use_sto=False): g, ge, glv, d = s.m(x, max_t, use_sto); return g, [], ge, [], glv, [], d

class Fuse(nn.Module):
    """mean / gate:两膝共享骨干 → 池化特征 f_v;全局 = 融合特征 → latent → 全局 POE 头;局部 = 各 f_v → 共享的局部 POE 头。"""
    def __init__(s, nout, mode):
        super().__init__(); base = POE_Baselie_Model(Aa.arch, nout, backbone_model=True); s.trunk = base.model; C = base.num_feat
        s.mode = mode; s.avg = nn.AdaptiveAvgPool2d(1)
        if mode == "gate": s.aV = nn.Linear(C, 128); s.aU = nn.Linear(C, 128); s.aw = nn.Linear(128, 1)
        s.latent = nn.Sequential(nn.Dropout(0.2), nn.Linear(C, C), nn.ReLU(True))
        nn.init.kaiming_normal_(s.latent[1].weight, mode="fan_out", nonlinearity="relu"); nn.init.zeros_(s.latent[1].bias)
        s.gh = Head(C, nout); s.lh = Head(C, nout)
    def forward(s, x, max_t=50, use_sto=False):
        B, V = x.shape[:2]; f = torch.flatten(s.avg(s.trunk(x.flatten(0, 1).repeat(1, 3, 1, 1))), 1).view(B, V, -1)
        z = f.mean(1) if s.mode == "mean" else (torch.softmax(s.aw(torch.tanh(s.aV(f)) * torch.sigmoid(s.aU(f))), 1) * f).sum(1)
        g, ge, glv = s.gh(s.latent(z), max_t, use_sto); L = [s.lh(f[:, v], max_t, use_sto) for v in range(V)]
        return g, [l[0] for l in L], ge, [l[1] for l in L], glv, [l[2] for l in L], None

class Zero4(nn.Module):
    """替代无标签的密度分支:输出全零 4 维。官方 POE_Baselie_Model 仍把它拼进 latent、Mammo_AGE 仍拼进每个局部头,但不再携带信息。"""
    def forward(s, x): return x.new_zeros(x.shape[0], 4)

def build():
    nout = 1 if Aa.reg else Aa.num_output
    m = Mammo_AGE(arch=Aa.arch, num_output=nout, nblock=Aa.nblock, hidden_size=Aa.hidden) if Aa.fusion == "tf" else NoTF(nout) if Aa.fusion == "concat" else Fuse(nout, Aa.fusion)
    if Aa.no_density:   # 无密度标签;置零使局部头只看单膝(该分支混入对侧膝信息)
        for mod in list(m.modules()):
            if isinstance(mod, POE_Baselie_Model): mod.predict_density = Zero4()
    c1 = ProbOrdiLoss(distance="JDistance", alpha_coeff=(0. if Aa.no_poe else 1e-5), beta_coeff=(0. if Aa.no_poe else 1e-4), margin=5,
                      main_loss_type=("reg" if Aa.reg else "cls"), criterion="l1")
    c2 = None if (Aa.reg or Aa.no_mvl) else MeanVarianceLoss(lambda_1=0.2, lambda_2=0.05, cumpet_ce_loss=False, start_age=0)
    return m.to(DEV), (c1, c2)

def compute_loss(pred, emb, lv, age, crit):
    c1, c2 = crit; _, _, _, pol = c1(pred, emb, lv, age, None, use_sto=use_sto)
    if c2 is None: return pol
    p_, a_ = (pred.view(-1, pred.shape[-1]), age.repeat(pred.shape[0])) if use_sto else (pred, age)
    return pol + c2(p_, a_)

def readouts(g, loc, age):
    """各头的期望年龄(已减 min_age):全局、两膝局部均值、各膝。"""
    f = (lambda o: cal_mae_acc_reg(o, age, use_sto)["pred_age"]) if Aa.reg else (lambda o: cal_mae_acc_cls(o, age, use_sto)["pred_age"])
    pg = np.asarray(f(g), float); pv = [np.asarray(f(l), float) for l in loc]
    return pg, (np.mean(pv, 0) if pv else pg), pv

def train_step(m, b, crit):
    """目标 = 年龄 − min_age − OFF(OFF 仅回归时为训练折平均年龄,免得输出层从 0 爬到 70 多岁)。"""
    img = b["img"].to(DEV, non_blocking=True); age = b["age"].to(DEV) - Aa.min_age - OFF
    g, loc, ge, le, glv, llv, _ = m(img, max_t=Aa.max_t, use_sto=use_sto)
    loss = compute_loss(g, ge, glv, age, crit)
    if loc: loss = loss + Aa.lam_local * sum(compute_loss(loc[i], le[i], llv[i], age, crit) for i in range(len(loc)))
    with torch.no_grad(): pg, pl, _ = readouts(g.detach(), [l.detach() for l in loc], age)
    return loss, (pl if Aa.primary == "local" else pg) + Aa.min_age + OFF

def predict(m, sq, par):
    m.eval(); out = {c: [] for c in ("pred_global", "pred_local", "pred_L", "pred_R")}; sh = Aa.min_age + OFF
    with torch.no_grad():
        for b in DataLoader(BagDS(sq, False, par), batch_size=Aa.bs, num_workers=8, pin_memory=True):
            img = b["img"].to(DEV, non_blocking=True); age = b["age"].to(DEV) - sh
            g, loc, *_ = m(img, max_t=Aa.max_t, use_sto=use_sto); pg, pl, pv = readouts(g, loc, age); n = len(pg)
            out["pred_global"] += list(pg + sh); out["pred_local"] += list(pl + sh)
            for sd in ("L", "R"):                      # 预测时视图顺序固定为 VS
                out[f"pred_{sd}"] += list(pv[VS.index(sd)] + sh) if (pv and sd in VS) else [np.nan] * n
    out = {c: np.array(v, dtype=float) for c, v in out.items()}
    out["pred_age"] = out["pred_local"] if Aa.primary == "local" else out["pred_global"]
    return out

# ---------------- 训练 ----------------
# 留出折盲化:日志只有内部验证 MAE;留出预测写进 blind/,模型定稿前不打开
ROOT = f"{OUT}/smoke" if Aa.smoke else OUT; BL = f"{ROOT}/blind"; os.makedirs(BL, exist_ok=True)
FO, FF, FH = f"{BL}/{Aa.tag}_oof.csv", f"{ROOT}/{Aa.tag}_folds.csv", f"{BL}/{Aa.tag}_holdout.csv"
for p in (FO, FF, FH):
    if os.path.exists(p): (os.remove(p) if Aa.smoke else os.rename(p, f"{p}.bak{int(time.time())}"))
t0 = time.time()
for k in ([5] if Aa.full else [int(x) for x in Aa.folds.split(",")]):
    if Aa.full:     # 最终模型:k = 5 只作种子,无留出
        tr_all = list(base); te = []
    else:
        tr_all = [s for s in base if fold_of[s] != k]; te = [s for s in base + extra if fold_of[s] == k]
    rng = random.Random(SPLIT + k); rng.shuffle(tr_all); nv = len(tr_all) // 10; va = tr_all[:nv]; tr = tr_all[nv:] + [s for s in extra if Aa.full or fold_of[s] != k]
    if Aa.target != "age": tr, va, te = [s for s in tr if s in lab], [s for s in va if s in lab], [s for s in te if s in lab]   # 划分后再剔除无目标值者
    if Aa.smoke: tr, va, te = tr[:Aa.smoke], va[:2 * Aa.bs], te[:2 * Aa.bs]
    OFF = float(np.mean([lab[s] for s in tr])) - Aa.min_age if Aa.reg else 0.0
    par = fit_norm(tr, k); m, crit = build()
    opt = torch.optim.Adam(m.parameters(), lr=Aa.lr, weight_decay=0.) if Aa.sched == "step" else torch.optim.AdamW(m.parameters(), lr=Aa.lr, weight_decay=Aa.wd)
    dl = DataLoader(BagDS(tr, True, par), batch_size=Aa.bs, shuffle=True, num_workers=8, drop_last=True, pin_memory=True)
    ne = 1 if Aa.smoke else Aa.epochs; ipe = len(dl); wu = Aa.warmup * ipe; tot = ne * ipe; it = 0
    best, best_ep, early, best_state = 1e9, 0, 0, None
    print(f"  fold {k}: 训练 {len(tr)} / 内部验证 {len(va)} / 留出 {len(te)} 人  归一化 {par}", flush=True)
    for ep in range(1, ne + 1):
        m.train(); tl, tm = [], []
        for b in dl:
            if Aa.sched == "step": lr = Aa.lr * (0.1 ** sum(ep >= ms for ms in [20, 40, 60, 80]))
            else: lr = Aa.lr * (it + 1) / wu if it < wu else Aa.lr * 0.5 * (1 + math.cos(math.pi * (it - wu) / max(tot - wu, 1)))
            for g_ in opt.param_groups: g_["lr"] = lr
            opt.zero_grad(); loss, pa = train_step(m, b, crit); loss.backward(); opt.step(); it += 1
            tl.append(loss.item()); tm.append(float(np.mean(np.abs(pa - b["age"].numpy()))))
        P = predict(m, va, par); vmae = float(np.mean(np.abs(P["pred_age"] - np.array([lab[s] for s in va]))))
        if vmae < best: best, best_ep, early = vmae, ep, 0; best_state = {k_: v.detach().clone() for k_, v in m.state_dict().items()}
        else: early += 1
        if Aa.smoke or ep % 5 == 0 or ep <= 2 or early > Aa.patience:
            print(f"  fold {k} ep {ep:3d} lr {lr:.1e}  train loss {np.mean(tl):.3f}  train MAE {np.mean(tm):.3f}  val MAE {vmae:.3f}  best {best:.3f} (ep {best_ep})  {(time.time() - t0) / 60:.1f} min", flush=True)
        if early > Aa.patience and ep >= Aa.min_ep: break
    if Aa.full:
        m.load_state_dict(best_state); sdir = f"{Aa.save_dir}/smoke" if Aa.smoke else Aa.save_dir; os.makedirs(sdir, exist_ok=True); ck = f"{sdir}/{Aa.tag}_full.pt"
        torch.save({"state_dict": best_state, "arch": Aa.arch, "par": par, "args": vars(Aa), "views": VS, "min_age": Aa.min_age, "OFF": OFF,
                    "inner_val_mae": best, "best_ep": best_ep, "stop_ep": ep, "n_train": len(tr), "val_seqn": va,
                    "note": "knee age v2 final model: all participants (10% internal validation for early stopping only), no held-out; "
                            "prediction = mean of the two knee-level expected ages; use ensemble weights and bias-correction coefficients calibrated on the five-fold out-of-fold predictions"}, ck)
        # 回读核对(POE 推理有随机采样,两次固定同一种子)
        torch.manual_seed(12345); P1 = predict(m, va, par)
        m2, _ = build(); m2.load_state_dict(torch.load(ck, map_location=DEV)["state_dict"]); torch.manual_seed(12345); P2 = predict(m2, va, par)
        dmax = float(np.nanmax(np.abs(P1["pred_age"] - P2["pred_age"])))
        assert dmax < 0.01, f"存档回读后预测不一致:最大差 {dmax}"
        row = dict(fold="full", inner_val_mae=round(best, 4), best_ep=best_ep, stop_ep=ep, minutes=round((time.time() - t0) / 60, 1))
        with open(FF, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row), lineterminator="\n"); w.writeheader(); w.writerow(row)
        with open(f"{ROOT}/{Aa.tag}_innerval.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n"); w.writerow(["SEQN", "age", "pred_age"])
            for i, s in enumerate(va): w.writerow([s, lab[s], round(float(P2["pred_age"][i]), 4)])
        print(f"  full done  内部验证 {best:.3f}  最优轮 {best_ep}  停止轮 {ep}  存档 {ck}  回读最大差 {dmax:.1e}  累计 {(time.time() - t0) / 60:.1f} min", flush=True)
        del m, m2, best_state; torch.cuda.empty_cache(); continue
    m.load_state_dict(best_state); P = predict(m, te, par); AG = np.array([lab[s] for s in te])
    mae = lambda p: float(np.mean(np.abs(p - AG))) if not np.isnan(p).all() else float("nan")
    row = dict(fold=k, inner_val_mae=round(best, 4), best_ep=best_ep, stop_ep=ep, minutes=round((time.time() - t0) / 60, 1))
    hrow = dict(fold=k, holdout_mae=mae(P["pred_age"]), holdout_global=mae(P["pred_global"]), holdout_local=mae(P["pred_local"]),
                holdout_L=mae(P["pred_L"]), holdout_R=mae(P["pred_R"]), holdout_global_local=mae((P["pred_global"] + P["pred_local"]) / 2))
    print("  fold {fold} done  内部验证 {inner_val_mae:.3f}  最优轮 {best_ep}  停止轮 {stop_ep}  累计 {minutes} min(留出结果已写入 blind/,定稿前不看)".format(**row), flush=True)
    if Aa.save_dir and not Aa.smoke: os.makedirs(Aa.save_dir, exist_ok=True); torch.save(best_state, f"{Aa.save_dir}/{Aa.tag}_fold{k}.pt")
    new = not os.path.exists(FO)
    with open(FO, "a", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        tx = Aa.target != "age"   # 非年龄目标:另记原始与训练用(取整截断)目标
        if new: w.writerow(["SEQN", "fold", "age"] + (["target_raw", "target"] if tx else []) + ["pred_age", "pred_global", "pred_local", "pred_L", "pred_R"])
        for i, s in enumerate(te): w.writerow([s, k, AGE[s]] + ([raw_t[s], lab[s]] if tx else []) + [round(float(P[c][i]), 4) for c in ("pred_age", "pred_global", "pred_local", "pred_L", "pred_R")])
    for path, r in ((FF, row), (FH, hrow)):
        new = not os.path.exists(path)
        with open(path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(r), lineterminator="\n")
            if new: w.writeheader()
            w.writerow(r)
    del m, best_state; torch.cuda.empty_cache()
print(f"=== {Aa.tag} 完成  {(time.time() - t0) / 60:.1f} min ===", flush=True)
