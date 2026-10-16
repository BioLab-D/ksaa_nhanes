# 解释性分析 · 第 1 步(不用模型):4,875 人双膝图像按训练时预处理后池化到 112 × 112,供平均图与差值图使用;
#   另写出 heldout_seqn.txt(5 个全数据模型共同的内部验证集,未参与训练)。输入 ksaa_for_saliency.csv 来自 saliency_input.R。
# 跑法(仓库根目录下):python code/interpretability/saliency_images.py
import os, sys, csv, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saliency_common import OUT, MODELS, R, roi_file, pre, pool
torch.set_num_threads(2)
rows = list(csv.DictReader(open(os.path.join(OUT, "ksaa_for_saliency.csv"))))
IMG = np.full((len(rows), 2, R, R), np.nan, np.float16); t0 = time.time()
for i, r in enumerate(rows):
    for v, sd in enumerate("LR"):
        IMG[i, v] = pool(pre(roi_file(r["SEQN"], sd), sd)).numpy()
    if i % 1000 == 0: print(i, f"{time.time() - t0:.0f}s", flush=True)
np.save(os.path.join(OUT, "img_all.npy"), IMG)
V = {}
for k in ["resnet18", "resnet50", "efficientnetb0", "convnexttiny", "densenet121"]:
    ck = torch.load(os.path.join(MODELS, f"m5_full_{k}_full.pt"), map_location="cpu", weights_only=False)
    V[k] = set(map(str, ck["val_seqn"]))
assert all(v == V["resnet18"] for v in V.values())                  # 5 个模型的内部验证集相同
HO = sorted((r["SEQN"] for r in rows if r["SEQN"] in V["resnet18"]), key=int)
with open(os.path.join(OUT, "heldout_seqn.txt"), "w") as fh: fh.write("\n".join(HO))
print(len(HO), "held-out participants in the analysis cohort")   
with open(os.path.join(OUT, "img_all_seqn.txt"), "w") as fh: fh.write("\n".join(r["SEQN"] for r in rows))
print("missing views:", int(np.isnan(IMG[:, :, 0, 0].astype(np.float32)).sum()), f"done {time.time() - t0:.0f}s")
