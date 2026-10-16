"""旋转救援:对检测失败的片子试转 90°,重跑检测。

审核图显示旋转是失败的压倒性主因。原 QC 只标出 2 张"旋转",是因为
判据要求两个膝都被检出 —— 而旋转片检测器根本检不出任何框(训练集全是
正立片),直接掉进 none。真实旋转片数量远高于 QC 报告。

不能整批丢弃这些片:审核图里已出现全膝关节置换病例,系统性剔除会制造
与骨骼健康直接相关的选择偏倚。

两个方向都试,取"恰好一左一右且水平并排"且平均置信度更高的那个。
"""
import os, csv, numpy as np
from PIL import Image
from ultralytics import YOLO
SRC = "/home/data/knee_yolo"

fails = [r for r in csv.DictReader(open(f"{SRC}/qc.csv")) if r["status"] != "ok"]
print(f"待救援 {len(fails)} 张", flush=True)
m = YOLO(f"{SRC}/runs/det_v1/weights/best.pt")
os.makedirs(f"{SRC}/rot", exist_ok=True)

def try_one(img):
    rr = m.predict(img, imgsz=1024, conf=0.25, device=0, verbose=False)[0]
    b = rr.boxes
    if b is None or len(b) == 0: return None
    cls = b.cls.cpu().numpy().astype(int); xywh = b.xywh.cpu().numpy(); cf = b.conf.cpu().numpy()
    if (cls==0).sum() != 1 or (cls==1).sum() != 1: return None
    pR, pL = xywh[cls==0][0], xywh[cls==1][0]
    if abs(pR[0]-pL[0]) < abs(pR[1]-pL[1]): return None      # 仍是上下堆叠
    if pR[0] > pL[0]: return None                            # 右膝应在画面左侧
    return float(cf.mean()), cls, xywh, rr.orig_shape

out = []
saved = {"90":0, "270":0}
for i, r in enumerate(fails):
    p = f"{SRC}/images/{r['stem']}.png"
    im0 = Image.open(p)
    best = None
    for ang, tag in ((90,"90"), (270,"270")):
        im = im0.rotate(ang, expand=True)
        tmp = f"{SRC}/rot/{r['stem']}_{tag}.png"; im.save(tmp)
        got = try_one(tmp)
        if got and (best is None or got[0] > best[0]):
            best = got; best_tag = tag
        os.remove(tmp)
    if best:
        conf, cls, xywh, shp = best
        H, W = shp
        saved[best_tag] += 1
        for c, xy in zip(cls, xywh):
            out.append(dict(stem=r["stem"], rot=best_tag, side="R" if c==0 else "L",
                            cx=round(xy[0]/W,6), cy=round(xy[1]/H,6),
                            w=round(xy[2]/W,6), h=round(xy[3]/H,6), conf=round(conf,4)))
    if (i+1) % 60 == 0: print(f"  {i+1}/{len(fails)}  已救回 {len(out)//2}", flush=True)

n = len(out)//2
print(f"\n救回 {n} / {len(fails)}  ({n/len(fails):.1%})")
print(f"  顺时针90°: {saved['90']}   逆时针90°: {saved['270']}")
print(f"  仍失败: {len(fails)-n}")
with open(f"{SRC}/boxes_rot.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["stem","rot","side","cx","cy","w","h","conf"])
    w.writeheader(); w.writerows(out)
print(f"-> {SRC}/boxes_rot.csv")
