import os, csv, numpy as np
from PIL import Image
from ultralytics import YOLO
SRC = "/home/data/knee_yolo"

fails = [r for r in csv.DictReader(open(f"{SRC}/qc.csv")) if r["status"] != "ok"]
print(f"{len(fails)} films to rescue", flush=True)
m = YOLO(f"{SRC}/runs/det_v1/weights/best.pt")
os.makedirs(f"{SRC}/rot", exist_ok=True)

def try_one(img):
    rr = m.predict(img, imgsz=1024, conf=0.25, device=0, verbose=False)[0]
    b = rr.boxes
    if b is None or len(b) == 0: return None
    cls = b.cls.cpu().numpy().astype(int)
    xywh = b.xywh.cpu().numpy()
    cf = b.conf.cpu().numpy()
    if (cls==0).sum() != 1 or (cls==1).sum() != 1: return None
    pR, pL = xywh[cls==0][0], xywh[cls==1][0]
    if abs(pR[0]-pL[0]) < abs(pR[1]-pL[1]): return None      # still stacked vertically
    # the participant's right knee is on the image left
    if pR[0] > pL[0]:
        return None
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
    if (i+1) % 60 == 0: print(f"  {i+1}/{len(fails)}  rescued {len(out)//2}", flush=True)

n = len(out)//2
print(f"\nrescued {n} / {len(fails)}  ({n/len(fails):.1%})")
print(f"  clockwise 90 deg: {saved['90']}   counterclockwise 90 deg: {saved['270']}")
print(f"  still failed: {len(fails)-n}")
with open(f"{SRC}/boxes_rot.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["stem","rot","side","cx","cy","w","h","conf"])
    w.writeheader(); w.writerows(out)
print(f"-> {SRC}/boxes_rot.csv")
