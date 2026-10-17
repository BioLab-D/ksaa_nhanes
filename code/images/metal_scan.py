"""Scan all ROIs for metal implants.

Review of the rescued films showed many knee replacements, screws and plates. Excluding these films would remove
participants with severe OA who had access to surgery, so the scan checks whether metal is concentrated among them.

Metal forms a saturated plateau, whereas the brightest bone gives only scattered pixels. The statistic is
frac_top, the fraction of pixels >= 0.995 x the image maximum; no threshold is fixed in advance.
"""
import os, re, csv, glob, numpy as np
from multiprocessing import Pool
from PIL import Image
ROI, SRC = "/home/data/knee_roi", "/home/data/knee_yolo"
low = {r["stem"] for r in csv.DictReader(open(f"{SRC}/boxes_lowconf.csv"))}
rot = {r["stem"] for r in csv.DictReader(open(f"{SRC}/boxes_rot.csv"))}

def one(f):
    a = np.asarray(Image.open(f)).astype(np.float32)
    mx = a.max()
    if mx <= 0: return None
    return (os.path.basename(f)[:-4], float((a >= 0.995*mx).mean()), float(np.percentile(a, 99.9)/mx))

if __name__ == "__main__":
    fs = sorted(glob.glob(f"{ROI}/*.png"))
    with Pool(12) as p: res = [r for r in p.map(one, fs, chunksize=40) if r]
    with open(f"{SRC}/metal_scan.csv","w",newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["roi","frac_top","p999_ratio"])
        w.writerows(res)
    grp = {"lowconf":[], "rot":[], "main":[]}
    for name, ft, _ in res:
        stem = name.rsplit("_",1)[0]
        grp["lowconf" if stem in low else "rot" if stem in rot else "main"].append(ft)
    print(f"{'group':<10}{'n':>7}{'median':>12}{'p90':>12}{'p99':>12}{'frac >1e-3':>12}")
    for k,v in grp.items():
        v = np.array(v)
        print(f"{k:<10}{len(v):>7}{np.median(v):>12.2e}{np.percentile(v,90):>12.2e}"
              f"{np.percentile(v,99):>12.2e}{(v>1e-3).mean():>12.1%}")
