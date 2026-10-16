"""全样本金属植入物筛查。

动机:低置信度救援片的目视表里,明显有一批全膝置换/螺钉/钢板。如果金属只
富集在"救援"这一组,那先前把它们排除在外 = 系统性地剔掉了关节置换患者
—— 这个人群与重度 OA、以及"有条件做手术"(医疗可及性)双向相关,非随机缺失。

判据:金属在 X 光上是饱和平台,而骨的最亮处只是零星像素。因此用
frac_top = 该图中 >= 0.995*max 的像素占比。正常膝该值极小(孤立峰),
金属会形成成片高值。阈值不预设,先看两组分布再定。
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
        w = csv.writer(fh); w.writerow(["roi","frac_top","p999_ratio"]); w.writerows(res)
    grp = {"lowconf":[], "rot":[], "main":[]}
    for name, ft, _ in res:
        stem = name.rsplit("_",1)[0]
        grp["lowconf" if stem in low else "rot" if stem in rot else "main"].append(ft)
    print(f"{'组':<10}{'n':>7}{'中位':>12}{'p90':>12}{'p99':>12}{'>1e-3 占比':>12}")
    for k,v in grp.items():
        v = np.array(v)
        print(f"{k:<10}{len(v):>7}{np.median(v):>12.2e}{np.percentile(v,90):>12.2e}"
              f"{np.percentile(v,99):>12.2e}{(v>1e-3).mean():>12.1%}")
