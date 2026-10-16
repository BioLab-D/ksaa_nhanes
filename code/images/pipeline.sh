#!/bin/bash
# 膝片预处理流水线:转换 -> 检测 -> 裁 ROI。
# 各步可重入(已有产物跳过)。
PY=${PY:-python}
cd /home/data/knee_yolo/work || exit 1

echo "=== [1/3] 等 PNG 转换结束 $(date '+%F %T') ==="
while pgrep -f convert_all.py >/dev/null; do sleep 30; done
echo "  PNG: $(ls /home/data/knee_yolo/images/*.png | wc -l) / 4908"

echo "=== [2/3] 全库检测 $(date '+%F %T') ==="
$PY /tmp/detect_all.py 2>&1 | grep -viE "^$"

echo "=== [3/3] 裁 ROI $(date '+%F %T') ==="
$PY /tmp/crop_rois.py 2>&1 | grep -viE "^$"

echo "=== 完成 $(date '+%F %T') ==="
echo "ROI 数: $(ls /home/data/knee_roi 2>/dev/null | wc -l)"
du -sh /home/data/knee_roi 2>/dev/null
