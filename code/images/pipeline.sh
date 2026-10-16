#!/bin/bash
# Preprocessing: PNG conversion -> knee detection -> ROI cropping. Existing outputs are skipped, so the script can be rerun.
PY=${PY:-python}
cd /home/data/knee_yolo/work || exit 1

echo "=== [1/3] waiting for PNG conversion to finish $(date '+%F %T') ==="
while pgrep -f convert_all.py >/dev/null; do sleep 30; done
echo "  PNG: $(ls /home/data/knee_yolo/images/*.png | wc -l) / 4908"

echo "=== [2/3] detection on all films $(date '+%F %T') ==="
$PY /tmp/detect_all.py 2>&1 | grep -viE "^$"

echo "=== [3/3] ROI cropping $(date '+%F %T') ==="
$PY /tmp/crop_rois.py 2>&1 | grep -viE "^$"

echo "=== done $(date '+%F %T') ==="
echo "ROI count: $(ls /home/data/knee_roi 2>/dev/null | wc -l)"
du -sh /home/data/knee_roi 2>/dev/null
