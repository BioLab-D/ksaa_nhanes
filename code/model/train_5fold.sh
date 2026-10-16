#!/bin/bash
# 五折交叉验证训练 Mammo-AGE 多视图膝龄模型(N = 4 融合层,448 ROI),输出逐折内部验证 MAE 与 OOF 预测。
# 训练配方按膝片适配:膝片增强 + 逐图 P0.5–P99.5 浮点归一化 + AdamW 余弦退火(wd 1e-4)。
# 第 1 段为稿件所用的五骨干;第 2、3 段为消融,稿件未报告。用法:bash train_5fold.sh
export TORCH_HOME=/home/data/torch PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True N_VIEWS=2 DROP=0 NO_ALBUMENTATIONS_UPDATE=1
PY=${PY:-python}; cd /tmp/mammo_age_official; D=${D:-/home/data/knee_yolo/v2}; mkdir -p $D
F='grep --line-buffered -vE check_version|fetch_version|validated_self|UserWarning|data.=.fetch'
C="--bs 16 --seed 20260831 --folds 0,1,2,3,4 --nblock 4"
BASE=${BASE:-"--aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15"}
Q=$D/p4_queue.log
run(){ t=$1; shift; L=$D/$t.log
  [ -s $D/${t}_folds.csv ] && [ "$(($(wc -l < $D/${t}_folds.csv)-1))" -ge 5 ] && { echo "$(date '+%F %T')  跳过 $t(已满五折)" >> $Q; return; }
  echo "=== $t $(date '+%F %T')  $* ===" > $L
  $PY knee_v2.py --tag $t $C $BASE "$@" 2>&1 | $F >> $L
  echo "=== 结束 $(date '+%F %T') ===" >> $L
  if [ -s $D/${t}_folds.csv ]; then
    v=$(tail -n +2 $D/${t}_folds.csv | cut -d, -f2)
    echo "$(date '+%F %T')  $t  $(echo $v | tr '\n' ' ') | 均值 $(echo "$v" | awk '{s+=$1;q+=$1*$1;n++}END{printf "%.4f ± %.4f", s/n, sqrt((q-s*s/n)/(n-1))}')" >> $Q
  else echo "$(date '+%F %T')  $t 无产出(见 $t.log)" >> $Q; fi
}
echo "=== 照搬 Mammo-AGE 的比较(N=4,五折)开始 $(date '+%F %T')  BASE=$BASE ===" > $Q

# ── 1 五骨干(稿件所用)
echo "$(date '+%F %T')  --- 第 1 段:五骨干 ---" >> $Q
run m4_resnet18       --arch resnet18       --img 448
run m4_resnet50       --arch resnet50       --img 448
run m4_efficientnetb0 --arch efficientnet_b0 --img 448
run m4_convnexttiny   --arch convnext_tiny  --img 448
run m4_densenet121    --arch densenet121    --img 448

# ── 2 模块消融 @ ResNet-18(稿件未报告)
echo "$(date '+%F %T')  --- 第 2 段:模块消融 ---" >> $Q
run m4_ab_wo_citf --arch resnet18 --img 448 --nblock 0     # w/o CI-Transformer
run m4_ab_wo_poe  --arch resnet18 --img 448 --no_poe       # w/o POE
run m4_ab_wo_mvl  --arch resnet18 --img 448 --no_mvl       # w/o MVL
run m4_ab_wo_mt   --arch resnet18 --img 448 --no_density   # w/o Multi-task(无膝部辅助标签)

# ── 3 输入尺寸 @ ResNet-18(稿件未报告;448 见第 1 段)
echo "$(date '+%F %T')  --- 第 3 段:输入尺寸 ---" >> $Q
run m4_sz112 --arch resnet18 --img 112
run m4_sz224 --arch resnet18 --img 224
run m4_sz320 --arch resnet18 --img 320
run m4_sz512 --arch resnet18 --img 512
echo "=== 照搬 Mammo-AGE 的比较完成 $(date '+%F %T') ===" >> $Q
