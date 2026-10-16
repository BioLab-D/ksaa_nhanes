#!/bin/bash
# 用与五折相同的配方,五个骨干各用全部 4,875 人训练一个最终模型并存档(内部验证 10% 只作早停,无留出折),供新片推理。
# 稿件结果来自五折 OOF,不用本队列。部署:单模型预测 = 两膝局部头期望年龄均值;集成权重与偏倚校正系数用五折 OOF 标定的值,
#   不能用最终模型的样本内预测重新标定(会低估偏倚)。每个配置无存档时自动重试一次。用法:bash train_full.sh
export TORCH_HOME=/home/data/torch PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True N_VIEWS=2 DROP=0 NO_ALBUMENTATIONS_UPDATE=1
PY=${PY:-python}; cd /tmp/mammo_age_official; D=/home/data/knee_yolo/v2; W=$D/weights; mkdir -p $W
F='grep --line-buffered -vE check_version|fetch_version|validated_self|UserWarning|data.=.fetch'
C="--bs 16 --seed 20260831 --nblock 4 --img 448 --full --save_dir $W"
BASE="--aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15"
Q=$D/p5_queue.log
run(){ t=$1; shift; L=$D/$t.log
  if [ -s $W/${t}_full.pt ]; then echo "$(date '+%F %T')  跳过 $t(存档已在)" >> $Q; return; fi
  for try in 1 2; do
    echo "=== $t 第 $try 次 $(date '+%F %T')  $* ===" >> $L
    $PY knee_v2.py --tag $t $C $BASE "$@" 2>&1 | $F >> $L
    if [ -s $W/${t}_full.pt ] && grep -q "full done" $L; then
      echo "$(date '+%F %T')  $t  $(tail -n 1 $D/${t}_folds.csv | tr -d '\r')  存档 $(du -h $W/${t}_full.pt | cut -f1)" >> $Q; return
    fi
    echo "$(date '+%F %T')  $t 第 $try 次没有存档(见 $t.log)" >> $Q; sleep 90
  done
}
echo "=== v2 最终模型存档队列开始 $(date '+%F %T') ===" > $Q
run m5_full_resnet18       --arch resnet18
run m5_full_efficientnetb0 --arch efficientnet_b0
run m5_full_resnet50       --arch resnet50
run m5_full_densenet121    --arch densenet121
run m5_full_convnexttiny   --arch convnext_tiny
echo "=== 完成 $(date '+%F %T') ===" >> $Q
