#!/bin/sh
#
# プロセスを分けてテストデータを生成する
#
CMD="python tdmaker.py"

DATA_NAMES="readout ddb dfits antenna skychop weather misti cabin"

date
NCPU=`python -c "import multiprocessing as m; print(m.cpu_count() - 1);"`
echo $DATA_NAMES | xargs -P${NCPU} -n1 $CMD
if [ $? -ne 0 ]; then
    echo "失敗:${CMD}"
    exit 1
fi
date

PARAMS="
     --p0                1.0 \
     --etaf              1.0 \
     --T0                0.0 \
     --Qr                0.25 \
     --linyfc            0.0 \
     --lower_cabin_temp -273.15 \
     --linear_readout    inc \
     --all_grad          1 \
     --prefix            testdata_linear_inc
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS

PARAMS="
     --p0                1.0 \
     --etaf              1.0 \
     --T0                0.0 \
     --Qr                0.25 \
     --linyfc            0.0 \
     --lower_cabin_temp -273.15 \
     --linear_readout    dec \
     --all_grad          1 \
     --prefix            testdata_linear_dec
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
