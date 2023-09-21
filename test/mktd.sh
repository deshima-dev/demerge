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

#
# find R, skyで利用するダミーデータ
# calibrate_to_power()関数で計算すると時間とともに線形に増加するTsignalを出力するためのダミーデータ
# Tsignalは0から300まで増加する
#
PARAMS="
     --p0                1.0                \
     --etaf              1.0                \
     --T0                0.0                \
     --Qr                0.25               \
     --linyfc            0.0                \
     --lower_cabin_temp -273.15             \
     --linear_readout    inc                \
     --all_grad          1                  \
     --prefix            testdata_linear_inc
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
$CMD dfits   $PARAMS

#
# find R, skyで利用するダミーデータ
# calibrate_to_power()関数で計算すると時間とともに線形に減少するTsignalを出力するためのダミーデータ
# Tsignalは300から0まで減少する
#
PARAMS="
     --p0                1.0                \
     --etaf              1.0                \
     --T0                0.0                \
     --Qr                0.25               \
     --linyfc            0.0                \
     --lower_cabin_temp -273.15             \
     --linear_readout    dec                \
     --all_grad          1                  \
     --prefix            testdata_linear_dec
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
$CMD dfits   $PARAMS

#
# shuttle機能をテストするためのダミーデータ
# lonの値が時間とともに線形に0から180まで変化する
#
PARAMS="
     --p0                1.0                    \
     --etaf              1.0                    \
     --T0                0.0                    \
     --Qr                0.25                   \
     --linyfc            0.0                    \
     --lower_cabin_temp -273.15                 \
     --linear_antenna    1                      \
     --all_grad          1                      \
     --prefix            testdata_linear_antenna
"
$CMD antenna $PARAMS


date
