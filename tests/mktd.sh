#!/bin/sh

#
# 実行スクリプトを指定
#
CMD="python tdmaker.py"

#
# 開始時刻を出力
#
date --iso-8601=seconds

#
# デフォルトで利用するダミーデータ
#
PARAMS="
     --time 3 \
     --measure_time 4
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
$CMD misti   $PARAMS
$CMD skychop $PARAMS

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
$CMD misti   $PARAMS
$CMD skychop $PARAMS

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
$CMD misti   $PARAMS
$CMD skychop $PARAMS

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
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
$CMD misti   $PARAMS
$CMD skychop $PARAMS

#
# 環境測定時間をreadoutよりも短くする
# skychop_stateやstateを補間したときの振る舞いを確認するためのダミーデータ
#
# ここで
# --timeはREADOUT時間、
# --measure_timeはREADOUT以外の環境測定時間、
# --over_timeはmeasure_timeに加算される調整用の時間。
#
PARAMS="
     --p0                1.0                    \
     --etaf              1.0                    \
     --T0                0.0                    \
     --Qr                0.25                   \
     --linyfc            0.0                    \
     --lower_cabin_temp -273.15                 \
     --prefix            testdata_short_measure \
     --time              3                      \
     --measure_time      2                      \
"
$CMD readout $PARAMS
$CMD ddb     $PARAMS
$CMD cabin   $PARAMS
$CMD antenna $PARAMS
$CMD weather $PARAMS
$CMD misti   $PARAMS
$CMD skychop $PARAMS

# 終了時刻を出力
date --iso-8601=seconds
