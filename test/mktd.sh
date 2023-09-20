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

$CMD readout \
     --p0                1.0 \
     --etaf              1.0 \
     --T0                0.0 \
     --linyfc            0.0 \
     --lower_cabin_temp -273.15 \
     --linear_readout    1 \
     --prefix            testdata_linear

