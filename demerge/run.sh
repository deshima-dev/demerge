#!/bin/sh
# 一連のDESHIMA解析を実行する
#
# Filename: run.sh
# (C)2023 内藤システムズ
#
# 実行環境条件
# ------------
#  - このスクリプトはshで動かしてください
#  - /bin/shが利用できること
#  - xargsが利用できること
#  - コマンド「python」でPython3.9が動作すること
#
#
# 使用方法
# --------
# このスクリプトの第一引数にobsid(観測ID)を与えてください。
#  $ ./run.sh -d <<データ保存先へのパス>> <<obsid>>
#
#
# 使用例
# ------
#  $ ./run.sh -d data 20171103184836
#
#
# 指定可能なオプション
# --------------------
# 観測IDは必須な引数ですが、これ以外にも指定できるオプションがあります。
#  -c キャッシュディレクトリを指定
#  -g グラフディレクトリを指定
#  -d 観測データディレクトリの指定
#  -b DDBファイルの指定
#  -o 出力データディレクトリの指定
#  -m マージオプションの指定
#  -p プロット実行オプションの指定
#

DEMERGE="$(cd "$(dirname "${BASH_SOURCE:-$0}")" && pwd)"
DEFAULT_DDB="${DEMERGE}/ddb_20231123.fits.gz"
NCPU=`python -c "import multiprocessing as m; print(m.cpu_count() - 1);"`

# ==================
# 指定可能オプション
# ==================
#  -c キャッシュディレクトリを指定
#  -g グラフディレクトリを指定
#  -d 観測データディレクトリの指定
#  -b DDBファイルの指定
#  -o 出力データディレクトリの指定
#  -m マージオプションの指定
#  -p プロット実行オプションの指定
#
while getopts c:g:d:b:o:m:p: OPT
do
    case $OPT in
    "c") CACHE_DIR="${OPTARG}";;
    "g") GRAPH_DIR="${OPTARG}";;
    "d") DATA_DIR="${OPTARG}";;
    "b") DDB_FILE="${OPTARG}";;
    "o") OUT_DIR="${OPTARG}";;
    "m") MERGE_OPTS="${OPTARG}";;
    "p") PLOT="${OPTARG}";;
    esac
done
shift $((OPTIND - 1))

OBSID=$1
if [ -z $OBSID ]; then
    echo "観測IDを指定してください。" 1>&2
    exit 1
fi

# オプションの規定値を設定
if [ -z "$CACHE_DIR" ]; then
    CACHE_DIR="cache" # 一時ファイルの場所の規定値
fi
if [ -z "$GRAPH_DIR" ]; then
    GRAPH_DIR="graph" # 作成したグラフを格納する場所の規定値
fi
if [ -z "$DATA_DIR" ]; then
    DATA_DIR="data" # 観測データの場所の規定値
fi
if [ -z "$DDB_FILE" ]; then
    DDB_FILE="${DEFAULT_DDB}" # DDBファイルの既定値
fi
if [ -z "$OUT_DIR" ]; then
    OUT_DIR="${CACHE_DIR}" # 出力ディレクトリの規定値
fi
if [ -z "$MERGE_OPTS" ]; then
    MERGE_OPTS="" # マージオプションの規定値
fi

# キャッシュやグラフを格納するディレクトリを作成する
if [ ! -d ${CACHE_DIR}/${OBSID} ]; then
    mkdir -p ${CACHE_DIR}/${OBSID}
    if [ $? -ne 0 ]; then
    exit 1
    fi
fi
if [ ! -d ${GRAPH_DIR}/${OBSID} ]; then
    mkdir -p ${GRAPH_DIR}/${OBSID}
    if [ $? -ne 0 ]; then
    exit 1
    fi
fi
if [ ! -d ${OUT_DIR}/${OBSID} ]; then
    mkdir -p ${OUT_DIR}/${OBSID}
    if [ $? -ne 0 ]; then
    exit 1
    fi
fi

START_TIME=`/bin/date +%s`

# 圧縮でも非圧縮でもTODファイルを扱えるようにする
TOD_FILE=""
if [ -f "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.fits.gz" ]; then
    TOD_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.fits.gz"
else
    TOD_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.fits"
fi

make_divided_data                                     \
       "${DATA_DIR}/cosmos_${OBSID}/kids.list"        \
       "${DATA_DIR}/cosmos_${OBSID}/localsweep.sweep" \
       "${TOD_FILE}"                                  \
       "${CACHE_DIR}/${OBSID}"
if [ $? -ne 0 ]; then
    echo "失敗:make_divided_data"
    exit 1
fi

ls ${CACHE_DIR}/${OBSID}/*.pkl | xargs -P${NCPU} -n1 calc_resonance_params
if [ $? -ne 0 ]; then
    echo "失敗:calc_resonance_params"
    exit 1
fi

rm -f "${CACHE_DIR}/${OBSID}/reduced_${OBSID}.fits"
if [ $? -ne 0 ]; then
    exit 1
fi

make_reduced_fits \
       "${CACHE_DIR}/${OBSID}" \
       "${CACHE_DIR}/${OBSID}/reduced_${OBSID}.fits"
if [ $? -ne 0 ]; then
    echo "失敗:make_reduced_fits"
    exit 1
fi

# ファイルの存在を確認する
MISTI_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.misti"
CABIN_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.cabin"
SKYCHOP_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.skychopper.dat.xz"
if [ ! -f $MISTI_FILE ]; then
    MISTI_FILE=""
fi
if [ ! -f $CABIN_FILE ]; then
    CABIN_FILE=""
fi
if [ ! -f $SKYCHOP_FILE ]; then
    SKYCHOP_FILE="${DATA_DIR}/cosmos_${OBSID}/${OBSID}.skychop"
    if [ ! -f $SKYCHOP_FILE ]; then
        SKYCHOP_FILE=""
    fi
fi

#
# 引数(上から順に)
# ================
# DDBファイルへの相対パス
# Master-to-KID ID対応ファイルへの相対パス
# reduced readoutファイルへの相対パス
# obsファイルへの相対パス
# antennaファイルへの相対パス
# skychopファイルへの相対パス
# weatherファイルへの相対パス
# mistiファイルへの相対パス
# cabinファイルへの相対パス
# その他のmergeに関するオプション
# TODデータとAntennaログの時刻のずれの補正値(ms)
# 出力するZarrファイルへの相対パス
#
merge_to_dems                                                \
    --ddb     "${DDB_FILE}"                                  \
    --corresp "${DATA_DIR}/cosmos_${OBSID}/kid_corresp.json" \
    --readout "${CACHE_DIR}/${OBSID}/reduced_${OBSID}.fits"  \
    --obs     "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.obs"     \
    --antenna "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.ant"     \
    --skychop "${SKYCHOP_FILE}"                              \
    --weather "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.wea"     \
    --misti   "${MISTI_FILE}"                                \
    --cabin   "${CABIN_FILE}"                                \
    --offset_time_antenna 20                                 \
    ${MERGE_OPTS}                                            \
    "${OUT_DIR}/${OBSID}/dems_${OBSID}.zarr.zip"

if [ $? -ne 0 ]; then
    echo "失敗:merge_to_dems"
    exit 1
fi

#
# 説明
# ====
# plot_sweep.pyには2つの引数を渡す必要がある。
# 引数はそれぞれキャッシュファイル名とグラフを格納するディレクトリ名。
# グラフを格納するディレクトリ名はすべての同じ。
# この「キャッシュファイル名」と「グラフを格納するディレクトリ名」のペアの一覧を作り、
# 一旦シェル変数に格納する。
# これをxargsにパイプで渡す。
# xargsには各コマンドの引数が2個であることを示す「-n2」オプションをつける。
#
if [ -n "$PLOT" ]; then
    FILENAMES=""
    for FILENAME in `ls ${CACHE_DIR}/${OBSID}/*.pkl`
    do
        FILENAMES="${FILENAME} ${GRAPH_DIR}/${OBSID} ${FILENAMES}"
    done
    echo $FILENAMES | xargs -P${NCPU} -n2 plot
    if [ $? -ne 0 ]; then
        echo "失敗:plot"
        exit 1
    fi
fi

END_TIME=`/bin/date +%s`
RUN_TIME=`expr ${END_TIME} - ${START_TIME}`
echo "実行時間: ${RUN_TIME}秒"
