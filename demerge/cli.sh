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
#  $ ./run.sh -d data/cosmos 20171103184836
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
#

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
#
while getopts c:g:d:b:o:m: OPT
do
    case $OPT in
	"c") CACHE_DIR="${OPTARG}";;
	"g") GRAPH_DIR="${OPTARG}";;
	"d") DATA_DIR="${OPTARG}";;
	"b") DDB_FILE="${OPTARG}";;
	"o") OUT_DIR="${OPTARG}";;
    "m") MERGE_OPTS="${OPTARG}";;
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
    DATA_DIR="data/cosmos" # 観測データの場所の規定値
fi
if [ -z "$DDB_FILE" ]; then
    DDB_FILE="data/ddb/ddb_20180619.fits.gz" # DDBファイルの規定値
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

START_TIME=`date +%s`

make_divided_data                                     \
       "${DATA_DIR}/cosmos_${OBSID}/kids.list"        \
       "${DATA_DIR}/cosmos_${OBSID}/localsweep.sweep" \
       "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.fits.gz" \
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

#
# 引数(上から順に)
# ================
# caldb fitsファイルへの相対パス
# reduced fitsファイルへの相対パス
# obsファイルへの相対パス
# antファイルへの相対パス
# skychopファイルへの相対パス
# weaファイルへの相対パス
# mistiファイルへの相対パス
# cabinファイルへの相対パス
# 出力するZarrファイルへの相対パス
#
merge_to_dems                                                \
    --ddb     "${DDB_FILE}"                                  \
    --readout "${CACHE_DIR}/${OBSID}/reduced_${OBSID}.fits"  \
    --obs     "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.obs"     \
    --antenna "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.ant"     \
    --skychop "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.skychop" \
    --weather "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.wea"     \
    --misti   "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.misti"   \
    --cabin   "${DATA_DIR}/cosmos_${OBSID}/${OBSID}.cabin"   \
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

END_TIME=`date +%s`
RUN_TIME=`expr ${END_TIME} - ${START_TIME}`

echo "実行時間: ${RUN_TIME}秒"
