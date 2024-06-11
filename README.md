# de:merge

[![Release](https://img.shields.io/pypi/v/demerge?label=Release&color=cornflowerblue&style=flat-square)](https://pypi.org/project/demerge/)
[![Python](https://img.shields.io/pypi/pyversions/demerge?label=Python&color=cornflowerblue&style=flat-square)](https://pypi.org/project/demerge/)
[![Downloads](https://img.shields.io/pypi/dm/demerge?label=Downloads&color=cornflowerblue&style=flat-square)](https://pepy.tech/project/demerge)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.10015892-cornflowerblue?style=flat-square)](https://doi.org/10.5281/zenodo.10015892)
[![Tests](https://img.shields.io/github/actions/workflow/status/deshima-dev/demerge/tests.yaml?label=Tests&style=flat-square)](https://github.com/deshima-dev/demerge/actions)

DESHIMA merge code for observed datasets

MergeToDfits()クラスの代わりとなるmerge_to_dems()関数を実装し、DEMSオブジェクトを生成できるようにしました。今回の更新では、付属のPythonパッケージ`demerge`をPython環境にインストールすることで、同名のコマンド`demerge`が利用可能となります（`demerge/run.sh`のラッパースクリプトです）。これを実行するとmerge_to_dems()関数が実行され、解析結果がまとまったZarrファイルが生成されます。DEMSの構造については[deshima-dev/dems](https://github.com/deshima-dev/dems)をご覧ください。

## 動作環境

- CPython: 3.9-3.12
- Dependent packages: [pyproject.toml](https://github.com/deshima-dev/demerge/blob/v2.13.0/pyproject.toml)をご確認ください

## インストール

### 一般ユーザ

pipで[PyPI](https://pypi.org)からインストールしてください。

```shell
pip install demerge==2.13.0
```

インストール後、`demerge`コマンドが利用可能になります。

### 開発者

以下のコマンドで実行環境にダウンロード（クローン）してください。

```shell
# テストデータを含める場合
git clone --recursive https://github.com/deshima-dev/demerge.git
cd demerge && git checkout v2.13.0

# 最小構成でダウンロードする場合
git clone --depth=1 https://github.com/deshima-dev/demerge.git
cd demerge && git checkout v2.13.0
```

続けて、`demerge`パッケージのPython環境へのインストールを行います。

```shell
cd /path/to/demerge

# pipでインストールする場合
pip install -e .

# Poetryでインストールする場合
poetry install
```

インストール後、`demerge`コマンドが利用可能になります。

## 解析の実行例

以下のコマンドでOBSID=20171103184436の解析を行います。観測データはdataに格納されているものとします。

```shell
demerge -d data 20171103184436
```

以下のコマンドでテスト用のデータをダウンロードすることができます。

```shell
deshima-rawdata download 20171103184436 -d data -e
```

どのデータがダウンロードできるかは、データリストをご確認ください。

```shell
deshima-rawdata list
```

### demergeコマンドの引数

以下の引数を指定すると、データの取得ディレクトリや保存先などを変更することができます。

引数 | 説明 | デフォルト
--- | --- | ---
`-d` | 観測データディレクトリの指定 | `data`
`-c` | キャッシュディレクトリを指定 | `cache`
`-g` | グラフディレクトリを指定 | `graph`
`-b` | DDBファイルを指定 | `dmerge/ddb_20231123.fits.gz`（パッケージに同梱）
`-i` | Master-to-KID ID対応ファイルを指定 | なし（なしの場合、DDBファイル中の対応関係を使用）
`-o` | 出力データディレクトリを指定 | `cache`
`-m` | マージオプションを指定（例：`-m "--coordinate radec"`） | なし

### 観測データの保存先の構造

```
data
├── cosmos_20171103184436
│   ├── 20171103184436_info.txt
│   ├── 20171103184436.ant
│   ├── 20171103184436.cabin
│   ├── 20171103184436.fits.gz
│   ├── 20171103184436.misti
│   ├── 20171103184436.obs
│   ├── 20171103184436.roomchop
│   ├── 20171103184436.skychop
│   ├── 20171103184436.tsky
│   ├── 20171103184436.wea
│   ├── kids.list
│   ├── localsweep_info.txt
│   ├── localsweep.sweep
│   └── pretune
│       ├── kids.list
│       ├── localsweep_info.txt
│       └── localsweep.sweep

((( 続く )))
```

### 解析結果の保存場所

`demerge`コマンドを実行すると、デフォルトでは直下にcacheとgraphというディレクトリが作成されます。これらのディレクトリの中にはさらにobsidのディレクトリが作成され、そこに解析結果が保存されます。最終結果のZarrファイルもcache/dems_((OBSID))/((OBSID)).zarr.zipというパスに保存されます。

以下にOBSID=20171103184436を解析した時のcache内のファイル構造を示します。
```
cache
├── 20171103184436
│   ├── dems_20171103184436.zarr.zip <--- これがmerge_to_dems.pyで生成されるファイル
│   ├── kid00000.pkl
│   ├── kid00001.pkl
│   ├── kid00002.pkl
│   ├── kid00003.pkl
│   ├── kid00004.pkl
│   ├── kid00005.pkl

(((中略)))

│   ├── kid00057.pkl
│   ├── kid00058.pkl
│   ├── kid00059.pkl
│   ├── kid00060.pkl
│   ├── kid00061.pkl
│   ├── kid00062.pkl
│   └── reduced_20171103184436.fits

(((続く)))
```

以下にOBSID=20171103184436を解析した時のgraph内のファイル構造を示します。

```
graph
├── 20171103184436
│   ├── sweep_kid00000.png
│   ├── sweep_kid00001.png
│   ├── sweep_kid00002.png
│   ├── sweep_kid00003.png
│   ├── sweep_kid00004.png
│   ├── sweep_kid00005.png

(((中略)))

│   ├── sweep_kid00059.png
│   ├── sweep_kid00060.png
│   ├── sweep_kid00061.png
│   ├── sweep_kid00062.png
│   ├── tod_kid00000.png
│   ├── tod_kid00001.png
│   ├── tod_kid00002.png
│   ├── tod_kid00003.png
│   ├── tod_kid00004.png
│   ├── tod_kid00005.png

(((中略)))

│   ├── tod_kid00059.png
│   ├── tod_kid00060.png
│   ├── tod_kid00061.png
│   └── tod_kid00062.png

(((続く)))
```

## テストの実行方法

```shell
cd /path/to/demerge/test
./mktd.sh                # テストに必要なダミーデータを生成します
python -m unittest       # merge_to_dems関数およびmerge_functionモジュールのテストを実行します
```

テストを実行するにはダミーデータとテストデータが必要になります。ダミーデータはmktd.shスクリプトを実行することで生成されます。テストデータはdata/cosmos_20171103184436/の中のものを使いますので、上記の手順であらかじめダウンロードしてください。

## merge_to_dems.pyについて

merge_to_dems()関数が定義されています。この関数は以下の8つのファイルからDEMSオブジェクトを生成します。

 - DDBファイル(.fits.gz)
 - Master-to-KID ID対応ファイル(.json)
 - obsファイル(.obs)
 - antennaファイル(.ant)
 - skychopファイル(.skychop)
 - weatherファイル(.wea)
 - mistiファイル(.misti)
 - cabinファイル(.cabin)
 - reduced readoutファイル(.fits)

reduced readoutファイルはdemerge/run.shを実行することによって生成されます。DDBファイルはdemergeのリポジトリに含まれているキャリブレーションデータです。その他のファイルは観測とともに得られるデータです。

merge_to_dems()関数は以下の必須引数とオプション引数をとることができます。

### 必須引数

引数名 | 型 | 説明
--- | --- | ---
`filename` | 文字列 | 出力ファイルへのパスを指定して下さい(.zarr.zip)
`--ddb` | 文字列 | DDBファイルへのパスを指定して下さい(.fits.gz)
`--corresp` | 文字列 | Master-to-KID ID対応ファイルへのパスを指定して下さい(.json)
`--obs` | 文字列 | obsファイルへのパスを指定して下さい(.obs)
`--antenna` | 文字列 | antennaファイルへのパスを指定して下さい(.antenna)
`--readout` | 文字列 | reduced readoutファイルへのパスを指定して下さい(.fits)
`--skychop` | 文字列 | skychopファイルへのパスを指定して下さい(.skychop)
`--weather` | 文字列 | weatherファイルへのパスを指定して下さい(.weather)
`--misti` | 文字列 | mistiファイルへのパスを指定して下さい(.misti)
`--cabin` | 文字列 | cabinファイルへのパスを指定して下さい(.cabin)

### オプション引数

引数名 | 型 | 既定値 | 説明
--- | --- | --- |---
`--coordinate` | 文字列 | `'azel'` | 座標系(azel/radec)を文字列で指定します
`--loadtype` | 文字列 | `'fshift'` | 読み込むデータを文字列で指定します
`--findR` | フラグ | なし | 指定するとFindR, Skyを実行します
`--ch` | 整数 | `0` | findRに利用するチャネルを整数で指定します
`--Rth` | 実数 | `280.0` | R閾値を実数で指定します
`--skyth` | 実数 | `150.0` | sky閾値を実数で指定します
`--cutnum` | 整数 | `1` | findRでのカット数を整数で指定します
`--still` | フラグ | なし | 指定するとstill観測用の解析を行います
`--period` | 整数 | `2` | still観測の1/2周期(秒)を整数で指定します
`--shuttle` | フラグ | なし | 指定するとshuttle観測用の解析を行います
`--lon_min_off` | 実数 | `0.0` | shuttle観測時のOFFにするlongitudeの最小値を実数で指定します
`--lon_max_off` | 実数 | `0.0` | shuttle観測時のOFFにするlongitudeの最大値を実数で指定します
`--lon_min_on` | 実数 | `0.0` | shuttle観測時のONにするlongitudeの最小値を実数で指定します
`--lon_max_on` | 実数 | `0.0` | shuttle観測時のONにするlongitudeの最大値を実数で指定します
`--offset_time_antenna` | 整数 | `0` | TODデータとAntennaログの時刻のずれの補正値(ms)を整数で指定します
`--debug` | フラグ | なし | 指定すると全ての引数の値をログとして表示します

`--loadtype`はTsignal(空の輝度温度)とfshift(MKIDのdf/f)を指定できます。demerge/run.shでは上記オプション引数を文字列として`-m "--coordinate radec --findR --debug"`のように指定することができます。

## tdmaker.pyについて

tdmaker.pyにはTestDataMakerクラスが定義されています。このクラスを利用するとmerge_to_dems()関数のテストに必要なダミーデータを生成することができます。以下のファイルを個別に生成することができます。

ファイルの種類 | 周期(秒) | 拡張子
--- | --- | ---
ddb | - | .fits.gz
antenna | 0.1 | .ant
readout | 0.00625 | .fits
skychop | 0.001 | .skychop
weather | 10 | .wea
misti | 0.1 | .misti
cabin | 60 | .cabin
dfits | - | .fits.gz

dfitsファイル(.fits.gz) (dfitsファイルはdeshima1.0のdemergeで生成されていたファイルです。これらのファイルの測定開始時刻は同一時刻に設定されます。同じ時刻から始まり各ファイル毎の周期でデータが記録されます。

TestDataMakerクラスはオプションを指定することで様々なデータを生成することができます。以下にtdmaker.pyで指定できる引数を示します。

引数名 | 型 | 既定値 | 説明
--- | --- | --- | ---
`data_name` | 文字列 | - | 個別にデータを生成する場合はデータ名を指定します。指定しない場合は全てのデータが生成されます。
`--time` | 整数 | `3` | 観測時間(分)を整数で指定して下さい
`--p0` | 実数 | `1.0` | p0をfloatで指定して下さい
`--etaf` | 実数 | `0.5` | etafをfloatで指定して下さい
`--T0` | 実数 | `1.0` | T0をfloatで指定して下さい
`--Qr` | 実数 | `1.1` | Qrをfloatで指定して下さい
`--linyfc` | 実数 | `0.25` | linyfcをfloatで指定して下さい
`--linear_readout` | 文字列 | - | readoutの値を線形に変化させる場合はinc/decのいずれかを指定して下さい
`--linear_antenna` | 真理値 | `False` | antennaのlonを線形に変化させる場合はTrueを指定して下さい
`--all_grad` | 真理値 | `False` | すべてのSCAN状態をGRADにする場合はTrueを指定して下さい
`--lower_cabin_temp` | 実数 | `15.0` | MainCabinの温度(degC)をfloatで指定して下さい
`--prefix` | 文字列 | `'testdata'` | 生成されるファイル名のprefixを指定して下さい
`--measure_time` | 整数 | - | 環境測定時間(分)を整数で指定して下さい

`data_name`を指定すると指定されたデータだけが生成されます。指定しない場合は全ての種類のデータが生成されます。`--time`は観測時間(分)を指定します。`--measure_time`を指定しなければ観測時間と環境測定時間は同じになります。観測時間とはreadoutの時間であり、環境測定時間はreadout以外のcabinやweatherなどの測定時間になります。通常`--time`だけを指定した場合、または`--time`も`--measure_time`も省略した場合は環境測定時間が観測時間よりも少し長くなるように生成されます。このようなデータではreadoutの時刻に合わせて環境測定時間を補間する場合に外挿処理は発生しません。外挿処理が発生した場合の振る舞いを調べる場合に`--measure_time`オプションで`--time`よりも短い時間を指定します。

`--p0`, `--etaf`, `--T0`, `--Qr`, `--linyfc`はmerge_function.pyで定義されているTlos_model()関数とcalibrate_to_power()関数で計算に利用される係数を設定します。mktd.shでは計算値のTsignalが定数になるように設定したり、0Kから300Kまで3分間で線形に変化するように係数を設定しています。

計算結果のTsignalを線形に変化させるためには`--linear_readout`にincまたはdecの文字列を指定します。incはTsignalが単調増加し、decは単調減少します。 同様にantennaのlongitudeも線形に変化させることができます。`--linear_antenna`にTrueを設定するとantennaのlongitudeが線形に増加します。

`--all_grad`をTrueにするとSCAN状態が全てGRADのantennaデータができます。省略(Falseを設定すると)すると観測時間の半分でGRADからONに切り替わるデータとなります。`--lower_cabin_temp`はMainCabinの温度を指定します。`--prefix`には生成するデータの接頭語を指定します。用途に合わせてわかりやすい名前を付けます。

TestDataMakerクラスの詳細な使い方はmktd.shスクリプトをご覧下さい。
