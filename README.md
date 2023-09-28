# dmerge
DFITSをASTE生データからマージするスクリプト

# 動作環境
 - python 3.9.10

| モジュール   | バージョン|
|--------------|-----------|
| numpy        | 1.24.4    |
| pandas       | 2.0.3     |
| astropy      | 5.2.2     |
| pyreadline   | 2.1       |
| matplotlib   | 3.7.2     |
| scipy        | 1.10.1    |
| lmfit        | 1.2.1     |
| sympy        | 1.12      |
| pyyaml       | 6.0       |
| dems         | 0.4.0     |

参考:https://github.com/deshima-dev/dems

# 実行例

```
$ pwd
/home/hoge/dmerge
$ ./run.sh -d data/deshima2.0 20171103184436
```

# テストの実行方法

```
$ pwd
/home/hoge/dmerge
$ cd test
$ ./mktd.sh                    # テストに必要なダミーデータを生成します
$ python test_merge_to_dems.py # merge_to_dems()関数のテストを実行します
```
