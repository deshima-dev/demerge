# dmerge
DFITSをASTE生データからマージするスクリプト

# 動作環境
 - Linux Kernel 5.4.0-153-generic
 - python 3.8.10

| モジュール   | バージョン|
|--------------|-----------|
| numpy        | 1.24.4    |
| pandas       | 2.0.3     |
| astropy      | 5.2.2     |
| pyreadline   | 2.1       |
| matplotlib   | 3.7.2     |
| scipy        | 1.10.1    |
| lmfit        | 1.2.1     |
| ad           | 1.3.2     |
| sympy        | 1.12      |
| pyyaml       | 6.0       |
| dems         | 0.2.0     |

# demsのインストール方法

```
$ pip install -q "dems==0.2.0"
```

参考:!https://github.com/deshima-dev/dems

# 実行例

```
$ ./run.sh -d data/deshima1.0 20171110114116
```
