# measure_halfwidth.py

Sol'Exを用いて撮影した画像から、簡易的に半値幅を求めグラフ化するスクリプト。

## 原理

- 白色LEDの様に波長にむらのない連続した光源が望ましい。
- 任意のナローバンドフィルターを Sol'Ex にセットし一番輝度の高い場所が中央になるようフレーミングし fits 形式でキャプチャする。
- SharpCap の Live Stacking などで十分な光量となるようにキャプチャし保存する。画像はRAWで画像する。ストレッチなど加工は行わない。
- 撮影した画像の左右中心をY方向に1ピクセルずつ幅 `--width` の平均輝度を求める。
- Y座標 x `--angstrom-per-pixel`をY座標0からの相対波長とする。
- 最大輝度のピクセルを100%の透過率としその50%のとなるピクセルの最小波長と最大波長の幅が半値幅。

## Usage example

```bash
$ python measure_halfwidth.py input.fits -o output.csv -w 2 -a 0.052 -c "Ha half Width"
```

## Options

```python
input.fits            入力FITSファイル（2次元画像）
--out, -o PATH        出力CSVファイル名（省略時は input.csv）
--width, -w N         中心から左右Nピクセルずつを平均（0なら中心列のみ）
--angstrom-per-pixel, -a F  1ピクセルあたりの波長幅（Å, float）。指定すると波長スケールでPNGプロットも出力
--caption, -c TEXT    プロットのタイトル（キャプション）を任意の文字列で指定
```

## 必要なパッケージ
- astropy
- numpy
- matplotlib

### Install

```bash
$ pip install -r requirements.txt
```

## 動作の確認状況

- SVBONY 7nm SIIフィルターでのみ確認している。光源や白黒モノクロカメラも特定の種類のみでの確認となる。
- 他の環境や対象で正しく計測されない可能性はある。
