measure_halfwidth.py

Usage example:

    python measure_halfwidth.py input.fits -o output.csv -w 2

Options:
    input.fits            入力FITSファイル（2次元画像）
    --out, -o PATH        出力CSVファイル名（省略時は input_profile.csv）
    --width, -w N         中心から左右Nピクセルずつを平均（0なら中心列のみ）
    --angstrom-per-pixel, -a F  1ピクセルあたりの波長幅（Å, float）。指定すると波長スケールでPNGプロットも出力
    --caption, -c TEXT    プロットのタイトル（キャプション）を任意の文字列で指定

例:
    python measure_halfwidth.py input.fits -o output.csv -w 2 -a 0.5 -c "My Spectrum"

This will read `input.fits`, average columns center-2..center+2 and write Y,intensity to `output.csv`.

Wavelength plotting:

    python measure_halfwidth.py input.fits -o output.csv -w 2 -a 0.5

ここで `-a 0.5` は1ピクセルあたり0.5 Å を意味します。指定すると、CSV に加えて `output_profile.png` という波長 (Å) 対 輝度 のプロットが保存されます。

Dependencies:
- astropy
- numpy

Install:

    pip install -r requirements.txt

Note: CSV uses header `y,intensity`. Y is 0-indexed row number.
