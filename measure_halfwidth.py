#!/usr/bin/env python3
"""
measure_halfwidth.py

FITS画像の中央の水平座標（列）から垂直方向に輝度を測定し、CSVと（任意で）波長スケールのプロットを保存するスクリプト。

Usage:
    python measure_halfwidth.py INPUT.fits [options]

Options:
    INPUT.fits                      入力FITSファイル（最初に見つかった2Dイメージ拡張を使用）
    --out, -o PATH                  出力CSVファイルパス（省略時は INPUT.csv）
    --width, -w N                   中心から左右Nピクセルずつを平均（0なら中心列のみ、デフォルト0）
    --angstrom-per-pixel, -a F      1ピクセルあたりの波長幅（Å, float）。指定すると波長スケールでPNGプロットを出力
    --caption, -c TEXT              プロットのタイトル（省略可）

CSV:
    1列目: y (行番号, 0-index)
    2列目: intensity (測定輝度, float)

Output PNG (when -a provided):
    CSVのベース名に ".png" を付けて保存します。
    プロットには以下の注釈を自動で追加します:
      - Peak（最大輝度）の水平線とピークマーカー
      - Half max（半値）の水平線
      - 半値での左右波長位置に垂直線
      - 左右垂直線間に両矢尻の矢印でFWHMを表示
      - Peak と Half max の数値注記（指数表記）

Example:
    python measure_halfwidth.py input.fits -o output.csv -w 2 -a 0.5 -c "My Spectrum"

Notes:
    - FWHMはピークの左右で半値に達する位置を線形補間して求めます。
    - 入力FITSに2D画像が見つからない場合はエラーになります。

"""
import argparse
from astropy.io import fits
import numpy as np
import csv
import os
import matplotlib.pyplot as plt


def read_image_from_fits(path):
    with fits.open(path) as hdul:
        # find first HDU with 2D data
        for h in hdul:
            data = h.data
            if data is None:
                continue
            if data.ndim == 2:
                return data.astype(float)
        raise ValueError("No 2D image data found in FITS file")


def vertical_profile(image, width=1):
    h, w = image.shape
    center_x = w // 2
    # determine columns to average
    half = max(0, int(width))
    cols = np.arange(center_x - half, center_x + half + 1)
    # clip to valid range
    cols = cols[(cols >= 0) & (cols < w)]
    profile = np.nanmean(image[:, cols], axis=1)
    return profile


def save_profile_csv(y, profile, outpath):
    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["y", "intensity"])
        for yi, val in zip(y, profile):
            writer.writerow([int(yi), float(val)])


def main():
    p = argparse.ArgumentParser(description="Extract vertical profile from center column of a FITS image and save to CSV")
    p.add_argument("fits_file", help="Input FITS path")
    p.add_argument("--out", "-o", help="Output CSV path")
    p.add_argument("--width", "-w", type=int, default=0, help="Half-width in pixels around center to average (0 means center column only). Use e.g. 2 to average center +/-2 -> 5 columns")
    p.add_argument("--angstrom-per-pixel", "-a", type=float, default=None, help="Wavelength per pixel in Angstrom (float). If provided, a PNG plot (wavelength vs intensity) will be saved.")
    p.add_argument("--caption", "-c", type=str, default=None, help="Optional caption/title to add to the plot")
    args = p.parse_args()

    if not os.path.exists(args.fits_file):
        print(f"FITS file not found: {args.fits_file}")
        raise SystemExit(1)

    image = read_image_from_fits(args.fits_file)
    # width param: if user specified 0 -> single column. We'll convert to half-width used in function
    half_width = args.width
    profile = vertical_profile(image, width=half_width)
    y = np.arange(len(profile))

    out = args.out
    if not out:
        base = os.path.splitext(os.path.basename(args.fits_file))[0]
        out = base + ".csv"

    save_profile_csv(y, profile, out)
    print(f"Saved profile CSV to {out}")

    # If angstrom per pixel provided, create a plot (wavelength on x-axis)
    if args.angstrom_per_pixel is not None:
        ang = float(args.angstrom_per_pixel)
        wavelengths = y * ang
        base = os.path.splitext(os.path.basename(out))[0]
        png_out = base + ".png"
        try:
            fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
            ax.plot(wavelengths, profile, lw=1)
            ax.set_xlabel("Relative Wavelength (Å)")
            ax.set_ylabel("Intensity")

            # set caption/title if provided
            if args.caption:
                ax.set_title(args.caption)

            # Peak and FWHM calculation
            max_val = np.nanmax(profile)
            half_val = max_val / 2.0
            # find indices where profile crosses half_val
            # we will find leftmost and rightmost crossing around the peak
            peak_idx = int(np.nanargmax(profile))
            # left side
            left_idx = None
            for i in range(peak_idx, 0, -1):
                if profile[i] <= half_val and profile[i-1] > half_val:
                    # linear interpolation for better precision
                    x1, x2 = wavelengths[i], wavelengths[i-1]
                    y1, y2 = profile[i], profile[i-1]
                    t = (half_val - y1) / (y2 - y1)
                    left_wl = x1 + t * (x2 - x1)
                    left_idx = left_wl
                    break
            # right side
            right_idx = None
            for i in range(peak_idx, len(profile)-1):
                if profile[i] >= half_val and profile[i+1] < half_val:
                    x1, x2 = wavelengths[i], wavelengths[i+1]
                    y1, y2 = profile[i], profile[i+1]
                    if (y2 - y1) == 0:
                        t = 0
                    else:
                        t = (half_val - y1) / (y2 - y1)
                    right_wl = x1 + t * (x2 - x1)
                    right_idx = right_wl
                    break

            # plot horizontal half-max line
            ax.hlines(half_val, wavelengths[0], wavelengths[-1], colors='C3', linestyles='--', label='Half max')
            # plot horizontal line at peak and mark peak
            peak_wavelength = wavelengths[peak_idx]
            ax.hlines(max_val, wavelengths[0], wavelengths[-1], colors='C4', linestyles='-.', label='Peak')
            ax.plot(peak_wavelength, max_val, marker='o', color='C4')
            # annotate numeric values for peak and half-max near right edge (scientific notation)
            text_x = wavelengths[-1]
            # place labels slightly left from right edge
            label_x = text_x - 0.02 * (wavelengths[-1] - wavelengths[0])
            ax.text(label_x, max_val, f"{max_val:.3e}", ha='right', va='bottom', color='C4', fontsize=9)
            ax.text(label_x, half_val, f"{half_val:.3e}", ha='right', va='bottom', color='C3', fontsize=9)

            # format y-axis (Intensity) to scientific notation similar to labels
            from matplotlib.ticker import ScalarFormatter
            class OOMFormatter(ScalarFormatter):
                def __init__(self, order=0, fmt="%1.1f", useOffset=True, useMathText=False):
                    self._order = order
                    ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
                def _set_order_of_magnitude(self, range):
                    self.orderOfMagnitude = self._order
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.3e}"))

            # plot vertical lines and arrow if both crossings found
            if left_idx is not None and right_idx is not None:
                ax.vlines([left_idx, right_idx], ymin=0, ymax=half_val, colors='C2', linestyles=':')
                # draw double-headed arrow between left_idx and right_idx at y = 0.9*half_val
                arrow_y = 0.9 * half_val
                ax.annotate('', xy=(right_idx, arrow_y), xytext=(left_idx, arrow_y),
                            arrowprops=dict(arrowstyle='<->', color='k'))
                fwhm_val = right_idx - left_idx
                ax.text(0.5*(left_idx+right_idx), arrow_y, f" FWHM = {fwhm_val:.3f} Å", ha='center', va='bottom')
            else:
                ax.text(0.5*(wavelengths[0]+wavelengths[-1]), 0.9*max_val, "FWHM not found", ha='center', va='bottom', color='r')

            ax.legend()
            # --- Show image strip from FITS above the plot ---
            try:
                # extract vertical strip from original image: center_x-10 .. center_x+9
                h_img, w_img = image.shape
                cx = w_img // 2
                left = max(0, cx - 10)
                right = min(w_img, cx + 10)  # exclusive
                strip = image[:, left:right]  # shape (ny, 20)
                # rotate counter-clockwise 90 deg
                strip_rot = np.rot90(strip, k=1)

                # Position inset exactly between data coordinates (0,0) and (wl_max,0)
                wl0 = wavelengths[0]
                wl_max = wavelengths[-1]
                # ensure transforms are up-to-date
                try:
                    fig.canvas.draw()
                except Exception:
                    pass
                try:
                    disp_lb = ax.transData.transform((0.0, 0.0))
                    disp_rb = ax.transData.transform((wl_max, 0.0))
                    figcoord_lb = fig.transFigure.inverted().transform(disp_lb)
                    figcoord_rb = fig.transFigure.inverted().transform(disp_rb)
                    left = figcoord_lb[0]
                    bottom = figcoord_lb[1]
                    width = figcoord_rb[0] - figcoord_lb[0]
                    # if width is non-positive for some reason, fallback to axis pos
                    if width <= 0:
                        pos = ax.get_position()
                        left = pos.x0
                        bottom = pos.y0
                        width = pos.width
                except Exception:
                    pos = ax.get_position()
                    left = pos.x0
                    bottom = pos.y0
                    width = pos.width

                # compute inset height so that it is exactly 20 display pixels tall
                fig_w, fig_h = fig.get_size_inches()
                dpi = fig.dpi
                desired_pixels = 20
                inset_height = desired_pixels / (dpi * fig_h)

                # create inset axes at the computed figure coordinates (left,bottom,width,height)
                axins = fig.add_axes([left, bottom, width, inset_height])

                # display the strip with imshow, extent mapped to wavelength axis
                axins.imshow(strip_rot, aspect='auto', cmap='gray', extent=[wl0, wl_max, 0, strip_rot.shape[0]], origin='lower', interpolation='nearest')
                axins.set_xlim(wl0, wl_max)
                axins.set_xticks([])
                axins.set_yticks([])
            except Exception as e:
                print(f"Failed to create inset image strip: {e}")

            ax.grid(alpha=0.3)
            # use bbox_inches='tight' to ensure labels are not cut off when constrained_layout is used
            fig.savefig(png_out, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved wavelength plot to {png_out}")
        except Exception as e:
            print(f"Failed to create plot: {e}")


if __name__ == '__main__':
    main()
