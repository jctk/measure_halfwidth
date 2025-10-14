#!/usr/bin/env python3
"""
generate_test_fits.py

簡易的なテストFITSを生成します。中心に垂直に強いトレース（ガウス）を置きます。
"""
from astropy.io import fits
import numpy as np

out = 'test_trace.fits'
ny = 200
nx = 300
image = np.random.normal(loc=0.0, scale=1.0, size=(ny, nx))

# add a vertical Gaussian trace at center_x
cx = nx // 2
ys = np.arange(ny)
amp = 100.0
sigma_y = 5.0
trace = amp * np.exp(-0.5 * ((ys - ny/2) / sigma_y) ** 2)
# add to center column and neighbors
for dx in [-1, 0, 1]:
    image[:, cx+dx] += trace * (0.8 if dx !=0 else 1.0)

hdu = fits.PrimaryHDU(image)
hdu.writeto(out, overwrite=True)
print('Wrote', out)
