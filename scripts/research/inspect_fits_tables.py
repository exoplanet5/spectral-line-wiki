"""Print headers and first rows of FITS binary tables passed as arguments."""
import sys
from astropy.io import fits
for p in sys.argv[1:]:
    print("=====", p)
    with fits.open(p) as h:
        for i, hdu in enumerate(h):
            print("--- HDU", i, type(hdu).__name__)
            print(repr(hdu.header))
            if hasattr(hdu, "columns") and hdu.columns is not None:
                print(hdu.columns)
                print(hdu.data[:3])
                print("nrows", len(hdu.data))
