import numpy as np
import xarray as xr


def compute(stack) -> xr.DataArray:
    """
    Green Forest Change detection using the Scene Classification Layer (SCL).
    Class 4 in SCL = Vegetation (green cover).
    Returns a binary mask: 4 where green cover is detected, 0 elsewhere.
    """
    scl = stack.sel(band="scl")
    return xr.where(np.round(scl) == 4, 4, 0)
