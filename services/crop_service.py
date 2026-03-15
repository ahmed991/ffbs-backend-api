import xarray as xr


WBI_THRESHOLD = 120
SCALE_FACTOR = 10000


def compute(stack) -> xr.DataArray:
    """
    Cotton crop detection using the Water Band Index (WBI).
    Based on a 10-band weighted formula using Sentinel-2 bands.
    Returns a binary mask: 1 = Cotton, 0 = Non-Cotton.
    """
    blue      = stack.sel(band='blue')      * SCALE_FACTOR
    green     = stack.sel(band='green')     * SCALE_FACTOR
    red       = stack.sel(band='red')       * SCALE_FACTOR
    rededge3  = stack.sel(band='rededge3')  * SCALE_FACTOR
    narrow_nir = stack.sel(band='nir08')    * SCALE_FACTOR
    nir       = stack.sel(band='nir')       * SCALE_FACTOR
    rededge2  = stack.sel(band='rededge2')  * SCALE_FACTOR
    swir1     = stack.sel(band='swir16')    * SCALE_FACTOR
    swir2     = stack.sel(band='swir22')    * SCALE_FACTOR

    wbi = (
        1.07  * blue
        - 0.68 * green
        - 0.24 * red
        + 0.17 * rededge3
        - 0.04 * narrow_nir
        - 0.39 * nir
        + 0.04 * rededge2
        + 0.36 * swir1
        - 0.01 * swir2
    )

    wbi = wbi.fillna(0)
    return xr.where(wbi > WBI_THRESHOLD, 1, 0)
