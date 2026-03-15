import xarray as xr


def classify_ndvi(ndvi: xr.DataArray) -> xr.DataArray:
    """
    Classify NDVI into 5 soil fertility categories:
        1 = Very Low  (NDVI < 0.2)
        2 = Low       (0.2 <= NDVI < 0.4)
        3 = Moderate  (0.4 <= NDVI < 0.6)
        4 = High      (0.6 <= NDVI < 0.8)
        5 = Very High (NDVI >= 0.8)
    """
    return xr.where(ndvi < 0.2, 1,
           xr.where(ndvi < 0.4, 2,
           xr.where(ndvi < 0.6, 3,
           xr.where(ndvi < 0.8, 4, 5))))


def compute(stack) -> xr.DataArray:
    nir = stack.sel(band="nir")
    red = stack.sel(band="red")
    ndvi = (nir - red) / (nir + red + 1e-6)
    return classify_ndvi(ndvi)
