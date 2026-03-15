import numpy as np
import xarray as xr


def compute(indicator: str, stack) -> xr.DataArray:
    if indicator == "NDVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        return (nir - red) / (nir + red + 1e-6)

    elif indicator == "NDWI":
        nir = stack.sel(band="nir")
        green = stack.sel(band="green")
        return (nir - green) / (nir + green + 1e-6)

    elif indicator == "PVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        return 1.5 * ((nir - 0.5 * red) / np.sqrt(1.25))

    elif indicator == "NDMI":
        nir = stack.sel(band="nir")
        swir = stack.sel(band="swir16")
        return (nir - swir) / (nir + swir + 1e-6)

    elif indicator == "EVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        blue = stack.sel(band="blue")
        return 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)

    elif indicator == "MSI":
        return stack.sel(band="swir16") / stack.sel(band="nir")

    elif indicator == "SAVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        L = 0.5
        return ((nir - red) / (nir + red + L)) * (1 + L)

    else:
        raise ValueError(f"Indicator {indicator} is not handled by vegetation_service.")
