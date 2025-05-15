from pystac_client import Client
import geopandas as gpd
import stackstac
import dask.diagnostics as dd
import rasterio
from rasterio.transform import from_bounds
import os
import numpy as np


def save_index_array(index_array, time_array, bounds, output_dir="results", indicator="INDEX"):
    import os
    os.makedirs(output_dir, exist_ok=True)

    for i, time_val in enumerate(time_array):
        arr = index_array[i, :, :].astype("float32")  # shape: (x, y)

        # Define raster transform from bounds
        transform = from_bounds(
            bounds[0], bounds[1], bounds[2], bounds[3],
            arr.shape[1], arr.shape[0]
        )

        profile = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": "float32",
            "crs": "EPSG:3857",
            "transform": transform
        }

        timestamp = np.datetime_as_string(time_val, unit='s').replace(':', '-')
        filename = f"{output_dir}/{indicator}_{timestamp}.tif"

        with rasterio.open(filename, "w", **profile) as dst:
            dst.write(arr, 1)

        print(f"Saved {filename}")


#TODO : Class based implementations for these functions, proper modulation and separate file systems to handle requests for Different data providers

def historical_viewer(params):
    # Step 1: Parse GeoJSON and extract bounds
    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # [xmin, ymin, xmax, ymax]

    # Step 2: Connect to STAC catalog
    catalog = Client.open("https://earth-search.aws.element84.com/v1")

    # Step 3: Query for Sentinel-2 L2A imagery
    query = {
        "eo:cloud_cover": {
            "gte": 0,
            "lte": 30
        }
    }

    result = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query=query,
        limit=50
    )

    items = list(result.get_items())
    if not items:
        return {"message": "No imagery found for the selected parameters."}

    # Step 4: Extract thumbnail URLs and bounding boxes
    thumbnails = []
    for item in items:
        asset = item.assets.get("thumbnail")
        if asset and item.bbox:
            thumbnails.append({
                "datetime": item.datetime.isoformat(),
                "thumbnail_url": asset.href,
                "id": item.id,
                "bbox": item.bbox  # [minLon, minLat, maxLon, maxLat]
            })

    return {
        "count": len(thumbnails),
        "thumbnails": thumbnails  # Limit to 20 for preview
    }


# TODO: Add a dedicated function for handling STAC queries from the server
# TODO split functions for Sentinel 2A, landast and planet
# TODO separate functions for each indicator, class based handling
def process_indicator(params):
    print(params)
    print(params.resample)
    sensor = params.satellite_sensor
    SENSOR_COLLECTION_MAP = {
        "sentinel-2": "sentinel-2-l2a",
        "sentinel-1": "sentinel-1-grd",
        "landsat": "landsat-c2-l2",
        "naip": "naip",
        "cop-dem-30": "cop-dem-glo-30",
        "cop-dem-90": "cop-dem-glo-90"
    }
    collection = SENSOR_COLLECTION_MAP.get(sensor)

    indicator = params.indicator.upper()
    cloud_cover = params.cloud_cover
    resample = params.resample

    INDICATOR_BAND_MAP = {
        "NDVI": ["nir", "red"],
        "NDWI": ["nir", "green"],
        "PVI": ["nir", "red"],
        "LAI": ["red", "nir", "swir16"],
        "NDMI": ["nir", "swir16"],
        "EVI": ["nir", "red", "blue"]
    }

    bands_required = INDICATOR_BAND_MAP.get(indicator)
    if not bands_required:
        return {"error": f"Indicator {indicator} is not supported."}

    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # [xmin, ymin, xmax, ymax]

    # Step 2: Connect to STAC
    catalog = Client.open("https://earth-search.aws.element84.com/v1")
    query = {
        "eo:cloud_cover": {
            "gte": 0,
            "lte": cloud_cover
        }
    }

    result = catalog.search(
        collections=[collection],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query=query,
    ).item_collection()

    items = list(result)
    if not items:
        return {"message": "No imagery found for the given parameters."}

    # Step 3: Stack with the correct bands
    stack = stackstac.stack(
        items=items,
        epsg=3857,
        assets=bands_required,
        bounds_latlon=bounds,
        resolution=10  # 10m resolution
    )

    # Step 4: Resample by time
    stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

    # Step 5: Compute the indicator
    if indicator == "NDVI":
        print("yes here")
        red = stack.sel(band="red").values
        nir = stack.sel(band="nir").values
        index = (nir - red) / (nir + red + 1e-6)


    elif indicator == "NDWI":
        green = stack.sel(band="green")
        nir = stack.sel(band="nir")
        index = (green - nir) / (green + nir + 1e-6)

    elif indicator == "PVI":
        red = stack.sel(band="red")
        nir = stack.sel(band="nir")
        index = 1.5 * ((nir - 0.5 * red) / np.sqrt(1.25))  # simplified example

    elif indicator == "NDMI":
        nir = stack.sel(band="nir")
        swir = stack.sel(band="swir16")
        index = (nir - swir) / (nir + swir + 1e-6)

    elif indicator == "EVI":
        print('yes evi')
        nir = stack.sel(band="nir").values
        red = stack.sel(band="red").values
        blue = stack.sel(band="blue").values
        index = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)



    else:
        return {"error": f"Indicator logic for {indicator} not implemented."}

# Assuming `index` is a (time, band=1, x, y) numpy array
# and `stack` is the xarray object you used for band selection
    print(index.shape)
    save_index_array(index, stack.time.values, bounds, indicator=indicator)

    return {
        "message": f"{indicator} index computed.",
        # "time_steps": [str(t.values) for t in index.time[:5]],  # preview
        "shape": index.shape
    }

  



