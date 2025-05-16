from pystac_client import Client
import geopandas as gpd
import stackstac
import dask.diagnostics as dd
import rasterio
from rasterio.transform import from_bounds
import os
import numpy as np
from urllib.parse import quote
from PIL import Image


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
# TODO add a separate endpoint for green forest change.
def process_indicator(params):
    print(params)
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
        "EVI": ["nir", "red", "blue"],
        "MSI":["nir", "swir16"],
        "SAVI":["nir", "red"],
        "SCL":["scl"]

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

    # Step 5: Compute the indicator
    if indicator == "NDVI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        print("yes here")
        red = stack.sel(band="red").values
        nir = stack.sel(band="nir").values
        index = (nir - red) / (nir + red + 1e-6)


    elif indicator == "NDWI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        green = stack.sel(band="green").values
        nir = stack.sel(band="nir").values
        index = (green - nir) / (green + nir + 1e-6)

    elif indicator == "PVI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        red = stack.sel(band="red").values
        nir = stack.sel(band="nir").values
        index = 1.5 * ((nir - 0.5 * red) / np.sqrt(1.25))  # simplified example

    elif indicator == "NDMI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        nir = stack.sel(band="nir").values
        swir = stack.sel(band="swir16").values
        index = (nir - swir) / (nir + swir + 1e-6)

    elif indicator == "EVI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        print('yes evi')
        nir = stack.sel(band="nir").values
        red = stack.sel(band="red").values
        blue = stack.sel(band="blue").values
        index = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)

    elif indicator == "MSI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        nir = stack.sel(band="nir").values
        swir = stack.sel(band="swir16").values  # swir16 is typically Band 11
        index = swir / nir 

    elif indicator == "SAVI":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        print('yes savi')
        nir = stack.sel(band="nir").values
        red = stack.sel(band="red").values
        L = 0.5
        index = ((nir - red) / (nir + red + L)) * (1 + L)
    elif indicator == "SCL":
        stack = stack.resample(time=resample).median("time", keep_attrs=True).compute()

        scl = stack.sel(band = "scl").values
        index = np.where(np.round(scl)==4,4,0)
        


    else:
        return {"error": f"Indicator logic for {indicator} not implemented."}

# Assuming `index` is a (time, band=1, x, y) numpy array
# and `stack` is the xarray object you used for band selection
    print(index.shape)

# Save GeoTIFFs and return list of filenames
    saved_files = []
    os.makedirs("results", exist_ok=True)

    for i, time_val in enumerate(stack.time.values):
        arr = index[i, :, :].astype("float32")
        transform = from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], arr.shape[1], arr.shape[0])
        profile = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": "float32",
            "crs": "EPSG:3857",
            "transform": transform
        }

        timestamp = np.datetime_as_string(time_val, unit='s').replace(":", "-")
        tif_filename = f"{indicator}_{timestamp}.tif"
        tif_path = f"results/{tif_filename}"
        png_path = tif_path.replace(".tif", ".png")

        with rasterio.open(tif_path, "w", **profile) as dst:
            dst.write(arr, 1)

        # Generate PNG
        arr = np.nan_to_num(arr, nan=-1)
        arr_norm = (arr + 1) / 2  # for NDVI range -1 to 1
        arr_scaled = (arr_norm * 255).clip(0, 255).astype("uint8")
        Image.fromarray(arr_scaled).save(png_path)

        saved_files.append({
            "timestamp": timestamp,
            "tif_url": f"http://localhost:8000/raster/{quote(tif_filename.replace('.tif', ''))}/tif",
            "png_url": f"http://localhost:8000/raster/{quote(tif_filename.replace('.tif', ''))}",
            "bounds": list(bounds)
        })

    return {
        "message": f"{indicator} index computed.",
        "products": saved_files
    }

