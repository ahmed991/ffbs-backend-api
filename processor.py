from pystac_client import Client
import geopandas as gpd
import stackstac
import xarray as xr
import numpy as np
import pandas as pd
import dask
import pyproj

def process_request(params):
    # Step 1: Load geometry and get bounding box
    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # xmin, ymin, xmax, ymax
    # Step 2: Load STAC items
    client = Client.open("https://earth-search.aws.element84.com/v1")
    search = client.search(
        collections=["sentinel-2-l2a"],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query={"eo:cloud_cover": {"lt": 30}},
        limit=100
    )

    items = list(search.get_items())
    print(f"number of items {len(items)}")
    if not items:
        return {"message": "No imagery found."}
    bounds_1 = [bounds[1],bounds[0],bounds[3],bounds[2]]

    # Step 3: Load assets into xarray with stackstac
    stack = stackstac.stack(items=items,epsg=4326, bounds_latlon=bounds_1)
    print(stack)
    data_lazy = stack.sel(band=["red","nir"])

    stack_resampled = data_lazy.resample(time = "MS").median("time", keep_attrs = True)
    print(stack_resampled)

    # Step 4: Compute NDVI or other index
    if params.index == "NDVI":
        nir = stack_resampled.sel(band="nir")
        red = stack_resampled.sel(band="red")

        index = (nir - red) / (nir + red + 1e-6)
        
    else:
        return {"message": f"Index {params.index} not implemented."}

    index.name = params.index
    return {
        "message": "index request is being completed"
    }
    # # Step 5: Group by aggregation
    # if params.aggregation == "monthly":
    #     grouped = index.groupby("time.month").median(dim="time")
    # elif params.aggregation == "weekly":
    #     week = index.time.dt.strftime("%Y-W%U")
    #     grouped = index.groupby(week).median(dim="time")
    # elif params.aggregation == "daily":
    #     grouped = index
    # else:
    #     return {"message": f"Unknown aggregation: {params.aggregation}"}

    # # Step 6: Return metadata and shape of results
    # return {
    #     "index": params.index,
    #     "aggregation": params.aggregation,
    #     "result_shape": list(grouped.shape),
    #     "time_steps": list(map(str, grouped.coords.get("time", grouped.coords.get("month", [])).values[:5]))
    # }
