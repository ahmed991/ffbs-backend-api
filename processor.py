from pystac_client import Client
import geopandas as gpd
import stackstac
import dask.diagnostics as dd

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

    # Step 4: Extract thumbnail URLs
    thumbnails = []
    for item in items:
        asset = item.assets.get("thumbnail")
        if asset:
            thumbnails.append({
                "datetime": item.datetime.isoformat(),
                "thumbnail_url": asset.href,
                "id": item.id
            })

    return {
        "count": len(thumbnails),
        "thumbnails": thumbnails  # limit to 20 for preview
    }

# TODO: Add a dedicated function for handling STAC queries from the server
def process_indicator(params):
    print(params)

    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # [xmin, ymin, xmax, ymax]
    bounds_1 = [bounds[1],bounds[0],bounds[3],bounds[2]]
    print(bounds_1)
    # Step 2: Connect to STAC catalog
    catalog = Client.open("https://earth-search.aws.element84.com/v1")
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
    ).item_collection()
    stack = stackstac.stack(items=result,epsg=4326,
                        bounds_latlon=bounds_1,
                        )
    red = stack.sel(band = 'red')
    nir = stack.sel(band = 'nir')
    ndvi = (nir-red)/(nir+red)
    with dd.ProgressBar():
        data = ndvi.compute()

    print(data)



