from pystac_client import Client
from datetime import datetime
import geopandas as gpd
import shapely.geometry

def process_request(params):
    """Handles loading Sentinel-2 data from STAC and prepares for index calculation."""
    # Step 1: Extract bounding box from GeoJSON
    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # xmin, ymin, xmax, ymax

    # Step 2: Connect to STAC
    client = Client.open("https://earth-search.aws.element84.com/v1")

    # Step 3: Search for Sentinel-2 imagery in the region and time range
    search = client.search(
        collections=["sentinel-2-l2a"],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query={"eo:cloud_cover": {"lt": 30}},  # Optional: Filter low-cloud imagery
        limit=100
    )

    items = list(search.get_items())

    if not items:
        return {"message": "No Sentinel-2 imagery found for the given parameters."}

    # Step 4: Extract useful metadata or asset URLs
    asset_list = []
    for item in items:
        assets = item.to_dict()["assets"]
        # You may later filter by index (NDVI = B08, B04)
        asset_list.append({
            "id": item.id,
            "date": item.datetime.strftime("%Y-%m-%d"),
            "assets": {
                "B04": assets.get("B04", {}).get("href", ""),
                "B08": assets.get("B08", {}).get("href", "")
            }
        })

    # Step 5: Return metadata, you will later add NDVI calculation
    return {
        "index_requested": params.index,
        "aggregation": params.aggregation,
        "results_found": len(asset_list),
        "data": asset_list[:5]  # Return first 5 for preview
    }
