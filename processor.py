from pystac_client import Client
import geopandas as gpd

def historical_viewer(params):
    # Step 1: Parse GeoJSON and extract bounding box
    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # [xmin, ymin, xmax, ymax]

    # Step 2: Connect to STAC and search items
    client = Client.open("https://earth-search.aws.element84.com/v1")
    search = client.search(
        collections=["sentinel-2-l2a"],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query={"eo:cloud_cover": {"lt": 20}},
        limit=50
    )

    items = list(search.get_items())
    if not items:
        return {"message": "No imagery found for the selected area and time range."}

    # Step 3: Collect thumbnails and dates
    thumbnails = []
    for item in items:
        thumb_asset = item.assets.get("thumbnail")
        if thumb_asset:
            thumbnails.append({
                "datetime": item.datetime.isoformat(),
                "thumbnail_url": thumb_asset.href,
                "id": item.id
            })

    return {
        "count": len(thumbnails),
        "thumbnails": thumbnails[:20]  # limit to 20 entries for preview
    }
