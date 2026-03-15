from pystac_client import Client
import geopandas as gpd

STAC_URL = "https://earth-search.aws.element84.com/v1"

SENSOR_COLLECTION_MAP = {
    "sentinel-2": "sentinel-2-l2a",
    "sentinel-1": "sentinel-1-grd",
    "landsat": "landsat-c2-l2",
    "naip": "naip",
    "cop-dem-30": "cop-dem-glo-30",
    "cop-dem-90": "cop-dem-glo-90"
}

INDICATOR_BAND_MAP = {
    "NDVI": ["nir", "red"],
    "NDWI": ["nir", "green"],
    "PVI": ["nir", "red"],
    "LAI": ["red", "nir", "swir16"],
    "NDMI": ["nir", "swir16"],
    "EVI": ["nir", "red", "blue"],
    "MSI": ["nir", "swir16"],
    "SAVI": ["nir", "red"],
    "SCL": ["scl"],
    "SFM": ["nir", "red"],
    "COTTON": [
        "blue", "green", "red",
        "rededge1", "rededge2", "rededge3",
        "nir", "nir08", "swir16", "swir22"
    ]
}

INDICATOR_ALIAS_MAP = {
    "Green Forest Change": "SCL",
    "GREEN FOREST CHANGE": "SCL",
    "Soil Fertility Map": "SFM",
    "SOIL FERTILITY MAP": "SFM",
    "Main Crop Identification": "COTTON",
    "MAIN CROP IDENTIFICATION": "COTTON",
}


def resolve_indicator(indicator: str) -> str:
    return INDICATOR_ALIAS_MAP.get(indicator, indicator.upper())


def get_bounds(geojson: dict):
    geometry = gpd.GeoDataFrame.from_features(geojson["features"])
    return geometry.total_bounds


def search_stac(collection: str, bounds, start_date: str, end_date: str, cloud_cover: float):
    catalog = Client.open(STAC_URL)
    query = {"eo:cloud_cover": {"gte": 0, "lte": cloud_cover}}
    result = catalog.search(
        collections=[collection],
        bbox=list(bounds),
        datetime=f"{start_date}/{end_date}",
        query=query
    ).item_collection()
    return list(result)


def historical_viewer(params):
    bounds = get_bounds(params.geojson)
    catalog = Client.open(STAC_URL)
    query = {"eo:cloud_cover": {"gte": 0, "lte": 30}}

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

    thumbnails = []
    for item in items:
        asset = item.assets.get("thumbnail")
        if asset and item.bbox:
            thumbnails.append({
                "datetime": item.datetime.isoformat(),
                "thumbnail_url": asset.href,
                "id": item.id,
                "bbox": item.bbox
            })

    return {
        "count": len(thumbnails),
        "thumbnails": thumbnails
    }
