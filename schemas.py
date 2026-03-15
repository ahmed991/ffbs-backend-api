from pydantic import BaseModel
from typing import Literal, Optional, Dict


class RequestParams(BaseModel):
    geojson: Dict
    start_date: str
    end_date: str
    satellite_sensor: Literal["sentinel-2", "sentinel-1", "landsat", "naip", "cop-dem-30", "cop-dem-90"]
    indicator: Literal["NDVI", "NDWI", "PVI", "LAI", "NDMI", "EVI", "SAVI", "MSI", "Green Forest Change", "Soil Fertility Map", "Main Crop Identification"]
    cloud_cover: Optional[float] = 100
    resample: Optional[str] = "MS"


class ViewerParams(BaseModel):
    geojson: dict
    start_date: str
    end_date: str
