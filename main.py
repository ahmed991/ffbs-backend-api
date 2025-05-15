from fastapi import FastAPI
from pydantic import BaseModel
from typing import Literal
from processor import  process_indicator,historical_viewer  # assume both exist
from typing import Literal, Optional, Dict


app = FastAPI()

class RequestParams(BaseModel):
    geojson: Dict
    start_date: str
    end_date: str
    satellite_sensor: Literal["sentinel-2", "sentinel-1", "landsat", "naip", "cop-dem-30", "cop-dem-90"]
    indicator: Literal["NDVI", "NDWI", "PVI", "LAI", "NDMI", "EVI"]
    cloud_cover: Optional[float] = 100  # default: allow any cloud cover
    resample: Optional[str] = "MS"  # default: weekly (1W, 1D, 1M)

class ViewerParams(BaseModel):
    geojson: dict
    start_date: str
    end_date: str

@app.post("/compute-index")
async def compute_index(params: RequestParams):
    result = process_indicator(params)
    return {"status": "success", "result": result}

@app.post("/historical-viewer")
async def get_thumbnails(params: ViewerParams):
    result = historical_viewer(params)
    return {"status": "success", "thumbnails": result}
