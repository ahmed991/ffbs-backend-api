from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from typing import Literal, Optional, Dict
from processor import process_indicator, historical_viewer
import os
import rasterio
import numpy as np
from PIL import Image
from urllib.parse import quote
from mangum import Mangum

app = FastAPI()
handler = Mangum(app)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

INDICATOR_RANGE = {
    "NDVI": (-1.0, 1.0),
    "NDWI": (-1.0, 1.0),
    "NDMI": (-1.0, 1.0),
    "EVI": (0.0, 3.0),
    "PVI": (0.0, 2.0),
    "LAI": (0.0, 6.0),
}

@app.get("/raster/{filename}")
def serve_raster_png(filename: str):
    png_path = f"results/{filename}.png"
    if os.path.exists(png_path):
        return FileResponse(png_path, media_type="image/png")
    else:
        return {"error": "PNG file not found."}

@app.get("/raster/{filename}/tif")
def serve_raster_tif(filename: str):
    tif_path = f"results/{filename}.tif"
    return FileResponse(tif_path, media_type="image/tiff")

@app.get("/raster/{filename}/bounds")
def get_raster_bounds(filename: str):
    tif_path = f"results/{filename}.tif"
    with rasterio.open(tif_path) as src:
        bounds = src.bounds
    return {
        "coordinates": [
            [bounds.left, bounds.top],
            [bounds.right, bounds.top],
            [bounds.right, bounds.bottom],
            [bounds.left, bounds.bottom]
        ]
    }

class RequestParams(BaseModel):
    geojson: Dict
    start_date: str
    end_date: str
    satellite_sensor: Literal["sentinel-2", "sentinel-1", "landsat", "naip", "cop-dem-30", "cop-dem-90"]
    indicator: Literal["NDVI", "NDWI", "PVI", "LAI", "NDMI", "EVI", "SAVI", "MSI", "SCL", "SFM"]
    cloud_cover: Optional[float] = 100
    resample: Optional[str] = "MS"

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
