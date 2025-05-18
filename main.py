from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse
from typing import Literal
from processor import  process_indicator,historical_viewer  # assume both exist
from typing import Literal, Optional, Dict
import rasterio
import numpy as np
import os
from PIL import Image
from fastapi.middleware.cors import CORSMiddleware
from mangum import Mangum
# TODO, add file clearner by date
app = FastAPI()
handler = Mangum(app)
# 🚀 Allow requests from your Vite frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://3.121.112.193:5173"],  # frontend origin
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
INDICATOR_RANGE = {
    "NDVI": (-1.0, 1.0),
    "NDWI": (-1.0, 1.0),
    "NDMI": (-1.0, 1.0),
    "EVI": (0.0, 3.0),
    "PVI": (0.0, 2.0),     # you can adjust this
    "LAI": (0.0, 6.0)
}
@app.get("/raster/{filename}")
@app.get("/raster/{filename}")
def serve_raster_png(filename: str):
    tif_path = f"results/{filename}.tif"
    png_path = f"results/{filename}.png"

    if not os.path.exists(png_path):
        with rasterio.open(tif_path) as src:
            arr = src.read(1)
            arr = np.nan_to_num(arr, nan=-9999)

            # Extract indicator name from filename
            indicator = filename.split("_")[0].upper()
            value_range = INDICATOR_RANGE.get(indicator, (arr.min(), arr.max()))

            min_val, max_val = value_range
            arr_norm = (arr - min_val) / (max_val - min_val + 1e-6)
            arr_scaled = (arr_norm * 255).clip(0, 255).astype("uint8")

            im = Image.fromarray(arr_scaled)
            im.save(png_path)

    return FileResponse(png_path, media_type="image/png")


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
            [bounds.left, bounds.top],     # top-left
            [bounds.right, bounds.top],    # top-right
            [bounds.right, bounds.bottom], # bottom-right
            [bounds.left, bounds.bottom]   # bottom-left
        ]
    }


class RequestParams(BaseModel):
    geojson: Dict
    start_date: str
    end_date: str
    satellite_sensor: Literal["sentinel-2", "sentinel-1", "landsat", "naip", "cop-dem-30", "cop-dem-90"]
    indicator: Literal["NDVI", "NDWI", "PVI", "LAI", "NDMI", "EVI","SAVI","MSI","SCL"]
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
