from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from mangum import Mangum
from concurrent.futures import ThreadPoolExecutor
import asyncio

from schemas import RequestParams, ViewerParams
from processor import process_indicator
from services.stac_service import historical_viewer
from services.raster_service import get_png_path, get_tif_path, get_raster_bounds

app = FastAPI()
handler = Mangum(app)
_executor = ThreadPoolExecutor()

app.mount("/static", StaticFiles(directory="static"), name="static")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.post("/compute-index")
async def compute_index(params: RequestParams):
    loop = asyncio.get_event_loop()
    result = await loop.run_in_executor(_executor, process_indicator, params)
    return {"status": "success", "result": result}


@app.post("/historical-viewer")
async def get_thumbnails(params: ViewerParams):
    loop = asyncio.get_event_loop()
    result = await loop.run_in_executor(_executor, historical_viewer, params)
    return {"status": "success", "result": result}


@app.get("/raster/{filename}")
def serve_raster_png(filename: str):
    path = get_png_path(filename)
    if not path:
        raise HTTPException(status_code=404, detail="PNG file not found.")
    return FileResponse(path, media_type="image/png")


@app.get("/raster/{filename}/tif")
def serve_raster_tif(filename: str):
    path = get_tif_path(filename)
    if not path:
        raise HTTPException(status_code=404, detail="TIF file not found.")
    return FileResponse(path, media_type="image/tiff")


@app.get("/raster/{filename}/bounds")
def serve_raster_bounds(filename: str):
    bounds = get_raster_bounds(filename)
    if not bounds:
        raise HTTPException(status_code=404, detail="TIF file not found.")
    return bounds
