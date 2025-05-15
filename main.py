from fastapi import FastAPI
from pydantic import BaseModel
from typing import Literal
from processor import  process_indicator,historical_viewer  # assume both exist

app = FastAPI()

class RequestParams(BaseModel):
    geojson: dict
    start_date: str
    end_date: str

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
