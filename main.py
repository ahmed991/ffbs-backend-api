from fastapi import FastAPI
from pydantic import BaseModel
from typing import Literal
from processor import process_request

app = FastAPI()

class RequestParams(BaseModel):
    geojson: dict
    product: Literal["sentinel-2"]
    index: Literal["NDVI", "NDWI", "PVI"]
    aggregation: Literal["daily", "weekly", "monthly"]
    start_date: str
    end_date: str

@app.post("/compute-index")
async def compute_index(params: RequestParams):
    result = process_request(params)
    return {"status": "success", "result": result}
