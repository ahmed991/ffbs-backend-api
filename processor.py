import stackstac

from services import stac_service, vegetation_service, soil_fertility_service, forest_service, crop_service, raster_service

VEGETATION_INDICATORS = {"NDVI", "NDWI", "PVI", "NDMI", "EVI", "MSI", "SAVI", "LAI"}


def process_indicator(params):
    indicator = stac_service.resolve_indicator(params.indicator)

    bands_required = stac_service.INDICATOR_BAND_MAP.get(indicator)
    if not bands_required:
        return {"error": f"Indicator {indicator} is not supported."}

    collection = stac_service.SENSOR_COLLECTION_MAP.get(params.satellite_sensor)
    bounds = stac_service.get_bounds(params.geojson)

    items = stac_service.search_stac(collection, bounds, params.start_date, params.end_date, params.cloud_cover)
    if not items:
        return {"message": "No imagery found for the given parameters."}

    stack = stackstac.stack(
        items=items,
        epsg=3857,
        assets=bands_required,
        bounds_latlon=bounds,
        resolution=10
    ).resample(time=params.resample).median("time", keep_attrs=True).compute()

    if indicator in VEGETATION_INDICATORS:
        index = vegetation_service.compute(indicator, stack)
    elif indicator == "SFM":
        index = soil_fertility_service.compute(stack)
    elif indicator == "SCL":
        index = forest_service.compute(stack)
    elif indicator == "COTTON":
        index = crop_service.compute(stack)
    else:
        return {"error": f"Indicator logic for {indicator} not implemented."}

    saved_files = raster_service.save_index_outputs(index, stack.time.values, bounds, indicator)

    return {
        "message": f"{indicator} index computed.",
        "products": saved_files
    }
