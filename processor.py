import matplotlib
matplotlib.use('Agg')  # Non-GUI backend suitable for servers
from pystac_client import Client
import geopandas as gpd
import stackstac
import rasterio
from rasterio.transform import from_bounds
import os
import numpy as np
from urllib.parse import quote
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import xarray as xr





from pystac_client import Client
import geopandas as gpd
import stackstac
import dask.diagnostics as dd
import rasterio
from rasterio.transform import from_bounds
import os
import numpy as np
from urllib.parse import quote
from PIL import Image

import numpy as np
import xarray as xr


def save_legend_image(output_path, indicator, colormap_used, value_range=None):
    import matplotlib.pyplot as plt
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize, ListedColormap

    fig, ax = plt.subplots(figsize=(4, 0.6))
    ax.set_title(indicator, fontsize=8)

    if indicator == "SFM":
        cmap = ListedColormap(["#a50026", "#f46d43", "#fdae61", "#a6d96a", "#1a9850"])
        bounds = [1, 2, 3, 4, 5, 6]
        norm = Normalize(vmin=1, vmax=5)
        cb = ColorbarBase(ax, cmap=cmap, norm=norm, boundaries=bounds, orientation="horizontal", ticks=[1, 2, 3, 4, 5])
        cb.ax.set_xticklabels(["1", "2", "3", "4", "5"])
    elif indicator == "SCL":
        cmap = ListedColormap(["green"])
        norm = Normalize(vmin=0, vmax=1)
        cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal", ticks=[0, 1])
        cb.ax.set_xticklabels(["0", "Green Cover"])
    else:
        cmap = plt.get_cmap(colormap_used)
        vmin, vmax = value_range if value_range else (0, 1)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal")
        cb.set_label(f"{indicator} ({vmin} to {vmax})", fontsize=7)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', transparent=True)
    plt.close()


def classify_ndvi(ndvi):
    """
    Classify NDVI or SFM values into fertility categories.

    Classes:
        1 = Very Low (NDVI < 0.2)
        2 = Low      (0.2 ≤ NDVI < 0.4)
        3 = Moderate (0.4 ≤ NDVI < 0.6)
        4 = High     (0.6 ≤ NDVI < 0.8)
        5 = Very High(NDVI ≥ 0.8)

    Parameters:
        ndvi: xarray.DataArray or numpy.ndarray

    Returns:
        Classified values with same shape/type
    """

    if isinstance(ndvi, xr.DataArray):
        return xr.where(ndvi < 0.2, 1,
               xr.where(ndvi < 0.4, 2,
               xr.where(ndvi < 0.6, 3,
               xr.where(ndvi < 0.8, 4, 5))))
    else:
        return np.where(ndvi < 0.2, 1,
               np.where(ndvi < 0.4, 2,
               np.where(ndvi < 0.6, 3,
               np.where(ndvi < 0.8, 4, 5))))

def save_index_array(index_array, time_array, bounds, output_dir="results", indicator="INDEX"):
    import os
    os.makedirs(output_dir, exist_ok=True)

    for i, time_val in enumerate(time_array):
        arr = index_array[i, :, :].astype("float32")  # shape: (x, y)

        # Define raster transform from bounds
        transform = from_bounds(
            bounds[0], bounds[1], bounds[2], bounds[3],
            arr.shape[1], arr.shape[0]
        )

        profile = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": "float32",
            "crs": "EPSG:3857",
            "transform": transform
        }

        timestamp = np.datetime_as_string(time_val, unit='s').replace(':', '-')
        filename = f"{output_dir}/{indicator}_{timestamp}.tif"

        with rasterio.open(filename, "w", **profile) as dst:
            dst.write(arr, 1)

        print(f"Saved {filename}")


#TODO : Class based implementations for these functions, proper modulation and separate file systems to handle requests for Different data providers

def historical_viewer(params):
    # Step 1: Parse GeoJSON and extract bounds
    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds  # [xmin, ymin, xmax, ymax]

    # Step 2: Connect to STAC catalog
    catalog = Client.open("https://earth-search.aws.element84.com/v1")

    # Step 3: Query for Sentinel-2 L2A imagery
    query = {
        "eo:cloud_cover": {
            "gte": 0,
            "lte": 30
        }
    }

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

    # Step 4: Extract thumbnail URLs and bounding boxes
    thumbnails = []
    for item in items:
        asset = item.assets.get("thumbnail")
        if asset and item.bbox:
            thumbnails.append({
                "datetime": item.datetime.isoformat(),
                "thumbnail_url": asset.href,
                "id": item.id,
                "bbox": item.bbox  # [minLon, minLat, maxLon, maxLat]
            })

    return {
        "count": len(thumbnails),
        "thumbnails": thumbnails  # Limit to 20 for preview
    }


# TODO: Add a dedicated function for handling STAC queries from the server
# TODO split functions for Sentinel 2A, landast and planet
# TODO separate functions for each indicator, class based handling
# TODO add a separate endpoint for green forest change.
# TODO add a separate endpoint for Soil firtility index
def process_indicator(params):
    from pystac_client import Client
    import geopandas as gpd
    import stackstac
    import rasterio
    from rasterio.transform import from_bounds
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import xarray as xr
    from urllib.parse import quote

    def classify_ndvi(ndvi):
        return xr.where(ndvi < 0.2, 1,
               xr.where(ndvi < 0.4, 2,
               xr.where(ndvi < 0.6, 3,
               xr.where(ndvi < 0.8, 4, 5))))

    def plot_tif_as_png(tif_path, png_path, indicator):
        with rasterio.open(tif_path) as src:
            data = src.read(1)
            fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
            ax.axis('off')
            colormap_used = None

            if indicator == "SFM":
                cmap = ListedColormap(["#a50026", "#f46d43", "#fdae61", "#a6d96a", "#1a9850"])
                colormap_used = "SFM_CUSTOM"
                ax.imshow(data, cmap=cmap, vmin=1, vmax=5)

            elif indicator == "SCL":
                mask = np.where(data == 4, 1, np.nan)
                cmap = ListedColormap(["green"])
                colormap_used = "SCL_GREEN"
                ax.imshow(mask, cmap=cmap)

            else:
                ranges = {
                    "NDVI": (-1, 1), "NDWI": (-1, 1), "NDMI": (-1, 1),
                    "EVI": (0, 3), "PVI": (0, 2), "LAI": (0, 6)
                }
                vmin, vmax = ranges.get(indicator, (np.nanmin(data), np.nanmax(data)))
                cmap = "RdYlGn"
                colormap_used = cmap
                ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax)

            plt.tight_layout()
            
            # Corrected save parameters for web compatibility:
            plt.savefig(
                png_path,
                dpi=300,
                bbox_inches='tight',
                pad_inches=0.05,
                format='png',
                facecolor='white',  # White background instead of transparent
                transparent=False   # Explicitly disable transparency
            )
            
            plt.close()
            return colormap_used

    def save_legend_image(output_path, indicator, colormap_used, value_range=None):
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize

        fig, ax = plt.subplots(figsize=(4, 0.6))
        ax.set_title(indicator, fontsize=8)

        if indicator == "SFM":
            cmap = ListedColormap(["#a50026", "#f46d43", "#fdae61", "#a6d96a", "#1a9850"])
            norm = Normalize(vmin=1, vmax=5)
            cb = ColorbarBase(ax, cmap=cmap, norm=norm, boundaries=[1, 2, 3, 4, 5, 6],
                              orientation="horizontal", ticks=[1, 2, 3, 4, 5])
            cb.ax.set_xticklabels(["1", "2", "3", "4", "5"])
        elif indicator == "SCL":
            cmap = ListedColormap(["green"])
            norm = Normalize(vmin=0, vmax=1)
            cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal", ticks=[0, 1])
            cb.ax.set_xticklabels(["0", "Green"])
        else:
            cmap = plt.get_cmap(colormap_used)
            vmin, vmax = value_range if value_range else (0, 1)
            norm = Normalize(vmin=vmin, vmax=vmax)
            cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal")
            cb.set_label(f"{indicator} ({vmin} to {vmax})", fontsize=7)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', transparent=True)
        plt.close()

    # Begin main logic
    sensor = params.satellite_sensor
    SENSOR_COLLECTION_MAP = {
        "sentinel-2": "sentinel-2-l2a",
        "sentinel-1": "sentinel-1-grd",
        "landsat": "landsat-c2-l2",
        "naip": "naip",
        "cop-dem-30": "cop-dem-glo-30",
        "cop-dem-90": "cop-dem-glo-90"
    }
    collection = SENSOR_COLLECTION_MAP.get(sensor)


    indicator = params.indicator
    print(indicator)
    if (indicator == "Green Forest Change" or indicator =="GREEN FOREST CHANGE"):
        indicator = "SCL"
    elif (indicator == "Soil Fertility Map" or indicator == "SOIL FERTILITY MAP"):
        indicator = "SFM"
    elif (indicator == "Main Crop Identification" or indicator =="MAIN CROP IDENTIFICATION"):
        indicator = "COTTON"

    print(indicator)

    indicator = indicator.upper()
    cloud_cover = params.cloud_cover
    resample = params.resample


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
        "COTTON":[
            'blue',        # B02
            'green',       # B03
            'red',         # B04
            'rededge1',    # B05
            'rededge2',    # B06
            'rededge3',    # B07
            'nir',         # B08
            'nir08',       # B8A (Narrow NIR)
            'swir16',      # B11
            'swir22'       # B12
        ]
    }

    bands_required = INDICATOR_BAND_MAP.get(indicator)
    if not bands_required:
        return {"error": f"Indicator {indicator} is not supported."}

    geometry = gpd.GeoDataFrame.from_features(params.geojson["features"])
    bounds = geometry.total_bounds

    catalog = Client.open("https://earth-search.aws.element84.com/v1")
    query = {"eo:cloud_cover": {"gte": 0, "lte": cloud_cover}}

    result = catalog.search(
        collections=[collection],
        bbox=list(bounds),
        datetime=f"{params.start_date}/{params.end_date}",
        query=query
    ).item_collection()

    items = list(result)
    if not items:
        return {"message": "No imagery found for the given parameters."}

    stack = stackstac.stack(
        items=items,
        epsg=3857,
        assets=bands_required,
        bounds_latlon=bounds,
        resolution=10
    ).resample(time=resample).median("time", keep_attrs=True).compute()

    # Compute index
    if indicator == "NDVI":
        index = (stack.sel(band="nir") - stack.sel(band="red")) / (stack.sel(band="nir") + stack.sel(band="red") + 1e-6)
    elif indicator == "NDWI":
        index = (stack.sel(band="nir") - stack.sel(band="green")) / (stack.sel(band="nir") + stack.sel(band="green") + 1e-6)
    elif indicator == "PVI":
        red = stack.sel(band="red")
        nir = stack.sel(band="nir")
        index = 1.5 * ((nir - 0.5 * red) / np.sqrt(1.25))
    elif indicator == "NDMI":
        nir = stack.sel(band="nir")
        swir = stack.sel(band="swir16")
        index = (nir - swir) / (nir + swir + 1e-6)
    elif indicator == "EVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        blue = stack.sel(band="blue")
        index = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
    elif indicator == "MSI":
        index = stack.sel(band="swir16") / stack.sel(band="nir")
    elif indicator == "SAVI":
        nir = stack.sel(band="nir")
        red = stack.sel(band="red")
        L = 0.5
        index = ((nir - red) / (nir + red + L)) * (1 + L)
    elif indicator == "SCL":
        scl = stack.sel(band="scl")
        index = xr.where(np.round(scl) == 4, 4, 0)
    elif indicator == "SFM":
        ndvi = (stack.sel(band="nir") - stack.sel(band="red")) / (stack.sel(band="nir") + stack.sel(band="red") + 1e-6)
        index = classify_ndvi(ndvi)
    elif indicator == "COTTON":
        # Scale factor to match original Sentinel-2 0–10000 range
        scale_factor = 10000

        # Extract needed bands and scale them
        blue = stack.sel(band='blue') * scale_factor
        green = stack.sel(band='green') * scale_factor
        red = stack.sel(band='red') * scale_factor
        rededge3 = stack.sel(band='rededge3') * scale_factor
        narrow_nir = stack.sel(band='nir08') * scale_factor
        nir = stack.sel(band='nir') * scale_factor
        rededge2 = stack.sel(band='rededge2') * scale_factor
        swir1 = stack.sel(band='swir16') * scale_factor
        swir2 = stack.sel(band='swir22') * scale_factor

        # === Step 8: Calculate WBI ===
        wbi = (
            1.07 * blue
            - 0.68 * green
            - 0.24 * red
            + 0.17 * rededge3
            - 0.04 * narrow_nir
            - 0.39 * nir
            + 0.04 * rededge2
            + 0.36 * swir1
            - 0.01 * swir2
        )

        # Clean up NaN or Inf
        wbi = np.nan_to_num(wbi)

        # === Step 9: Threshold WBI to Detect Cotton ===
        # Based on the paper, optimal WBI threshold is typically 100–220
        threshold = 120
        cotton_mask = xr.where(wbi > threshold, 1, 0).astype('uint8')
  # 1 = Cotton, 0 = Non-Cotton
        cotton_pixels = np.sum(cotton_mask[0] == 1)

        # Each pixel is 100 m²
        area_m2 = cotton_pixels * 100

        # Convert to hectares
        area_ha = area_m2 / 10000
        # For output consistency, assign cotton_mask to index
        index = cotton_mask
    
    else:
        return {"error": f"Indicator logic for {indicator} not implemented."}

    os.makedirs("results", exist_ok=True)
    saved_files = []

    for i, time_val in enumerate(stack.time.values):
        arr = index.sel(time=time_val).values.astype("float32")
        transform = from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], arr.shape[1], arr.shape[0])
        profile = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": "float32",
            "crs": "EPSG:3857",
            "transform": transform
        }
    
        timestamp = np.datetime_as_string(time_val, unit='s').replace(":", "-")
        tif_filename = f"{indicator}_{timestamp}.tif"
        tif_path = f"results/{tif_filename}"
        png_path = tif_path.replace(".tif", ".png")
        legend_path = tif_path.replace(".tif", "_legend.png")

        with rasterio.open(tif_path, "w", **profile) as dst:
            dst.write(arr, 1)

        colormap_used = plot_tif_as_png(tif_path, png_path, indicator)

        if indicator in ["NDVI", "NDWI", "NDMI", "EVI", "PVI", "LAI"]:
            value_range = {
                "NDVI": (-1, 1), "NDWI": (-1, 1), "NDMI": (-1, 1),
                "EVI": (0, 3), "PVI": (0, 2), "LAI": (0, 6)
            }.get(indicator)
        else:
            value_range = None

        save_legend_image(legend_path, indicator, colormap_used, value_range)
        # Save legend image
        print(f"http://3.121.112.193:8000/raster/{quote(legend_path.split('/')[-1].replace('.png', ''))}")
        saved_files.append({
            "timestamp": timestamp,

            "tif_url": f"http://3.121.112.193:8000/raster/{quote(tif_filename.replace('.tif', ''))}/tif",
            "png_url": f"http://3.121.112.193:8000/raster/{quote(tif_filename.replace('.tif', ''))}",
            "legend_url": f"http://3.121.112.193:8000/raster/{quote(legend_path.split('/')[-1].replace('.png', ''))}",

            "bounds": list(bounds),
            "colormap_used": colormap_used

        })
    if indicator == "COTTON":
        return {
            {
        "message": f"{indicator} index computed.",
        "Cotton Actual Area (ha)": area_ha,
        "products": saved_files
    }
        }
    else:

        return {
            "message": f"{indicator} index computed.",
            "products": saved_files
        }

