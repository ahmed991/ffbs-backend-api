import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from rasterio.transform import from_bounds
from urllib.parse import quote

SERVER_URL = "http://3.121.112.193:8000"
RESULTS_DIR = "results"

INDICATOR_VALUE_RANGES = {
    "NDVI": (-1, 1),
    "NDWI": (-1, 1),
    "NDMI": (-1, 1),
    "EVI": (0, 3),
    "PVI": (0, 2),
    "LAI": (0, 6),
    "COTTON": (0, 1),
}


def save_geotiff(arr: np.ndarray, bounds, tif_path: str):
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
    with rasterio.open(tif_path, "w", **profile) as dst:
        dst.write(arr, 1)


def plot_tif_as_png(tif_path: str, png_path: str, indicator: str) -> str:
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
    plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.05,
                format='png', facecolor='white', transparent=False)
    plt.close()
    return colormap_used


def save_legend_image(output_path: str, indicator: str, colormap_used: str, value_range=None):
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


def build_file_entry(indicator: str, tif_filename: str, legend_path: str, bounds, colormap_used: str, arr: np.ndarray = None, timestamp: str = ""):
    entry = {
        "timestamp": timestamp,
        "tif_url": f"{SERVER_URL}/raster/{quote(tif_filename.replace('.tif', ''))}/tif",
        "png_url": f"{SERVER_URL}/raster/{quote(tif_filename.replace('.tif', ''))}",
        "legend_url": f"{SERVER_URL}/raster/{quote(legend_path.split('/')[-1].replace('.png', ''))}",
        "bounds": list(bounds),
        "colormap_used": colormap_used,
    }

    if indicator == "COTTON" and arr is not None:
        cotton_pixels = (arr == 1).sum()
        area_m2 = cotton_pixels * 100
        area_ha = area_m2 / 10000
        entry["cotton_area_ha"] = round(float(area_ha), 2)

    return entry


def get_raster_bounds(filename: str):
    tif_path = f"{RESULTS_DIR}/{filename}.tif"
    if not os.path.exists(tif_path):
        return None
    with rasterio.open(tif_path) as src:
        b = src.bounds
    return {
        "coordinates": [
            [b.left, b.top],
            [b.right, b.top],
            [b.right, b.bottom],
            [b.left, b.bottom]
        ]
    }


def get_png_path(filename: str):
    path = f"{RESULTS_DIR}/{filename}.png"
    return path if os.path.exists(path) else None


def get_tif_path(filename: str):
    path = f"{RESULTS_DIR}/{filename}.tif"
    return path if os.path.exists(path) else None


def save_index_outputs(index, stack_time_values, bounds, indicator: str):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    saved_files = []

    for time_val in stack_time_values:
        arr = index.sel(time=time_val).values.astype("float32")
        timestamp = np.datetime_as_string(time_val, unit='s').replace(":", "-")

        tif_filename = f"{indicator}_{timestamp}.tif"
        tif_path = f"{RESULTS_DIR}/{tif_filename}"
        png_path = tif_path.replace(".tif", ".png")
        legend_path = tif_path.replace(".tif", "_legend.png")

        save_geotiff(arr, bounds, tif_path)
        colormap_used = plot_tif_as_png(tif_path, png_path, indicator)

        value_range = INDICATOR_VALUE_RANGES.get(indicator)
        save_legend_image(legend_path, indicator, colormap_used, value_range)

        entry = build_file_entry(indicator, tif_filename, legend_path, bounds, colormap_used, arr, timestamp)
        saved_files.append(entry)

    return saved_files
