
#install packages
import geopandas as gpd
import rasterio
import os
import pandas as pd
import matplotlib.pyplot as  plt
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from shapely.geometry import box
import rasterio.mask
import fiona
from rasterio.plot import show
import numpy as np
from scipy.spatial import distance

#shapefiles
shps= r'C:\Users\Vicente\repo\gapcontagion\shapefiles' 
shp_files = [file for file in os.listdir(shps) if file.endswith('.shp')]

raster_resolution = 1
outfolder = r'C:\Users\Vicente\repo\gapcontagion\disturbances' 
if not os.path.exists(outfolder): 
    os.makedirs(outfolder)


for date in shp_files:
    print("rasterizing", os.path.join(shps,date) )
    gap_polygons= gpd.read_file(os.path.join(shps,date))
    xmin,ymin,xmax,ymax=gap_polygons.total_bounds
    print("read shp and got bounds successful") 
    width = int((xmax - xmin) / raster_resolution)
    height = int((ymax - ymin) / raster_resolution)
    transform = from_origin(xmin, ymax, raster_resolution, raster_resolution)
    mask = geometry_mask(gap_polygons.geometry, transform=transform, out_shape=(height, width), invert=True)
    mask = mask.astype(np.uint8) 
    plt.figure(figsize=(8, 8))
    plt.imshow(mask, extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='gray')
    plt.title(f"Mask for {date}")
    plt.show()
    out_path = os.path.join(outfolder, os.path.splitext(date)[0] + '.tif')
    with rasterio.open(out_path, 'w', driver='GTiff', height=height, width=width, count=1, dtype='uint8', crs=gap_polygons.crs, transform=transform) as dst:
        dst.write(mask.astype(rasterio.uint8), 1)
    print("successfully rasterized")
 
disturbances_raster = [file for file in os.listdir(outfolder) if file.endswith('.tif')]
print(disturbances_raster)


with fiona.open(os.path.join(shps,"boundary.shp"), "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
outfolder2 =r'C:\Users\P_pol\repo\bci_canopies_disturbances\disturbances_cropped' 
if not os.path.exists(outfolder2): 
    os.makedirs(outfolder2)

#crop the rasters to same extent
for date in disturbances_raster:
    path = os.path.join(outfolder, date)
    print("working with path:", path)
    with rasterio.open(path) as crocodile:     
        print("numpy array read")
        outimage, outtransform = rasterio.mask.mask(crocodile, shapes, crop=True)
        print("slay")
        out_meta = crocodile.meta.copy()
        out_meta.update({
            'height': outimage.shape[1],
            'width': outimage.shape[2],
            'transform': outtransform
        })
        newpath=os.path.join(outfolder2, date)
        with rasterio.open(newpath, 'w', **out_meta) as dest:
            dest.write(outimage)
        print(f"raster overriden and cropped to {path}")

#read and plot sanity check
distance_dir =r'C:\Users\P_pol\repo\bci_canopies_disturbances\distances' 
if not os.path.exists(distance_dir): 
    os.makedirs(distance_dir)
disturbances= [file for file in os.listdir(outfolder2) if file.endswith('.tif')]
search_radius=100
for date in disturbances:
    date=disturbances[0]
    path = os.path.join(outfolder2, date)
    with rasterio.open(path) as src:
        data= src.read(1)
        row_indices, col_indices = np.where(data== 1)
        distances = np.zeros_like(data, dtype=float)
        for row in range(distances.shape[0]):
            for col in range(distances.shape[1]):
                pixel_coordinates = np.array([[row, col]])
                pixel_coords_2023 = np.column_stack((row_indices, col_indices))

                valid_indices = np.where(
                    (np.abs(row - row_indices) <= search_radius) & 
                    (np.abs(col - col_indices) <= search_radius)
                )
                distances_to_ones = distance.cdist(pixel_coordinates, pixel_coords_2023[valid_indices])
                min_distance = np.min(distances_to_ones)
                distances[row, col] = min_distance
        new_distances_raster_path = os.path.join(outfolder2, 'distances_to_ones.tif')
        with rasterio.open(new_distances_raster_path, 'w', **src.profile) as dst:
            dst.write(distances, 1)
        print("Distances raster created:", new_distances_raster_path)


import rasterstat

def mymean(x):
    return np.ma.mean(x)

zonal_stats("tests/data/polygons.shp",
    "tests/data/slope.tif",
    stats="count",
    add_stats={'mymean':mymean})



#crop to just island pictures

#so just run cdist in non NA pixels = if  no is.na(pixel_coords_2023 NA the continue 

#we need rasterstats


#memory issues arg

import os as that









allbounds = []

#sanity check: checking shapes, bounds <3
for date in disturbances_raster:
    path = os.path.join(outfolder2, date)
    print("working with path:", path)
    with rasterio.open(path) as crocodile:
        data = crocodile.read()
        shape = data.shape
        bounds = crocodile.bounds
    allbounds.append(bounds)
    print(date, shape, bounds)


#to create boundary for crop
df_bounds = pd.DataFrame(allbounds, columns=['xmin', 'ymin', 'xmax', 'ymax'])
box_xmin, box_ymin, box_xmax, box_ymax = df_bounds["xmin"].max(), df_bounds["ymin"].max(), df_bounds["xmax"].min(), df_bounds["ymax"].min()
crop_box = box(box_xmin, box_ymin, box_xmax, box_ymax)
print(crop_box.wkt)
df = gpd.GeoDataFrame({"id":1,"geometry":[crop_box]})
df.to_file(os.path.join(shps,"boundary.shp"))

shp_2020=r'C:\Users\P_pol\Downloads\drive-download-20230812T214513Z-001\gaps18to20sp.shp'
#rasterize 18-20
rasterize_2020= gpd.read_file(shp_2020)
rasterize_2020.crs
xmin, ymin, xmax, ymax= rasterize_2020.total_bounds


width = int((xmax - xmin) / raster_resolution)
height = int((ymax - ymin) / raster_resolution)

transform = from_origin(xmin, ymax, raster_resolution, raster_resolution)
mask = geometry_mask(rasterize_2020.geometry, transform=transform, out_shape=(height, width), invert=True)

out= r'C:\Users\P_pol\Downloads\rasterize18_20.tif'
with rasterio.open(out, 'w', driver='GTiff', height=height, width=width, count=1, dtype='uint8', crs=rasterize_2020.crs, transform=transform) as dst:
    dst.write(mask.astype(rasterio.uint8), 1)