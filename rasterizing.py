
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

#shapefiles
shps= r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\shapefiles' 
shp_files = [file for file in os.listdir(shps) if file.endswith('.shp')]

raster_resolution = 1
outfolder = r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\outshapefiles' 
if not os.path.exists(outfolder): 
    os.makedirs(outfolder)

shp_2019 = r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\shapefiles\gaps18to19spnew.shp'
shp_2020
shp_2021
shp_2022
shp_2023

for 
transform = from_origin(xmin, ymax, raster_resolution, raster_resolution)
mask = geometry_mask(rasterize_2020.geometry, transform=transform, out_shape=(height, width), invert=True)

out= r'C:\Users\P_pol\Downloads\rasterize18_20.tif'
with rasterio.open(out, 'w', driver='GTiff', height=height, width=width, count=1, dtype='uint8', crs=rasterize_2020.crs, transform=transform) as dst:
    dst.write(mask.astype(rasterio.uint8), 1)



for date in shp_files:
    date=shp_files[1]
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
        dst.write(mask, 1)
    print("successfully rasterized")

width = int((xmax - xmin) / raster_resolution)
height = int((ymax - ymin) / raster_resolution)

transform = from_origin(xmin, ymax, raster_resolution, raster_resolution)
mask = geometry_mask(rasterize_2020.geometry, transform=transform, out_shape=(height, width), invert=True)

out= r'C:\Users\P_pol\Downloads\rasterize18_20.tif'
with rasterio.open(out, 'w', driver='GTiff', height=height, width=width, count=1, dtype='uint8', crs=rasterize_2020.crs, transform=transform) as dst:
    dst.write(mask.astype(rasterio.uint8), 1)

disturbances_rasterized= r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\outshapefiles' 
disturbances_raster = [file for file in os.listdir(disturbances_rasterized) if file.endswith('.tif')]
print(disturbances_raster)

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
with fiona.open(os.path.join(shps,"boundary.shp"), "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
outfolder2 = r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\croppedrasters' 
if not os.path.exists(outfolder2): 
    os.makedirs(outfolder2)

#crop the rasters to same extent
for date in disturbances_raster:
    path = os.path.join(disturbances_rasterized, date)
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
disturbances_dir= r'C:\Users\mitchellme\OneDrive - Smithsonian Institution\Documents\croppedrasters' 
disturbances= [file for file in os.listdir(disturbances_dir) if file.endswith('.tif')]
for date in disturbances:
    path = os.path.join(disturbances_dir, date)
    print("working with path:", path)
    with rasterio.open(path) as src:
        plt.figure(figsize=(8, 8))
        show(src, title=date)
        plt.show()






















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