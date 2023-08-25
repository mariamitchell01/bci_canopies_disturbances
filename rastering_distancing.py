#### Getting distance .tif to the nearest gap #####

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
from scipy.ndimage import distance_transform_edt

#### read in disturbance shapefiles ####
shps= r'/Volumes/LaCie/stri_thesis/processed_data/shapefiles_disturbances' 
aux_shps = r'/Volumes/LaCie/stri_thesis/processed_data/aux_shps'
shp_files = [file for file in os.listdir(shps) if file.endswith('.shp')]
path= os.path.join(shps,shp_files[1])

#### create outfolder for disturbance rasters #####
raster_resolution = 1
outfolder = r'/Volumes/LaCie/stri_thesis/processed_data/rasters_disturbances' 
if not os.path.exists(outfolder): 
    os.makedirs(outfolder)

##### rasterize the disturbance shapefiles #####
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
 
 ##### read in the disturbance rasters ###
rasters_disturbances = [file for file in os.listdir(outfolder) if file.endswith('.tif')]
print(rasters_disturbances)


#### create folder for cropped disturbance reasters #####
with fiona.open(os.path.join(aux_shps,"boundary.shp"), "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
outfolder2 =r'/Volumes/LaCie/stri_thesis/processed_data/rasters_disturbances_cropped' 
if not os.path.exists(outfolder2): 
    os.makedirs(outfolder2)

#### crop the rasters to same extent ######
for date in rasters_disturbances:
    path = os.path.join(outfolder, date)
    print("working with path:", path)
    with rasterio.open(path) as crocodile:     
        print("numpy array read")
        outimage, outtransform = rasterio.mask.mask(crocodile, shapes, crop=True)
        print("finished")
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

##### create folder for distance-to-gap rasters ####
distance_dir =r'/Volumes/LaCie/stri_thesis/processed_data/raster_distances' 
if not os.path.exists(distance_dir): 
    os.makedirs(distance_dir)

##### overwrite previous raster_disturbances #####
raster_disturbances = [file for file in os.listdir(outfolder2) if file.endswith('.tif')]

##### use distance_transform _edt for finding the distance to gaps ###
for date in rasters_disturbances:
    path = os.path.join(outfolder2, date)
    with rasterio.open(path) as src:
        data = src.read(1)
        distances = distance_transform_edt(data == 0).astype(np.float32)
        profile = src.profile
        profile.update(dtype=rasterio.float32)
        new_name = "distances_" + date
        new_distances_raster_path = os.path.join(distance_dir, new_name)
        with rasterio.open(new_distances_raster_path, 'w', **profile) as dst:
            dst.write(distances, 1)
        print("Distances raster created:", new_distances_raster_path)



