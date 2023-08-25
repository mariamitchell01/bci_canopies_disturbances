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
from shapely.geometry import MultiPolygon

from scipy.spatial import distance
from scipy.ndimage import distance_transform_edt

#### read in disturbance shapefiles ####
shps= r'C:\Users\P_pol\repo\bci_canopies_disturbances\shapefiles_disturbances' 
aux_shps = r'/Volumes/LaCie/stri_thesis/processed_data/aux_shps'
path_bci = r'C:\Users\P_pol\repo\bci_canopies_disturbances\aux_shps\BCI_Outline_Minus25.shp'
distance_dir =r'C:\Users\P_pol\repo\bci_canopies_disturbances\distances' 
outfolder2=r'C:\Users\P_pol\repo\bci_canopies_disturbances\disturbances_cropped' 
outfolder = r'C:\Users\P_pol\repo\bci_canopies_disturbances\disturbances' 

shp_files = [file for file in os.listdir(shps) if file.endswith('.shp')]

if not os.path.exists(distance_dir): 
    os.makedirs(distance_dir)
if not os.path.exists(outfolder2): 
    os.makedirs(outfolder2)  
if not os.path.exists(distance_dir):
    os.makedirs(distance_dir)
if not os.path.exists(outfolder): 
    os.makedirs(outfolder)

#### create outfolder for disturbance rasters #####
raster_resolution = 1


##### rasterize the disturbance shapefiles #####
for date in shp_files:
    date= shp_files[0]
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
    # i have largest_polygon which is polygon i want everything out of that polygon to be an NA
    out_image, out_transform= rasterio.mask.mask(mask,[largest_polygon], crop=True)
    plt.figure(figsize=(8, 8))
    plt.imshow(out_image[0], extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='gray')
    plt.title(f"Mask for {date}")
    plt.show()
    out_path = os.path.join(outfolder, os.path.splitext(date)[0] + '.tif')
    with rasterio.open(out_path, 'w', driver='GTiff', height=height, width=width, count=1, dtype='uint8', crs=gap_polygons.crs, transform=transform) as dst:
        dst.write(mask.astype(rasterio.uint8), 1)
    print("successfully rasterized")
 
 ##### read in the disturbance rasters ###
rasters_disturbances = [file for file in os.listdir(outfolder) if file.endswith('.tif')]
print(rasters_disturbances)


#### read bci shapefile #### 
shapes = gpd.read_file(path_bci)                                #reads geopandas data frame containin shapes ion ['geometry'] column
shapes['geometry'] = shapes['geometry'].to_crs(epsg=32617)      #reproject to desiered CRS
if isinstance(shapes, MultiPolygon):                            #if it is a multipolygon
                multi_polygon = shapes
                polygons = list(multi_polygon.geoms)            #list of polygons that compose the multipolygon
                largest_polygon = max(polygons, key=lambda polygon: polygon.area)

largest_polygon.bounds
for date in rasters_disturbances:
    path = os.path.join(outfolder, date)
    print("working with path:", path)
    with rasterio.open(path) as crocodile:     
        print("numpy array read")
        data = crocodile.read(1)
        outimage, outtransform = rasterio.mask.mask(crocodile,[largest_polygon], crop=True, nodata=-9999)  #!the brackets are important very important because it needs to be iterable
        plt.imshow(outimage[0], cmap='viridis')  # You can change the cmap as needed
        plt.show()                                #take out as needed
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

##### overwrite previous raster_disturbances #####
raster_disturbances = [file for file in os.listdir(outfolder2) if file.endswith('.tif')]

##### use distance_transform _edt for finding the distance to gaps ###
for date in rasters_disturbances:
    date = rasters_disturbances[0]
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


