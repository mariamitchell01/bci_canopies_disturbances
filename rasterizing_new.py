
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

#shapefiles
shps= r'C:\Users\P_pol\repo\bci_canopies_disturbances\shapefiles_disturbances' 
shp_files = [file for file in os.listdir(shps) if file.endswith('.shp')]
path= os.path.join(shps,shp_files[1])

raster_resolution = 1
outfolder = r'C:\Users\Vicente\repo\gapcontagion\disturbances' 
if not os.path.exists(outfolder): 
    os.makedirs(outfolder)

#with rasterio.open(path) as src:
#if src values equal to na then 0
#    data = src.read(1)


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
#get shape bounds
xmin,ymin,xmax,ymax=gap_polygons.total_bounds

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
        s
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

#distance_transform _edt is the answer
disturbances
for date in disturbances:
    path = os.path.join(outfolder2, date)
    with rasterio.open(path) as src:
        data = src.read(1)
        distances = distance_transform_edt(data == 0).astype(np.float32)
        profile = src.profile
        profile.update(dtype=rasterio.float32)
        new_distances_raster_path = os.path.join(distance_dir, date)
        with rasterio.open(new_distances_raster_path, 'w', **profile) as dst:
            dst.write(distances, 1)
        print("Distances raster created:", new_distances_raster_path)



#delta of the distances yeart to year, archaic method and chat gpt loop version

distance_dir =r'C:\Users\P_pol\repo\bci_canopies_disturbances\distances' 
distances_files= [file for file in os.listdir(distance_dir) if file.endswith('.tif')]

with rasterio.open(os.path.join(distance_dir,distances_files[0])) as src:
    data18= src.read(1)

with rasterio.open(os.path.join(distance_dir,distances_files[1])) as src:
    data19= src.read(1)
    
with rasterio.open(os.path.join(distance_dir,distances_files[2])) as src:
    data20= src.read(1)

with rasterio.open(os.path.join(distance_dir,distances_files[3])) as src:
    data21= src.read(1)
    
with rasterio.open(os.path.join(distance_dir,distances_files[4])) as src:
    data22= src.read(1)

# Calculate the difference
delta1 = data18 - data19
delta2= data19-data20
delta3=data20-data21
delta4=data21-data22

dirout=r'C:\Users\P_pol\repo\bci_canopies_disturbances\distancess_delta'
if not os.path.exists(dirout):
    os.makedirs(dirout)



distance_dir = r'C:\Users\P_pol\repo\bci_canopies_disturbances\distances'
output_dir = r'C:\Users\P_pol\repo\bci_canopies_disturbances\distancess_delta'

distances_files = [file for file in os.listdir(distance_dir) if file.endswith('.tif')]
distances_files.sort()

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#writing the loop
for i in range(len(distances_files) - 1):
    with rasterio.open(os.path.join(distance_dir, distances_files[i])) as src1, \
         rasterio.open(os.path.join(distance_dir, distances_files[i + 1])) as src2:       
        data1 = src1.read(1)
        data2 = src2.read(1)
        delta = data1 - data2
        output_path = os.path.join(output_dir, f'delta_{i}.tif')
        profile = src1.profile
        profile.update(dtype=rasterio.uint16)  
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(delta, 1)
        print("Delta data raster saved:", output_path)


