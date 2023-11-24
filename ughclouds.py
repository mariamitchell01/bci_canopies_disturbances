#packages
import rasterio
import shapely
from shapely.geometry import box
from rasterio.mask import mask
import numpy as np
import matplotlib.pyplot as plt
from rasterio.features import shapes
from shapely.geometry import shape
import geopandas as gpd
from rasterio.plot import show
import cv2
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import pandas as pd
import os
from shapely.geometry import box as box1

def tile_ortho(sub, rows, columns, buffer, output_folder):
    with rasterio.open(sub) as src:
        bounds = src.bounds
        xmin, ymin, xmax, ymax = bounds
        tile_height = (ymax - ymin) / rows
        tile_width = (xmax - xmin) / columns
        xmins = np.arange(xmin, (xmax - tile_width + 1), tile_width)
        xmaxs = np.arange((xmin + tile_width), xmax + 1, tile_width)
        ymins = np.arange(ymin, (ymax - tile_height + 1), tile_height)
        ymaxs = np.arange((ymin + tile_height), ymax + 1, tile_height)
        X, Y = np.meshgrid(xmins, ymins)
        Xmax, Ymax = np.meshgrid(xmaxs, ymaxs)
        gridInfo = pd.DataFrame({
            'xmin': X.flatten(),
            'ymin': Y.flatten(),
            'xmax': Xmax.flatten(),
            'ymax': Ymax.flatten(),
        })
        print(gridInfo)
    with rasterio.open(sub) as src:
        for idx, row in gridInfo.iterrows():
            geom2 = box(row['xmin'] - buffer, row['ymin'] - buffer, row['xmax'] + buffer, row['ymax'] + buffer)
            out_image, out_transform = rasterio.mask.mask(src, [geom2], crop=True)

            out_meta = src.meta
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform})
            output_filename = f"output_raster_{idx}.tif"
            filename = os.path.join(output_folder, output_filename)
            with rasterio.open(filename, "w", **out_meta) as dest:
                dest.write(out_image)

clouded=r"/Volumes/LaCie/stri_thesis/raw_data/2019_06_27_BCI_Dipteryx_18cm.tif"
whole_island_shp= r"/Volumes/LaCie/stri_thesis/processed_data/aux_shps/BCI_Outline_Minus25.shp"
if not os.path.exists(r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/clouded_tiles"):
    os.makedirs(r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/clouded_tiles")


tile_ortho(clouded, 5, 5, 0, r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/clouded_tiles")
list_of_clouded_tiles= os.listdir(r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/clouded_tiles")

final_gdf= gpd.GeoDataFrame()
list_of_dfs = []

for tile in list_of_clouded_tiles:
    clouded= os.path.join(r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/clouded_tiles", tile)
    with rasterio.open(clouded) as src:
        data2=src.read()
        transform=src.transform
        data2=data2.transpose(1,2,0)
    #calculate the standard deviation of the clouded image
    std_deviati = np.std(data2, axis=-1).astype(np.int8)
    mask_dv=np.where(std_deviati<40,1,0)

    #calculate the mean of the clouded image
    mean = np.mean(data2, axis=-1)
    mask_mean=np.where(mean>175,1,0)
    
    blue = data2[:,:,2]
    p90 = np.percentile(blue, 80)
    blue_filtered = blue[blue > p90]

    if blue_filtered.size > 0:
        blue_upper_10 = blue_filtered.min()
    # Continue with your code
    else:
        blue_upper_10 = 0 

    mask_blue=np.where(blue>blue_upper_10,2,0)

    albedo= data2[:,:,0]+data2[:,:,1]+data2[:,:,2]
    mask_albedo=np.where(albedo>195,2,0)
    clouds1= mask_dv+mask_mean+mask_blue+mask_albedo

    maskfinal=np.where(clouds1>1,0,1)
    kernel = np.ones((8, 8), np.uint8)  
    mask_dilated = cv2.dilate(maskfinal.astype(np.uint8), kernel, iterations=3)
    contours, _ = cv2.findContours(mask_dilated, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    contour_vis = np.zeros_like(mask_dilated)
    cv2.drawContours(contour_vis, contours, -1, (255, 255, 255), 2)

    smooth_mask = np.zeros_like(mask_dilated)
    cv2.drawContours(smooth_mask, contours, -1, 1, thickness=cv2.FILLED)

    results = list(shapes(smooth_mask.astype(np.uint8), transform=transform))
    geoms = [shape(result[0]) for result in results]
    
    # Create a GeoDataFrame with the polygons
    gdf = gpd.GeoDataFrame(geometry=geoms)
    gdf['area'] = gdf.area
    gdf = gdf.sort_values(by=['area'], ascending=False)
    gdf = gdf.iloc[1:]
    gdf.crs = "EPSG:32617"

    # Convert to a Pandas DataFrame
    df_to_append = pd.DataFrame({
    'geometry': gdf['geometry'],
    'area': gdf['area']
    })
    list_of_dfs.append(df_to_append)
final_df = pd.concat(list_of_dfs, ignore_index=True)

geometry = final_df['geometry']
crs = "EPSG:32617"  # Update with your desired coordinate reference system
gdf_final = gpd.GeoDataFrame(final_df, geometry=geometry, crs=crs)
gdf_final = gdf_final.sort_values(by=['area'], ascending=False)
percentile_90 = gdf_final['area'].quantile(0.95)
gdf_final_filtered = gdf_final[gdf_final['area'] > percentile_90]

whole_island = gpd.GeoDataFrame(whole_island, geometry=whole_island.geometry, crs="EPSG:32617")

final_gdf2 = gpd.sjoin(gdf_final_filtered, whole_island, how="inner", op='intersects')
final_gdf3 = final_gdf2.drop(columns='index_right')

final_gdf3['geometry'] = final_gdf3.buffer(25) #i want 25 meter buffer
union_geom = unary_union(final_gdf3['geometry'])

union_gdf = gpd.GeoDataFrame(geometry=[union_geom])
union_gdf.crs = "EPSG:32617"

union_gdf.to_file(r"/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/cloud_masked2019.shp")