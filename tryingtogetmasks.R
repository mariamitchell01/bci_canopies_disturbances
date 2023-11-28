#cropping disturbances
library(raster)
library(sf)
library(terra)
library(sp)


bci <- st_read("/Volumes/LaCie/stri_thesis/processed_data/aux_shps/BCI_Outline_Minus25.shp")
bci <- st_transform(bci, crs = 32617)
bciraster <- rasterize(bci, mask2018, field = 1)
bciextent <- extent(624080, 629751, 1009717, 1014917)

#importing rasters
#2018
mask2018 <- raster("/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/CHM_2018_QAQC.tif")
nodata2018 <- is.na(mask2018)
#2019
mask2019 <- raster("/Volumes/LaCie/stri_thesis/processed_data/cloudmasks/raster_cloudmask2019.tif")
mask2019 <- raster::crop(mask2019, bci)
mask2019 <- setExtent(mask2019, bciextent)
mask2019 <- raster::resample(mask2019, bciraster, method = "bilinear")
mask2019[mask2019 > 0] <- 1

#2020
mask2020
#2021
mask2021
#2022
mask2022
#2023 (no clouds)


#cropping the disturbances