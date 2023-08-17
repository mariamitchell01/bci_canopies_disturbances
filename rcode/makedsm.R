install.packages("groundhog")
groundhog::groundhog.library("lidR","2020-11-01")
groundhog::groundhog.library("raster", "2020-09-01")
groundhog::groundhog.library("sp","2020-09-01")
groundhog::groundhog.library("rgdal","2020-09-01")


#read data

# Read polygon buffer 25 m inland from lake
buffer <- rgdal::readOGR("D:/TroubleshootingV/Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  

cat22at <- lidR::catalog("D:/TroubleshootingV/tiles/BCI22Tiles_alignedto21Trim/")
cat23f <-lidR::catalog("D:/TroubleshootingV/tiles/BCI23FTiles_alignedto22Trim/")

#make dsms
  
  
dsm22 <- lidR::grid_canopy(cat22at,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.01))

dsm23 <- lidR::grid_canopy(cat23f,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01))

#save dsms

raster::writeRaster(dsm22, file = "D:/TroubleshootingV/dsms/DSM_2022.tif")
raster::writeRaster(dsm23, file = "D:/TroubleshootingV/dsms/DSMF_2023.tif")