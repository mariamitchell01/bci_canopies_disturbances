
#install packages

groundhog::groundhog.library("lidR","2020-11-01")
groundhog::groundhog.library("raster", "2020-09-01")
groundhog::groundhog.library("sp","2020-09-01")
groundhog::groundhog.library("rgdal","2020-09-01")

# read dsms
dsm20 <- raster::raster("D:/TroubleshootingV/dsms/DSM_2021_corrected.tif")
dsm21 <- raster::raster("D:/TroubleshootingV/dsms/DSM_2021_corrected.tif")
dsm22 <- raster::raster("D:/TroubleshootingV/dsms/DSM_2022.tif")
dsm23 <- raster::raster("D:/TroubleshootingV/dsms/DSMF_2023.tif")

#read forest age
age <- rgdal::readOGR("D:/TroubleshootingV/Data_Ancillary/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
age$AgeClass <- "Other"
age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
ageUse <- age[!(age$AgeClass=="Other"),]

#read polygon buffer25m inland from lake
buffer <- rgdal::readOGR("D:/TroubleshootingV/Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  

# Read in BCI DEM (from 2009 lidar data)
dem <- raster::raster("D:/TroubleshootingV/dsms/LidarDEM_BCI.tif")
dem <- raster::crop(dem, raster::extent(ageUse))
dem <- raster::resample(dem,dsm22)

#remove raster areas outside BCI perimeter
dsm20 <- raster::mask(dsm20, buffer)
dsm21 <- raster::mask(dsm21, buffer)
dsm22 <- raster::mask(dsm22, buffer)
dsm23 <- raster::mask(dsm23, buffer)

#remove raster areas in clearings
dsm20 <- raster::mask(dsm20, ageUse)
dsm21 <- raster::mask(dsm21, ageUse)
dsm22 <- raster::mask(dsm22, ageUse)
dsm23 <- raster::mask(dsm23, ageUse)

#crop to ensure each raster has same extent
dsm20 <- raster::crop(dsm20, raster::extent(ageUse))
dsm21 <- raster::crop(dsm21, raster::extent(ageUse))
dsm22 <- raster::crop(dsm22, raster::extent(ageUse))
dsm23 <- raster::crop(dsm23, raster::extent(ageUse))

#subtract ground elevation
chm20 <- dsm20-dem
chm21 <- dsm21-dem
chm22 <- dsm22-dem
chm23 <- dsm23-dem

#make sure all years have the same extent
chm21 <- raster::crop(chm21, raster::extent(chm21))
chm22 <- raster::crop(chm22, raster::extent(chm22))
chm23 <- raster::crop(chm23, raster::extent(chm22))

#calculate the change in canopy height
d22to21 <- chm22-chm21
d23to22 <- chm23-chm22

#save rasters of canopy height change
raster::writeRaster(d22to21, file="D:/TroubleshootingV/chms/dCHM22to21.tif",  overwrite=TRUE)
raster::writeRaster(d23to22, file= "D:/TroubleshootingV/chms/dCHMF23to22.tif",  overwrite=TRUE)

#install forest gaps and define function
install.packages("ForestGapR")

getForestGaps <- function (chm_layer, threshold = 10, size = c(1, 10^4)) 
{
  chm_layer[chm_layer > threshold] <- NA
  chm_layer[chm_layer <= threshold] <- 1
  gaps <- raster::clump(chm_layer, directions = 4, gap = FALSE)
  rcl <- raster::freq(gaps)
  rcl[, 2] <- rcl[, 2] * raster::res(chm_layer)[1]^2
  rcl <- cbind(rcl[, 1], rcl)
  z <- raster::reclassify(gaps, rcl = rcl, right = NA)
  z[is.na(gaps)] <- NA
  gaps[z > size[2]] <- NA
  gaps[z < size[1]] <- NA
  gaps <- raster::clump(gaps, directions = 4, gap = FALSE)
  names(gaps) <- "gaps"
  return(gaps)
}


# Define gap height threshold, min gap size, and max gap size
gapSzMin <- 25
gapSzMax <- 10^6
gapHtThresh <- -5

# Identify gaps  
gaps22to21 <- getForestGaps(d22to21,
                            threshold = gapHtThresh ,
                            size=c(gapSzMin,gapSzMax))
gaps23to22 <- getForestGaps(d23to22,
                            threshold = gapHtThresh ,
                            size=c(gapSzMin,gapSzMax))
# Create a Spatial Polygon Data Frame object, where each polygon is a gap
gaps22to21sp <- ForestGapR::GapSPDF(gaps22to21)

gaps23to22sp <- ForestGapR::GapSPDF(gaps23to22)

#write shapefiles

rgdal::writeOGR(gaps22to21sp,
                dsn = "gaps22to21_shapefile",
                layer = "gaps22to21sp", 
                driver = "ESRI Shapefile")


rgdal::writeOGR(gaps23to22sp,
                dsn = "gapsF23to22_shapefile",
                layer = "gapsF23to22sp", 
                driver = "ESRI Shapefile")






