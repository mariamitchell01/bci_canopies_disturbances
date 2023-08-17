# Define new canopy disturbances
#### Packages ####
  
  # Install the following packages, if needed
  
   install.packages("ForestGapR") # version 0.0.3 used
   install.packages("raster") # version 3.3-13 used
   install.packages("rgdal") # version 1.5-16 used
   install.packages("sp") # version 1.4-2 used

#### READ RAW DATA ####

# Read in digital surface model rasters:
  dsm09 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2009.tif")
  #dsm15 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2015_corrected.tif")
  dsm18 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2018_corrected.tif")
  dsm19 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2019_corrected.tif")
  dsm20 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2020_corrected.tif")
  dsm21 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/DSM_2021_corrected.tif")
  
# Read in forest age (used to exclude recent clearings)
  age <- rgdal::readOGR("C:/Users/mitchellme/Documents/Data_Ancillary/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  age$AgeClass <- "Other"
  age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
  age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
  ageUse <- age[!(age$AgeClass=="Other"),]

# Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("C:/Users/mitchellme/Documents/Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
  
# Read in BCI DEM (from 2009 lidar data)
  dem <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/LidarDEM_BCI.tif")
  dem <- raster::crop(dem, raster::extent(ageUse))
  dem <- raster::resample(dem,dsm19)
 
#### PROCESS AND SAVE CANOPY HEIGHT DATA ####
  
# Remove raster areas outside BCI perimeter (exclude within 25 m of lake)
 dsm09 <- raster::mask(dsm09, buffer)  
# dsm15 <- raster::mask(dsm15, buffer)  
  dsm18 <- raster::mask(dsm18, buffer)
  dsm19 <- raster::mask(dsm19, buffer)
  dsm20 <- raster::mask(dsm20, buffer)
  dsm21 <- raster::mask(dsm21, buffer)

# Remove raster areas in clearings
  dsm09 <- raster::mask(dsm09, ageUse)  
#  dsm15 <- raster::mask(dsm15, ageUse)  
   dsm18 <- raster::mask(dsm18, ageUse)
   dsm19 <- raster::mask(dsm19, ageUse)
   dsm20 <- raster::mask(dsm20, ageUse)
   dsm21 <- raster::mask(dsm21, ageUse)

# Crop to ensure each raster has same extent
  dsm09 <- raster::crop(dsm09, raster::extent(ageUse))  
#  dsm15 <- raster::crop(dsm15, raster::extent(ageUse))  
   dsm18 <- raster::crop(dsm18, raster::extent(ageUse)) 
   dsm19 <- raster::crop(dsm19, raster::extent(ageUse))
   dsm20 <- raster::crop(dsm20, raster::extent(ageUse))
   dsm21 <- raster::crop(dsm21, raster::extent(ageUse))
  
# Subtract ground elevation
  chm09 <- dsm09-dem
#  chm15 <- dsm15-dem
   chm18 <- dsm18-dem
   chm19 <- dsm19-dem  
   chm20 <- dsm20-dem
   chm21 <- dsm21-dem
  
# Make sure all years have the same extent
#  chm15 <- raster::crop(chm15, raster::extent(chm09))
   chm18 <- raster::crop(chm18, raster::extent(chm09))
   chm19 <- raster::crop(chm19, raster::extent(chm09))  
   chm20 <- raster::crop(chm20, raster::extent(chm09))
   chm21 <- raster::crop(chm21, raster::extent(chm09))
  
# Load cloud and QAQC masks for 2015-2020
#  cloudMask15 <- raster::raster("Data_QAQC/CloudMasks/CloudMask_2015.tif")
  cloudMask18 <- raster::raster("Data_QAQC/CloudMasks/CloudMask_2018.tif")
  cloudMask19 <- raster::raster("C:/Users/mitchellme/Documents/ArcGIS/Projects/cloudmasks/Cloudmask_2019.tif")
 # cloudMask20 <- raster::raster("Data_QAQC/CloudMasks/CloudMask_2020.tif")
  cloudMask21 <- raster::raster("C:/Users/mitchellme/Documents/ArcGIS/Projects/cloudmasks/newCloudmask_2021.tif")
  
  #qaqcMask15 <- raster::raster("Data_QAQC/QAQCMask_2015.tif")
  #qaqcMask18 <- raster::raster("Data_QAQC/QAQCMask_2018.tif")
  qaqcMask19 <- raster::raster("C:/Users/mitchellme/Documents/ArcGIS/Projects/cloudmasks/QAQCMask_2019.tif")
  #qaqcMask20 <- raster::raster("Data_QAQC/QAQCMask_2020.tif")
  qaqcMask21 <- raster::raster("C:/Users/mitchellme/Documents/ArcGIS/Projects/cloudmasks/QAQCMask_2021.tif")
    
  # Resample to extent and resolution of CHMs
    #cloudMask15 <- raster::resample(cloudMask15, chm15)
   # cloudMask18 <- raster::resample(cloudMask18, chm18)
    cloudMask19 <- raster::resample(cloudMask19, chm19)
  #  cloudMask20 <- raster::resample(cloudMask20, chm20)
    cloudMask21 <- raster::resample(cloudMask21, chm21)
    
    qaqcMask15 <- raster::resample(qaqcMask15, chm15)
    #qaqcMask18 <- raster::resample(qaqcMask18, chm18)
    qaqcMask19 <- raster::resample(qaqcMask19, chm19)
    #qaqcMask20 <- raster::resample(qaqcMask20, chm20)
    qaqcMask21 <- raster::resample(qaqcMask21, chm21)
  # What % of pixels have clouds or fail QAQC?
    # Clouds
     # 100*length(which(!is.na(raster::values(cloudMask15)) & raster::values(cloudMask15)<=0.99))/length(which(!is.na(raster::values(cloudMask15))))
     # 100*length(which(!is.na(raster::values(cloudMask18)) & raster::values(cloudMask18)<=0.99))/length(which(!is.na(raster::values(cloudMask18))))
      100*length(which(!is.na(raster::values(cloudMask19)) & raster::values(cloudMask19)<=0.99))/length(which(!is.na(raster::values(cloudMask19))))
      #100*length(which(!is.na(raster::values(cloudMask20)) & raster::values(cloudMask20)<=0.99))/length(which(!is.na(raster::values(cloudMask20))))
      100*length(which(!is.na(raster::values(cloudMask21)) & raster::values(cloudMask21)<=0.99))/length(which(!is.na(raster::values(cloudMask21))))
    # QAQC
      #100*length(which(!is.na(raster::values(cloudMask15)) & raster::values(qaqcMask15)<=0.99))/length(which(!is.na(raster::values(cloudMask15))))
      #100*length(which(!is.na(raster::values(cloudMask18)) & raster::values(qaqcMask18)<=0.99))/length(which(!is.na(raster::values(cloudMask18))))
      100*length(which(!is.na(raster::values(cloudMask19)) & raster::values(qaqcMask19)<=0.99))/length(which(!is.na(raster::values(cloudMask19))))
     # 100*length(which(!is.na(raster::values(cloudMask20)) & raster::values(qaqcMask20)<=0.99))/length(which(!is.na(raster::values(cloudMask20))))
      100*length(which(!is.na(raster::values(cloudMask21)) & raster::values(qaqcMask21)<=0.99))/length(which(!is.na(raster::values(cloudMask21))))
      
  # Remove masked pixels
   # chm15[!(cloudMask15>0.99)] <- NA
   # chm18[!(cloudMask18>0.99)] <- NA
    chm19[!(cloudMask19>0.99)] <- NA
   # chm20[!(cloudMask20>0.99)] <- NA
    chm21[!(cloudMask21>0.99)] <- NA
    
  #  chm15[!(qaqcMask15>0.99)] <- NA
   # chm18[!(qaqcMask18>0.99)] <- NA
    chm19[!(qaqcMask19>0.99)] <- NA
   # chm20[!(qaqcMask20>0.99)] <- NA
    chm21[!(qaqcMask21>0.99)] <- NA
  
# Do 2015 canopy height correction based on 2009 and 2018 values
    
    vals09 <- raster::values(chm09)
    vals15 <- raster::values(chm15)
    vals18 <- raster::values(chm18)
    
    toChange <- which((vals15-vals09) > 5 & (vals18-vals15) < -1)
    newVals <- vals09[toChange]
    chm15c <- chm15
    raster::values(chm15c)[toChange] <- newVals
       
#  Save
   # raster::writeRaster(chm15, "Data_HeightRasters/CHM_2015_QAQC_wBias.tif")
  #  raster::writeRaster(chm15c, "Data_HeightRasters/CHM_2015_QAQC.tif")
  #  raster::writeRaster(chm18, "Data_HeightRasters/CHM_2018_QAQC.tif")
    raster::writeRaster(chm19, "C:/Users/mitchellme/Documents/Data_HeightRasters/CHM_2019_QAQC.tif", overwrite = TRUE)
   # raster::writeRaster(chm20, "Data_HeightRasters/CHM_2020_QAQC.tif")
    raster::writeRaster(chm21, "C:/Users/mitchellme/Documents/Data_HeightRasters/CHM_2021_QAQC.tif", overwrite = TRUE)
    
    chm18 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/CHM_2018_QAQC.tif")
    chm20 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/CHM_2020_QAQC.tif")
  
#### MAKE AND SAVE CANOPY HEIGHT CHANGE RASTERS ####  
  
  # Calculate the change in canopy height for each interval  
#    d15to18 <- chm18-chm15c
#    d18to20 <- chm20-chm18
     d18to19 <- chm19-chm18
     d19to20 <- chm20-chm19
     d20to21 <- chm21-chm20
    
  # Make additional rasters where values below an initial height threshold are
  # omitted 
#    d15to18tall <- d15to18
#    d18to20tall <- d18to20
     d18to19tall <- d18to19
    d19to20tall <- d19to20
     d20to21tall <- d20to21
    
  # Mask out areas that are initially < 10 m in height and decrease between two years
    shortThresh <- 10
    
    # Find proportion of cells that are < 10 m
#    100*length(which(raster::values(chm15c)<shortThresh & !is.na(raster::values(chm15c))))/length(which(!is.na(raster::values(chm15c))))
     100*length(which(raster::values(chm18)<shortThresh & !is.na(raster::values(chm18))))/length(which(!is.na(raster::values(chm18))))
    100*length(which(raster::values(chm19)<shortThresh & !is.na(raster::values(chm19))))/length(which(!is.na(raster::values(chm19))))
     100*length(which(raster::values(chm20)<shortThresh & !is.na(raster::values(chm20))))/length(which(!is.na(raster::values(chm20))))
     100*length(which(raster::values(chm21)<shortThresh & !is.na(raster::values(chm21))))/length(which(!is.na(raster::values(chm21))))
    
    
    # 2015 - 2018
    #  short15 <- rep(0, length(raster::values(chm15c)))
    #  short15[raster::values(chm15c)<shortThresh & !is.na(raster::values(chm15c))] <- 1
    #  d15to18tall@data@values[short15==1] <- NA
        
    # 2018 - 2020
    #  short18 <- rep(0, length(raster::values(chm18)))
    #  short18[raster::values(chm18)<shortThresh & !is.na(raster::values(chm18))] <- 1
    #    d18to20tall@data@values[short18==1] <- NA  
     
    # 2018 - 2019
       short18 <- rep(0, length(raster::values(chm18)))
       short18[raster::values(chm18)<shortThresh & !is.na(raster::values(chm18))] <- 1
        d18to19tall@data@values[short18==1] <- NA 
     
    # 2019 - 2020
      short19 <- rep(0, length(raster::values(chm19)))
     short19[raster::values(chm19)<shortThresh & !is.na(raster::values(chm19))] <- 1
        d19to20tall@data@values[short19==1] <- NA
     
    # 2020 - 2021
      short20 <- rep(0, length(raster::values(chm20)))
      short20[raster::values(chm20)<shortThresh & !is.na(raster::values(chm20))] <- 1
        d20to21tall@data@values[short20==1] <- NA
    
      
  # Save rasters of canopy height change, with and without height mask
     
     #raster::writeRaster(d15to18, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM15to18.tif")
     #raster::writeRaster(d18to20, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM18to20.tif")
     raster::writeRaster(d18to19, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM18to19.tif", overwrite=TRUE)
     raster::writeRaster(d19to20, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM19to20.tif", overwrite=TRUE)
     raster::writeRaster(d20to21, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM20to21.tif",  overwrite=TRUE)
     
     #raster::writeRaster(d15to18tall, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM15to18tall.tif")
     #raster::writeRaster(d18to20tall, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM18to20tall.tif")
     raster::writeRaster(d18to19tall, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM18to19tall.tif",  overwrite=TRUE)
     raster::writeRaster(d19to20tall, file="C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM19to20tall.tif",  overwrite=TRUE)
     raster::writeRaster(d20to21tall, file= "C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM20to21tall.tif",  overwrite=TRUE)
     
     #Read rasters
     
     d18to19 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM18to19.tif")
     d19to20 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM19to20.tif")
     d20to21 <- raster::raster("C:/Users/mitchellme/Documents/Data_HeightRasters/dCHM20to21.tif")
     
    
        
#### BINARY NEW GAPS #### 
  
  # Use ForestGapR package to delineate new gaps
    
    # Define new function (modified from original ForestGapR function) that does not
     # include diagonal pixels in the same gap
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
      
      # 2015 - 2018
      
      # Identify gaps  
      gaps15to18 <- getForestGaps(d15to18,
                                  threshold = gapHtThresh ,
                                  size=c(gapSzMin,gapSzMax))
      
      # Create a Spatial Polygon Data Frame object, where each polygon is a gap
      gaps15to18sp <- ForestGapR::GapSPDF(gaps15to18)
      
      # Calculate the area and perimeter from each gap object
      gaps15to18sp@data$area <- NA
      gaps15to18sp@data$perimeter <- NA
      for(i in 1:length(gaps15to18sp)){
        gaps15to18sp[gaps15to18sp$gap_id==i,"area"] <- raster::area(gaps15to18sp[gaps15to18sp$gap_id==i,])
        
        perims_j <- c()
        for(j in 1:length(gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons)){
          
          coordsj <- gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
          
          lengths_k <- c()
          for(k in 2:dim(coordsj)[1]){
            lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
            
          }
          perims_j[j] <- sum(lengths_k)
          
          if(gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
            perims_j[j] <- 0
          }
          
        }
        gaps15to18sp[gaps15to18sp$gap_id==i,"perimeter"] <- sum(perims_j)
      }
      

      
  
    # 2018 - 2020
        
      # Identify gaps  
      gaps18to20 <- getForestGaps(d18to20,
                                  threshold = gapHtThresh ,
                                  size=c(gapSzMin,gapSzMax))
      
        # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps18to20sp <- ForestGapR::GapSPDF(gaps18to20)
        
        # Calculate the area and perimeter from each gap object
        gaps18to20sp@data$area <- NA
        gaps18to20sp@data$perimeter <- NA
        for(i in 1:length(gaps18to20sp)){
          gaps18to20sp[gaps18to20sp$gap_id==i,"area"] <- raster::area(gaps18to20sp[gaps18to20sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps18to20sp[gaps18to20sp$gap_id==i,"perimeter"] <- sum(perims_j)
        }
        
        # 2018 - 2019
        
        # Identify gaps  
        gaps18to19 <- getForestGaps(d18to19,
                                    threshold = gapHtThresh ,
                                    size=c(gapSzMin,gapSzMax))
        
        # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps18to19sp <- ForestGapR::GapSPDF(gaps18to19)
        
        # Calculate the area and perimeter from each gap object
        gaps18to19sp@data$area <- NA
        gaps18to19sp@data$perimeter <- NA
        for(i in 1:length(gaps18to19sp)){
          gaps18to19sp[gaps18to19sp$gap_id==i,"area"] <- raster::area(gaps18to19sp[gaps18to19sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps18to19sp[gaps18to19sp$gap_id==i,"perimeter"] <- sum(perims_j)
        }
        
        # 2019 - 2020
        
        # Identify gaps  
        gaps19to20 <- getForestGaps(d19to20,
                                    threshold = gapHtThresh ,
                                    size=c(gapSzMin,gapSzMax))
        
        # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps19to20sp <- ForestGapR::GapSPDF(gaps19to20)
        
        # Calculate the area and perimeter from each gap object
        gaps19to20sp@data$area <- NA
        gaps19to20sp@data$perimeter <- NA
        for(i in 1:length(gaps19to20sp)){
          gaps19to20sp[gaps19to20sp$gap_id==i,"area"] <- raster::area(gaps19to20sp[gaps19to20sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps19to20sp[gaps19to20sp$gap_id==i,"perimeter"] <- sum(perims_j)
        }
        
        # 2020 - 2021
        
        # Identify gaps  
        gaps20to21 <- getForestGaps(d20to21,
                                    threshold = gapHtThresh ,
                                    size=c(gapSzMin,gapSzMax))
        
        # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps20to21sp <- ForestGapR::GapSPDF(gaps20to21)
        
        # Calculate the area and perimeter from each gap object
        gaps20to21sp@data$area <- NA
        gaps20to21sp@data$perimeter <- NA
        for(i in 1:length(gaps20to21sp)){
          gaps20to21sp[gaps20to21sp$gap_id==i,"area"] <- raster::area(gaps20to21sp[gaps20to21sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps20to21sp[gaps20to21sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps20to21sp[gaps20to21sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps20to21sp[gaps20to21sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps20to21sp[gaps20to21sp$gap_id==i,"perimeter"] <- sum(perims_j)
        }
        
        
    
  # Save rasters of new gap pixels
       #raster::writeRaster(gaps15to18, file="Data_HeightRasters/newGaps15to18.tif")
       #raster::writeRaster(gaps18to20, file="Data_HeightRasters/newGaps18to20.tif")
       raster::writeRaster(gaps18to19, file="C:/Users/mitchellme/Documents/Data_HeightRasters/newGaps18to19.tif", overwrite=TRUE)
       raster::writeRaster(gaps19to20, file="C:/Users/mitchellme/Documents/Data_HeightRasters/newGaps19to20.tif", overwrite=TRUE)
       raster::writeRaster(gaps20to21, file="C:/Users/mitchellme/Documents/Data_HeightRasters/newGaps20to21.tif", overwrite=TRUE)
      

#### ASSIGN GAP POLYGONS TO BOOTSTRAP GROUPS ####
  blocks <- read.csv("C:/Users/mitchellme/Documents/Data_Ancillary/bootstrapBlocks.csv")
       
  #gaps15to18sp$block <- NA 
  #gaps18to20sp$block <- NA
  gaps18to19sp$block <- NA
  gaps19to20sp$block <- NA
  gaps20to21sp$block <- NA
  
  for(i in 1:dim(blocks)[1]){
    gaps18to19sp@data[gaps18to19sp@data[,1]>=blocks$xmin[i] & gaps18to19sp@data[,1] < blocks$xmax[i]
                      & gaps18to19sp@data[,2]>=blocks$ymin[i] & gaps18to19sp@data[,2] < blocks$ymax[i],"block"] <- i
   
    gaps19to20sp@data[gaps19to20sp$X1>=blocks$xmin[i] & gaps19to20sp$X1 < blocks$xmax[i]
                      & gaps19to20sp$X2>=blocks$ymin[i] & gaps19to20sp$X2 < blocks$ymax[i],"block"] <- i
     
    gaps20to21sp@data[gaps20to21sp$X1>=blocks$xmin[i] & gaps20to21sp$X1 < blocks$xmax[i]
                       & gaps20to21sp$X2>=blocks$ymin[i] & gaps20to21sp$X2 < blocks$ymax[i],"block"] <- i
  }  
    
    
#### CALCULATE MEAN HT CHANGE IN EACH GAP #####
  
  #gaps15to18sp$htDrop <- NA
  #valsGap <- raster::values(gaps15to18)
  #valsHt <- raster::values(d15to18)
  
  #for(i in 1:length(gaps15to18sp)){
    #gapCells <- which(valsGap==gaps15to18sp$gap_id[i] & !is.na(valsGap))
   # htCells <- valsHt[gapCells]
   # gaps15to18sp$htDrop[i] <- mean(htCells)
  #}
  
  gaps18to19sp$htDrop <- NA
  valsGap <- raster::values(gaps18to19)
  valsHt <- raster::values(d18to19)
  
  for(i in 1:length(gaps18to19sp)){
    gapCells <- which(valsGap==gaps18to19sp$gap_id[i] & !is.na(valsGap))
    htCells <- valsHt[gapCells]
    gaps18to19sp$htDrop[i] <- mean(htCells)
  }
  
  gaps19to20sp$htDrop <- NA
  valsGap <- raster::values(gaps19to20)
  valsHt <- raster::values(d19to20)
  
  for(i in 1:length(gaps19to20sp)){
    gapCells <- which(valsGap==gaps19to20sp$gap_id[i] & !is.na(valsGap))
    htCells <- valsHt[gapCells]
    gaps19to20sp$htDrop[i] <- mean(htCells)
  }
  
  gaps20to21sp$htDrop <- NA
  valsGap <- raster::values(gaps20to21)
  valsHt <- raster::values(d20to21)
  
  for(i in 1:length(gaps20to21sp)){
    gapCells <- which(valsGap==gaps20to21sp$gap_id[i] & !is.na(valsGap))
    htCells <- valsHt[gapCells]
    gaps20to21sp$htDrop[i] <- mean(htCells)
  }
  
#### SAVE GAP POLYGONS ####
    
    # Save gap shapefiles
   
    rgdal::writeOGR(gaps18to19sp,
                    dsn = "gaps18to19_shapefilenew",
                    layer = "gaps18to19spnew", 
                    driver = "ESRI Shapefile")
    
    rgdal::writeOGR(gaps19to20sp,
                    dsn = "gaps19to20_shapefilenew",
                    layer = "gaps19to20spnew", 
                    driver = "ESRI Shapefile")
    rgdal::writeOGR(gaps20to21sp,
                    dsn = "gaps20to21_shapefilenew",
                    layer = "gaps20to21spnew", 
                    driver = "ESRI Shapefile")
  
    rgdal::writeOGR(gaps15to18sp,
                    dsn = "gaps15to18_shapefile",
                    layer = "gaps15to18sp", 
                    driver = "ESRI Shapefile")
    
    rgdal::writeOGR(gaps18to20sp,
                    dsn = "gaps18to20_shapefile",
                    layer = "gaps18to20sp", 
                    driver = "ESRI Shapefile")
    
                    
    