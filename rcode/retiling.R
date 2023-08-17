
#original author: kc cushman
#edited by mia mitchell & vicente vasquez

# Although not a full solution -- to still access the
the packages that kc used at the time, you can use the groundhog package
# Some functions nonetheless that are depracted are replaced
but some comparable functions have not be updated in lidR that I am aware of


# Install packages

install.packages("groundhog")
groundhog::groundhog.library("lidR","2020-11-01")
groundhog::groundhog.library("raster", "2020-09-01")
groundhog::groundhog.library("sp","2020-09-01")
groundhog::groundhog.library("rgdal","2020-09-01")


# Read in the catalogs

cat18 <- lidR::catalog("/pointclouds/2021_06_29_BCI_Dipteryx.las")
cat19 <- lidR::catalog("/pointclouds/2021_06_29_BCI_Dipteryx.las")
cat20 <- lidR::catalog("/pointclouds/2021_06_29_BCI_Dipteryx.las")
cat21 <- lidR::catalog("/pointclouds/2021_06_29_BCI_Dipteryx.las")
cat22 <- lidR::catalog("/pointclouds/2022_07_21_medium_1.las")
cat23 <- lidR::catalog("/pointclouds/BCI_whole_2023_06_19_medium.las")


#### Create grid for new .las tiles

xmin <- 623400
xmax <- 630100
ymin <- 1009700
ymax <- 1015200

tileSz <- 150
xmins <- seq(xmin,(xmax-tileSz),tileSz)
xmaxs <- seq((xmin+tileSz),xmax,tileSz)
ymins <- seq(ymin,(ymax-tileSz),tileSz)
ymaxs <- seq((ymin+tileSz),ymax,tileSz)

gridInfo <- data.frame(xmin=rep(xmins, length(ymins)),
xmax=rep(xmaxs, length(ymins)),
ymin=rep(ymins, each=length(xmins)),
ymax=rep(ymaxs, each=length(xmins)))

# Only keep grid cells that actually overlay BCI, using BCI soils shapefile

soils <- rgdal::readOGR("Data_Ancillary/BCI_Soils/BCI_Soils.shp")
soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

gridInfo$Use <- NA

for(i in 1:dim(gridInfo)[1]){
poly_i <- as(raster::extent(as.numeric(gridInfo[i,1:4])),
'SpatialPolygons')

sp::proj4string(poly_i) <- sp::proj4string(soils)
test_i <- sp::over(x = poly_i, y = soils)
if(is.na(test_i$ID[1])){
gridInfo$Use[i] <- F}
if(!is.na(test_i$ID[1])){
gridInfo$Use[i] <- T}
}

# Only keep tiles that overlap BCI soils shapefile

gridInfo <- gridInfo[gridInfo$Use==T,]
gridInfo$ID <- 1:dim(gridInfo)[1]

write.csv(gridInfo, row.names = F, file = ("/Data_Ancillary/gridInfo.csv"))

# Retiling

overlap <- 30

#Retile 2023 photogrammetry point cloud with overlap


# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat23f, 
                                 xleft = gridInfo$xmin[i] - overlap,
                                 ybottom = gridInfo$ymin[i] - overlap,
                                 xright = gridInfo$xmax[i] + overlap,
                                 ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("/BCI23Tiles/BCI23_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat23f, 
                                 xleft = gridInfo$xmin[i] - overlap,
                                 ybottom = gridInfo$ymin[i] - overlap,
                                 xright = gridInfo$xmax[i] + overlap,
                                 ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI23Tiles_dec/BCI23d_",i,".laz"))
  }
  
  print(i)
}

# Retile 2022 photogrammetry point cloud with overlap 

# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat22, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI22Tiles/BCI22_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat22, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI22Tiles_dec/BCI22d_",i,".laz"))
  }
  
  print(i)
}

# Retile 2021 photogrammetry point cloud with overlap 

# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat21, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI22Tiles/BCI21_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat21, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI21Tiles_dec/BCI21d_",i,".laz"))
  }
  
  print(i)
}

# Retile 2020 photogrammetry point cloud with overlap 

# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat20, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI20Tiles/BCI20_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat20, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI20Tiles_dec/BCI20d_",i,".laz"))
  }
  
  print(i)
}


# Retile 2019 photogrammetry point cloud with overlap 

# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat19, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI19Tiles/BCI19_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat19, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI19Tiles_dec/BCI19d_",i,".laz"))
  }
  
  print(i)
}

# Retile 2018 photogrammetry point cloud with overlap 

# full-resolution
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat18, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI19Tiles/BCI18_",i,".laz"))
  }
  
  print(i)
}

#decimated
for(i in 1:dim(gridInfo)[1]){
  data <- lidR::clip_rectangle(cat18, 
                               xleft = gridInfo$xmin[i] - overlap,
                               ybottom = gridInfo$ymin[i] - overlap,
                               xright = gridInfo$xmax[i] + overlap,
                               ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("D:/TroubleshootingV/tiles/BCI18Tiles_dec/BCI18d_",i,".laz"))
  }
  
  print(i)
}


#### Rename transition matrix files ####

# renames output from running do22to21.bat before running apply2022transformation.bat

# 2022

for(i in 1:dim(gridInfo)[1]){
  file.rename(from = list.files("D:/TroubleshootingV/tiles/BCI22Tiles_dec/",
                                full.names = T, pattern = paste0("BCI22d_",i,"_REG")),
              to = paste0("D:/TroubleshootingV/tiles/BCI22Tiles_dec/BCI22mat_",i,".txt"))
}

# renames output from running do23to22.bat before running apply2023transformation.bat
# 2023

for(i in 1:dim(gridInfo)[1]){
  file.rename(from = list.files("D:/TroubleshootingV/tiles/BCI23FTiles_dec/",
                                full.names = T, pattern = paste0("BCI23fd_",i,"_REG")),
              to = paste0("D:/TroubleshootingV/tiles/BCI23FTiles_dec/BCI23fmat_",i,".txt"))
}


#### Remove overlap from aligned point clouds ####

# THESE DATA FOLDERS ARE ARCHIVED


# 2022
for(i in 1:dim(gridInfo)[1]){
  
  data <- lidR::clip_rectangle(las = lidR::readLAS("D:/TroubleshootingV/tiles/BCI22Tiles_alignedto21Full/BCI22af_",i,".las"),
                               xleft=gridInfo$xmin[i],
                               xright=gridInfo$xmax[i],
                               ybottom=gridInfo$ymin[i],
                               ytop=gridInfo$ymax[i])
  if(length(data@data$X)>0){
    lidR::writeLAS(las = data,
                   file = paste0("D:/TroubleshootingV/tiles/BCI22Tiles_alignedto21Trim/BCI22at_",i,".laz"))
  }
}

#2023

for (i in 1:dim(gridInfo)[1]) {
  data <- lidR::clip_rectangle(
    las = lidR::readLAS(paste0("D:/TroubleshootingV/tiles/BCI23FTiles_alignedto22Full/BCI23faf_", as.character(i), ".las")),
    xleft = gridInfo$xmin[i],
    xright = gridInfo$xmax[i],
    ybottom = gridInfo$ymin[i],
    ytop = gridInfo$ymax[i]
  )
  if (length(data@data$X) > 0) {
    lidR::writeLAS(
      las = data,
      file = paste0("D:/TroubleshootingV/tiles/BCI23FTiles_alignedto22Trim/BCI23fat_", as.character(i), ".laz")
    )
  }
}
