install.packages("rgeos")

#import gap shapefiles
gaps2021 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps20to21_shapefilenew/gaps20to21spnew.shp")
gaps2020 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps19to20_shapefilenew/gaps19to20spnew.shp")
gaps2019 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps18to19_shapefilenew/gaps18to19spnew.shp")


gaps2021$area > 40
# Filter gaps based on area 25-35
gaps2019_25to35 <- gaps2019[gaps2019$area > 25 & gaps2019$area < 35,]

gaps2019_35to45 <- gaps2019[gaps2019$area > 35 & gaps2019$area < 45,]

gaps2019_45to55 <- gaps2019[gaps2019$area > 45 & gaps2019$area < 55,]

gaps2019_grt55 <- gaps2019[gaps2019$area > 55,]

#write shapefiles

sf::st_write(gaps2019_25to35, dsn = "C:/Users/mitchellme/OneDrive - Smithsonian Institution/Desktop/Processing",layer = "gaps2019_25to35", driver = "ESRI Shapefile")

sf::st_write(gaps2019_35to45, dsn = "C:/Users/mitchellme/OneDrive - Smithsonian Institution/Desktop/Processing",layer = "gaps2019_35to45", driver = "ESRI Shapefile")

sf::st_write(gaps2019_45to55, dsn = "C:/Users/mitchellme/OneDrive - Smithsonian Institution/Desktop/Processing",layer = "gaps2019_45to55", driver = "ESRI Shapefile")

sf::st_write(gaps2019_grt55, dsn = "C:/Users/mitchellme/OneDrive - Smithsonian Institution/Desktop/Processing",layer = "gaps2019_grt55", driver = "ESRI Shapefile")