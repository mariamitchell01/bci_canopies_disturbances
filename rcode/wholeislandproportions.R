#Mia Mitchell


#import gap shapefiles
gaps2021 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps20to21_shapefilenew/gaps20to21spnew.shp")
gaps2020 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps19to20_shapefilenew/gaps19to20spnew.shp")
gaps2019 <- sf::st_read("D:/Data_HeightRasters/Data_HeightRasters/Data_GapShapefiles/gaps18to19_shapefilenew/gaps18to19spnew.shp")


#import bci shapefile
bcishape <- sf::st_read("D:/Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
bcishape <- sf::st_transform(bcishape,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
bciarea <- sf::st_area(bcishape)

#calculate yearly disturbance rates

gaps2021_area <- gaps2021$area
gaps2021_sum <- sum(gaps2021_area)
prop2021 <- gaps2021_sum / bciarea
print(prop2021)

gaps2020_area <- gaps2020$area
gaps2020_sum <- sum(gaps2020_area)
prop2020 <- gaps2020_sum / bciarea
print(prop2020)

gaps2019_area <- gaps2019$area
gaps2019_sum <- sum(gaps2019_area)
prop2019 <- gaps2019_sum / bciarea
print(prop2019)
