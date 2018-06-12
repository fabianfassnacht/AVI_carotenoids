library(rgdal)
library(raster)
library(maptools)
###Einlesen###
hyp<-stack("U:/EnMap/AA_BaumartenErkennung/Hyperspektralszenen/Karlsruhe/100820_Karlsruhe_02_rad_atm_geo_masked.bsq")
shp<-readOGR("U:/EnMap/AA_BaumartenErkennung/ReferenceDatasets/Karlsruhe/Samples_Baumarten/Final_Samples_Shapes_60",layer="Samples_60.shp")
###Projektion, Koordinatensystem anpassen###
#shp<-spTransform(shp,CRS=CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
###Plot###
plot(hyp,10,col=gray.colors(100))
plot(shp,col="yellow",add=T)
###Spektren extrahieren###
ana<-extract(hyp,shp)
write.csv(ana,file="U:/EnMap/AA_BaumartenErkennung/Hyperspektralszenen/Karlsruhe/training_ka.csv")
###Testplot###
plot(c(1:125),ana[1,],type="l", col="red")