pkgs<-c("rgdal","caret","raster","foreign", "randomForest")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("F:/KIT_Forschung/11_AG_Angular_Index/1_new_runs")

# load data
# read full table
spectra <- read.table("corrected_data_ag.csv", sep="\t", head=T)
head(spectra)
#spectra <- spectra[1:90,]


# extract target parameter total chlorophyll and carotenoid content

species <- spectra[,1]
chl <- spectra[,4]
car <- spectra[,5]
refl_o <- spectra[,6:43]
wavel <- seq(410, 780, 10)

species

###############
###############
### Set dataset
###############
###############

#### maple

car_m <- car[1:38]
chl_m <- chl[1:38]
refl_m <- refl_o[1:38,]
refl_m <- refl_m * 100

#### chestnut

car_c <- car[39:56]
chl_c <- chl[39:56]
refl_c <- refl_o[39:56,]
refl_c <- refl_c*100

#### beech

car_b <- car[57:90]
chl_b <- chl[57:90]
refl_b <- refl_o[57:90,]
refl_b <- refl_b*100




#### plot CHL vs. CAR


setwd("F:/KIT_Forschung/11_AG_Angular_Index/99_Paper/6_chl_vs_car_plots")
pdf(file="chl_car_plots_fit.pdf", width=9, height=3)
par(mfrow=c(1,3))
#plot(car_m, chl_m, main="Maple leaves", ylab="tot. chl. content [mg/m²]", xlab="tot. car. content [mg/m²")
#plot(car_c, chl_c, main="Chestnut leaves", ylab="tot. chl. content [mg/m²]", xlab="tot. car. content [mg/m²")
#plot(car_b, chl_b, main="Beech leaves", ylab="tot. chl. content [mg/m²]", xlab="tot. car. content [mg/m²")

plot(car_m, chl_m, lty=2, ylab="", xlab="", cex=0.8, main="Maple leaves")
mtext("tot. car. content [mg/m²]",side=1,col="black",line=2.5, cex=1)
mtext("tot. chl. content [mg/m²]",side=2,col="black",line=2.5, cex=1)
r2 <- (cor(chl_m, car_m))^2
r2_r <- round(r2, digits=2)
mtext(side=1, line=-1.0, text=paste("r²: ", r2_r, "  "), padj=-0.2, adj=1)
#fit line
f1 <- lm(chl_m ~ car_m)
cf1 <- coef(f1)
newenv <- seq(0,180, (180/500))
newenv2 <- cf1[1]+cf1[2]*newenv
lines(newenv, newenv2, lty=2, col="blue")

plot(car_c, chl_c, lty=2, ylab="", xlab="", cex=0.8, main="Chestnut leaves")
mtext("tot. car. content [mg/m²]",side=1,col="black",line=2.5, cex=1)
mtext("tot. chl. content [mg/m²]",side=2,col="black",line=2.5, cex=1)
r2 <- (cor(chl_c, car_c))^2
r2_r <- round(r2, digits=2)
mtext(side=1, line=-1.0, text=paste("r²: ", r2_r, "  "), padj=-0.2, adj=1)
#fit line
f2 <- lm(chl_c ~ car_c)
cf2 <- coef(f2)
newenv <- seq(0,180, (180/500))
newenv2 <- cf2[1]+cf2[2]*newenv
lines(newenv, newenv2, lty=2, col="blue")


plot(car_b, chl_b, lty=2, ylab="", xlab="", cex=0.8, main="Beech leaves")
mtext("tot. car. content [mg/m²]",side=1,col="black",line=2.5, cex=1)
mtext("tot. chl. content [mg/m²]",side=2,col="black",line=2.5, cex=1)
r2 <- (cor(chl_b, car_b))^2
r2_r <- round(r2, digits=2)
mtext(side=1, line=-1.0, text=paste("r²: ", r2_r, "  "), padj=-0.2, adj=1)
#fit line
f3 <- lm(chl_b ~ car_b)
cf3 <- coef(f3)
newenv <- seq(0,180, (180/500))
newenv2 <- cf3[1]+cf3[2]*newenv
lines(newenv, newenv2, lty=2, col="blue")


dev.off()


