pkgs<-c("rgdal","caret","raster","foreign", "randomForest", "mgcv", "stats", "pspline")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

#setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs")
setwd("/media/fabi/KINGSTON/11_AG_Angular_Index/1_new_runs")

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



## all data

refl_car <- refl_o*100

#################### playing around ##############################

# best index 1:   1   13  15
# best index 2:   5   13  14

a <- 1
b <- 13
c <- 15

par(mfrow=c(2,2))
refl_o1 <- refl_o * 100


####
#### Graph 1: AVI
####

# calculate Indices for all data and for species data separately
AVI <- AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])
AVI_m <- AVI2(wavel[a],wavel[b],wavel[c],refl_m[,a],refl_m[,b],refl_m[,c])
AVI_c <- AVI2(wavel[a],wavel[b],wavel[c],refl_c[,a],refl_c[,b],refl_c[,c])
AVI_b <- AVI2(wavel[a],wavel[b],wavel[c],refl_b[,a],refl_b[,b],refl_b[,c])

# indices gitelson paper
git <- (1/refl_o[,11] - 1/refl_o[,33])*refl_o[,35]
git_m <-(1/refl_m[,11] - 1/refl_m[,33])*refl_m[,35]
git_c <-(1/refl_c[,11] - 1/refl_c[,33])*refl_c[,35]
git_b <-(1/refl_b[,11] - 1/refl_b[,33])*refl_b[,35]

# chappele index
chap <- refl_o[,36] / refl_o[,10]
chap_m <- refl_m[,36] / refl_m[,10]
chap_c <- refl_c[,36] / refl_c[,10]
chap_b <- refl_b[,36] / refl_b[,10]

# merged index mCRI + AVI

git_AVI <- range01(git) + range01(AVI)
git_AVI_m <- range01(git_m) + range01(AVI_m)
git_AVI_c <- range01(git_c) + range01(AVI_c)
git_AVI_b <- range01(git_b) + range01(AVI_b)

# merged index Chappele + AVI

chap_AVI <- range01(chap) + range01(AVI)
chap_AVI_m <- range01(chap_m) + range01(AVI_m)
chap_AVI_c <- range01(chap_c) + range01(AVI_c)
chap_AVI_b <- range01(chap_b) + range01(AVI_b)

#

#################
########### AVI functions
#################
dev.off()

############################################

AVI1 <- function(wave1, wave2, wave3, band1, band2, band3) {
  

  v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
  v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
  v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
  
  
  AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
  return(AVI1)
  
}

AVI2 <- function(wave1, wave2, wave3, band1, band2, band3) {
  
  
  v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
  v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
  v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
  
  
  AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
  return(AVI2)
  
}

AVI3 <- function(wave1, wave2, wave3, band1, band2, band3) {
  
  
  v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
  v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
  v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
  
  
  AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
  return(AVI3)
  
}

# scale function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# calculate lm and rmse
rmse_lm <- function(resp, pred) {
  
  gpqu <- lm(resp~pred)
  sumgpqu<-summary(gpqu)
  rmse <- sqrt(sum(sumgpqu$residuals^2) / length(pred))
  return(rmse)
}

#################################
####### function for ploting:

min_car <- refl_b[1,1:19]
med_car <- refl_b[18,1:19]
max_car <- refl_b[33,1:19]
minr <- min(min_car, med_car, max_car)
maxr <- max(min_car, med_car, max_car)

maxr <- maxr + 5


setwd("/media/fabi/KINGSTON/11_AG_Angular_Index/99_Paper/12_plot_found_index")


#pdf(file="found_AVI_old.pdf", width=3, height=19)
tiff(file="found_AVI.tif", width=3, height=19, units="in", compression="lzw", res=600)

plot(wavel[1:19], min_car, ylim=c(500,maxr), xlim=c(400,600), col="white", xlab="Wavelength [nm]", ylab = "Reflectance [10000 = 100%]", cex=1.7)
par(new=T)
lines(wavel[1:19], min_car, ylim=c(500,maxr))
points(wavel[1:19], min_car, ylim=c(500,maxr), pch=15, cex=0.7)
lines(wavel[1:19], med_car, add=T, ann=F, axes=F, ylim=c(0,maxr), col="darkgreen")
points(wavel[1:19], med_car, ylim=c(500,maxr), pch="+", cex=0.7)
lines(wavel[1:19], max_car, add=T, ann=F, axes=F, ylim=c(0,maxr), col="darkorange")
points(wavel[1:19], max_car, ylim=c(500,maxr), pch=0, cex=0.7, col="grey30")


# plot triangles
# for car_min
par(new=T)
bnds <- c(1,13,15)
bnds2 <- c(1,15)
yval <- min_car[bnds]
wvl <- wavel[bnds]
yval2 <- min_car[bnds2]
wvl2 <- wavel[bnds2]
lines(wvl, yval, col="red")
lines(wvl2,yval2, col="red")


# for car_med
yval <- med_car[bnds]
wvl <- wavel[bnds]
yval2 <- med_car[bnds2]
wvl2 <- wavel[bnds2]
lines(wvl, yval, col="red")
lines(wvl2,yval2, col="red")


# for car_max
yval <- max_car[bnds]
wvl <- wavel[bnds]
yval2 <- max_car[bnds2]
wvl2 <- wavel[bnds2]
lines(wvl, yval, col="red")
lines(wvl2,yval2, col="red")


dev.off()

