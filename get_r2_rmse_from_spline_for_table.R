pkgs<-c("rgdal","caret","raster","foreign", "randomForest")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("D:/KIT_Forschung/11_AG_Angular_Index/99_Paper/5_in_data")

# load data
# read full table
spectra <- read.table("corrected_data_ag.csv", sep="\t", head=T)

head(spectra)
spectra <- spectra[1:90,]


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

########## plot + retrieve rmse and r²


######### chlorophyll

rmser2(chl, AVI, indexname="AVI")
rmser2(chl_m, AVI_m, indexname="AVI")
rmser2(chl_b, AVI_b, indexname="AVI")
rmser2(chl_c, AVI_c, indexname="AVI")

rmser2(chl, git, indexname="mCRI")
rmser2(chl_m, git_m, indexname="mCRI")
rmser2(chl_b, git_b, indexname="mCRI")
rmser2(chl_c, git_c, indexname="mCRI")

rmser2(chl, chap, indexname="CI")
rmser2(chl_m, chap_m, indexname="CI")
rmser2(chl_b, chap_b, indexname="CI")
rmser2(chl_c, chap_c, indexname="CI")

rmser2(chl, chap_AVI, indexname="merged AVI-CI")
rmser2(chl_m, chap_AVI_m, indexname="merged AVI-CI")
rmser2(chl_b, chap_AVI_b, indexname="merged AVI-CI")
rmser2(chl_c, chap_AVI_c, indexname="merged AVI-CI")

rmser2(chl, git_AVI, indexname="merged AVI-mCRI")
rmser2(chl_m, git_AVI_m, indexname="merged AVI-mCRI")
rmser2(chl_b, git_AVI_b, indexname="merged AVI-mCRI")
rmser2(chl_c, git_AVI_c, indexname="merged AVI-mCRI")


######### carotenoid

rmser2(car, AVI, indexname="AVI")
rmser2(car_m, AVI_m, indexname="AVI")
rmser2(car_b, AVI_b, indexname="AVI")
rmser2(car_c, AVI_c, indexname="AVI")

rmser2(car, git, indexname="mCRI")
rmser2(car_m, git_m, indexname="mCRI")
rmser2(car_b, git_b, indexname="mCRI")
rmser2(car_c, git_c, indexname="mCRI")

rmser2(car, chap, indexname="CI")
rmser2(car_m, chap_m, indexname="CI")
rmser2(car_b, chap_b, indexname="CI")
rmser2(car_c, chap_c, indexname="CI")

rmser2(car, chap_AVI, indexname="merged AVI-CI")
rmser2(car_m, chap_AVI_m, indexname="merged AVI-CI")
rmser2(car_b, chap_AVI_b, indexname="merged AVI-CI")
rmser2(car_c, chap_AVI_c, indexname="merged AVI-CI")

rmser2(car, git_AVI, indexname="merged AVI-mCRI")
rmser2(car_m, git_AVI_m, indexname="merged AVI-mCRI")
rmser2(car_b, git_AVI_b, indexname="merged AVI-mCRI")
rmser2(car_c, git_AVI_c, indexname="merged AVI-mCRI")

#### plot and retrieve rmse + r²

rmser2 <- function(car, index, indexname){

plot(car, index, ylab=indexname, xlab="tot. car. content [mg/m²]")
ind.spl <- smooth.spline(car, index, cv=T, nknots=5)
lines(ind.spl, col="blue")

ind.spl2 <- smooth.spline(index,car, cv=T, nknots=5)
rmse <- sqrt(sum((car-fitted(ind.spl2))^2) / length(car))
rmse_r <- round(rmse, digits=2)
r2 <- (cor(fitted(ind.spl2), car))^2
r2_r <- round(r2, digits=2)
mtext(side=1, line=-1.0, text=paste("RMSE: ", rmse_r, "  r²: ", r2_r, "  "), padj=-0.2, adj=1)

}


############################################
#### needed basic functions
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

