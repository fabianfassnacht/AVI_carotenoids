pkgs<-c("rgdal","caret","raster","foreign", "randomForest")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("/home/fabi/KIT/11_AG_Angular_Index/99_Paper/5_in_data")

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



## all data

refl_car <- refl_o*100

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/99_Paper/1_AVI_search/NDVI/systematic_NDVI")
path = "NDVI_lin_r2.txt"  
sink(path)

# define iteration variables
k <- 1
j <- 2
h <- 3

# start loop
for(k in 1:36) { 
  
  j1 = k + 1
  
  for (j in j1:37) {
    
    h1 = j + 1 
    for (h in h1:38) {  

      #calculate AVIs
      
      wave1 <- wavel[k]
      wave2 <- wavel[j]
      wave3 <- wavel[h]
      
      band1 <- refl_car[,k]
      band2 <- refl_car[,j]
      band3 <- refl_car[,h]
      
#       AVI1 <- band1/band2
#       AVI2 <- band1/band3
#       AVI3 <- band2/band3
      
#       v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
#       v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
#       v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
# 
# 
#       AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
#       AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
#       AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      AVI1 <- (band1-band2)/(band1+band2)
      AVI2 <- (band2-band1)/(band2+band1)
      
      AVI1[is.na(AVI1)] <- 0
      AVI2[is.na(AVI2)] <- 0
      #AVI3[is.na(AVI3)] <- 0
      
      #AVI1 <- range01(git)+range01(AVI1)
      #AVI2 <- range01(git)+range01(AVI2)
      #AVI3 <- range01(git)+range01(AVI3)
      
      c1 <- lm(car~AVI1)
      c1_sum <- summary(c1)
      rmse1 <- sqrt(sum(c1_sum$residuals^2) / length(car))
      
      c2 <- lm(car~AVI2)
      c2_sum <- summary(c2)
      rmse2 <- sqrt(sum(c2_sum$residuals^2) / length(car))
      
      #c3 <- lm(car~AVI3)
      #c3_sum <- summary(c3)
      #rmse3 <- sqrt(sum(c3_sum$residuals^2) / length(car))

      th <- 0.70
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      #if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      #if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      #if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
#       
#       if (rmse1 < th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
#       if (rmse1 < th) print(AVI_index, attach=TRUE)
#       if (rmse1 < th) print(rmse1, attach=TRUE)
#       
#       if (rmse2 < th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
#       if (rmse2 < th) print(AVI_index, attach=TRUE)
#       if (rmse2 < th) print(rmse2, attach=TRUE)
#       
#       if (rmse3 < th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
#       if (rmse3 < th) print(AVI_index, attach=TRUE)
#       if (rmse3 < th) print(rmse3, attach=TRUE)

    }
  }
}


sink()


#################################################
#################################################


## test best found results

k <- 35
j <- 36
i <- 37

###
### on beech data

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[i]

band1 <- refl_b[,k]
band2 <- refl_b[,j]
band3 <- refl_b[,i]


# v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
# v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
# v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
# 
# 
# AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
# AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
# AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

AVI1 <- (band1-band2)/(band1+band2)
AVI2 <- (band2-band1)/(band2+band1)

c1 <- lm(car_b~AVI1)
c1_sum <- summary(c1)
rmse1 <- sqrt(sum(c1_sum$residuals^2) / length(car))

c2 <- lm(car_b~AVI2)
c2_sum <- summary(c2)
rmse2 <- sqrt(sum(c2_sum$residuals^2) / length(car))
rmse2
c2_sum$adj.r.squared

c3 <- lm(car_b~AVI3)
c3_sum <- summary(c3)
rmse3 <- sqrt(sum(c3_sum$residuals^2) / length(car))

c1_sum
c2_sum
c3_sum

###
### on maple data

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[i]

band1 <- refl_m[,k]
band2 <- refl_m[,j]
band3 <- refl_m[,i]


# v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
# v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
# v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
# 
# 
# AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
# AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
# AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

AVI1 <- (band1-band2)/(band1+band2)
AVI2 <- (band2-band1)/(band2+band1)

c1 <- lm(car_m~AVI1)
c1_sum <- summary(c1)
rmse1 <- sqrt(sum(c1_sum$residuals^2) / length(car))

c2 <- lm(car_m~AVI2)
c2_sum <- summary(c2)
rmse2 <- sqrt(sum(c2_sum$residuals^2) / length(car))
rmse2
c1_sum$adj.r.squared

c3 <- lm(car_m~AVI3)
c3_sum <- summary(c3)
rmse3 <- sqrt(sum(c3_sum$residuals^2) / length(car))

c1_sum
c2_sum
c3_sum

###
### on chestnut data

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[i]

band1 <- refl_c[,k]
band2 <- refl_c[,j]
band3 <- refl_c[,i]


# v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
# v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
# v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
# 
# 
# AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
# AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
# AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

AVI1 <- (band1-band2)/(band1+band2)
AVI2 <- (band2-band1)/(band2+band1)

c1 <- lm(car_c~AVI1)
c1_sum <- summary(c1)
rmse1 <- sqrt(sum(c1_sum$residuals^2) / length(car))

c2 <- lm(car_c~AVI2)
c2_sum <- summary(c2)
rmse2 <- sqrt(sum(c2_sum$residuals^2) / length(car))
rmse2
c1_sum$adj.r.squared

c3 <- lm(car_c~AVI3)
c3_sum <- summary(c3)
rmse3 <- sqrt(sum(c3_sum$residuals^2) / length(car))

c1_sum
c2_sum
c3_sum




#################
########### AVI functions
#################


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


