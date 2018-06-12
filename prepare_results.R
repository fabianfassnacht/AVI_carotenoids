pkgs<-c("rgdal","caret","raster","foreign", "randomForest")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs")

# load data
# read full table
spectra <- read.table("all_data_prepared_tp.txt", sep="\t", head=T)

head(spectra)
spectra <- spectra[1:106,]


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

car_m <- car[1:47]
chl_m <- chl[1:47]
refl_m <- refl_o[1:47,]
refl_m <- refl_m * 100

#### chestnut

car_c <- car[48:69]
chl_c <- chl[48:69]
refl_c <- refl_o[48:69,]
refl_c <- refl_c*100

#### beech

car_b <- car[70:106]
chl_b <- chl[70:106]
refl_b <- refl_o[70:106,]
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

#found indices
plot(chl, AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c]))
plot(car, AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c]))

# indices gitelson paper
git <- 1/refl_o[,11] - 1/refl_o[,15]
git2 <- 100/refl_o[,11] - 100/refl_o[,30]
chap <- refl_o[,36] / refl_o[,10]
gamon <- (refl_o[,13] - refl_o[,17]) / (refl_o[,13] + refl_o[,17])
datt <- 0.0049*(refl_o[,27]/(refl_o[,15]*refl_o[,31])^0.748)

plot(car, git)
plot(car, git2)
plot(car, chap)
plot(car, gamon)
plot(car, datt)

# merged index 

plot(car, (chap+(AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])*120)))
plot(car, range01(chap)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c]))))

# prepare results

t <-AVI3(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])

# calculate lm between index and carotenoids

chapavi <- range01(chap)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])))
gpqu <- lm(car~chapavi)
sumgpqu<-summary(gpqu)
rmse <- sqrt(sum(sumgpqu$residuals^2) / 106)
rmse

plot(car, chapavi)

##### species dependent results ######

chap_m <- refl_m[,36] / refl_m[,10]
chap_c <- refl_c[,36] / refl_c[,10]
chap_b <- refl_b[,36] / refl_b[,10]

plot(car_m, (chap_m+(AVI2(wavel[a],wavel[b],wavel[c],refl_m[,a],refl_m[,b],refl_m[,c])*120)))
plot(car_c, (chap_c+(AVI2(wavel[a],wavel[b],wavel[c],refl_c[,a],refl_c[,b],refl_c[,c])*120)))
plot(car_b, (chap_b+(AVI2(wavel[a],wavel[b],wavel[c],refl_b[,a],refl_b[,b],refl_b[,c])*120)))

cor(car_m, (chap_m+(AVI2(wavel[a],wavel[b],wavel[c],refl_m[,a],refl_m[,b],refl_m[,c])*120)))
cor(car_c, (chap_c+(AVI2(wavel[a],wavel[b],wavel[c],refl_c[,a],refl_c[,b],refl_c[,c])*120)))
cor(car_b, (chap_b+(AVI2(wavel[a],wavel[b],wavel[c],refl_b[,a],refl_b[,b],refl_b[,c])*120)))

plot(car_m, (range01(chap_m)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_m[,a],refl_m[,b],refl_m[,c]))))
plot(car_c, (range01(chap_c)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_c[,a],refl_c[,b],refl_c[,c])))))
plot(car_b, (range01(chap_b)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_b[,a],refl_b[,b],refl_b[,c])))))
     
cor(car_m, (range01(chap_m)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_m[,a],refl_m[,b],refl_m[,c])))))
cor(car_c, (range01(chap_c)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_c[,a],refl_c[,b],refl_c[,c])))))
cor(car_b, (range01(chap_b)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_b[,a],refl_b[,b],refl_b[,c])))))
     
     
     ############# make range of carotenoid values smaller     
     
     par(mfrow=c(3,1))
     
     t <-  car < 140 & car > 70
     
     refl_o2 <- refl_o1[t,]
     chl_o2 <- chl[t]
     car_o2 <- car[t]
     
     chap <- refl_o2[,36] / refl_o2[,10]
     
     plot(chl_o2, AVI2(wavel[a],wavel[b],wavel[c],refl_o2[,a],refl_o2[,b],refl_o2[,c]))
     plot(car_o2, AVI2(wavel[a],wavel[b],wavel[c],refl_o2[,a],refl_o2[,b],refl_o2[,c]))
     plot(car_o2, chl_o2)
     plot(car_o2, (chap+(AVI2(wavel[a],wavel[b],wavel[c],refl_o2[,a],refl_o2[,b],refl_o2[,c])*100)))
     




     ##############################################
     ######## create graph spectral curve #########
     
     tiff(filename = "Found_indices_graphically.tif",
          width = 30, height = 20, units = "cm", pointsize = 18,
          compression = "lzw", bg = "white", res = 600)
     
     # plot of sample with chl=126.9
     test <- data.matrix(refl_o[75,])
     test <- t(test)
     
     test
     test2
     # plot of sample with chl=449.95
     test2 <- data.matrix(refl_o[105,])
     test2 <- t(test2)
     
     test3<-rbind(test, test2)
     plot (test3, xlab = "wavelength [no unit]", ylab="reflectance [10000 = 100%]")
     
     
     #abline(a=0,b=0, v=4)
     #abline(a=0,b=0, v=31)
     #abline(a=0,b=0, v=32)
     mtext("car = 43.51", at = c(12,5000))
     
     #abline(a=0,b=0, v=42)
     #abline(a=0,b=0, v=69)
     #abline(a=0,b=0, v=70)
     mtext("car = 137.18", at = c(51, 5000))
     
     abline(a=0,b=0, v=1, col = 4)
     abline(a=0,b=0, v=13, col = 4)
     abline(a=0,b=0, v=15, col = 4)
     abline(a=0,b=0, v=36, col = 5)
     abline(a=0,b=0, v=10, col = 5)
     
     abline(a=0,b=0, v=39, col = 4)
     abline(a=0,b=0, v=51, col = 4)
     abline(a=0,b=0, v=53, col = 4)
     abline(a=0,b=0, v=48, col = 5)
     abline(a=0,b=0, v=74, col = 5)
     
     str(test2)
     
     dev.off()
     
     
     ###



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
     
     # scale function
     range01 <- function(x){(x-min(x))/(max(x)-min(x))}