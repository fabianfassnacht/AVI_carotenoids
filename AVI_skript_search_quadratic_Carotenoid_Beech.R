pkgs<-c("rgdal","caret","raster","foreign")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("/home/fabi/AG_Angular_Index/Input")

# load data
# read full table
beech<-read.table("Beech-reflectance.txt", header=TRUE, sep="\t")
#maple<-read.table("Maple-reflectance.txt", header=TRUE, sep="\t")
#birch1<-read.table("Birch_1.txt", header=TRUE, sep="\t")
#birch2<-read.table("Birch_2.txt", header=TRUE, sep="\t")

head(beech)
#head(maple)
#head(birch1)

# extract target parameter total chlorophyll



# define used wavelengths
wavel <- seq(410, 780, 10)
#head(birch1)

#dummy <- wavel
#dummy <- dummy - 395
#dummy


# save reflectances of each wavelength in variable refl
refl1 <- beech[,5:42]
refl2 <- maple[,6:43]
refl3 <- birch1[,dummy]
refl4 <- birch2[,dummy]

refl_comb <- rbind(refl1, refl2, refl3, refl4)


refl <- refl_comb*100
# create array with wavelength
#str(refl)
#wavel <- seq(410, 780, 10)


###############
###############
### Set dataset
###############
###############

# #### maple
# 
# chl <- chl[1:
# refl <- refl[1:

#### beech

beech_chl <- beech[,4]
beech_refl <- beech[,5:42]

# #### birch1
# 
# chl <- chl[100:149]
# refl <- refl[100:149,]
# 
# ### birch2
# 
# chl <- chl[150:197]
# refl <- refl[150:197,]


# open textfile to store results
setwd("/home/fabi/AG_Angular_Index/Out_Car")
path = "car_bbeech.txt"  
sink(path)

# define iteration variables
k <- 1
j <- 13
h <- 15

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
      
      band1 <- beech_refl[,k]*100
      band2 <- beech_refl[,j]*100
      band3 <- beech_refl[,h]*100
      

      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

      c1 <- lm(beech_chl~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(beech_chl~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(beech_chl~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)

      
      if (c1_sum$adj.r.squared > 0.89) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > 0.89) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > 0.89) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > 0.89) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > 0.89) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > 0.89) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > 0.89) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > 0.89) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > 0.89) print(c3_sum$adj.r.squared, attach=TRUE)
      
      
        
    }
  }
}


sink()


#################################################
#################################################


## test best found results

###
### 1


k <- 16
j <- 17
h <- 34

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[h]

band1 <- beech_refl[,k]*100
band2 <- beech_refl[,j]*100
band3 <- beech_refl[,h]*100


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(beech_chl~AVI1+I(AVI1^2))
c1_sum <- summary(c1)

c2 <- lm(beech_chl~AVI2+I(AVI2^2))
c2_sum <- summary(c2)

c3 <- lm(beech_chl~AVI3+I(AVI3^2))
c3_sum <- summary(c3)

#AVI2 <- AVI2-(refl[,38]^-1)
# fit quadratic model
c2qu <- lm(beech_chl~AVI3+I(AVI3^2))
sumc2qu<-summary(c2qu)
sumc2qu

plot(AVI3, beech_chl)


# calculate rmse values
summary(sumc2qu$residuals)

rmse <- sqrt(sum(sumc2qu$residuals^2) / length(beech_chl))
rmse



##############################################
######## create graph fitted values  #########

tiff(filename = "Beech_16_17_34.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot chl vs. index
plot(AVI2, beech_chl)
# add curve to plot
cq = coef(c2qu)
newenv = seq(min(AVI2, na.rm=TRUE), max(AVI2, na.rm=TRUE), by = (max(AVI2, na.rm=TRUE) - min(AVI2, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')

dev.off()

##############################################
######## create graph fitted values  #########

tiff(filename = "Birch2_9_10_14.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot chl vs. index
plot(AVI3, chl)
# add curve to plot
cq = coef(c2qu)
newenv = seq(min(AVI3, na.rm=TRUE), max(AVI3, na.rm=TRUE), by = (max(AVI3, na.rm=TRUE) - min(AVI3, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')

dev.off()
##############################################
######## create graph spectral curve #########

tiff(filename = "Found_index_graphically.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot of sample with chl=126.9
test <- data.matrix(refl[5,])
test <- t(test)

test
test2
# plot of sample with chl=449.95
test2 <- data.matrix(refl[25,])
test2 <- t(test2)

test3<-rbind(test, test2)
plot (test3, xlab = "wavelength [no unit]", ylab="reflectance [10000 = 100%]")


abline(a=0,b=0, v=1)
abline(a=0,b=0, v=31)
abline(a=0,b=0, v=33)
mtext("chl-tot = 126.9", at = c(12,5000))

abline(a=0,b=0, v=39)
abline(a=0,b=0, v=69)
abline(a=0,b=0, v=71)
mtext("chl-tot = 449.95", at = c(51, 5000))

str(test2)

dev.off()

###################################
###################################
#### add additional information ###
#green peak
gp <- AVI2 + (40/refl[,13])
#NIR
gp <- AVI2 + (1/refl[,38])
#both
gp <- AVI2 + (40/refl[,13]) + (1/refl[,38])

plot(gp,chl)

mgp <- lm(chl~gp)
summgp <- summary(mgp)
summgp

gpqu <- lm(chl~gp+I(gp^2))
sumgpqu<-summary(gpqu)
sumgpqu

# add curve to plot
cq = coef(sumgpqu)
newenv = seq(min(gp, na.rm=TRUE), max(gp, na.rm=TRUE), by = (max(gp, na.rm=TRUE) - min(gp, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')

rmse <- sqrt(sum(sumgpqu$residuals^2) / 197)
rmse
##############################################
#################################################

## test best found results

###
### 2


wave1 <- wavel[4]
wave2 <- wavel[31]
wave3 <- wavel[32]

band1 <- refl[,4]
band2 <- refl[,31]
band3 <- refl[,32]


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1a <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2a <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3a <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(chl~AVI1)
c1_sum <- summary(c1)

c2 <- lm(chl~AVI2)
c2_sum <- summary(c2)

c3 <- lm(chl~AVI3)
c3_sum <- summary(c3)

c1_sum
c2_sum
c3_sum

#AVI3 <- AVI3-(refl[,38]^-1)
# fit quadratic model
c3qu <- lm(chl~AVI3+I(AVI3^2))
sumc2qu<-summary(c2qu)
sumc2qu


# fit quadratic model
c3qu <- lm(chl~AVI3+I(AVI3^2))
sumc3qu<-summary(c3qu)
sumc3qu

# calculate rmse values
summary(sumc3qu$residuals)

rmse <- sqrt(sum(sumc3qu$residuals^2) / nrs)
rmse



##############################################
######## create graph fitted values  #########

tiff(filename = "fit_index_graphically.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot chl vs. index
plot(AVI3, chl)
# add curve to plot
cq = coef(c3qu)
newenv = seq(min(AVI3, na.rm=TRUE), max(AVI3, na.rm=TRUE), by = (max(AVI3, na.rm=TRUE) - min(AVI3, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')

dev.off()


##############################################
######## create graph spectral curve #########

tiff(filename = "Found_indices_graphically.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot of sample with chl=126.9
test <- data.matrix(refl[5,])
test <- t(test)

test
test2
# plot of sample with chl=449.95
test2 <- data.matrix(refl[25,])
test2 <- t(test2)

test3<-rbind(test, test2)
plot (test3, xlab = "wavelength [no unit]", ylab="reflectance [10000 = 100%]")


abline(a=0,b=0, v=4)
abline(a=0,b=0, v=31)
abline(a=0,b=0, v=32)
mtext("chl-tot = 126.9", at = c(12,5000))

abline(a=0,b=0, v=42)
abline(a=0,b=0, v=69)
abline(a=0,b=0, v=70)
mtext("chl-tot = 449.95", at = c(51, 5000))

abline(a=0,b=0, v=9, col = 4)
abline(a=0,b=0, v=10, col = 4)
abline(a=0,b=0, v=14, col = 4)

abline(a=0,b=0, v=47, col = 4)
abline(a=0,b=0, v=48, col = 4)
abline(a=0,b=0, v=52, col = 4)

str(test2)

dev.off()


###



##############################################
######## calculate earlier presented VIs #####

###### (R800-R680) / (R800+R680)  ## 780 instead of 800

wavel[28]

band1 <- refl[,38] 
band2 <- refl[,28]

vi <- (band1-band2)/(band1+band2)

mvi <- lm(chl~vi)
summary(mvi)

mviqu <- lm(chl~vi+I(vi^2))
summviqu<-summary(mviqu)
summviqu

plot(vi, chl)

# calculate rmse values
rmse <- sqrt(sum(mviqu$residuals^2) / 197)
rmse

###### (R800-R700) / (R800+R700) ## 780 instead of 800

wavel[28]

band1 <- refl[,38] 
band2 <- refl[,30]

vi1 <- (band1-band2)/(band1+band2)
#vi1 <- vi1-(refl[,38]^-1)

mvi1 <- lm(chl~vi1)
summary(mvi1)

plot(vi1, chl)

mvi1qu <- lm(chl~vi1+I(vi1^2))
summvi1qu<-summary(mvi1qu)
summvi1qu

# calculate rmse values
rmse <- sqrt(sum(mvi1qu$residuals^2) / 48)
rmse


###### R800 / R680   ## 780 instead of 800

wavel[28]

band1 <- refl[,38] 
band2 <- refl[,28]

vi2 <- band1/band2

mvi2 <- lm(chl~vi2)
summary(mvi2)

mvi2qu <- lm(chl~vi2+I(vi2^2))
summvi2qu<-summary(mvi2qu)
summvi2qu

plot(vi2, chl)

# calculate rmse values
rmse <- sqrt(sum(mviqu2$residuals^2) / 48)
rmse

###### (R750/R720) ^-1

band1 <- refl[,38] 
band2 <- refl[,32]

vi3 <- (band1/band2)^(-1)

mvi3 <- lm(chl~vi3)
summary(mvi3)

mvi3qu <- lm(chl~vi3+I(vi3^2))
summvi3qu<-summary(mvi3qu)
summvi3qu

plot(vi3,chl)

# calculate rmse values
rmse <- sqrt(sum(mvi3qu$residuals^2) / 48)
rmse
 