pkgs<-c("rgdal","caret","raster","foreign")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("/home/fabi/AG_Angular_Index/Input")

# load data
# read full table
comb<-read.table("carotenoid_prepared.csv", header=TRUE, sep="\t")
comb

comb_car <- comb[,1]
comb_refl <- comb[,2:39]

# extract target parameter total chlorophyll



# define used wavelengths
wavel <- seq(410, 780, 10)
#head(birch1)

# dummy <- wavel
# dummy <- dummy - 395
# dummy
# 
# 
# # save reflectances of each wavelength in variable refl
# refl1 <- beech[,5:42]
# refl2 <- maple[,6:43]
# refl3 <- birch1[,dummy]
# refl4 <- birch2[,dummy]
# 
# refl_comb <- rbind(refl1, refl2, refl3, refl4)
# 
# 
# refl <- refl_comb*100
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
#maple_chl <- maple[,4]
#maple_refl <- maple[,6:43]

#### beech

#beech_chl <- beech[,4]
#beech_refl <- beech[,5:42]

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
path = "car_combined.txt"  
sink(path)

# define iteration variables
k <- 1
j <- 13
h <- 14

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
      
      band1 <- comb_refl[,k]*100
      band2 <- comb_refl[,j]*100
      band3 <- comb_refl[,h]*100
      

      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

      c1 <- lm(comb_car~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(comb_car~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(comb_car~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)

      
      if (c1_sum$adj.r.squared > 0.8) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > 0.8) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > 0.8) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > 0.8) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > 0.8) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > 0.8) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > 0.8) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > 0.8) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > 0.8) print(c3_sum$adj.r.squared, attach=TRUE)
      
      
        
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
h <-34

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[h]

band1 <- comb_refl[,k]*100
band2 <- comb_refl[,j]*100
band3 <- comb_refl[,h]*100


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(comb_car~AVI1+I(AVI1^2))
c1_sum <- summary(c1)

c2 <- lm(comb_car~AVI2+I(AVI2^2))
c2_sum <- summary(c2)

c3 <- lm(comb_car~AVI3+I(AVI3^2))
c3_sum <- summary(c3)

#AVI2 <- AVI2-(refl[,38]^-1)
# fit quadratic model

c2qu <- lm(comb_car~AVI3+I(AVI3^2))
sumc2qu<-summary(c2qu)
sumc2qu

plot(AVI3, comb_car)


# calculate rmse values
summary(sumc2qu$residuals)

rmse <- sqrt(sum(sumc2qu$residuals^2) / length(comb_car))
rmse



##############################################
######## create graph fitted values  #########

tiff(filename = "Car_comb_1_13_14.tif",
     width = 30, height = 20, units = "cm", pointsize = 18,
     compression = "lzw", bg = "white", res = 600)

# plot chl vs. index
plot(AVI2, comb_car)
# add curve to plot
cq = coef(c2qu)
newenv = seq(min(AVI2, na.rm=TRUE), max(AVI2, na.rm=TRUE), by = (max(AVI2, na.rm=TRUE) - min(AVI2, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')

dev.off()

