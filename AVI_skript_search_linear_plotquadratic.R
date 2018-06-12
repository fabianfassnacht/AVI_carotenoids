pkgs<-c("rgdal","caret","raster","foreign")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("/home/fabi/AG_Angular_Index/Input")

# load data
# read full table
beech<-read.table("Beech-reflectance.txt", header=TRUE, sep="\t")
maple<-read.table("Maple-reflectance.txt", header=TRUE, sep="\t")
birch1<-read.table("Birch_1.txt", header=TRUE, sep="\t")
birch2<-read.table("Birch_2.txt", header=TRUE, sep="\t")

head(beech)
head(maple)
head(birch1)

# extract target parameter total chlorophyll
chl <- beech[,1:3]
chl2 <- maple[,1:3]
chl3 <- birch1[,2:4]
chl4 <- birch2[,2:4]
str(chl)
str(chl4)

chl_comb <- rbind(chl, chl2, chl3, chl4)
chl <- chl_comb[,3]
str(chl)
# define used wavelengths
wavel <- seq(410, 780, 10)
head(birch1)

dummy <- wavel
dummy <- dummy - 395
dummy

refl3

# save reflectances of each wavelength in variable refl
refl <- beech[,5:42]
refl2 <- maple[,6:43]
refl3 <- birch1[,dummy]
refl4 <- birch2[,dummy]

refl_comb <- rbind(refl, refl2, refl3, refl4)

refl <- refl_comb*100
# create array with wavelength
str(refl)
wavel <- seq(410, 780, 10)


# open textfile to store results

path = "best_subset_regression_all.txt"  
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
      
      band1 <- refl[,k]
      band2 <- refl[,j]
      band3 <- refl[,h]
      

      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

      c1 <- lm(chl~AVI1)
      c1_sum <- summary(c1)
      
      c2 <- lm(chl~AVI2)
      c2_sum <- summary(c2)
      
      c3 <- lm(chl~AVI3)
      c3_sum <- summary(c3)

      
      if (c1_sum$adj.r.squared > 0.921) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > 0.921) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > 0.921) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > 0.921) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > 0.921) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > 0.921) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > 0.921) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > 0.921) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > 0.921) print(c3_sum$adj.r.squared, attach=TRUE)
        
    }
  }
}


sink()


#################################################
#################################################


## test best found results

wave1 <- wavel[1]
wave2 <- wavel[31]
wave3 <- wavel[33]

band1 <- refl[,1]
band2 <- refl[,31]
band3 <- refl[,33]


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(chl~AVI1)
c1_sum <- summary(c1)

c2 <- lm(chl~AVI2)
c2_sum <- summary(c2)

c3 <- lm(chl~AVI3)
c3_sum <- summary(c3)

c1_sum
c2_sum
c3_sum


# fit quadratic model
c2qu <- lm(chl~AVI2+I(AVI2^2))
sumc2qu<-summary(c2qu)
sumc2qu

# plot chl vs. index
plot(AVI2, chl)
# add curve to plot
cq = coef(c2qu)
newenv = seq(min(AVI2, na.rm=TRUE), max(AVI2, na.rm=TRUE), by = (max(AVI2, na.rm=TRUE) - min(AVI2, na.rm=TRUE))/500)
sp.quad = cq[1] + cq[2]*newenv +cq[3]*newenv^2
lines(newenv,sp.quad, col='red')


 