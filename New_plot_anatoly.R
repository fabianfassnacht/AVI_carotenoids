
#Set working directory and load libraries

setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs")

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

#################### playing around ##############################

# best index 1:   1   13  15
# best index 2:   5   13  14

a <- 1
b <- 13
c <- 15

#par(mfrow=c(2,2))

refl_o1 <- refl_o * 100

chap <- refl_o[,36] / refl_o[,10]
git <- (1/refl_o[,11] - 1/refl_o[,33])*refl_o[,35]
test <-range01(git)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])))
test2 <-range01(chap)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])))
AVI <- AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])
plot(car, test)



#################
#################
#### preparing NE plot merged mCRI-AVI
#################


## calculate linear model and quadratic model between merged index and car
## + calculate rmse
#test<-range01(test)
c1 <- lm(car~test+I(test^2))
c2 <- lm(car~test)
sumc1 <- summary(c1)
rmse1 <- sqrt(sum(sumc1$residuals^2) / length(car))
rmse1
sumc2 <- summary(c2)
rmse2 <- sqrt(sum(sumc2$residuals^2) / length(car))
rmse2

## plot best fit functions
#plot(car, test)
cq = coef(c1)
newenv = seq(min(test, na.rm=TRUE), max(test, na.rm=TRUE), by = (max(test, na.rm=TRUE) - min(test, na.rm=TRUE))/500)
sp.quad = cq[1]+cq[2]*newenv + cq[3]*(newenv^2) 
#lines(sp.quad,newenv, col='red')


cq2 = coef(c2)
newenv = seq(min(test, na.rm=TRUE), max(test, na.rm=TRUE), by = (max(test, na.rm=TRUE) - min(test, na.rm=TRUE))/500)
sp.lin = cq2[1] + cq2[2]*newenv
#lines(sp.lin,newenv, col='blue')

#dev.off()


#####
#####
## plot NE for merged mCRI-AVI
#####
#####

## found quadratic best fit model = 16.641 + 53.260 * mergedAVI + 6.063 * mergedAVI^2
## first derivative = 53.260 + 2 * 6.063 * mergedAVI

vls <- seq(0, max(car), ((max(car)/500)))
vls
t2 <- cq[2]-2*cq[3]*vls
t2
t3 <- rmse1/t2
plot(vls, t3, ylim = c(0, 0.1), xlab="total car", ylab="NE car", col="white")
lines(vls, t3, ylim=c(0,1), xlim=c(min(car), max(car)))

## found quadratic best fit model = 12.98 + 64.07 * mergedAVI 
## first derivative = 64.07

t4 = rmse2/(2*cq2[2])
#abline(a=t4, b= 0)

t4

#################
#################
#### preparing NE plot merged chap-AVI
#################


## calculate linear model and quadratic model between merged index and car
## + calculate rmse
test2<-range01(test2)
c1c <- lm(car~test2+I(test2^2))
c2c <- lm(car~test2)
sumc1c <- summary(c1c)
rmse1c <- sqrt(sum(sumc1c$residuals^2) / length(car))
rmse1c
sumc2c <- summary(c2c)
rmse2c <- sqrt(sum(sumc2c$residuals^2) / length(car))
rmse2c

## plot best fit functions
#plot(car, test2)
cqc = coef(c1c)
newenvc = seq(min(test2, na.rm=TRUE), max(test2, na.rm=TRUE), by = (max(test2, na.rm=TRUE) - min(test2, na.rm=TRUE))/500)
sp.quadc = cqc[1]+cqc[2]*newenv + cqc[3]*(newenvc^2) 
#lines(sp.quadc,newenvc, col='red')


cq2c = coef(c2c)
newenvc = seq(min(test2, na.rm=TRUE), max(test2, na.rm=TRUE), by = (max(test2, na.rm=TRUE) - min(test2, na.rm=TRUE))/500)
sp.linc = cq2c[1] + cq2c[2]*newenvc
#lines(sp.linc,newenvc, col='blue')

#dev.off()


#####
#####
## plot NE for merged chap-AVI
#####
#####

## found quadratic best fit model = 7.920 + 60.219 * chap-AVI + 6.711 * chap-AVI^2
## first derivative = 60.219 + 2 * 6.711 * mergedAVI

vlsc <- seq(0, max(car), ((max(car)/500)))
vlsc
t2c <- cqc[2]+(2*cqc[3]*vlsc)
t2c
t3c <- rmse1c/t2c
par(new=T)
#plot(vls, t3c, ylim = c(0, 0.2), axes=F, xlab="", ylab="", col="red")
lines(vlsc, t3c, col="red")

## found quadratic best fit model = 3.242 + 72.564 * mergedAVI 
## first derivative = 64.07

t4c = rmse2c/(2*cq2c[2])
#abline(a=t4c, b= 0, col="red")


#################
#################
#### preparing NE plot mCRI
#################


## calculate linear model and quadratic model between merged mCRI and car
## + calculate rmse
#git<-range01(git)
c1a <- lm(car~git+I(git^2))
c2a <- lm(car~git)
sumc1a <- summary(c1a)
rmse1a <- sqrt(sum(sumc1a$residuals^2) / length(car))
rmse1a
sumc2a <- summary(c2a)
rmse2a <- sqrt(sum(sumc2a$residuals^2) / length(car))
rmse2a

## plot best fit functions
#plot(car, git)
cqa = coef(c1a)
newenva = seq(min(git, na.rm=TRUE), max(git, na.rm=TRUE), by = (max(git, na.rm=TRUE) - min(git, na.rm=TRUE))/500)
sp.quada = cqa[1]+cqa[2]*newenva + cqa[3]*(newenva^2) 
#lines(sp.quada,newenva, col='red')


cq2a = coef(c2a)
sp.lina = cq2a[1] + cq2a[2]*newenva
#lines(sp.lina,newenva, col='red')

#dev.off()


#####
#####
## plot NE for mCRI index
#####
#####

## found quadratic best fit model = 19.408 + 1.828 * mCRI + 1.205 * mCRI ^ 2
## first derivative = 1.828 + 2 * 1.205 * mCRI

vlsa <- seq(0, max(car), ((max(car)/500)))
vlsa
t2a <- cqa[2]+(2 *cqa[3]*vlsa)
t2a
t3a <- rmse1a/t2a
par(new=T)
plot(vlsa, t3a, ylim = c(0, 0.2), axes=F, xlab="", ylab="", col="green")
lines(vlsa, t3a, col="green")


## found quadratic best fit model = -3.377 + 13.211 * mergedAVI 
## first derivative = 13.211

t4a = rmse2a/(2*cq2a[2])
#abline(a=t4a, b=0, col="green")



#################
#################
#### preparing NE plot AVI
#################


## calculate linear model and quadratic model between merged mCRI and car
## + calculate rmse
#AVI<-range01(AVI)
c1b <- lm(car~AVI+I(AVI^2))
c2b <- lm(car~AVI)
sumc1b <- summary(c1b)
rmse1b <- sqrt(sum(sumc1b$residuals^2) / length(car))
rmse1b
sumc2b <- summary(c2b)
rmse2b <- sqrt(sum(sumc2b$residuals^2) / length(car))
rmse2b

## plot best fit functions
#plot(car, AVI)
cqb = coef(c1b)
newenvb = seq(min(AVI, na.rm=TRUE), max(AVI, na.rm=TRUE), by = (max(AVI, na.rm=TRUE) - min(AVI, na.rm=TRUE))/500)
sp.quadb = cqb[1]+cqb[2]*newenvb + cqb[3]*(newenvb^2) 
#lines(sp.quadb,newenvb, col='red')


cq2b = coef(c2b)
sp.linb = cq2b[1] + cq2b[2]*newenvb
#lines(sp.linb,newenvb, col='red')

#dev.off()


#####
#####
## plot NE for AVI
#####
#####

## found quadratic best fit model = 34.14 + 2576.8 * AVI + -17300.24 * AVI^2
## first derivative = 2576.8 - 2 * 17300.24 * AVI

vlsb <- seq(0, max(car), ((max(car)/500)))
vlsb
t2b <- cqb[2] +(2*cqb[3]*vlsb)
t2b
t3b <- rmse1b/t2b
par(new=T)
#plot(vlsb, t3b, ylim = c(0, 0.2), axes=F, xlab="", ylab="", col="darkgreen")
lines(vlsb, t3b, col="blue")

## found quadratic best fit model = -3.377 + 13.211 * mergedAVI 
## first derivative = 13.211

t4b = rmse2b/(2*cq2b[2])
#abline(a=t4b, b=0, col="darkgreen")




#################
#################
#### preparing NE plot chap
#################


## calculate linear model and quadratic model between merged mCRI and car
## + calculate rmse
chap<-range01(chap)
c1d <- lm(car~chap+I(chap^2))
c2d <- lm(car~chap)
sumc1d <- summary(c1d)
rmse1d <- sqrt(sum(sumc1d$residuals^2) / length(car))
rmse1d
sumc2d <- summary(c2d)
rmse2d <- sqrt(sum(sumc2d$residuals^2) / length(car))
rmse2d

## plot best fit functions
#plot(car, chap)
cqd = coef(c1d)
newenvd = seq(min(chap, na.rm=TRUE), max(chap, na.rm=TRUE), by = (max(chap, na.rm=TRUE) - min(chap, na.rm=TRUE))/500)
sp.quadd = cqd[1]+cqd[2]*newenvd + cqd[3]*(newenvd^2) 
#lines(sp.quadb,newenvd, col='red')


cq2d = coef(c2d)
sp.lind = cq2d[1] + cq2d[2]*newenvd
#lines(sp.lind,newenvd, col='red')

#dev.off()


#####
#####
## plot NE for chap
#####
#####

## found quadratic best fit model = 34.14 + 2576.8 * AVI + -17300.24 * AVI^2
## first derivative = 2576.8 - 2 * 17300.24 * AVI
## car = 2576.8 - 2 * 17300.24 * AVI



vlsd <- seq(0, max(car), ((max(car)/89)))
vlsd
t2d <- cqd[2]+(2*cqd[3]*vlsd)
t2d
t3d <- rmse1d/t2d
par(new=T)

#dev.off()
plot(vlsd, t3d, ylim = c(0, 0.2), col="white")
lines(vlsd, t3d, col="grey80")

## found quadratic best fit model = -3.377 + 13.211 * mergedAVI 
## first derivative = 13.211
t4d = rmse2d/(2*cq2d[2])
#abline(a=t4d, b=0, col="grey80")





##########################################

# scale function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

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
