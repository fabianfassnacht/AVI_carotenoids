
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

AVI <- AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])
chap <- refl_o[,36] / refl_o[,10]
git <- (1/refl_o[,11] - 1/refl_o[,33])*refl_o[,35]
AVI_mCRI <-range01(git)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])))
AVI_chap <- range01(chap) + range01(AVI)

##########
### make plots
##########

pdf(file="EN_plot_quad_fun.pdf", width=7, height=6)

pigpl_EN_qu(car, git, add=F, col="black", lty=1)
pigpl_EN_qu(car, AVI, add=T, col="black", lty=2, lwd=2)
pigpl_EN_qu(car, chap, add=T, col="black", lty=3)
pigpl_EN_qu(car, AVI_mCRI, add=T, col="red", lty=1)
pigpl_EN_qu(car, AVI_chap, add=T, col="red", lty=3, lwd=2)
legend(x="topright", legend=c("mCRI", "AVI", "CI", "AVI_mCRI", "AVI_CI"), lwd=c(1,2,1,1,2), lty=c(1,2,3,1,3), col=c("black","black", "black","red", "red"), title="", cex=0.7, bty="n")

dev.off()


pdf(file="EN_plot_lin_fun.pdf", width=7, height=6)

pigpl_EN_lin(car, git, add=F, col="black", lty=1)
pigpl_EN_lin(car, AVI, add=T, col="black", lty=2, lwd=2)
pigpl_EN_lin(car, chap, add=T, col="black", lty=3)
pigpl_EN_lin(car, AVI_mCRI, add=T, col="red", lty=1)
pigpl_EN_lin(car, AVI_chap, add=T, col="red", lty=3, lwd=2)
legend(x="topright", legend=c("mCRI", "AVI", "CI", "AVI_mCRI", "AVI_CI"), lwd=c(1,2,1,1,2), lty=c(1,2,3,1,3), col=c("black","black", "black","red", "red"), title="", cex=0.7, bty="n")

dev.off()


pdf(file="EN_plot_exp_fun.pdf", width=7, height=6)

pigpl_EN_exp(car, git, add=F, col="black", lty=1)
pigpl_EN_exp(car, AVI, add=T, col="black", lty=2, lwd=2)
pigpl_EN_exp(car, chap, add=T, col="black", lty=3)
pigpl_EN_exp(car, AVI_mCRI, add=T, col="red", lty=1)
pigpl_EN_exp(car, AVI_chap, add=T, col="red", lty=3, lwd=2)
legend(x="topright", legend=c("mCRI", "AVI", "CI", "AVI_mCRI", "AVI_CI"), lwd=c(1,2,1,1,2), lty=c(1,2,3,1,3), col=c("black","black", "black","red", "red"), title="", cex=0.7, bty="n")
dev.off()



pdf(file="EN_plot_spline_fun.pdf", width=7, height=6)

spl_NE(car, git, add=F, col="black", lty=1)
spl_NE(car, AVI, add=T, col="black", lty=2, lwd=2)
spl_NE(car, chap, add=T, col="black", lty=3)
spl_NE(car, AVI_mCRI, add=T, col="red", lty=1)
spl_NE(car, AVI_chap, add=T, col="red", lty=3, lwd=2)
legend(x="topright", legend=c("mCRI", "AVI", "CI", "AVI_mCRI", "AVI_CI"), lwd=c(1,2,1,1,2), lty=c(1,2,3,1,3), col=c("black","black", "black","red", "red"), title="", cex=0.7, bty="n")

dev.off()



#############
### EN plot with quadratic function
#############


pigpl_EN_qu <- function(car, index, add=T, ...) {

  c1 <- lm(index~car+I(car^2))
  sumc1 <- summary(c1)
  rmse <- sqrt(sum(sumc1$residuals^2) / length(car))
  cq = coef(c1)

  newenv = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/(length(car)-1))
  newenv2 = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/500)


  if (add==TRUE) {
    dv1 <- cq[2] + 2 * cq[3] * newenv2
    EN_car <- rmse/dv1
    lines(newenv2, EN_car, ...) 
  }

  else {
    plot(newenv, car, ylim=c(0,(max(car)+5)), col="white", ylab="NE car", xlab="total car", main="NE for quadratic fit function")
    dv1 <- cq[2] + 2 * cq[3] * newenv2
    EN_car <- rmse/dv1
    lines(newenv2, EN_car, ...) 
  }
}

#############
### EN plot with linear function
#############

pigpl_EN_lin <- function(car, index, add=T, ...) {
  
  c1 <- lm(index~car)
  sumc1 <- summary(c1)
  rmse <- sqrt(sum(sumc1$residuals^2) / length(car))
  cq = coef(c1)
  
  newenv = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/(length(car)-1))
  newenv2 = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/500)
  
  
  if (add==TRUE) {
    dv1 <- cq[2]
    EN_car <- rmse/dv1
    abline(a= EN_car, b=0, ...) 
  }
  
  else {
    plot(newenv, car, ylim=c(0,(max(car)+5)), col="white", ylab="NE car", xlab="total car", main="NE for linear fit function")
    dv1 <- cq[2] 
    EN_car <- rmse/dv1
    abline(a= EN_car, b=0, ...) 
  }
}


#############
### EN plot with exponential function
#############

pigpl_EN_exp <- function(car, index, add=T, ...) {
  
  c1 <- nls(index ~ exp(a + b * car), start = list(a = 0, b = 0))
  sumc1 <- summary(c1)
  rmse <- sqrt(sum(sumc1$residuals^2) / length(car))
  cq = coef(c1)

  
  newenv = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/(length(car)-1))
  newenv2 = seq(min(car, na.rm=TRUE), max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/500)
  
  
  if (add==TRUE) {
    
    
    dv1 <- cq[2]*exp(cq[1] + cq[2]*newenv2)
    EN_car <- rmse/dv1
    lines(newenv2, EN_car, ...)
  }
  
  else {
    plot(newenv, car, ylim=c(0,(max(car)+5)), col="white", ylab="NE car", xlab="total car", main="NE for exponential fit function")
    dv1 <- cq[2]*exp(cq[1] + cq[2]*newenv2)
    EN_car <- rmse/dv1
    lines(newenv2, EN_car, ...) 
  }
}

#############
### EN plot with spline function
#############

spl_NE <- function(car, index, filename, add=TRUE, ...) {
  
  ind.spl <- smooth.spline(car, index, cv=T)
  (ind.spl)
  
  rmse <- sqrt(sum((index-fitted(ind.spl))^2) / length(car))
  rmse
  
  newenv = seq(0, max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/(500))
  newenv2 = seq(0, max(car, na.rm=TRUE), by = (max(car, na.rm=TRUE) - min(car, na.rm=TRUE))/(79))
  
  spl_1dv<-predict(ind.spl,newenv, deriv=1)
  
  EN<-rmse/spl_1dv$y
  
  if (add==TRUE) {
    
    lines(newenv, EN, ...) }
  
  else {
    
    plot(car, newenv2, col="white")
    lines(newenv, EN, ...) }
  
}



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
