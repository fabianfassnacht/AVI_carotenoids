pkgs<-c("rgdal","caret","raster","foreign", "randomForest", "mgcv", "stats", "pspline")
lapply(pkgs,require, character.only=T)


#Set working directory and load libraries

setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs")

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

##########
##########
# create graphs

pigPlot(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, AVI, AVI_m, AVI_c, AVI_b, filename="1_AVI.pdf", indexname="AVI")
pigPlot(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, git, git_m, AVI_c, AVI_b, filename="2_mCRI.pdf", indexname="mCRI")
pigPlot(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, chap, chap_m, chap_c, chap_b, filename="3_Chappele.pdf", indexname="CI")
pigPlot(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, git_AVI, git_AVI_m, git_AVI_c, git_AVI_b, filename="4_merged_AVI_mCRI.pdf", indexname="Merged AVI mCRI")
pigPlot(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, chap_AVI, chap_AVI_m, chap_AVI_c, chap_AVI_b, filename="5_merged_AVI_CI.pdf", indexname="Merged AVI CI")

pigPlot2(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, filename="6_chl_var.pdf")

pigPlot3(car, chl, AVI, species, filename="7_AVI_all_in_one.pdf", indexname="AVI")
pigPlot3(car, chl, git, species, filename="8_mCRI_all_in_one.pdf", indexname="mCRI")
pigPlot3(car, chl, chap, species, filename="9_chap_all_in_one.pdf", indexname="CI")
pigPlot3(car, chl, chap_AVI, species, filename="10_merged_chap_AVI_all_in_one.pdf", indexname="merged AVI CI")
pigPlot3(car, chl, git_AVI, species, filename="11_merged_mCRI_AVI_all_in_one.pdf", indexname="merged AVI mCRI")

setwd("D:/KIT_Forschung/11_AG_Angular_Index/99_Paper/2_Graphs_fitted_function")
pigPlot4(car, AVI, species, filename="2_AVI_species_pch.pdf", indexname="AVI")
pigPlot4(car, git, species, filename="3_mCRI_species_pch.pdf", indexname="mCRI")
pigPlot4(car, chap, species, filename="4_CI_species_pch.pdf", indexname="CI")
pigPlot4(car, chap_AVI, species, filename="5_AVI_CI_species_pch.pdf", indexname="merged AVI CI")
pigPlot4(car, git_AVI, species, filename="6_AVI_mCRI_species_pch.pdf", indexname="merged mCRI AVI")

setwd("D:/KIT_Forschung/11_AG_Angular_Index/99_Paper/2_Graphs_fitted_function")
pigPlot5(car, AVI, species, filename="2_AVI_species_pch_nspl.pdf", indexname="AVI")
pigPlot5(car, git, species, filename="3_mCRI_species_pch_nspl.pdf", indexname="mCRI")
pigPlot5(car, chap, species, filename="4_CI_species_pch_nspl.pdf", indexname="CI")
pigPlot5(car, chap_AVI, species, filename="5_AVI_CI_species_pch_nspl.pdf", indexname="merged AVI CI")
pigPlot5(car, git_AVI, species, filename="6_AVI_mCRI_species_pch_nspl.pdf", indexname="merged mCRI AVI")



dev.off()
###################
##### calculate RMSE values (linear / quadratic)

rmse_lm(car,AVI)
rmse_lm(chl,AVI)
rmse_lm(car_m,AVI_m)
rmse_lm(chl_m,AVI_m)
rmse_lm(car_c,AVI_c)
rmse_lm(chl_c,AVI_c)
rmse_lm(car_b,AVI_b)
rmse_lm(chl_b,AVI_b)

rmse_lm(car,git)
rmse_lm(chl,git)
rmse_lm(car_m,git_m)
rmse_lm(chl_m,git_m)
rmse_lm(car_c,git_c)
rmse_lm(chl_c,git_c)
rmse_lm(car_b,git_b)
rmse_lm(chl_b,git_b)

rmse_lm(car,chap)
rmse_lm(chl,chap)
rmse_lm(car_m,chap_m)
rmse_lm(chl_m,chap_m)
rmse_lm(car_c,chap_c)
rmse_lm(chl_c,chap_c)
rmse_lm(car_b,chap_b)
rmse_lm(chl_b,chap_b)

rmse_lm(car,chap_AVI)
rmse_lm(chl,chap_AVI)
rmse_lm(car_m,chap_AVI_m)
rmse_lm(chl_m,chap_AVI_m)
rmse_lm(car_c,chap_AVI_c)
rmse_lm(chl_c,chap_AVI_c)
rmse_lm(car_b,chap_AVI_b)
rmse_lm(chl_b,chap_AVI_b)

rmse_lm(car,git_AVI)
rmse_lm(chl,git_AVI)
rmse_lm(car_m,git_AVI_m)
rmse_lm(chl_m,git_AVI_m)
rmse_lm(car_c,git_AVI_c)
rmse_lm(chl_c,git_AVI_c)
rmse_lm(car_b,git_AVI_b)
rmse_lm(chl_b,git_AVI_b)



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

##########################
##########################
## pig plot 1
##########################


pigPlot <- function(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, index, index_m, index_c, index_b, filename, indexname) {

#plot indices in 4 panel layout

setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs/out_pdfs")
pdf(file=filename, width=8, height=8)

par(mfrow=c(2,2))


######### all data

# set for the right side
#plot reflectance values
plot(car, index, axes=F, lty=2,  xlim=range(car),ann=FALSE, cex=0.8, pch=1)
axis(3, xlim=range(car),col="black")
mtext("total car content [mg/m²]",side=3,col="black",line=2.5, cex=0.8)

#prepare second plot
par(new=T)
mtext(indexname,side=2,line=2.5, cex=0.8)
axis(2, col="black",col.axis="black")
box()
par(new=T)

# plot GA score values
plot(chl, index, xlim=range(chl), axes=F, main="all data", ann=FALSE, pch=3, cex=0.6, col="grey70")
mtext("total chl content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
axis(1, xlim=range(chl), col="black",col.axis="black")
legend(x="topleft", legend=c("chl", "car"), pch=c(3,1), col=c("grey70","black"), title="all data")

########### maple data

# set for the right side
#plot reflectance values
plot(car_m, index_m, axes=F, lty=2,  xlim=range(car_m),ann=FALSE, cex=0.8, pch=1)
axis(3, xlim=range(car_m),col="black")
mtext("total car content [mg/m²]",side=3,col="black",line=2.5, cex=0.8)

#prepare second plot
par(new=T)
mtext(indexname,side=2,line=2.5, cex=0.8)
axis(2, col="black",col.axis="black")
box()
par(new=T)

# plot GA score values
plot(chl_m, index_m, xlim=range(chl_m), axes=F, main="maple data", ann=FALSE, pch=3, cex=0.6, col="grey70")
mtext("total chl content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
axis(1, xlim=range(chl_m), col="black",col.axis="black")
legend(x="topleft", legend=c("chl", "car"), pch=c(3,1), col=c("grey70","black"),title="maple data")

########### chestnut data

# set for the right side
#plot reflectance values
plot(car_c, index_c, axes=F, lty=2,  xlim=range(car_c),ann=FALSE, cex=0.8, pch=1)
axis(3, xlim=range(car_c),col="black")
mtext("total car content [mg/m²]",side=3,col="black",line=2.5, cex=0.8)

#prepare second plot
par(new=T)
mtext(indexname,side=2,line=2.5, cex=0.8)
axis(2, col="black",col.axis="black")
box()
par(new=T)

# plot GA score values
plot(chl_c, index_c, xlim=range(chl_c), axes=F, main="chestnut data", ann=FALSE, pch=3, cex=0.6, col="grey70")
mtext("total chl content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
axis(1, xlim=range(chl_c), col="black",col.axis="black")
legend(x="topleft", legend=c("chl", "car"), pch=c(3,1), col=c("grey70","black"), title="chestnut data")


########### beech data

# set for the right side
#plot reflectance values
plot(car_b, index_b, axes=F, lty=2, xlim=range(car_b),ann=FALSE, cex=0.8, pch=1)
axis(3, xlim=range(car_b),col="black")
mtext("total car content [mg/m²]",side=3,col="black",line=2.5, cex=0.8)

#prepare second plot
par(new=T)
mtext(indexname,side=2,line=2.5, cex=0.8)
axis(2, col="black",col.axis="black")
box()
par(new=T)

# plot GA score values
plot(chl_b, index_b, xlim=range(chl_b), axes=F, main="beech data", ann=FALSE, pch=3, cex=0.6, col="grey70")
mtext("total chl content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
axis(1, xlim=range(chl_b), col="black",col.axis="black")
legend(x="topleft", legend=c("chl", "car"), pch=c(3,1), col=c("grey70","black"), title="beech data")

#### close device
dev.off()

}


##########################
##########################
## pig plot 2
##########################



pigPlot2 <- function(car, car_m, car_c, car_b, chl, chl_m, chl_c, chl_b, filename=NULL) {
  
  #plot indices in 4 panel layout
  
  if(!is.null(filename)) {
    setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs/out_pdfs")
    pdf(file=filename, width=8, height=8)
  }
  par(mfrow=c(2,2))
  
  ######### all data
  
  # set for the right side
  #plot reflectance values
  plot(car, chl, axes=F, lty=2, ann=F, cex=0.8, pch=1)

  
  #prepare second plot
  par(new=T)
  mtext("total chl content [mg/m²]",side=2,line=2.5, cex=0.8, col="black")
  mtext("all data",side=3, cex=1, line = 0.5)
  axis(2, col="black",col.axis="black")
  box()
  par(new=T)
  
  # plot GA score values
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  axis(1, col="black",col.axis="black")
  
  ########### maple data
  
  # set for the right side
  #plot reflectance values
  
  plot(car_m, chl_m, axes=F, lty=2, ann=F, cex=0.8, pch=1)
  
  
  #prepare second plot
  par(new=T)
  mtext("total chl content [mg/m²]",side=2,line=2.5, cex=0.8)
  mtext("maple data",side=3, cex=1, line = 0.5)
  axis(2, col="black",col.axis="black")
  box()
  par(new=T)
  
  # plot GA score values
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  axis(1, col="black",col.axis="black")
  
  ########### chestnut data
  
  # set for the right side
  #plot reflectance values
  plot(car_c, chl_c, axes=F, lty=2, ann=F, cex=0.8, pch=1)
  
  
  #prepare second plot
  par(new=T)
  mtext("total chl content [mg/m²]",side=2,line=2.5, cex=0.8)
  mtext("chestnut data",side=3, cex=1, line = 0.5)
  axis(2, col="black",col.axis="black")
  box()
  par(new=T)
  
  # plot GA score values
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  axis(1, col="black",col.axis="black")
  
  
  ########### beech data
  
  # set for the right side
  #plot reflectance values
  plot(car_b, chl_b, axes=F, lty=2, ann=F, cex=0.8, pch=1)
  
  
  #prepare second plot
  par(new=T)
  mtext("total chl content [mg/m²]",side=2,line=2.5, cex=0.8)
  mtext("beech data",side=3, cex=1, line = 0.5)
  axis(2, col="black",col.axis="black")
  box()
  par(new=T)
  
  # plot GA score values
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  axis(1, col="black",col.axis="black")
  
  #### close device
  if(!is.null(filename)) {
    dev.off()
  }
}


##########################
##########################
## pig plot 3
##########################


pigPlot3 <- function(car, chl, index, species, filename, indexname) {
  
  colX1 = c("darkmagenta","chartreuse3","darkorange")[as.factor(species)]
  colX2 = c("grey70","grey80","grey90")[as.factor(species)] 
  #plot indices in 4 panel layout
  
  setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs/out_pdfs")
  pdf(file=filename, width=5.5, height=5.5)
  
  #par(mfrow=c(2,2))

  ######### all data
  
  # set for the right side
  #plot reflectance values
  plot(car, index, axes=F, lty=2,  xlim=range(car),ann=FALSE, cex=0.8, pch=1, col=colX1)
  axis(3, xlim=range(car),col="black")
  mtext("total car content [mg/m²]",side=3,col="black",line=2.5, cex=0.8)
  
  #prepare second plot
  par(new=T)
  mtext(indexname,side=2,line=2.5, cex=0.8)
  axis(2, col="black",col.axis="black")
  box()
  par(new=T)
  
  # plot GA score values
  plot(chl, index, xlim=range(chl), axes=F, main="all data", ann=FALSE, pch=3, cex=0.6, col=colX2)
  mtext("total chl content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  axis(1, xlim=range(chl), col="black",col.axis="black")
  legend(x="topleft", legend=c("chl m", "car m","chl c", "car c","chl b", "car b"), pch=c(3,1,3,1,3,1), col=c("grey70","darkmagenta","grey80", "chartreuse3","grey90","darkorange"))
  
  
  dev.off()

}


##########################
##########################
## pig plot 4
##########################


pigPlot4 <- function(car, index, species, filename, indexname) {
  
  pchX1 = c(1,2,3)[as.factor(species)]

  #setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs/out_pdfs")
  pdf(file=filename, width=5.5, height=5.5)

  ######### all data
  
  # set for the right side
  #plot reflectance values
  plot(car, index, axes=F, lty=2,  xlim=range(car),ann=FALSE, cex=0.8, pch=pchX1)
  axis(1, xlim=range(car),col="black")
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  
  #prepare second plot
  par(new=T)
  mtext(indexname,side=2,line=2.5, cex=0.8)
  axis(2, col="black",col.axis="black")
  box()
  #par(new=T)
  
  legend(x="topleft", legend=c( "car m","car c", "car b"), pch=c(1,2,3))
  
  ind.spl <- smooth.spline(car, index, cv=T, nknots=5)
  lines(ind.spl, col="blue")
  
  ind.spl2 <- smooth.spline(index,car, cv=T, nknots=5)
  rmse <- sqrt(sum((car-fitted(ind.spl2))^2) / length(car))
  rmse_r <- round(rmse, digits=2)
  r2 <- (cor(fitted(ind.spl2), car))^2
  r2_r <- round(r2, digits=2)
  mtext(side=1, line=-1.0, text=paste("RMSE: ", rmse_r, "  r²: ", r2_r, "  "), padj=-0.2, adj=1)

  dev.off()


}



##########################
##########################
## pig plot 5
##########################




pigPlot5 <- function(car, index, species, filename, indexname) {
  
  pchX1 = c(1,2,3)[as.factor(species)]
  
  #setwd("D:/KIT_Forschung/11_AG_Angular_Index/1_new_runs/out_pdfs")
  pdf(file=filename, width=5.5, height=5.5)
  
  ######### all data
  
  # set for the right side
  #plot reflectance values
  plot(car, index, axes=F, lty=2,  xlim=range(car),ann=FALSE, cex=0.8, pch=pchX1)
  axis(1, xlim=range(car),col="black")
  mtext("total car content [mg/m²]",side=1,col="black",line=2.5, cex=0.8)
  
  #prepare second plot
  par(new=T)
  mtext(indexname,side=2,line=2.5, cex=0.8)
  axis(2, col="black",col.axis="black")
  box()
  #par(new=T)
  
  legend(x="topleft", legend=c( "car b","car c", "car m"), pch=c(1,2,3))
 
  ind.spl2 <- sm.spline(index, car, df=3)
  
  
  
  lines(sm.spline(car, index, df=3), lty=2, col = "blue")

  newenv = seq(min(index, na.rm=TRUE), max(index, na.rm=TRUE), by = (max(index, na.rm=TRUE) - min(index, na.rm=TRUE))/500)

  y_pred<-predict(ind.spl2, index)
  #y_pred <- ind.spl2$ysmth
  #y_pred <- y_pred[,1]
  rmse <- sqrt(sum((car-y_pred)^2) / length(car))
  rmse_r <- round(rmse, digits=2)
  r2 <- (cor(y_pred, car))^2
  r2_r <- round(r2, digits=2)
  mtext(side=1, line=-1.0, text=paste("RMSE: ", rmse_r, "  r²: ", r2_r, "  "), padj=-0.2, adj=1)
  
  dev.off()
  
  
}
