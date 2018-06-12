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

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "rmse.txt"  
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
      

      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

      c1 <- lm(car~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      rmse1 <- sqrt(sum(c1_sum$residuals^2) / 106)
      rmse2 <- sqrt(sum(c2_sum$residuals^2) / 106)
      rmse3 <- sqrt(sum(c3_sum$residuals^2) / 106)
      

      
      AVI_index <-paste("##### ##### AVI1: ", k, " ", j, " ", h)
      print(AVI_index, attach=TRUE)
      print(rmse1, attach=TRUE)
      
      AVI_index <-paste("##### ##### AVI2: ", k, " ", j, " ", h)
      print(AVI_index, attach=TRUE)
      print(rmse2, attach=TRUE)
      
      AVI_index <-paste("##### ##### AVI3: ", k, " ", j, " ", h)
      print(AVI_index, attach=TRUE)
      print(rmse3, attach=TRUE)

        
    }
  }
}


sink()


#################################################
#################################################


## test best found results

k <- 1
j <- 13
i <- 15

###
### on beech data

wave1 <- wavel[k]
wave2 <- wavel[j]
wave3 <- wavel[i]

band1 <- refl_b[,k]
band2 <- refl_b[,j]
band3 <- refl_b[,i]


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(car_b~AVI1+I(AVI1^2))
c1_sum <- summary(c1)

c2 <- lm(car_b~AVI2+I(AVI2^2))
c2_sum <- summary(c2)

c3 <- lm(car_b~AVI3+I(AVI3^2))
c3_sum <- summary(c3)

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


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(car_m~AVI1+I(AVI1^2))
c1_sum <- summary(c1)

c2 <- lm(car_m~AVI2+I(AVI2^2))
c2_sum <- summary(c2)

c3 <- lm(car_m~AVI3+I(AVI3^2))
c3_sum <- summary(c3)

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


v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)


AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

c1 <- lm(car_c~AVI1+I(AVI1^2))
c1_sum <- summary(c1)

c2 <- lm(car_c~AVI2+I(AVI2^2))
c2_sum <- summary(c2)

c3 <- lm(car_c~AVI3+I(AVI3^2))
c3_sum <- summary(c3)

c1_sum
c2_sum
c3_sum


##################
##################
## chestnut

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "car_chest.txt"  
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
      
      band1 <- refl_c[,k]
      band2 <- refl_c[,j]
      band3 <- refl_c[,h]
      
      
      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
      
      
      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      c1 <- lm(car_c~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car_c~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car_c~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      th <- 0.7
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
      
    }
  }
}


sink()

##################
##################
## beech

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "car_beech.txt"  
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
      
      band1 <- refl_b[,k]
      band2 <- refl_b[,j]
      band3 <- refl_b[,h]
      
      
      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
      
      
      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      c1 <- lm(car_b~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car_b~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car_b~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      th <- 0.85
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
      
    }
  }
}


sink()



##################
##################
## maple

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "car_maple.txt"  
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
      
      band1 <- refl_m[,k]
      band2 <- refl_m[,j]
      band3 <- refl_m[,h]
      
      
      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
      
      
      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      c1 <- lm(car_m~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car_m~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car_m~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      th <- 0.6
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
      
    }
  }
}


sink()

###################
##################
## beech + chestnut

car_bc <- c(car_c, car_b) 
refl_bc <- rbind(refl_c, refl_b)


# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "car_beech_chest.txt"  
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
      
      band1 <- refl_bc[,k]
      band2 <- refl_bc[,j]
      band3 <- refl_bc[,h]
      
      
      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
      
      
      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      c1 <- lm(car_bc~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car_bc~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car_bc~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      th <- 0.85
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
      
    }
  }
}


sink()

##################
## low car

t <-  car < 140 & car > 100
plot(car[t], chl[t])
cor(car[t], chl[t])

car_sm <- car[t]
refl_sm <- refl_car[t,]

# open textfile to store results
setwd("/home/fabi/KIT/11_AG_Angular_Index/1_new_runs/out")
path = "car_sm_90_140.txt"  
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
      
      band1 <- refl_sm[,k]
      band2 <- refl_sm[,j]
      band3 <- refl_sm[,h]
      
      
      v_a <- sqrt((wave1 - wave2)^2+(band1 - band2)^2)
      v_b <- sqrt((wave2 - wave3)^2+(band2 - band3)^2)
      v_c <- sqrt((wave1 - wave3)^2+(band1 - band3)^2)
      
      
      AVI1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
      AVI2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
      AVI3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))
      
      c1 <- lm(car_sm~AVI1+I(AVI1^2))
      c1_sum <- summary(c1)
      
      c2 <- lm(car_sm~AVI2+I(AVI2^2))
      c2_sum <- summary(c2)
      
      c3 <- lm(car_sm~AVI3+I(AVI3^2))
      c3_sum <- summary(c3)
      
      th <- 0.4
      
      if (c1_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI1: ", k, " ", j, " ", h)
      if (c1_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c1_sum$adj.r.squared > th) print(c1_sum$adj.r.squared, attach=TRUE)
      
      if (c2_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI2: ", k, " ", j, " ", h)
      if (c2_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c2_sum$adj.r.squared > th) print(c2_sum$adj.r.squared, attach=TRUE)
      
      if (c3_sum$adj.r.squared > th) AVI_index <- paste("##### ##### AVI3: ", k, " ", j, " ", h)
      if (c3_sum$adj.r.squared > th) print(AVI_index, attach=TRUE)
      if (c3_sum$adj.r.squared > th) print(c3_sum$adj.r.squared, attach=TRUE)
      
    }
  }
}


sink()



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
# scale function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# calculate lm between index and carotenoids

chapavi <- range01(chap)+range01((AVI2(wavel[a],wavel[b],wavel[c],refl_o1[,a],refl_o1[,b],refl_o1[,c])))
gpqu <- lm(car~chapavi)
sumgpqu<-summary(gpqu)
rmse <- sqrt(sum(sumgpqu$residuals^2) / 106)
rmse

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


