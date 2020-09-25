##Median statistics for dry season models

##calculating Median Statistics for each model

#1) for  a model, calculate the R2 value per station, then take the median

#read in actual rainfall
rm(list=ls())


#set wd
setwd("C:/Users/ksanf/Documents/SD/Actual_Rain")
observation=read.csv("actual_DRY.csv")
observation=observation[,-1]


#read in LOOCV reconstruction files and make into list
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Reconstruction")
files = list.files(pattern="*.csv")
reconstruction.list =lapply(files, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))
reconstruction.list=lapply(reconstruction.list, function(x) x[,-1])

val=NULL
median.r2=NULL
r2=vector()

for (a in 1:length(reconstruction.list)){
  for (i in 1:ncol(observation)) {
    mod=lm(observation[,i]~reconstruction.list[[a]][,i])
    r2[i]<-summary(mod)$r.squared
  }
  val=median(r2)
  median.r2[a]=val
}

plot(median.r2, main = "Median R2 Values", ylab="Median R2", ylim=c(.21,.34))
#calc mod.2015 first
abline(h=median.r2.2015, lty=2, col="blue")

hist(median.r2, main="Median R2 Per Model", xlab="Median R2", xlim=c(.21,.34), breaks=15)
points(median.r2.2015, 0, pch=16,  col= "blue", cex=2)


##add 2015 model
mod.2015=read.csv("C:/Users/ksanf/Documents/SD/2015_tests/2015_LOOCV_reconstruction_DRY.csv")

mod.2015=mod.2015[,-1]

r2.2015=vector()

for (i in 1:ncol(mod.2015)) {
  mod=lm(observation[,i]~mod.2015[,i])
  r2.2015[i]<-summary(mod)$r.squared
}
val=median(r2.2015)
median.r2.2015=val

##################################################################################################

#2) for a model, calculate the RMSE value per station, then take the median

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

val=NULL
median.rmse=NULL
rmses=vector()

for (a in 1:length(reconstruction.list)){
  for (i in 1:ncol(observation)) {
    rmse=RMSE(reconstruction.list[[a]][,i], observation[,i])
    rmses[i]=rmse
  }
  val=median(rmses)
  median.rmse[a]=val
}

head(reconstruction.list[[1]])
plot(median.rmse, main = "Median RMSE Values", ylab="Median RMSE", ylim=c(121.5, 135))

length(median.rmse)
abline(h=median.rmse.2015, lty=2, col="blue")

hist(median.rmse, main="Median RMSE Per Model", xlab="Median RMSE", xlim=c(121.5,135),breaks=19)
points(median.rmse.2015, 0, pch=16,  col= "blue", cex=2)





##add 2015

for (i in 1:ncol(mod.2015)) {
  rmse=RMSE(mod.2015[,i], observation[,i])
  rmses[i]=rmse
}
val=median(rmses)
median.rmse.2015=val


##RE


MSE = function(m, o){
  mean((m - o)^2)
}


##need to calculate MSE in order to calculate RE
#code to calculate MSEs for every station

mses=NULL
val=NULL
median.re=NULL
REs=vector()

for (a in 1:length(reconstruction.list)){
  for (i in 1:ncol(observation)) {
    mse=MSE(reconstruction.list[[a]][,i], observation[,i])
    mses[i]=mse
    RE=1-(mses[i]/MSE(mean(observation[,i]), observation[,i]))
    REs[i]=RE
  }
  val=median(REs)
  median.re[a]=val
}

length(median.re)

plot(median.re, main = "Median RE Values", ylab="Median RE", ylim=c(.17,.33))
abline(h=median.re.2015, lty=2, col= "blue")

hist(median.re, main="Median RE Per Model", xlab="Median RE", xlim=c(.17,.36))
points(median.re.2015, 0, pch=16,  col= "blue", cex=2)

max(median.re)



plot(median.rmse, median.r2)

plot(median.re, median.rmse)



for (i in 1:ncol(observation)) {
  mse=MSE(test.data[,i], observation[,i])
  mses[i]=mse
}




REs=NULL

for (k in 1:length(mses)){
  RE=1-(mses[k]/MSE(mean(observation[,k]), observation[,k]))
  REs[k]=RE
}

median(REs)  


#calculate RE for 2015


for (i in 1:ncol(observation)) {
  mse=MSE(mod.2015[,i], observation[,i])
  mses[i]=mse
  RE=1-(mses[i]/MSE(mean(observation[,i]), observation[,i]))
  REs[i]=RE
}

val=median(REs)
median.re.2015=val






