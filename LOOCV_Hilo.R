##Hilo Reconstruction LOOCV

### plot actual with range of estimates for a single station #HILO
rm(list=ls())

setwd("C:/Users/ksanf/Documents/Statistical_DS/LOOCV/Reconstruction")


filenames=list.files(path="C:/Users/ksanf/Documents/Statistical_DS/LOOCV/Reconstruction", full.names=TRUE)
datalist= lapply(filenames, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))


df=data.frame(matrix(1:28))
names(df)="year"

for (i in 1:length(datalist)){
  data=datalist[[i]][,c(1,112)]
  df=merge(df, data, by.x="year", by.y=1)
}



df=t(df)
df=df[-1,]
dim(df)

time=1980:2007
setwd("C:/Users/ksanf/Documents/SD/Actual_Rain")
observation=read.csv("actual_wet.csv", check.names = F)

#remove year column
observation=observation[,-1]


actual=observation[,"87"]


##read in 2015 reconstruction
mod.2015=read.csv("C:/Users/ksanf/Documents/SD/2015_tests/2015_LOOCV_reconstruction_WET.csv", check.names = F)
mod.2015=mod.2015[,-1]

mod.2015=mod.2015[,"87"]

boxplot(df,  names=time, las=2, cex.axis=.7, main= "HILO") 
lines(actual, col="blue")
lines(mod.2015, col="red")




which( colnames(datalist[[1]])=="87" )


