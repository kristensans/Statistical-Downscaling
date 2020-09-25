##HONOLULU Reconstruction LOOCV

### plot actual with range of estimates for a single station #Hono Airport
rm(list=ls())

setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/Reconstruction")


filenames=list.files(path="C:/Users/ksanf/Documents/SD/LOOCV_WET/Reconstruction", full.names=TRUE)
datalist= lapply(filenames, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))


df=data.frame(matrix(1:28))
names(df)="year"

for (i in 1:length(datalist)){
  data=datalist[[i]][,c(1,531)]
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


actual=observation[,"703"]


##read in 2015 reconstruction
mod.2015=read.csv("C:/Users/ksanf/Documents/SD/2015_tests/2015_LOOCV_reconstruction_WET.csv", check.names = F)
mod.2015=mod.2015[,-1]

mod.2015=mod.2015[,"703"]




boxplot(df,  names=time, las=2, cex.axis=.7, main= "HONOLULU") 
lines(actual, col="blue")
lines(mod.2015, col="red")


colnames(datalist[[1]])=stations


which( colnames(datalist[[1]])=="703" )

which(colnames(Xfit)=="703")

hono2015=Xfit[,111]

