###################################################################################
#Rainfall reconstrution using every possible model combination for PCs 1-4 (WET SEASON)
###################################################################################
##clear everything
rm(list=ls())

#load packages
library(stringr)


#load in dredge tables; navigate to folder where dredge tables were saved

setwd("C:/Users/ksanf/Documents/Test_SD/Dredge_Tables")

dredge_W1<-read.csv("wet_PC1_dredge.csv")
dredge_W2<-read.csv("wet_PC2_dredge.csv")
dredge_W3<-read.csv("wet_PC3_dredge.csv")
dredge_W4<-read.csv("wet_PC4_dredge.csv")

#put dredge tables in a list

dredge.list<-list(dredge_W1, dredge_W2, dredge_W3, dredge_W4)

#load in data frames containing predictor/predictand information

setwd("C:/Users/ksanf/Documents/Test_SD/DFs")

wet_df1<-read.csv("wet_PC1.csv")
wet_df2<-read.csv("wet_PC2.csv")
wet_df3<-read.csv("wet_PC3.csv")
wet_df4<-read.csv("wet_PC4.csv")

#put dataframes in a list

df.list<-list(wet_df1, wet_df2, wet_df3, wet_df4)

##shorten dredge lists
#for each dredge table, remove all models where delta AIC is greater than 2. 
#remove all columns that contain other information

dredge.list<-lapply(dredge.list, function(x) subset(x, subset=delta<2))
dredge.list<-lapply(dredge.list, function(x) x[,2:(ncol(x)-6)])

#for every dredge table, set the coefficient in the cell to the column (variable) name
for (j in 1:length(dredge.list)){
  for (i in 1:ncol(dredge.list[[j]])){
    dredge.list[[j]][!is.na(dredge.list[[j]][,i]), i ]<-names(dredge.list[[j]])[i]
  }
}

#add a column to each drege table that lists the included variables, skips NA variables
#separte the variable names with a "+" sign
for (j in 1:length(dredge.list)){
  dredge.list[[j]]$varlist <- apply(dredge.list[[j]], 1, function(x) 
    paste(str_trim(x[!is.na(x)]), collapse=" + "))
}

#save column of listed variable names as new objects

varlist1<-dredge.list[[1]]$varlist
varlist2<-dredge.list[[2]]$varlist
varlist3<-dredge.list[[3]]$varlist
varlist4<-dredge.list[[4]]$varlist

#save variable lists to file

write.csv(varlist1, file="C:/Users/ksanf/Documents/Test_SD/VarLists_wet/varlist_d1.csv", row.names = FALSE)
write.csv(varlist2, file="C:/Users/ksanf/Documents/Test_SD/Varlists_wet/varlist_d2.csv", row.names = FALSE)
write.csv(varlist3, file="C:/Users/ksanf/Documents/Test_SD/Varlists_wet/varlist_d3.csv", row.names = FALSE)
write.csv(varlist4, file="C:/Users/ksanf/Documents/Test_SD/Varlists_wet/varlist_d4.csv", row.names = FALSE)


# for every list of variables for each dredge table, make a linear model using varlist
#as the predictors and "y" of the corresponding PC as the predictand (from df)
#output the model fits for each model

for (k in 1:length(dredge.list)){
  fits_df=data.frame(matrix(ncol=nrow(dredge.list[[k]]), nrow=28))
  j=0
  for (i in dredge.list[[k]]$varlist){
    j=j+1
    model <- lm(paste("y ~", i), data=df.list[[k]])
    fits=(model$fit)
    fits.df<-data.frame(fits)
    fits_df[,j]<-fits.df[,1]
    colnames(fits_df)[j]<-paste0("PC", k, "_model_",  j)
  }
  write.csv(fits_df, file=paste0("C:/Users/ksanf/Documents/Test_SD/Fits_for_Reconstruction/wet_fits_PC", k, ".csv"))
}




###########################################################

#RAINFALL RECONSTRUCTION

######################################################################

##Xmatrix was created at beginning of script "PCA_Comp_Analysis_Dredging"; 
#the following code is a continuation from that script
Xmean<-apply(Xmatrix,2,mean)

Xvar<-apply(Xmatrix,2,var)

#set wokrking directory to folder where fits for linear models are saved;
#remove first column (record number column)
setwd("C:/Users/ksanf/Documents/Test_SD/Fits_for_Reconstruction")
fits_PC1=read.csv("wet_fits_PC1.csv")
fits_PC1=fits_PC1[2:ncol(fits_PC1)]

fits_PC2=read.csv("wet_fits_PC2.csv")
fits_PC2=fits_PC2[2:ncol(fits_PC2)]

fits_PC3=read.csv("wet_fits_PC3.csv")
fits_PC3=fits_PC3[2:ncol(fits_PC3)]

fits_PC4=read.csv("wet_fits_PC4.csv")
fits_PC4=fits_PC4[2:ncol(fits_PC4)]


#set pc truncation number
ntrunc=4
#matrix multiplication with eigenvectors
Xanotrunc<-pc[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])


#reconstruction for every combination of PC1, PC2, Pc3, PC4 models
#loops through list of models for each PC
#saves fully reconstructed dataset for every model combination:
##naming convention ABCD (A = model number of pc1, B = model number of PC2, C= model number of PC3, D= model number of PC4)
for (a in 1:ncol(fits_PC1)){
  for (b in 1:ncol(fits_PC2)){
    for (c in 1:ncol(fits_PC3)){
      for (d in 1:ncol(fits_PC4)){
        pcfit<-cbind(fits_PC1[a],fits_PC2[b],fits_PC3[c],fits_PC4[d])
        pcfit=as.matrix(pcfit)
        write.csv(pcfit, file=paste0("C:/Users/ksanf/Documents/Test_SD/Fits_All_Combos_Wet/fits", a,b,c,d, ".csv"))
        Xanofit<-pcfit[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
        Xfit  <-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        Xtrunc<-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        i=1
        while (i<=length(Xvar)) {
          Xfit[,i]<-Xanofit[,i]*sqrt(Xvar[i])+Xmean[i]
          Xtrunc[,i]<-Xanotrunc[,i]*sqrt(Xvar[i])+Xmean[i]
          i=i+1
        }
        write.csv(Xtrunc, file=paste0("C:/Users/ksanf/Documents/Test_SD/Actual_Rain/Truncated_Precip_Wet.csv"))
        write.csv(Xfit, file=paste0("C:/Users/ksanf/Documents/Test_SD/Fitted_Precip_All_Combos_Wet/xfit", a,b,c,d, ".csv"))
        
      }
    }
  }
}


###################################################################
##PLOTTING
#plot for a single station (change number for different station)

plot(Xmatrix[,17], type = "line", col= "blue", main="Station 17, All Years")
lines(Xtrunc[,17], type="line", col="red")
lines(Xfit[,17], type = "line", col="green")
legend(20, 600, legend=c("actual", "truncated", "projected"),
       col=c("blue", "red",  "green"), lty=1, cex=0.8)


##################################################################
#LOOCV
###############################################################

##clear everything
rm(list=ls())

#repeat steps to calculate fits - using LOOCV

setwd("C:/Users/ksanf/Documents/Test_SD/DFs")
wet_df1<-read.csv("wet_PC1.csv")
wet_df2<-read.csv("wet_PC2.csv")
wet_df3<-read.csv("wet_PC3.csv")
wet_df4<-read.csv("wet_PC4.csv")

df.list<-list(wet_df1, wet_df2, wet_df3, wet_df4)

#read in varlist for each PC
setwd("C:/Users/ksanf/Documents/Test_SD/VarLists_wet")
varlist1<-read.csv("varlist_d1.csv")
varlist2<-read.csv("varlist_d2.csv")
varlist3<-read.csv("varlist_d3.csv")
varlist4<-read.csv("varlist_d4.csv")

var.lists<-list(varlist1, varlist2, varlist3, varlist4)


##USE LOOCV TO CALCULATE FITS
#do PCs one at a time
#PC1
#calculate and save fits for PC1, all models
for (a in 1:nrow(var.lists[[1]])){
  fit=NULL
  for (i in 1:nrow(df.list[[1]])){ #1:28 years
    datause = df.list[[1]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[1]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[1]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC1/PC1_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC2
#calculate and save fits for PC2, all models
for (a in 1:nrow(var.lists[[2]])){
  fit=NULL
  for (i in 1:nrow(df.list[[2]])){ #1:28 years
    datause = df.list[[2]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[2]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[2]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC2/PC2_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC3
#calculate and save fits for PC3, all models
for (a in 1:nrow(var.lists[[3]])){
  fit=NULL
  for (i in 1:nrow(df.list[[3]])){ #1:28 years
    datause = df.list[[3]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[3]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[3]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC3/PC3_predict_mod_", a, ".csv"), row.names=TRUE)
}


#PC4
#calculate and save fits for PC4, all models
for (a in 1:nrow(var.lists[[4]])){
  fit=NULL
  for (i in 1:nrow(df.list[[4]])){ #1:28 years
    datause = df.list[[4]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[4]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[4]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC4/PC4_predict_mod_", a, ".csv"), row.names=TRUE)
}

#########################################################################
#USE OUTPUT FROM PREVIOUS STEP TO RECONSTRUCT RAINFALL 
#caluclate rmse and mae within the same loop


##Calculating RMSE
#m is for model (fitted) values, o is for observed (true) values.
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

##calculating MAE
#actual and predicted are vectors
library(ModelMetrics)


##read in fits as a list
#PC1
setwd("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC1")

LOOCV_fits1 = list.files(pattern="*.csv")
LOOCV_fits1 = lapply(LOOCV_fits1, read.csv)

#PC2
setwd("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC2")
LOOCV_fits2 = list.files(pattern="*.csv")
LOOCV_fits2 = lapply(LOOCV_fits2, read.csv)


#PC#
setwd("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC3")
LOOCV_fits3 = list.files(pattern="*.csv")
LOOCV_fits3 = lapply(LOOCV_fits3, read.csv)


#PC4
setwd("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Predict_PC4")
LOOCV_fits4 = list.files(pattern="*.csv")
LOOCV_fits4 = lapply(LOOCV_fits4, read.csv)



#reconstruction loop using LOOCV_fits
#need pc, eigenvector, and Xmatrix objects from PCA_Comp_Analysis_Dredging code
#calculate RMSE and mae in same code
ntrunc=4
Xanotrunc<-pc[1,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
Xmean<-apply(Xmatrix,2,mean)
Xvar<-apply(Xmatrix,2,var)
Xfit_all=NULL
vals=data.frame(matrix(nrow=864, ncol=2))
vals2=data.frame(matrix(nrow=864, ncol=2))
r2=data.frame(matrix(nrow=903, ncol=2))
j=0
stations=vector()
stations=colnames(Xmatrix)

for (a in 1:length(LOOCV_fits1)){
  for (b in 1:length(LOOCV_fits2)){
    for (c in 1:length(LOOCV_fits3)){
      for (d in 1:length(LOOCV_fits4)){
        for (n in 1:28){
          pcfit<-cbind(LOOCV_fits1[[a]][n,2],LOOCV_fits2[[b]][n,2],LOOCV_fits3[[c]][n,2],LOOCV_fits4[[d]][n,2])
          pcfit=as.matrix(pcfit)
          Xanofit<-pcfit[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
          Xfit  <-matrix(0, 1,ncol(Xmatrix))
          Xtrunc<-matrix(0, 1,ncol(Xmatrix))
          i=1
          while (i<=length(Xvar)) {
            Xfit[,i]<-Xanofit[,i]*sqrt(Xvar[i])+Xmean[i]
            Xtrunc[,i]<-Xanotrunc[,i]*sqrt(Xvar[i])+Xmean[i]
            i=i+1
          }
          Xfit_all <- rbind(Xfit_all, Xfit)
        }
        colnames(Xfit_all)=stations
        write.csv(Xfit_all, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/Reconstruction/Reconstruction_", a, "_", b, "_", c, "_", d, ".csv"), row.names=TRUE)
        rmsevals=RMSE(Xfit_all, Xmatrix)
        MAEvals=mae(Xmatrix, Xfit_all)
        rmse=data.frame(matrix(ncol=2, nrow=1))
        MAE=data.frame(matrix(ncol=2, nrow=1))
        
        rmse=c(paste0(a,"_", b, "_", c, "_", d), rmsevals)
        MAE=c(paste0(a,"_", b, "_", c, "_", d), MAEvals)
        #
        
        for (i in 1:ncol(Xmatrix)) {
          mod=lm(Xmatrix[,i]~Xfit_all[,i])
          r2[i,2]<-summary(mod)$r.squared
          
        }
        
        colnames(r2)<-c("Station", paste0(a, "_", b, "_", c, "_", d))
        r2[,1]<-stations
        #
        write.csv(r2, file=paste0("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/LOOCV_R2/R2_", a, "_", b, "_", c, "_", d, ".csv"), row.names=F)
        
        Xfit_all=NULL
        
        j=j+1
        vals[j,]=rmse
        vals2[j,]=MAE
        colnames(vals)<-c("Combination", "RMSE")
        colnames(vals2)<-c("Combination", "MAE")
        write.csv(vals, file="C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/RMSE_MAE/rmse.csv", row.names=FALSE)
        write.csv(vals2, file="C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/RMSE_MAE/mae.csv", row.names=FALSE)
        
      }
    }
  }
}

#############################################################################
#calculate median r2 per station and map

#########################################################################

#set wd
rm(list=ls())
setwd("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/LOOCV_R2")

multmerge=function(mypath) {
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist= lapply(filenames, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))
  Reduce(function(x,y) merge(x,y), datalist)}

#merge all r2 per station per model

mymergeddata = multmerge("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/LOOCV_R2")
write.csv(mymergeddata, file="C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/merged_r2.csv", row.names=FALSE)
mymergeddata=read.csv("C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/merged_r2.csv", check.names= FALSE)
#calculate median per station
mymergeddata$median=apply(mymergeddata[,2:865],1, median, na.rm = TRUE)
#take only median and save
wet_station_r2=mymergeddata[,c(1,866)]
write.csv(wet_station_r2, file= "C:/Users/ksanf/Documents/Test_SD/LOOCV_WET/wet_r2_per_station.csv", row.names = FALSE)

#scatterplot and histogram of median
plot(mymergeddata$median, xlab = "Station", ylab = "R^2", main = "Median R^2 Value per Station")
hist(mymergeddata$median, xlab = "r-squared", main= "Median R^2 Value per Station", xlim=c(0,1), breaks=10)
