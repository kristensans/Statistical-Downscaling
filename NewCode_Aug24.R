###This code was last updated on Aug 24, 2020. It is updated to use the new regridded
#variables skt and air 2m, since it was discovered on August 19th that the variables dowloaded had a 
#different grid size (192 vs 144). Oliver regridded these 2 variables to match the 144 grid size.
#directories are also edited to output all files to a new folder called SD.

#this code includes: PCA on rainfall, composite analysis,projection index calculation, multicollinearity analysis
#with correlation matrices, dredging, calculating model fits for competing models, resonstructing rainfall using every
#model combo, LOOCV, and calculating median statistics (r^2, rmse, mae) per station
#NOTE: median statitsics per model can be found in different script  

##clear everything
rm(list=ls())
#######################################################
#all created functions remain the same
#######################################################

#load packages
require(ncdf4)
require(Matrix)
library(car)
library(stats)
library(arm)
library(MuMIn)
library(boot)


#set working directory
setwd("C:/Users/ksanf/Documents/SD/Predictor_Variables")

################################################################
#create functions for data processing
################################################################

load_rain_data<-function(datapath,filename,season='wet',typ='CSV',index='year'){
  # input parameter:
  # datapath: string with the local directory path
  # filename: string with the file name (file type CSV expected)
  # season: 'wet' or 'dry'
  # index: a string that identifies the column with the time information
  # used as an index system for the samples in the other columns
  print(paste("load data from CSV file ",filename))
  print (paste("local directory path: ",datapath))
  filename<-paste(datapath,filename,sep='/')
  df<-as.data.frame(read.csv(file=filename,header = TRUE, stringsAsFactors = FALSE, check.names = F))
  if (index %in% names(df)){
    print (paste("index column '",index,"' found in data table. Ok!",sep=''))
  }
  return(df)
}


subsample_index<-function(X,i1,i2){
  # subsamples from data frame the rows using the assigned index
  # column (time index, years) 
  return(X[i1:i2,])
}


remove_NA<-function(X){
  #removes every column in the data that contains a value of NA
  #X: data frame
  return(X[,colSums(is.na(X)) ==0])
}



# TODO write this function up correctly!!!
load_lonlat<-function(file_loc,island='BI',season='wet'){
  ###attach lat/long
  # 2019-06-08 OET: updated the file name to make it adjustable 
  # for the local system environment, see top of script 
  # for the variable definition
  # file_loc
  #read in the data
  
  # TODO, this part should be done right at the top when we select the stations for data analysis.
  # and simplify the process....
  locations<-read.csv(file=file_loc, header = TRUE, stringsAsFactors = FALSE)
  locations[,c(1,3:4)]
  
  
  #need wet and dry transposed
  WET<-read.csv(file=file_merge_wet, header = TRUE, stringsAsFactors = FALSE)
  DRY<-read.csv(file=file_merge_dry, header = TRUE, stringsAsFactors = FALSE)
  
  
  #merge with wet and dry
  locations_wet<-merge(locations, WET, by="SKN")
  locations_dry<-merge(locations, DRY, by="SKN")
  
  #keep only 1980 - 2007
  locations_wet2<-locations_wet[,c("SKN", "Lat_DD","Lon_DD","ElevFT","X1980","X1981","X1982","X1983","X1984","X1985","X1986","X1987","X1988","X1989", "X1990","X1991","X1992","X1993","X1994","X1995", "X1996", "X1997","X1998","X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006", "X2007")]
  locations_dry2<-locations_dry[,c("SKN", "Lat_DD","Lon_DD","ElevFT","X1980","X1981","X1982","X1983","X1984","X1985","X1986","X1987","X1988","X1989", "X1990","X1991","X1992","X1993","X1994","X1995", "X1996", "X1997","X1998","X1999", "X2000", "X2001", "X2002", "X2003", "X2004", "X2005", "X2006", "X2007")]
  
  #remove NAS
  locations_wet_nonas<-na.omit(locations_wet2)
  locations_dry_nonas<-na.omit(locations_dry2)
  
  #only keep location data
  locations_wet3<-locations_wet_nonas[,c("SKN","Lat_DD","Lon_DD","ElevFT")]
  locations_dry3<-locations_dry_nonas[,c("SKN","Lat_DD","Lon_DD","ElevFT")]
  
}


# NetCDF data functions 


load_climate_data<-function(localpath,filename,varname,dimnames=c('lon','lat','time'),
                            missing_val=9E9){
  ncfile<-nc_open(paste(localpath,filename,sep='/'))
  print(paste("load climate data from netcdf file ",filename,sep=''))
  print(paste("local directory: ",localpath,sep=''))
  print(paste("selected variable: ",varname,sep=''))
  variable<-ncvar_get(ncfile, varname)
  lon<-ncvar_get(ncfile, "lon")
  lat<-ncvar_get(ncfile, "lat")
  time<-ncvar_get(ncfile,'time')
  time_units<-(ncatt_get(ncfile,'time','units'))$value
  var_units<-(ncatt_get(ncfile,varname,'units'))$value
  ###set values to NA
  data<-replace(variable,(abs(variable)>missing_val),NA)
  ncdata<-list('data'=data,'lon'=lon,'lat'=lat,'time'=time,
               'time_units'=time_units,'var_units'=var_units)
  print(time_units)
  #nc_close(paste(localpath,filename,sep=''))
  return(ncdata)
}


select_climate_region<-function(ncdata, region='HI'){
  # input parameter:
  # ncdata: object of type list 
  # It must be consistent with the data returned from function load_climate_data
  # region: string identifier for the selected region
  # 
  # the index range must be adjusted for different climate data sets
  # that have other than the NCEP reanalysis 2.5 x 2.5 resolution
  # NOTE: latitude in ncep reanalysis are ordered from north to south 
  if (region=='HI') {
    iwest<-73
    ieast<-97
    inorth<-21
    isouth<-41
  }
  if (region=='NPAC') {
    iwest<-53
    ieast<-117
    inorth<-11
    isouth<-41
  }
  ncdata_out<-list("data"=ncdata$data[iwest:ieast,inorth:isouth,],
                   'lon'=ncdata$lon[iwest:ieast],
                   'lat'=ncdata$lat[inorth:isouth],
                   'time'=ncdata$time,
                   'var_units'=ncdata$var_units,'time_units'=ncdata$time_units)
  return(ncdata_out)
}



# create seasonal mean from monthly mean input 

clim_to_seasonal<-function(X,season='wet'){
  # function calculates the wet or dry season from monthly mean data
  # input parameter:
  # X: 3-dim data array with lon,lat,time dimensions
  # Returns:
  # S: 3-dim data array with lon,lat,time dimension (time dim reduced to length iyr)
  # Important: 
  # (1) The first time step in the netCDF data is the month January.
  # (2) data dimensions must be lon,lat,time
  nlon=dim(X)[1]
  nlat=dim(X)[2]
  ntime=dim(X)[3]
  iyr<-floor(ntime/12)-1 
  S<-array(0,dim=c(nlon,nlat,iyr))
  if (season=='wet'){
    #wet season (Nov - Apr)
    for (a in 1:iyr) {
      S[,,a]<-rowMeans(X[,,(12*(a-1)+11):(12*a+4)],dims=2)
    }
  }
  if (season=='dry'){
    #dry season (May - Oct)
    for (a in 1:iyr) {
      S[,,a]<-rowMeans(X[,,(12*a+5):(12*a+10)],dims=2)
    }
  }
  print(paste("calculated seasonal mean from monthly mean input"))
  print(paste("Season: ",season,sep=''))
  print(paste("Number of seasons: ",iyr),sep='')
  print(paste("Dimension of returned array (nlon,nlat,ntime): ",nlon,nlat,iyr))
  return(S)
}

composite<-function(y,x,xcrit,typ='>',outfile=NULL,is.sample=FALSE) {
  dim.y<-dim(y)
  inum<-length(x)
  if(dim.y[2]!=inum) {
    print("Error: dimension mismatch in matrix row dimension and vector x length")} 
  else{
    if (typ=='>') {
      print("composite analysis for y conditional upon x > xcrit")
      iscomp<-(x>xcrit)
    }
    if (typ=='<') {
      print("composite analysis for y conditional upon x < xcrit")
      iscomp<-(x<xcrit)
    }
    if(sum(iscomp,na.rm=TRUE)>0) {
      ncomp<-sum(iscomp)
      buffer<-y[,iscomp] 
      ymean<-apply(buffer,1,mean)
      ymedian<-apply(buffer,1,median)
      ysd<-apply(buffer,1,sd)
      ymin<-apply(buffer,1,min)
      ymax<-apply(buffer,1,max)
      if (is.sample) {
        print("ok")
        res<-list(comp.mean=ymean,comp.median=ymedian,comp.sd=ysd,comp.min=ymin,comp.max=ymax,n=ncomp,fraction=ncomp/inum,sample=buffer) 
      } else {
        res<-list(comp.mean=ymean,comp.median=ymedian,comp.sd=ysd,comp.min=ymin,comp.max=ymax,n=ncomp,fraction=ncomp/inum)
      }    
    } else {
      print("composite criteria never fulfilled. Empty result")
      res<-(-1)
    }
  }
  
  if (!is.null(outfile)) {
    
    write.table(round(cbind(res$comp.mean,res$comp.median,res$comp.sd,res$comp.min,res$comp.max),4),file=outfile)
    print(paste("wrote results to file ",outfile,sep=''))
  }
  return (res)
}



my.proj<- function(field=NA,pattern=NA) {
  # normalize p
  # History:
  
  #20190515 The mask file was not implemented yet
  #         and I removed the code related to this
  #         (This proj-function was used in connection
  #         with statistical downscaling without 
  #         a mask array in the past.
  #         The experimental version (proj2exp.R)
  #         has the mask array included in the calculations.)
  #
  #20130101 added a mask file to exclude
  #         certain areas in projection
  #         (i.e. the pattern that is interpreted as
  #          a unit-length vector is zero in
  #          the masked zones (mask is 0 or 1)
  #20110601 limit to valid data points      
  x1<-as.vector(field)
  p1<-as.vector(pattern)
  isvalidx1<-!is.na(x1)
  isvalidp1<-!is.na(p1)
  isvalid<-(isvalidx1 & isvalidp1)
  if (sum(isvalid)>0){
    pbuffer<-p1[isvalid]
    xbuffer<-x1[isvalid]
    pnorm<-sqrt(t(pbuffer)%*%pbuffer)
    pbuffer<-pbuffer/pnorm
    result<-(t(xbuffer)%*%pbuffer)
  } else {
    result<-NA
  }
  print (paste("proj(): check mag of pattern",pnorm,sum(isvalidx1),sum(isvalidp1)))
  result
}



###############################################################
# Main part of the script
###############################################################

season<-'dry'
start_year=1980
end_year=2007


# load the rainfall station data

#if (season=='wet'){
# file_rain=local.file_rain_wet
#}
#if (season=='dry'){
# file_rain=local.file_rain_dry
#}

#Xall<-load_rain_data(paste(local.data_path,local.dir_rainfall,sep=''),file_rain)

#load rainfall station data with my data
if (season=='wet'){
  file_rain="RF_wet_transposed.xlsx.csv"
}
if (season=='dry'){
  file_rain="RF_dry_transposed.csv"
}

Xall<-load_rain_data(paste("C:/Users/ksanf/Documents/SD/rainfall_data"),file_rain)




# get the time index range for the selected years

i1<-which(Xall['year']==start_year)[1]
i2<-which(Xall['year']==end_year)[1]

# select rainfall by years 
X<-subsample_index(Xall,i1,i2)

#Remove NAs
X_nonas<-X[ , colSums(is.na(X)) == 0]

#####make data frames matrices ##remove year column first

ncol=dim(X_nonas)[2]
Xmatrix<-as.matrix(X_nonas[,2:ncol])


######################################################################
#run PCA 
######################################################################

pr<-prcomp(Xmatrix,center=T,scale=T)


# eigenvectors of the PCA (weight factors for the stations) 
eigenvectors<-pr$rotation

notScaled=eigenvectors[,1:10]

write.csv(notScaled, file = "C:/Users/ksanf/Documents/SD/rainfall_data/eigenvectors_dry.csv")

######################################################################
#OET 2020-01-09 
#OET use function scale 
######################################################################

Xano<-scale(Xmatrix,center=TRUE,scale=TRUE)


# for the reconstruction of precipitation anomalies 
# we need later the rescaling and re-centering of the
# 

pc<-Xano%*%eigenvectors  # columns have for each eigenvector the PC time series
# now use full eigenvectors for testing reconstruction 
# below to make sure we can 100% reconstruct Xano
write.csv(pc, file = "C:/Users/ksanf/Documents/SD/rainfall_data/pc_dry.csv")


test<-pc%*%t(eigenvectors)

######################################################################
# Part two: composite analysis of NetCDF climate data
######################################################################


# load the 3-d field data to form a composite analysis


##open NetCDF file
#ncdir<-paste(local.data_path,local.dir_reanalysis,'/',sep='/')
ncdir<-paste("C:/Users/ksanf/Documents/SD/Predictor_Variables")

##ALL VARIABLE LIST
ncfilelist<-c("omega500.mon.mean.nc", "shum_x_vwnd.700.mon.mean.nc", "hgt1000.mon.mean.nc", "slp.mon.mean.nc", "air.1000-500.mon.mean.nc", "shum_x_vwnd.925.mon.mean.nc", "pottmp.1000-500.mon.mean.nc", "shum_x_uwnd.925.mon.mean.nc", "hgt500.mon.mean.nc", "shum_x_uwnd.700.mon.mean.nc", "pottmp.1000-850.mon.mean.nc", "skt.mon.mean.regridded.nc", "air.2m.mon.mean.regridded.nc", "pwtr.mon.mean.nc", "shum925.mon.mean.nc", "shum700.mon.mean.nc")
ncvarlist<-c("omega", "shum", "hgt", "slp", "air", "shum", "pottmp", "shum", "hgt", "shum", "pottmp", "skt", "air", "pr_wtr", "shum", "shum")
labelslist<-c("omega", "merid_transport700", "hgt1000", "slp", "airdiff", "merid_trans925", "pottmp1000_500", "zonal_trans_925", "hgt500", "zonal_trans_700", "pottmp1000_850", "skt", "air2m", "pr_wtr", "shum925", "shum700")


# for finding the right years in the seasonal data set I need the first for 
# of the first seasonal data (calculated here below from monthly mean NetCDF data)
nc_start_wet<-1949 # for wet season
nc_start_dry<-1948 # for dry season


# region that is selected from the (global) NetCDF data
region<-'HI'

# PC time series used to form composites and used as predictand in regression model
pcamode=4

# choose the number of lowest ncomp and highest ncomp years to be used in composites
ncomp<-4


##############################################################
# get the first and last year from the rainfall data  
# to find the matching years in the climate data
# (could be improved to make it less error prone)
##############################################################
hi_rain_start<-X$year[1]
hi_rain_end<-X$year[length(X$year)]
# get time index for the netcdf data to use 
if (season=='wet') {
  nc_start<-nc_start_wet
}
if (season=='dry') {
  nc_start<-nc_start_dry
}
itime1<-hi_rain_start-nc_start+1
itime2<-hi_rain_end-nc_start+1
iselect<-itime1:itime2

##############################################################
# predictand time series
##############################################################
y<-pc[,pcamode]
iyrs<-length(y)
##############################################################
# array to store the projection index time series
# for all climate variables (organized in columns)
##############################################################

iclimvar<-length(ncvarlist)
xproj<-array(0,dim=c(iyrs,iclimvar))

##############################################################
# loop over climate variables
##############################################################
library(sp)
library(rgdal)
hawaii <- readOGR("C:/Users/ksanf/Documents/Statistical_DS/state_map_info/dems/coast", "coast_geo_shp")


for (ivar in 1:iclimvar){
  ncfile<-ncfilelist[ivar]
  ncvar<-ncvarlist[ivar]
  labels<-labelslist[ivar]
  clim_all<-load_climate_data(ncdir,ncfile,ncvar)
  clim_reg<-select_climate_region(clim_all,region=region)
  xclim<-clim_to_seasonal(clim_reg$data,season=season)
  
  # 3-dim array with lon,lat, time dimensions
  # convert into 2dim with (lon*lat, time dimensions for use 
  # with composite analysis function
  dim3d<-dim(xclim[,,iselect])
  xclim2<-array(xclim[,,iselect],dim=c(dim3d[1]*dim3d[2],dim3d[3]))
  
  ###################################################################
  # Composite analysis & regression model
  ###################################################################
  iyrs<-length(y)
  result<-sort(y, index.return=TRUE)
  # since we want to use the composite function
  # I have to create a rank index series
  # this will find the right years to form composites
  rank<-1:length(y) 
  rank[result$ix]<-rank
  low<-composite(xclim2,rank,ncomp+1,typ='<')
  high<-composite(xclim2,rank,iyrs-ncomp,typ='>')
  # reshape composites back to 2dim fields for plotting
  comp<-array(high$comp.mean-low$comp.mean,dim=c(dim3d[1],dim3d[2]))
  
  
  ##########################################################
  # Plotting the composite pattern
  ##########################################################
  # change lats from north to south to south to north order
  lon2<-clim_reg$lon
  lat2rev<-rev(clim_reg$lat)
  irev<-rev(1:length(lat2rev))
  filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode, sep=""),plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})
  
  imgpath2=file.path("C:", "Users","ksanf","Documents","SD","Composite_Plots", paste(labels, "_PC",pcamode,"_", season,"Comp_Plot.jpg",sep=""))
  jpeg(imgpath2, width = 500, height = 350)
  filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode, ",", season, sep=""),plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})

  
  p<-filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode,", ", season, sep=""))
  plot(hawaii,add=T)
  print(p)
  dev.off()
  
  #######################################################
  # Projection index calculation
  #######################################################
  for (i in 1:iyrs) {
    xproj[i,ivar]<-my.proj(xclim2[,i],comp)
  }
}


# y is the eigenvectors
# save dataframe for each of the 4 pcs


df<-data.frame('y'=y,'x'=xproj)
colnames(df)[2:17]<-labelslist

write.csv(df, "C:/Users/ksanf/Documents/SD/DFs/dry_PC4.csv", row.names=FALSE)


####################################################################
#review variable rank and multicollinearity
#################################################################

##loop for making univariable models
#one rsq table for pc and season
rsq<-data.frame(matrix(NA, nrow=length(ncvarlist), ncol=2))
colnames(rsq)[1]<-"Variable"
colnames(rsq)[2]<-paste("Rsq_PC", pcamode,"_",season,sep="")


for(k in 1:length(ncvarlist)){
  labels<-labelslist[k]
  
  rsq[k,1]<-labels
  mod1<-lm(df$y~df[[k+1]])
  rsq[k,2]<-summary(mod1)$r.squared
  
  imgpath<-file.path("C:", "Users","ksanf","Documents","SD","Residual_Plots", paste(labels, "_PC",pcamode,"_", season,"Res_Plot.png",sep=""))
  dpi=300
  png(file=imgpath,width=6.5*dpi,height=6.5*dpi,res=dpi)
  par(mfrow=c(2,2))
  p<-plot(mod1, main=paste("LM_PC", pcamode,season, labels,sep="_"))
  print(p)
  dev.off()
  
}

rsq<-rsq[order(rsq[[2]], decreasing=T),]
rsq
setwd("C:/Users/ksanf/Documents/SD/RSQ_Tables")
write.table(rsq,file=paste0("Rsq_PC",pcamode, season,".csv"),sep=",",col.names= T, row.names= F)



######################################

#Making Correlation Matrices
######################################
##new code
dev.off()
library(corrplot)
newdf<-subset(df, select = -c(y))
df.corr<-cor(newdf)
df.corr<-df.corr*df.corr

#The default value for mar is c(5.1, 4.1, 4.1, 2.1)
par(xpd=NA,oma=c(0,0,3,0))

corrplot(df.corr,method="shade", type="lower", tl.cex=.7,tl.offset=.1,cl.lim=c(0,1),col=colorRampPalette(c("blue", "white", "red"))(500),tl.col="black", number.cex = .5, addCoef.col = "black")
title(paste(season, "season","PC",pcamode, sep="_"), adj=.75, line=2)
###manually save corr plot here


#############################################################################
###make dredge tables  ##do this 1 pc at a time
#make model with all remaining predictors after eliminating based on multicollinearity
#wet season pc 1
res<-glm(y~omega + merid_transport700 + hgt1000 + merid_trans925 + zonal_trans_925 + zonal_trans_700  + pottmp1000_850 + air2m, df, family=gaussian, na.action = "na.fail")
#wet season pc 2
res<-glm(y~omega +merid_transport700 + airdiff +merid_trans925 +zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + skt +pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
#wet season pc 3
res<-glm(y~ omega + merid_transport700 + airdiff + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
#wet season pc4
res<-glm(y~omega+ merid_transport700 + +hgt1000+ airdiff+ merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + pr_wtr+ shum925 +shum700, df, family=gaussian, na.action = "na.fail")


#dry season pc1
res<-glm(y~omega + merid_transport700 + slp + merid_trans925 + pottmp1000_500 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 +air2m + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
#dry season pc2
res<-glm(y~omega + merid_transport700 +hgt1000 + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 + skt + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
#dry season pc3
res<-glm(y~omega + merid_transport700 + slp + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
#dry season pc4
res<-glm(y~omega + merid_transport700 + hgt1000 + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 + air2m  + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")


####################################################################################
#all possible models
dredge=dredge(res, rank = "AICc",extra="R^2")
write.csv(dredge,"C:/Users/ksanf/Documents/SD/Dredge_Tables/dry_PC4_dredge.csv", row.names = FALSE)


###################################################################################
#Model Combinations

###################################################################################
##clear everything
rm(list=ls())

setwd("C:/Users/ksanf/Documents/SD/Dredge_Tables")

dredge_W1<-read.csv("wet_PC1_dredge.csv")
dredge_W2<-read.csv("wet_PC2_dredge.csv")
dredge_W3<-read.csv("wet_PC3_dredge.csv")
dredge_W4<-read.csv("wet_PC4_dredge.csv")

dredge.list<-list(dredge_W1, dredge_W2, dredge_W3, dredge_W4)


setwd("C:/Users/ksanf/Documents/SD/DFs")
wet_df1<-read.csv("wet_PC1.csv")
wet_df2<-read.csv("wet_PC2.csv")
wet_df3<-read.csv("wet_PC3.csv")
wet_df4<-read.csv("wet_PC4.csv")


df.list<-list(wet_df1, wet_df2, wet_df3, wet_df4)

##shorten dredge lists

dredge.list<-lapply(dredge.list, function(x) subset(x, subset=delta<2))
dredge.list<-lapply(dredge.list, function(x) x[,2:(ncol(x)-6)])

nrow(dredge_W4)
nrow(dredge.list[[1]])
df.list[[2]]



##LOOP 1: sets coefficients = variable name
for (j in 1:length(dredge.list)){
  for (i in 1:ncol(dredge.list[[j]])){
    dredge.list[[j]][!is.na(dredge.list[[j]][,i]), i ]<-names(dredge.list[[j]])[i]
  }
}

library(stringr)

#LOOP 2: makes column with variable list
for (j in 1:length(dredge.list)){
  dredge.list[[j]]$varlist <- apply(dredge.list[[j]], 1, function(x) 
    paste(str_trim(x[!is.na(x)]), collapse=" + "))
}

varlist1<-dredge.list[[1]]$varlist
varlist2<-dredge.list[[2]]$varlist
varlist3<-dredge.list[[3]]$varlist
varlist4<-dredge.list[[4]]$varlist

write.csv(varlist1, file="C:/Users/ksanf/Documents/SD/VarLists_wet/varlist_d1.csv", row.names = FALSE)
write.csv(varlist2, file="C:/Users/ksanf/Documents/SD/Varlists_wet/varlist_d2.csv", row.names = FALSE)
write.csv(varlist3, file="C:/Users/ksanf/Documents/SD/Varlists_wet/varlist_d3.csv", row.names = FALSE)
write.csv(varlist4, file="C:/Users/ksanf/Documents/SD/Varlists_wet/varlist_d4.csv", row.names = FALSE)



#LOOP 3: Makes model for every variable list of every PC; PC dredge table and df must match; Output of loop = tables of model fits

fits.list=data.frame(matrix(ncol=1, nrow=28))
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
  fits.list<-cbind(fits.list, fits_df)
  write.csv(fits_df, file=paste0("C:/Users/ksanf/Documents/SD/Fits_for_Reconstruction/wet_fits_PC", k, ".csv"))
}

#write.csv(fits.list, file="C:/Users/ksanf/Documents/SD/Fits_for_Reconstruction/fits_list_wet.csv", row.names = FALSE)


###########################################################

#reconstruction code

######################################################################
#OET 2020-01-09 
#OET use function scale 
######################################################################
##clear everything
rm(list=ls())

##need to get these from beginning of code "xmatrix"
Xmean<-apply(Xmatrix,2,mean)

Xvar<-apply(Xmatrix,2,var)


setwd("C:/Users/ksanf/Documents/SD/Fits_for_Reconstruction")
fits_PC1=read.csv("wet_fits_PC1.csv")
fits_PC1=fits_PC1[2:ncol(fits_PC1)]

fits_PC2=read.csv("wet_fits_PC2.csv")
fits_PC2=fits_PC2[2:ncol(fits_PC2)]

fits_PC3=read.csv("wet_fits_PC3.csv")
fits_PC3=fits_PC3[2:ncol(fits_PC3)]

fits_PC4=read.csv("wet_fits_PC4.csv")
fits_PC4=fits_PC4[2:ncol(fits_PC4)]


##loop for making pc fit with every combination
ntrunc=4
Xanotrunc<-pc[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])

for (a in 1:ncol(fits_PC1)){
  for (b in 1:ncol(fits_PC2)){
    for (c in 1:ncol(fits_PC3)){
      for (d in 1:ncol(fits_PC4)){
        pcfit<-cbind(fits_PC1[a],fits_PC2[b],fits_PC3[c],fits_PC4[d])
        pcfit=as.matrix(pcfit)
        write.csv(pcfit, file=paste0("C:/Users/ksanf/Documents/SD/Fits_All_Combos_Wet/fits", a,b,c,d, ".csv"))
        Xanofit<-pcfit[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
        Xfit  <-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        Xtrunc<-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        i=1
        while (i<=length(Xvar)) {
          Xfit[,i]<-Xanofit[,i]*sqrt(Xvar[i])+Xmean[i]
          Xtrunc[,i]<-Xanotrunc[,i]*sqrt(Xvar[i])+Xmean[i]
          i=i+1
        }
        write.csv(Xtrunc, file=paste0("C:/Users/ksanf/Documents/SD/Trunc_Precip_Wet/xtrunc", a,b,c,d, ".csv"))
        write.csv(Xfit, file=paste0("C:/Users/ksanf/Documents/SD/Fitted_Precip_All_Combos_Wet/xfit", a,b,c,d, ".csv"))
        
      }
    }
  }
}


plot(Xmatrix[,17], type = "line", col= "blue", main="Station 17, All Years")
lines(Xtrunc[,17], type="line", col="red")
lines(Xfit[,17], type = "line", col="green")
legend(20, 600, legend=c("actual", "truncated", "projected"),
       col=c("blue", "red",  "green"), lty=1, cex=0.8)


#########################################################################
#DRY SEASON

#Making model Combinations - dry season

##clear everything
rm(list=ls())

setwd("C:/Users/ksanf/Documents/SD/Dredge_Tables")

dredge_W1<-read.csv("dry_PC1_dredge.csv")
dredge_W2<-read.csv("dry_PC2_dredge.csv")
dredge_W3<-read.csv("dry_PC3_dredge.csv")
dredge_W4<-read.csv("dry_PC4_dredge.csv")

dredge.list<-list(dredge_W1, dredge_W2, dredge_W3, dredge_W4)


setwd("C:/Users/ksanf/Documents/SD/DFs")
dry_df1<-read.csv("dry_PC1.csv")
dry_df2<-read.csv("dry_PC2.csv")
dry_df3<-read.csv("dry_PC3.csv")
dry_df4<-read.csv("dry_PC4.csv")


df.list<-list(dry_df1, dry_df2, dry_df3, dry_df4)

##shorten dredge lists

dredge.list<-lapply(dredge.list, function(x) subset(x, subset=delta<2))
dredge.list<-lapply(dredge.list, function(x) x[,2:(ncol(x)-6)])

nrow(df.list[[1]])

##LOOP 1: sets coefficients = variable name
for (j in 1:length(dredge.list)){
  for (i in 1:ncol(dredge.list[[j]])){
    dredge.list[[j]][!is.na(dredge.list[[j]][,i]), i ]<-names(dredge.list[[j]])[i]
  }
}

dredge.list[[1]]


library(stringr)

#LOOP 2: makes column with variable list
for (j in 1:length(dredge.list)){
  dredge.list[[j]]$varlist <- apply(dredge.list[[j]], 1, function(x) 
    paste(str_trim(x[!is.na(x)]), collapse=" + "))
}

varlist1<-dredge.list[[1]]$varlist
varlist2<-dredge.list[[2]]$varlist
varlist3<-dredge.list[[3]]$varlist
varlist4<-dredge.list[[4]]$varlist

write.csv(varlist1, file="C:/Users/ksanf/Documents/SD/Varlists_dry/varlist_d1.csv", row.names = FALSE)
write.csv(varlist2, file="C:/Users/ksanf/Documents/SD/Varlists_dry/varlist_d2.csv", row.names = FALSE)
write.csv(varlist3, file="C:/Users/ksanf/Documents/SD/Varlists_dry/varlist_d3.csv", row.names = FALSE)
write.csv(varlist4, file="C:/Users/ksanf/Documents/SD/Varlists_dry/varlist_d4.csv", row.names = FALSE)



#LOOP 3: Makes model for every variable list of every PC; PC dredge table and df must match; Output of loop = tables of model fits

fits.list=data.frame(matrix(ncol=1, nrow=28))
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
  fits.list<-cbind(fits.list, fits_df)
  write.csv(fits_df, file=paste0("C:/Users/ksanf/Documents/SD/Fits_for_Reconstruction/dry_fits_PC", k, ".csv"))
}

#write.csv(fits.list, file="C:/Users/ksanf/Documents/Statistical_DS/Fits_for_Reconstruction/fits_list_dry.csv", row.names = FALSE)

###########################################################

#reconstruction code

######################################################################
#OET 2020-01-09 
#OET use function scale 
######################################################################
##clear everything
rm(list=ls())

##Run PCA code for rainfall objects

##need to get these from beginning of code "xmatrix"
#get xmtarix of dry season from beginning of code

Xmean<-apply(Xmatrix,2,mean)

Xvar<-apply(Xmatrix,2,var)

setwd("C:/Users/ksanf/Documents/SD/Fits_for_Reconstruction")
fits_PC1=read.csv("dry_fits_PC1.csv")
fits_PC1=fits_PC1[2:ncol(fits_PC1)]

fits_PC2=read.csv("dry_fits_PC2.csv")
fits_PC2=fits_PC2[2:ncol(fits_PC2)]

fits_PC3=read.csv("dry_fits_PC3.csv")
fits_PC3=fits_PC3[2:ncol(fits_PC3)]

fits_PC4=read.csv("dry_fits_PC4.csv")
fits_PC4=fits_PC4[2:ncol(fits_PC4)]


##loop for making pc fit with every combination
ntrunc=4
Xanotrunc<-pc[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])

for (a in 1:ncol(fits_PC1)){
  for (b in 1:ncol(fits_PC2)){
    for (c in 1:ncol(fits_PC3)){
      for (d in 1:ncol(fits_PC4)){
        pcfit<-cbind(fits_PC1[a],fits_PC2[b],fits_PC3[c],fits_PC4[d])
        pcfit=as.matrix(pcfit)
        write.csv(pcfit, file=paste0("C:/Users/ksanf/Documents/SD/Fits_All_Combos_Dry/fits", a,b,c,d, ".csv"))
        Xanofit<-pcfit[,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
        Xfit  <-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        Xtrunc<-matrix(0, nrow(Xmatrix),ncol(Xmatrix))
        i=1
        while (i<=length(Xvar)) {
          Xfit[,i]<-Xanofit[,i]*sqrt(Xvar[i])+Xmean[i]
          Xtrunc[,i]<-Xanotrunc[,i]*sqrt(Xvar[i])+Xmean[i]
          i=i+1
        }
        write.csv(Xtrunc, file=paste0("C:/Users/ksanf/Documents/SD/Trunc_Precip_Dry/xtrunc", a,b,c,d, ".csv"))
        write.csv(Xfit, file=paste0("C:/Users/ksanf/Documents/SD/Fitted_Precip_All_Combos_Dry/xfit", a,b,c,d, ".csv"))
        
      }
    }
  }
}



plot(Xmatrix[,17], type = "line", col= "blue", main="Station 17, All Years")
lines(Xtrunc[,17], type="line", col="red")
lines(Xfit[,17], type = "line", col="green")
legend(20, 600, legend=c("actual", "truncated", "projected"),
       col=c("blue", "red",  "green"), lty=1, cex=0.8)

#######################
#LOOCV

##WET SEASON

##clear everything
rm(list=ls())

setwd("C:/Users/ksanf/Documents/SD/DFs")
wet_df1<-read.csv("wet_PC1.csv")
wet_df2<-read.csv("wet_PC2.csv")
wet_df3<-read.csv("wet_PC3.csv")
wet_df4<-read.csv("wet_PC4.csv")

df.list<-list(wet_df1, wet_df2, wet_df3, wet_df4)

#read in varlist for each PC
setwd("C:/Users/ksanf/Documents/SD/VarLists_wet")
varlist1<-read.csv("varlist_d1.csv")
varlist2<-read.csv("varlist_d2.csv")
varlist3<-read.csv("varlist_d3.csv")
varlist4<-read.csv("varlist_d4.csv")

var.lists<-list(varlist1, varlist2, varlist3, varlist4)


##LOOCV
#wet season
#do PCs one at a time
#PC1
for (a in 1:nrow(var.lists[[1]])){
  fit=NULL
  for (i in 1:nrow(df.list[[1]])){ #1:28 years
    datause = df.list[[1]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[1]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[1]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC1/PC1_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC2
for (a in 1:nrow(var.lists[[2]])){
  fit=NULL
  for (i in 1:nrow(df.list[[2]])){ #1:28 years
    datause = df.list[[2]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[2]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[2]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC2/PC2_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC3
for (a in 1:nrow(var.lists[[3]])){
  fit=NULL
  for (i in 1:nrow(df.list[[3]])){ #1:28 years
    datause = df.list[[3]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[3]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[3]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC3/PC3_predict_mod_", a, ".csv"), row.names=TRUE)
}


#PC4
for (a in 1:nrow(var.lists[[4]])){
  fit=NULL
  for (i in 1:nrow(df.list[[4]])){ #1:28 years
    datause = df.list[[4]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[4]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[4]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC4/PC4_predict_mod_", a, ".csv"), row.names=TRUE)
}



#####STEP 2: USE OUTPUT FROM STEP 1 TO RECONSTRUCT RAINFALL AND CALCULATE STATS


##Calculating RMSE
#m is for model (fitted) values, o is for observed (true) values.

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}


##calculating MAE
#actual and predicted are vectors
library(ModelMetrics)


#write.csv(Xmatrix, file="C:/Users/ksanf/Documents/Statistical_DS/actual_DRY.csv")

#mae(Xmatrix, Xfit_all)



##read in fits as a list?
#PC1
setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC1")

LOOCV_fits1 = list.files(pattern="*.csv")
LOOCV_fits1 = lapply(LOOCV_fits1, read.csv)

#PC2
setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC2")
LOOCV_fits2 = list.files(pattern="*.csv")
LOOCV_fits2 = lapply(LOOCV_fits2, read.csv)


#PC#
setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC3")
LOOCV_fits3 = list.files(pattern="*.csv")
LOOCV_fits3 = lapply(LOOCV_fits3, read.csv)


#PC4
setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/Predict_PC4")
LOOCV_fits4 = list.files(pattern="*.csv")
LOOCV_fits4 = lapply(LOOCV_fits4, read.csv)



#reconstruction loop
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
        write.csv(Xfit_all, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/Reconstruction/Reconstruction_", a, "_", b, "_", c, "_", d, ".csv"), row.names=TRUE)
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
        write.csv(r2, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_WET/LOOCV_R2/R2_", a, "_", b, "_", c, "_", d, ".csv"), row.names=F)
        
        Xfit_all=NULL
        
        j=j+1
        vals[j,]=rmse
        vals2[j,]=MAE
        colnames(vals)<-c("Combination", "RMSE")
        colnames(vals2)<-c("Combination", "MAE")
        write.csv(vals, file="C:/Users/ksanf/Documents/SD/LOOCV_WET/RMSE_MAE/rmse.csv", row.names=FALSE)
        write.csv(vals2, file="C:/Users/ksanf/Documents/SD/LOOCV_WET/RMSE_MAE/mae.csv", row.names=FALSE)
        
      }
    }
  }
}


#set wd
rm(list=ls())
setwd("C:/Users/ksanf/Documents/SD/LOOCV_WET/LOOCV_R2")

multmerge=function(mypath) {
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist= lapply(filenames, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))
  Reduce(function(x,y) merge(x,y), datalist)}


mymergeddata = multmerge("C:/Users/ksanf/Documents/SD/LOOCV_WET/LOOCV_R2")


write.csv(mymergeddata, file="C:/Users/ksanf/Documents/SD/LOOCV_WET/merged_r2.csv", row.names=FALSE)

mymergeddata=read.csv("C:/Users/ksanf/Documents/SD/LOOCV_WET/merged_r2.csv", check.names= FALSE)

mymergeddata$median=apply(mymergeddata[,2:865],1, median, na.rm = TRUE)

wet_station_r2=mymergeddata[,c(1,866)]

#attach location data for mapping
#dry_r2=cbind(locations_dry3[,2:3], dry_station_r2)
write.csv(wet_station_r2, file= "C:/Users/ksanf/Documents/SD/LOOCV_WET/wet_r2_per_station.csv", row.names = FALSE)


plot(mymergeddata$median, xlab = "Station", ylab = "R^2", main = "Median R^2 Value per Station")
hist(mymergeddata$median, xlab = "r-squared", main= "Median R^2 Value per Station", xlim=c(0,1), breaks=10)



####################################################################################
rmses=read.csv("C:/Users/ksanf/Documents/SD/LOOCV_WET/RMSE_MAE/rmse.csv")

plot(rmses$RMSE, main= "LOOCV RMSE Values", ylab= "RMSE")
abline(h=rmse2015, col= "blue", lty=2)

#find rmse2015 calculation in script "2015_LOOCV_code"


hist(rmses$RMSE, xlab= "RMSE", main = "LOOCV RMSE Values")
points(rmse2015, 0, pch=16,  col= "blue", cex=2)



mae=read.csv("C:/Users/ksanf/Documents/Statistical_DS/LOOCV_DRY/RMSE_MAE/mae.csv")

plot(mae$MAE, ylim=c(126,143), main= "LOOCV MAE Values")
abline(h=mae2015, lty=2, col="blue")





####2015 vals

write.csv(Xmatrix, "C:/Users/ksanf/Documents/SD/Actual_Rain/actual_wet.csv")




#####################################################################3
#DRY SEASON LOOCV

##clear everything
rm(list=ls())


setwd("C:/Users/ksanf/Documents/SD/DFs")
dry_df1<-read.csv("dry_PC1.csv")
dry_df2<-read.csv("dry_PC2.csv")
dry_df3<-read.csv("dry_PC3.csv")
dry_df4<-read.csv("dry_PC4.csv")

df.list<-list(dry_df1, dry_df2, dry_df3, dry_df4)

#read in varlist for each PC
setwd("C:/Users/ksanf/Documents/SD/VarLists_dry")
varlist1<-read.csv("varlist_d1.csv")
varlist2<-read.csv("varlist_d2.csv")
varlist3<-read.csv("varlist_d3.csv")
varlist4<-read.csv("varlist_d4.csv")

var.lists<-list(varlist1, varlist2, varlist3, varlist4)


#do PCs one at a time
#PC1
for (a in 1:nrow(var.lists[[1]])){
  fit=NULL
  for (i in 1:nrow(df.list[[1]])){ #1:28 years
    datause = df.list[[1]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[1]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[1]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC1/PC1_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC2
for (a in 1:nrow(var.lists[[2]])){
  fit=NULL
  for (i in 1:nrow(df.list[[2]])){ #1:28 years
    datause = df.list[[2]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[2]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[2]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC2/PC2_predict_mod_", a, ".csv"), row.names=TRUE)
}

#PC3
for (a in 1:nrow(var.lists[[3]])){
  fit=NULL
  for (i in 1:nrow(df.list[[3]])){ #1:28 years
    datause = df.list[[3]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[3]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[3]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC3/PC3_predict_mod_", a, ".csv"), row.names=TRUE)
}


#PC4
for (a in 1:nrow(var.lists[[4]])){
  fit=NULL
  for (i in 1:nrow(df.list[[4]])){ #1:28 years
    datause = df.list[[4]][-i,] #remove 1 year at a time
    model<-lm(paste("y ~", var.lists[[4]][a,]), data=datause) #27 years 
    y=predict(model, df.list[[4]][i,])
    fit <- rbind(fit, y)
  }
  row.names(fit) <- c(1:28)
  
  write.csv(fit, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC4/PC4_predict_mod_", a, ".csv"), row.names=TRUE)
}



#####STEP 2: USE OUTPUT FROM STEP 1 TO RECONSTRUCT RAINFALL AND CALCULATE STATS


##Calculating RMSE
#m is for model (fitted) values, o is for observed (true) values.

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}


##calculating MAE
#actual and predicted are vectors
library(ModelMetrics)


#write.csv(Xmatrix, file="C:/Users/ksanf/Documents/Statistical_DS/actual_DRY.csv")

#mae(Xmatrix, Xfit_all)



##read in fits as a list?
#PC1
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC1")

LOOCV_fits1 = list.files(pattern="*.csv")
LOOCV_fits1 = lapply(LOOCV_fits1, read.csv)

#PC2
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC2")
LOOCV_fits2 = list.files(pattern="*.csv")
LOOCV_fits2 = lapply(LOOCV_fits2, read.csv)


#PC#
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC3")
LOOCV_fits3 = list.files(pattern="*.csv")
LOOCV_fits3 = lapply(LOOCV_fits3, read.csv)


#PC4
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Predict_PC4")
LOOCV_fits4 = list.files(pattern="*.csv")
LOOCV_fits4 = lapply(LOOCV_fits4, read.csv)



#reconstruction loop

#need xmatrix and eigenvectors for dry season from code above

ntrunc=4
Xanotrunc<-pc[1,1:ntrunc]%*%t(eigenvectors[,1:ntrunc])
Xmean<-apply(Xmatrix,2,mean)
Xvar<-apply(Xmatrix,2,var)
Xfit_all=NULL
vals=data.frame(matrix(nrow=700, ncol=2))
vals2=data.frame(matrix(nrow=700, ncol=2))
r2=data.frame(matrix(nrow=851, ncol=2))
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
        write.csv(Xfit_all, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/Reconstruction/Reconstruction_", a, "_", b, "_", c, "_", d, ".csv"), row.names=TRUE)
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
        write.csv(r2, file=paste0("C:/Users/ksanf/Documents/SD/LOOCV_DRY/LOOCV_R2/R2_", a, "_", b, "_", c, "_", d, ".csv"), row.names=F)
        
        Xfit_all=NULL
        
        j=j+1
        vals[j,]=rmse
        vals2[j,]=MAE
        colnames(vals)<-c("Combination", "RMSE")
        colnames(vals2)<-c("Combination", "MAE")
        write.csv(vals, file="C:/Users/ksanf/Documents/SD/LOOCV_DRY/RMSE_MAE/rmse.csv", row.names=FALSE)
        write.csv(vals2, file="C:/Users/ksanf/Documents/SD/LOOCV_DRY/RMSE_MAE/mae.csv", row.names=FALSE)
        
      }
    }
  }
}


#set wd
rm(list=ls())
setwd("C:/Users/ksanf/Documents/SD/LOOCV_DRY/LOOCV_R2")

multmerge=function(mypath) {
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist= lapply(filenames, function (x) read.csv(file=x, header=TRUE, check.names=FALSE))
  Reduce(function(x,y) merge(x,y), datalist)}


mymergeddata = multmerge("C:/Users/ksanf/Documents/SD/LOOCV_DRY/LOOCV_R2")


write.csv(mymergeddata, file="C:/Users/ksanf/Documents/SD/LOOCV_DRY/merged_r2.csv", row.names=FALSE)

mymergeddata=read.csv("C:/Users/ksanf/Documents/SD/LOOCV_DRY/merged_r2.csv", check.names = FALSE)

mymergeddata$median=apply(mymergeddata[,2:701],1, median, na.rm = TRUE)

dry_station_r2=mymergeddata[,c(1,702)]

#attach location data for mapping
#dry_r2=cbind(locations_dry3[,2:3], dry_station_r2)
write.csv(dry_station_r2, file= "C:/Users/ksanf/Documents/SD/LOOCV_DRY/dry_r2_per_station.csv", row.names = FALSE)

plot(mymergeddata$median, xlab = "Station", ylab = "R^2", main = "Median R^2 Value per Station")
hist(mymergeddata$median, xlab = "r-squared", main= "Median R^2 Value per Station", xlim=c(0,1), breaks=10)


rmses=read.csv("C:/Users/ksanf/Documents/Statistical_DS/LOOCV_DRY/RMSE_MAE/rmse.csv")

plot(rmses$RMSE, main= "LOOCV RMSE Values", ylab= "RMSE", ylim=c(195,227))
abline(h=rmse2015, col= "blue", lty=2)

#find rmse2015 calculation in script "2015_LOOCV_code"


hist(rmses$RMSE, xlab= "RMSE", main = "LOOCV RMSE Values", xlim=c(195,230))
points(rmse2015, 0, pch=16,  col= "blue", cex=2)



mae=read.csv("C:/Users/ksanf/Documents/Statistical_DS/LOOCV_DRY/RMSE_MAE/mae.csv")

plot(mae$MAE, ylim=c(126,143), main= "LOOCV MAE Values")
abline(h=mae2015, lty=2, col="blue")





####2015 vals

write.csv(Xmatrix, "C:/Users/ksanf/Documents/SD/Actual_Rain/actual_wet.csv")

