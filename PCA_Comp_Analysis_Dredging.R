#this code includes: PCA on rainfall, composite analysis,projection index calculation, 
#and saves a df of predictand and predictor info each time the loop is run.
#to test multicollinearity, use this script in conjunction with script titled
#"Multicollinearity_and_Variable_Ranking". 
#To skip this step, continue on to "Dredging" incuded in this script

##clear everything
rm(list=ls())
#######################################################
# (1) import package libraries and define functions
#######################################################

# (1.1) load packages
require(ncdf4)
require(Matrix)
library(car)
library(stats)
library(arm)
library(MuMIn)
library(boot)
library(sp)
library(rgdal)


################################################################
# (1.2) create functions for data processing
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
  # OET use file.path to make it OS independent
  filename<-file.path(datapath,filename)
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
  #OET 2020-12-17 use file.path to make OS independent
  ncfile<-nc_open(file.path(localpath,filename))
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
  # 2020-12-17 OET: This function ONLY WORKS with 
  # NCEP2 reanalysis data with 2.5 by 2.5 degree
  # resolution and standard lat,long order
  # TODO: develop a better generalized method to 
  # find the right index ranges for extracting HI region

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
  # IMPORTANT: Function assumes monthly mean climate data starting
  # with January!!!
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
  # This functions takes vector x checks if values in x pass a threshold
  # xcrit and then selects the corresponding values from vector y
  # to calculate statistics.
  # vector y is a second vector with elements selected conditioned on
  # the check where in vector x the threshold is passed.
  # e.g. it could be a vector with dates, or an index counter.
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
  #20201217 OET: TODO check the warning message and 
  #         update the code accordingly
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
    #20201222 updated code, both are vectors
    # inner product of two vectors is returned as 
    # matrix (1x1), drop() removes dimension
    pnorm<-sqrt(drop(pbuffer%*%pbuffer))
    pbuffer<-pbuffer/pnorm
    result<-drop(xbuffer%*%pbuffer)
  } else {
    result<-NA
  }
  print (paste("proj(): check mag of pattern",pnorm,sum(isvalidx1),sum(isvalidp1)))
  result
}



###############################################################
# (2) Main part of the script
###############################################################

###############################################################
# (2.1) User defined parameters
###############################################################

 

###############################################################
# (2.1.1) User-controlled selections and options 
###############################################################
season<-'dry'
# years used for the calculations
start_year=1980
end_year=2007


###############################################################
# (2.1.2) User defined path names
###############################################################

# OET 2020-12-17 Use of file.path() 
# to make path names independent of the OS


datapath=file.path('C:','Users','timm','Documents','Github',
                   'Statistical-Downscaling_KS','Data')
shpdir=datapath # shape file(s)

outputdir=file.path('C:','Users','timm','Documents','R',
                    'SDHI_KS_test')
outputdir.train=file.path(outputdir,'results','training')
outputdir.figs= file.path(outputdir,'figures','training')


# NetCDF files (climate predictors)
#ncdir<-paste(local.data_path,local.dir_reanalysis,'/',sep='/')
#ncdir=file.path(datapath,'reanalysis')
ncdir<-datapath


###############################################################
# (2.1.3) Netcdf climate data (predictors)
###############################################################

##ALL VARIABLE LIST
#build a list of all possible variables to be evaluated
# 2020-12-17 OET: commented out the files that are currently missing
# 2020-12-17 OET: TODO: a better way would de be to put this into a csv
# file with netcdf file name, ncvar name, label name, and additional column 
# to tell script in the dredge part which variables to include/ exclude)

ncfilelist<-c("omega500.mon.mean.nc", "shum_x_vwnd.700.mon.mean.nc", 
              "hgt1000.mon.mean.nc", "slp.mon.mean.nc", 
              #OET "air.1000-500.mon.mean.nc", 
              "shum_x_vwnd.925.mon.mean.nc", 
              "pottmp.1000-500.mon.mean.nc", "shum_x_uwnd.925.mon.mean.nc", 
              "hgt500.mon.mean.nc", "shum_x_uwnd.700.mon.mean.nc", 
              "pottmp.1000-850.mon.mean.nc", "skt.mon.mean.regridded.nc", 
              "air.2m.mon.mean.regridded.nc", "pwtr.mon.mean.nc", 
              "shum925.mon.mean.nc", "shum700.mon.mean.nc")

ncvarlist<-c("omega", "shum", "hgt", "slp", 
             #OET "air", 
             "shum", "pottmp", "shum", 
             "hgt", "shum", "pottmp", "skt", 
             "air", "pr_wtr", "shum", "shum")
labelslist<-c("omega", "merid_transport700", "hgt1000", "slp", 
              #OET "airdiff", 
              "merid_trans925", "pottmp1000_500", "zonal_trans_925", 
              "hgt500", "zonal_trans_700", "pottmp1000_850", "skt", 
              "air2m", "pr_wtr", "shum925", "shum700")


# set starting years to calculate seasonal averages
nc_start_wet<-1949 # for wet season
nc_start_dry<-1948 # for dry season

# NCEP reanalysis data are organized in latitude from north to south
# (this must be taken into account in the plotting)
# set to TRUE if the latitude values are ordered from highest to lowest 
# value, else set FALSE

lat_north_south=TRUE


# select region from the NetCDF data
# ('HI' by default, but individual islands e.g. 'BI' could be used
# see function )
region<-'HI'

###############################################################
# (2.1.4) PCA analysis of rainfall (predictand)
###############################################################



# select the PCA mode, i.e. the PC time series
# to be used as predictand in the regression model
# (pcamode <= pcamax)

pcamode=3


###############################################################
# (2.1.5) Composite analysis (applied to PC time series)
###############################################################

# choose the number of lowest ncomp and highest ncomp 
# years to be used in anomaly composites
ncomp<-4

###############################################################
# Dredge process (fiding the best model using AIC)
###############################################################

# select top ranking models passing the delta-AIC check
# models from the dredge process to be passed through as
# equivalently skillful models (delta<=delta_AIC_max)

delta_AIC_max=2.0


###############################################################
# (2.2) Load the rainfall station data
###############################################################
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

#navigate to folder containing rainfall data files (wet and dry)
Xall<-load_rain_data(datapath,file_rain)

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

###############################################################
# (2.3) Run PCA on rainfall data (predictands) 
###############################################################
# run PCA on the rainfall data and retain the 
# eigenvectors of the first 10 PCs
# <OET 2021-01-17>
# The highest PC time series used to form composites 
# and used as predictand in regression model 
# should can be increased if needed.
# The number may depend on the season
# (we used North Rule of Thumb to find number of PCA modes
# to work with).
# You must adjust the code before the dredge-function call.
# If you want more PC modes define more GLM models 
# in dredge-function (see variable res).
#</OET>
pca_max=4

pr<-prcomp(Xmatrix,center=T,scale=T)

# eigenvectors of the PCA (weight factors for the stations) 
eigenvectors<-pr$rotation

# OET 2021-01-17 (not used elsewhere in this code)
#notScaled=eigenvectors[,1:pcamax]

###############################################################
#OET 2020-01-09 
#OET use function scale 
###############################################################
#scale rainfall data

Xano<-scale(Xmatrix,center=TRUE,scale=TRUE)

#pc is scaled rainfall data multiplied by the eigenvectors
#using truncated number of PCs, will only be 
#able to reconstruct Xano at best

pc<-Xano%*%eigenvectors  

###############################################################
#(2.4) Composite analysis of NetCDF climate data (predictors)
###############################################################


# load the 3-d field data to form a composite analysis


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


# predictand time series
# selects the PC time series used as predictand
print(paste("PCA mode chosen as predictand: ",pcamode))
y<-pc[,pcamode]

iyrs<-length(y)
if (pcamode>pca_max){
  print("!!Warning: higher PCA modes should be considered with caution")
  print(paste("Recommended maximum PCA modes to use:",pca_max))
  print("!!No call model for dredge function defined.")
  print("!!Code execution stopped.")
  stop()
}

# array to store the projection index time series
# for all climate variables (organized in columns)


iclimvar<-length(ncvarlist)
xproj<-array(0,dim=c(iyrs,iclimvar))

##############################################################
# Loop over climate variables
# Note: this loop must be performed 8 times:
# 2 different seasons and 4 different PCs.
# The pc number and season must be changed each time 
# in the code above. each time the loop is
# Run, a unique dataframe for the selected season 
# and pc must be saved.
##############################################################

#read in shapefile of hawaii state from local folder 
# to add outline to composite plots

#OET: problem reading the shp file / or problem usign the function?
#hawaii <- readOGR(datapath, layer="coast_geo_shp")

#this code loops through the list of variables and for each variable: calculates the composite anomaly 
#(high minus low), creates and saves a contour plot for each variable, performs the vector
#projection

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
  
  
  # Composite analysis & regression model
  
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
  
  
 
  # Plotting the composite pattern
  # 2020-12-17 OET: Specific to the NCEP reanalysis data
  # 
  # change lats from north to south to south to north order
  imgfile=file.path(outputdir.figs,
                    paste(labels, "_PC",pcamode,"_", 
                          season,"Comp_Plot.jpg",sep=""))
  
  lon2<-clim_reg$lon
  if (lat_north_south==TRUE){
    print("in")
    lat2rev<-rev(clim_reg$lat)
    irev<-rev(1:length(lat2rev))
    filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode, sep=""))
    # OET ,plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})
    # updated path
    jpeg(imgfile, width = 500, height = 350)
    filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode, ",", season, sep=""))
    #OET ,plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})
    
    p<-filled.contour(clim_reg$lon-360,lat2rev,comp[,irev],main=paste(labels, ", PC",pcamode,", ", season, sep=""))
    #OET plot(hawaii,add=T)
  }else {
    filled.contour(clim_reg$lon-360,clim_reg$lat,comp,main=paste(labels, ", PC",pcamode, sep=""))
    # OET ,plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})
    # updated path
    jpeg(imgfile, width = 500, height = 350)
    filled.contour(clim_reg$lon-360,clim_reg$lat,comp,main=paste(labels, ", PC",pcamode, ",", season, sep=""))
    #OET ,plot.axes={ axis(1); axis(2); plot(hawaii, add=T )})
    
    p<-filled.contour(clim_reg$lon-360,clim_reg$lat,comp,main=paste(labels, ", PC",pcamode,", ", season, sep=""))
    #OET plot(hawaii,add=T)
  }
  #print(p)
  dev.off()
  
  #######################################################
  # Projection index calculation
  #######################################################
  for (i in 1:iyrs) {
    xproj[i,ivar]<-my.proj(xclim2[,i],comp)
  }
}


# y is the eigenvectors, xproj is the predictor information
# save dataframe for each of the 4 pcs


df<-data.frame('y'=y,'x'=xproj)
# OET columns relabeling depends of length of labelslist
# which is the number of climate field variables  (predictors)
colnames(df)[2:(length(labelslist)+1)]<-labelslist

#change the season and pc number in the file name to match the selections used
#OET: use a variable name for the output directory
#     and the file name is created based on selected season 
#     and PCA mode
filename=paste(season,"_PC",pcamode,".csv",sep="")
print("write predictand and predictor data table:")
print(file.path(outputdir.train,filename))

write.csv(df,file.path(outputdir.train,filename), row.names=FALSE)



####################################################################################
#Dredge Code
#note: each res needs to be made for each season/pc combination. The following steps 
#must only be completed while the corresponding season and pc is selected
###################################################################################

#to make the dredge table, first make a model ("res") while the season and pc is selected 
#in above code and then run code for "All possible models" using the res
#############################################################################
###make dredge tables  ##do this 1 pc at a time
#make model with all remaining predictors after eliminating based on multicollinearity

##########################################################################
# OET: This section has to be done by the user who builds the model
# Advantage of using specific variable names
# + pre-selection can be done by users
# + adjustment for each PC mode possible
# + transparent and physical meaning attached with names of predictors
# Cons:
# - not a plug-and-play script !
# - modification must be done every time when one predictor is changed
##########################################################################

#<OET 2020-12-17>
# res<-glm(y~omega + merid_transport700 + hgt1000 + merid_trans925 + zonal_trans_925 + zonal_trans_700  + pottmp1000_850 + air2m, df, family=gaussian, na.action = "na.fail")
# #wet season pc 2
# res<-glm(y~omega +merid_transport700 + airdiff +merid_trans925 +zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + skt +pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
# #wet season pc 3
# res<-glm(y~ omega + merid_transport700 + airdiff + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
# #wet season pc4
# res<-glm(y~omega+ merid_transport700 + +hgt1000+ airdiff+ merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + pr_wtr+ shum925 +shum700, df, family=gaussian, na.action = "na.fail")
#wet season pc 1

# OET: make the selection dependent on the season and PC mode
if (season =='wet'){
  if (pcamode==1){
    #wet season pc 2
    res<-glm(y~omega + merid_transport700 + hgt1000 
         +merid_trans925 + zonal_trans_925
         + zonal_trans_700  + pottmp1000_850 + air2m, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==2){
    #wet season pc 2
    res<-glm(y~omega +merid_transport700 
           #+ airdiff 
           + merid_trans925 +zonal_trans_925 
         + hgt500 + zonal_trans_700 + pottmp1000_850 
         + skt +pr_wtr + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==3){
    #wet season pc 3
    res<-glm(y~ omega + merid_transport700  
           #+ airdiff 
           + zonal_trans_925 + hgt500 + zonal_trans_700 
         + pottmp1000_850 + air2m + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==4){
    #wet season pc4
    res<-glm(y~omega+ merid_transport700 + hgt1000 
           #+ airdiff 
           + merid_trans925 + zonal_trans_925 + 
           hgt500 + zonal_trans_700 + pottmp1000_850 + 
           air2m + pr_wtr+ shum925 +shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
   
} # end wet season

if (season=="dry"){
  # #dry season pc1
  # res<-glm(y~omega + merid_transport700 + slp + merid_trans925 + pottmp1000_500 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 +air2m + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
  # #dry season pc2
  # res<-glm(y~omega + merid_transport700 +hgt1000 + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 + skt + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
  # #dry season pc3
  # res<-glm(y~omega + merid_transport700 + slp + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700 + pottmp1000_850 + air2m + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")
  # #dry season pc4
  # res<-glm(y~omega + merid_transport700 + hgt1000 + airdiff + merid_trans925 + zonal_trans_925 + hgt500 + zonal_trans_700  + pottmp1000_850 + air2m  + pr_wtr + shum925 + shum700, df, family=gaussian, na.action = "na.fail")

  if (pcamode==1){
    #dry season pc1
    res<-glm(y~omega + merid_transport700 + slp 
         + merid_trans925 + pottmp1000_500 
         + zonal_trans_925 + hgt500 
         + zonal_trans_700  + pottmp1000_850 
         +air2m + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==2){
    # #dry season pc2
    res<-glm(y~omega + merid_transport700 
         +hgt1000 
         #+ airdiff 
         + merid_trans925 + zonal_trans_925 
         + hgt500 + zonal_trans_700  + pottmp1000_850 
         + skt + pr_wtr + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==3){
    #dry season pc3
    res<-glm(y~omega + merid_transport700 + slp 
         #+ airdiff 
         + merid_trans925 + zonal_trans_925 + hgt500 
         + zonal_trans_700 + pottmp1000_850 + air2m 
         + pr_wtr + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
  if (pcamode==4){
    # #dry season pc4
    res<-glm(y~omega + merid_transport700 + hgt1000 
         #+ airdiff 
         + merid_trans925 + zonal_trans_925 
         + hgt500 + zonal_trans_700  + pottmp1000_850 
         + air2m  + pr_wtr + shum925 + shum700, 
         df, family=gaussian, na.action = "na.fail")
  }
}
print(paste("Dredge model call for season ",
            season," PC #",pcamode,sep=""))
print(res$call)
#</OET>


####################################################################################
#all possible models
print("execute dredge-function ...")
dredge=dredge(res, rank = "AICc",extra="R^2")
print("... done.")
#save the dredge table to directory each time (change season name and pc number)
filename=paste(season,"_PC",pcamode,"_dredge.csv",sep="")
print("write predictand and predictor data table:")
# if sub-folder does not exist create
dhelp=file.path(outputdir.train,'dredge')
if (!dir.exists(dhelp)){
  dir.create(dhelp)
  print("created subfolder:")
  print(dhelp)
}
print(file.path(dhelp,filename))

ipass=dredge$delta<=delta_AIC_max
write.csv(dredge[ipass,],file.path(dhelp,filename), row.names = FALSE)



