##Continued from Script "PCA_Comp_Analysis_Dredging"
#this script assesses multicollinearity in order to eliminate variables 
#before the dredge process. Needs to be run 8 times for 2 seasons and 4 PCs
#using y (eigenvectors) and df (scaled predictor information) that is unique to each 
#season/PC.


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

##once this is performed for all season/PC combinations, variables that have .9 or higher 
#correlation with other variables (based on correlation matrices) are elimated based on ranking
#from rsquare tables. lower ranked vaiable displaying multicollinearity is removed.
