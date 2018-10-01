####################################################################################################################
# Try to analyze individual runs of STEPPS prediction
####################################################################################################################
library(rioja)
###############
#########################################################################################################################
#
#########################################################################################################################
library(rstan)
library(RColorBrewer)
library(fields)
setwd('~/workflow_stepps_prediction/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
output.loc <- 'output/'
#########################################################################################################################
source(paste(help.fun.loc,'pred_helper_funs.r',sep=''))
source(paste(help.fun.loc,'pred_plot_funs.r',sep=''))


fit <- read_stan_csv(paste(output.loc,'prediction_nd_reduced_domain_114_sites.csv',sep=''))
saveRDS(fit,paste(data.loc,'fit_data.RDS',sep=''))
fit <- readRDS(paste(data.loc,'fit_data.RDS',sep=''))
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
rm(fit)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
post_dat <- list(post = post,par_names = par_names)
rm(post)
saveRDS(post_dat,paste(data.loc,'post_data.RDS',sep=''))
post_dat <- readRDS(paste(data.loc,'post_data.RDS',sep=''))
r <- build_r_nb(post_dat=post_dat,N = 554,T=19,K=13)
rm(post_dat)
saveRDS(r,paste(data.loc,'r_large.RDS',sep=''))
r <- readRDS(paste(data.loc,'r_large.RDS',sep=''))
r_mean <- apply(r$r,c(1,2),mean)
rm(r)

load(paste(data.loc,'prediction_13_taxa_554_cells_83_knots_cal_pl_Ka_Kgamma_EPs_114_sites_full_pred.rdata',sep=''))
taxa <- colnames(y)








#realocate data of one time slice into 19 lists

r_new_format <- list()

for(i in 1:T){
  r_new_format[[i]] <- r$r[seq(i,N*T,T),,]  
}
# hopefully one list now stands for predictions for one time slice. 
rm(r)

#########################################################################################################################
#Harvard Forest coordinates
#########################################################################################################################
coord.harv.for <- data.frame(lon = -72.2,lat = 42.53)
#coord.harv.for <- data.frame(lon = -68.77,lat = 44.9)#Caribou Bog
#coord.harv.for <- data.frame(lon = -71.61,lat = 42.52)#

sputm <- SpatialPoints(coord.harv.for , proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
hf_us <- spgeo@coords

hf_index <- which.min(paldist2(coord.agg.final,hf_us, dist.method = 'euclidean'))

########################################################################################################################
#look at beech for Harvard Forest 
time_scale <-seq(150,1950,100) 



hf_site <- 
  lapply(1:K,function(taxon){
    sapply(1:length(r_new_format),function(x){
        r_new_format[[x]][hf_index,taxon,]
    })
})

#estimate quantiles for all sites
pdf(paste(plot.loc,'harvard_forest.pdf',sep=''),height=5,width = 15)
#par(mfrow=c(1,2))
layout(mat=matrix(ncol=3,c(1,1,2)))
sapply(1:K,function(taxon){  
  
  prediction <- hf_site[[taxon]]
  qt_hf <- apply(prediction,2,function(x) quantile(x, probs=c(0.025,0.16,0.5,0.83,0.975)))

#matplot(time_scale,t(qt_hf),type='l',xlim=c(2000,0),ylim=c(0,1))

  matplot(time_scale,t(qt_hf),type='n',xlim=c(2000,0),ylim=c(0,1),main = taxa[taxon])
  polygon(c(time_scale,rev(time_scale)),c(qt_hf['2.5%',],rev(qt_hf['97.5%',])),col='gray')
  polygon(c(time_scale,rev(time_scale)),c(qt_hf['16%',],rev(qt_hf['83%',])),col='lightblue',border='lightblue')
  lines(time_scale,qt_hf['50%',],lwd=2)
  abline(v = time_scale[time_scale%in%c(1250,150)],lty=2)

  #difference firs sample last sample
  diff_sample <- prediction[,time_scale==150] - prediction[,time_scale==1250]
  quantile_change <- quantile(diff_sample,probs=c(0.025,0.16,0.5,0.83,0.975))
  #histogram of differences for individual runs
  hist(diff_sample,main =taxa[taxon],xlim=c(-1,1))
  abline(v = quantile_change,lty=2,col=2:6,lwd=2)
  
  #probability of positive difference
  data.frame(probability = sum(diff_sample>0)/2000,quantile=t(round(quantile_change,2)))
})
dev.off()


