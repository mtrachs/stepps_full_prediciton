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
setwd('~/workflow_stepps_prediction/analysis')
help.fun.loc <- 'utils/'
data.loc <- '~/workflow_stepps_prediction/prediction/data/'
plot.loc <- 'plots/'
output.loc <- '~/workflow_stepps_prediction/prediction/output/'
#########################################################################################################################
source(paste(help.fun.loc,'pred_helper_funs.r',sep=''))
source(paste(help.fun.loc,'pred_plot_funs.r',sep=''))
source(paste(help.fun.loc,'difference_estimation.R',sep=''))


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

differences_final <- difference_bins(dataset = r_new_format,
                                     coordinates=coord.agg.final,
                                     quantiles = 0.5)


##################################################################################################################3
#
##################################################################################################################
load(paste(data.loc,'prediction_13_taxa_554_cells_83_knots_cal_pl_Ka_Kgamma_EPs_116_sites_full_pred_holocene.rdata',sep=''))
taxa <- colnames(y)
limits <- list(xlims=range(coord.agg.final$east),ylims=range(coord.agg.final$north))

#saveRDS(r_mean,paste(data.loc,'r_mean.RDS',sep=''))

breaks = c(0, 0.02, 0.1, 0.2, 0.3, 0.4, 0.6, 1)
breaks <- unique(c(-rev(breaks),breaks))
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })


plot_pred_maps(r_mean=differences_final,
               centers = coord.agg.final,
               taxa = taxa,
               t = seq(150,1850,100),
               N = N,
               K = K,
               T = 18,
               thresh = breaks,
               limits = limits,
               type = 'prop',
               suff = 'reduced_domain_diff_two_millennia',
               save_plots =TRUE ,
               fpath = plot.loc,
               height = 12,
               width = 12)

  