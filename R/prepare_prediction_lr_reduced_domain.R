#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)
library(dplyr)
library(DT)
library(neotoma)
library(sp)
library(fields)
library(rgdal)
library(abind)
library(rstan)

#-----------------------------------------------------------------------------------------------
setwd('~/git_upload/stepps_full_prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'

########################################################################################################################
#read the huge posterior file of the vegetation model
########################################################################################################################
fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/output/township_120knots_final.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('eta','rho') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])

eta <- est.e.bayes[grep('eta',names(est.e.bayes))]
rho <- est.e.bayes[grep('rho',names(est.e.bayes))]
########################################################################################################################
# load knot coords 
########################################################################################################################
load(paste(data.loc,'township_data_13_taxa_6796_cells_120_knots.rdata',sep=''))
#########################################################################################################################
#load raw vegetation data
load(paste(data.loc,'elicitation_neus_certainty_median_80_sites_only_abies_new_species.RData',sep=''))
#######################################################################################################################
#load pollen data and coordinates of pollen to make a new idx_files
#######################################################################################################################
y <- readRDS(paste(data.loc,'pollen_def.RDS',sep=''))
y <- round(y)
pollen_coords <- readRDS(paste(data.loc,'coordinates_def.RDS',sep='')) 
pollen_coords <- matrix(ncol=2,as.numeric(pollen_coords))
colnames(pollen_coords) <- c('lon','lat')
##########################################################################################################################
# remove pollen data not in domain
##########################################################################################################################
coord_remove <- (pollen_coords[,'lon'] <(-76.5) | pollen_coords[,'lat'] <40.25)  
rm_index <- which(coord_remove == TRUE)

rm_index_total <- sapply(rm_index,function(x) ((x-1)*19+1) : ((x)*19))

rm_index_total <- matrix(ncol=1, rm_index_total)
y_test <- y[-rm_index_total,]
pollen_coords_test <- pollen_coords[coord_remove==FALSE,]
y <- y_test
pollen_coords <- pollen_coords_test 

saveRDS(y,paste(data.loc,'pollen_reduced_domain.RDS',sep=''))
saveRDS(pollen_coords,paste(data.loc,'coordinates_reduced_domain.RDS',sep=''))

######################################################################################################################
# make low resolution coordinates
######################################################################################################################
coord.east.unique <- sort(unique(veg_coords[,1]))
coord.north.unique <- sort(unique(veg_coords[,2]))

coord.east.agg <- seq(coord.east.unique[2],coord.east.unique[length(coord.east.unique)],24000)
coord.north.agg <- seq(coord.north.unique[2],coord.north.unique[length(coord.north.unique)],24000)

coord.agg <- expand.grid(coord.east.agg,coord.north.agg)

veg_coords <- data.frame(veg_coords)
colnames(coord.agg) <- c('east','north')
coord.agg <- data.frame(coord.agg)

dist.coords <- rdist(coord.agg,veg_coords)
min.dist.coords <- apply(dist.coords,1,min)
coord.agg.final <- coord.agg[min.dist.coords==0,] 

########################################################################################################################
coords.neus <- matrix(ncol=2,c(rep(c(-76.5,-67),each =2),rep(c(40.25,47.5),2)))
coords.neus <- as.data.frame(coords.neus)
colnames(coords.neus) <- c('Lon','Lat')


veg_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                     '+init=epsg:4326', 
                     '+init=epsg:3175')

reconst_grid <- build_grid(veg_box, resolution = 24000, proj = '+init=epsg:3175')

coord.agg.final <- coord.agg.final[(coord.agg.final$east>veg_box[1]) & (coord.agg.final$north>veg_box[2]),]
N_cells <- nrow(coord.agg.final)



veg_agg_fake <- r[1:N_cells,]

#######################################################################################################################
#---------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------
calib_trans <- cbind(pollen_coords,as.data.frame(y[seq(1,nrow(y),19),]))
veg_trans <- cbind(coord.agg.final,veg_agg_fake)

colnames(calib_trans)[c(1:2,11:12)] <- c('lon','lat','Other conifer','Other hardwood')

veg_table <- to_stepps_shape(veg_trans,   '~ east + north',      '+init=epsg:3175')
pol_table <- to_stepps_shape(calib_trans, '~ lon + lat', '+init=epsg:4326')

target_taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
                      'Beech','Ash','Tamarack','Spruce','Pine','Oak','Hemlock'))

#change this there is a better location for this file!!
source('R/prep_input_modified.R')
neus_agg <- prep_input(veg = veg_table, 
                       pollen = pol_table, 
                       target_taxa = target_taxa,
                       grid   = reconst_grid,hood = 7e+05)

neus_agg$d <- neus_agg$d/10^6
N <- neus_agg$N_cells
N_cores <- neus_agg$N_cores
idx_cores <- neus_agg$idx_cores
d <- neus_agg$d




#--------------------------------------------------------------------------------------------
# have to load parameters from calibration model
#--------------------------------------------------------------------------------------------
  library(rstan)
  source(paste(help.fun.loc,'pred_helper_funs.r',sep='')) #make sure I use medians of a
  source(paste(help.fun.loc,'process_funs.r',sep='')) #make sure I use medians of a
  source('R/build_cal_main.r') # this is strange....
 
    run <- runs[[4]]
    kernel <- run$kernel
    num_a <- run$one_a
    one_psi <- run$one_psi
    handle <- strsplit(run$suff_fit,'_A')[[1]][1]
    
    #look at that again...
    fname = paste(data.loc,'cal_pl_Ka_Kgamma_EPs_modified_a_110.csv',sep='')
    fit <- read_stan_csv(fname)
    post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)

   


  param.names <-colnames(post[,1,]) #find parameter names of posterior
  param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
  param.e.bayes <- param.greek%in%c('gamma','phi') #only take eta and rho
  est.e.bayes <- colMeans(post[,1,param.e.bayes])
  gamma.all <- est.e.bayes[grep('gamma',names(est.e.bayes))]
  gamma <- gamma.all
  phi.all <- est.e.bayes[grep('phi',names(est.e.bayes))]
  phi <- phi.all

###################################################################################################################
# weight matrixes
###################################################################################################################
w <- build_weight_matrix(post = post,d = t(neus_agg$d),idx_cores = neus_agg$idx_cores,
                         N = neus_agg$N_cells,N_cores =  N_cores,run = run) #check N



sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = N_pot,
                            d_pot =  d_pot,run = run)


#-------------------------------------------------------------------------------------------
#have to rescale gamma
# simplest is to build the kernel and add the weight of the 8 closest locations...
# 
weights_pot <- build_w_pot(post = post, K=K,N_pot = N_pot,d_pot =  d_pot,run = run)
gamma.added <- (weights_pot[1,] + weights_pot[2,])/sum_w_pot
gamma.new <- gamma + gamma.added
gamma <- gamma.new

#################################################################################################################
# knot distance matrix
#################################################################################################################
#should probably only use knots within the domain

knot_coords <- as.data.frame(knot_coords)
colnames(knot_coords) <- c('east','north')

knot_coords_test <- knot_coords[((knot_coords$east>min(coord.agg.final$east))&knot_coords$north>min(coord.agg.final$north)),] 
knot_coords <- knot_coords_test

d_inter <- rdist(coord.agg.final,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6  
N_knots <- nrow(knot_coords)

#--------------------------------------------------------------------------------------------
res <- 3 
ages    = seq(150,1950,100)
T       = length(ages) 
lag     = unname(as.matrix(dist(matrix(ages), upper=TRUE)))

#-------------------------------------------------------------------------------------------
num_sites <-length(idx_cores)
#--------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
#store data in .dump file
stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores','lag', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
    ), 
file=paste(data.loc,
                 'prediction_',K,'_taxa_',N,'_cells_',N_knots,'_knots_',handle,'_',num_sites,
                 '_sites_full_pred_old_calib.dump',sep=""))


save(K, N,T, N_knots, N_cores,lag, 
       y, res,
       sum_w_pot,
       rho,eta,gamma,phi,
       idx_cores,
       d,d_knots,d_inter,w,
     coord.agg.final, 
file=paste(data.loc,
           'prediction_',K,'_taxa_',N,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'
           _sites_full_pred_old_calib.rdata',sep=""))


###############################################################################################################################
# test for the most moden site
###########################################################################################################################
y <- as.data.frame(y[seq(1,nrow(y),19),])
T <- 1
res <- 3

stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores','lag', 
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
  ), 
  file=paste(data.loc,'test_prediction_old_calibration.dump',sep=''))

save(K, N,T, N_knots, N_cores,lag, 
     y, res,
     sum_w_pot,
     rho,eta,gamma,phi,
     idx_cores,
     d,d_knots,d_inter,w,coord.agg.final,
     file=paste(data.loc,'test_prediction_old_calibration.rdata',sep=''))
