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


fit <- read_stan_csv(paste(output.loc,'prediction_nd_reduced_domain_holocene.csv',sep=''))
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
r <- build_r_nb(post_dat=post_dat,N = 554,T=20,K=13)
rm(post_dat)
saveRDS(r,paste(data.loc,'r_large.RDS',sep=''))
r <- readRDS(paste(data.loc,'r_large.RDS',sep=''))
r_mean <- apply(r$r,c(1,2),median)
rm(r)

load(paste(data.loc,'prediction_13_taxa_554_cells_83_knots_cal_pl_Ka_Kgamma_EPs_116_sites_full_pred_holocene.rdata',sep=''))
taxa <- colnames(y)
limits <- list(xlims=range(coord.agg.final$east),ylims=range(coord.agg.final$north))

#saveRDS(r_mean,paste(data.loc,'r_mean.RDS',sep=''))

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })


plot_pred_maps(r_mean=r_mean,
               centers = coord.agg.final,
               taxa = taxa,
               t = seq(500,10000,500),
               N = N,
               K = K,
               T = 20,
               thresh = breaks,
               limits = limits,
               type = 'prop',
               suff = 'reduced_domain_holocene_median',
               save_plots =TRUE ,
               fpath = plot.loc,
               height = 12,
               width = 12)
