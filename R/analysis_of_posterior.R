########################################################################################################################
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


fit <- read_stan_csv(paste(output.loc,'prediction_nd_reduced_domain_old_calibration.csv',sep=''))
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
post_dat <- list(post = post,par_names = par_names)

r <- build_r_nb(post_dat=post_dat,N = 554,T=19,K=13)
r_mean <- apply(r$r,c(1,2),median)
saveRDS(r$r,paste(output.loc,'predicted_vegetation.RDS',sep=''))


load(paste(data.loc,'test_prediction_old_calibration.rdata',sep=''))
limits <- list(xlims=range(coord.agg.final$east),ylims=range(coord.agg.final$north))
colours <- rev(brewer.pal(10,'RdYlBu'))
taxa <- colnames(y)


breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
 breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                     function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
# 
#r_mean_category <- cut(r_mean, breaks, include.lowest=TRUE, labels=FALSE)
#r_mean_category <- matrix(r_mean_category,ncol=13)


plot_pred_maps(r_mean=r_mean,
               centers = coord.agg.final,
               taxa = taxa,
               t = seq(150,1950,100),
               N = 554,
               K = K,
               T = 19,
               thresh = breaks,
               limits = limits,
               type = 'prop',
               suff = 'test_downcore',
               save_plots =TRUE ,
               fpath = plot.loc,
               height = 36,
               width = 36)
