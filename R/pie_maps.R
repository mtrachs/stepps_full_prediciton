#########################################################################################################################
# make pie plots of fossil data
#########################################################################################################################
library(sp)
library(maptools)
setwd('~/git_upload/stepps_full_prediction/')
data.loc <- 'data/'
plot.loc <- 'plots/' 
source('utils/dataPlotFuns.r')
source('utils/simDataFuns.r')

#################################################################################################################
#load vegetation coordinates
#################################################################################################################
load(paste(data.loc,'prediction_13_taxa_554_cells_83_knots_cal_pl_Ka_Kgamma_EPs_87
           _sites_full_pred_old_calib.rdata',sep=''))
rm(y)
centers   = coord.agg.final#centers_veg
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])
shift=30000
######################################################################################################################

#load pollen data
pollen.def <- readRDS(paste(data.loc,'pollen_reduced_domain.RDS',sep=''))
pollen_coordinates <- readRDS(paste(data.loc,'coordinates_reduced_domain.RDS',sep=''))
pollen_coordinates <- matrix(as.numeric(pollen_coordinates),ncol=2)
sputm <- SpatialPoints(pollen_coordinates, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords

new.order = order(colSums(pollen.def), decreasing=TRUE)
taxa <- colnames(pollen.def)

taxa = taxa[new.order]
pollen.def = pollen.def[,new.order]



col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
             "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",'cyan',"#CAB2D6")

years <- seq(150,1950,100) 

pdf(paste(plot.loc,'pie_maps_downcore.pdf',sep=''),height =15, width = 15)
sapply(1:19, function(x){
  pollen.use <- pollen.def[seq(x,nrow(pollen.def),19),]
  pollen = pollen.use
  pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
  colnames(pollen_props) = taxa
  
  centers <- pollen_coord_us
  colnames(centers) = c('x', 'y')
  
  centers <- centers[is.finite(rowSums(pollen_props)),]
  pollen_props <- pollen_props[is.finite(rowSums(pollen_props)),]
  
  par(mfrow=c(1,1))
  pieMap(proportions = pollen_props, 
         centers  = centers,
         restrict = FALSE,
         inputRestricted = FALSE,
         xlim   = c(xlo+shift, xhi-shift),
         ylim   = c(ylo+shift, yhi-shift),
         radius = 18000,
         scale  = 1,
         xlab   = 'x',
         ylab   = 'y', 
         add_legend=TRUE,
         main_title= paste(years[x], 'cal BP') ,
         col_list=col_list)
  
})
dev.off()



