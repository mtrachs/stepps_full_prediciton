difference_bins <- function(dataset,coordinates,quantiles){
  
  list_diffs <- #first keep all those grid points separated
    lapply(1:nrow(coordinates),function(zz){
      
      site <- # aggregate data of one site into a list 
        lapply(1:K,function(taxon){
          sapply(1:length(dataset),function(x){
            dataset[[x]][zz,taxon,]
          })
        })
    
      diffs <- 
        sapply(1:length(site),function(x){
          differences <- apply(site[[x]],1,diff)
          qts <- apply(differences,1,function(x) quantile(x,probs=quantiles))
        })
      
    })
  
  differences_final <- list_diffs[[1]] # estimate differences
  
  for(i in 2:length(list_diffs)){
    differences_final <- rbind(differences_final,list_diffs[[i]])
  }
 
    differences_final #output, matrix N *T-1 X K
} 
  


