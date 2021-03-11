#' reshuffle interactions from an edge list
#'
#' @param edge_list columns site,focal,neighbour,magnitue
#' @param scale reshuffle at the "site" or "network" scale
#' @param keep.diagonal flag, whether diagonal (intra-specific) terms
#' are kept fixed (TRUE) or randomized as well (FALSE)
#' @return dataframe with same columns, magnitude reshuffled
#' @export
#'
#' @examples
reshuffle_spatial_network <- function(edge_list, scale = "site", keep.diagonal = FALSE){
  
  if(scale == "network"){
    
    if(keep.diagonal){
      els <- edge_list[complete.cases(edge_list),]
      m <- sample(els$magnitude[which(els$focal != els$neighbour)])
      els$magnitude[which(els$focal != els$neighbour)] <- m
    }else{
      els <- edge_list[complete.cases(edge_list),]
      m <- sample(els$magnitude)
      els$magnitude <- m
    }# if-else keep diagonal
  }else{
    if(keep.diagonal){
      els <- NULL
      sites <- sort(unique(edge_list$site))
      for(i.site in 1:length(sites)){
        elss <- edge_list[edge_list$site == sites[i.site],]
        elss <- elss[complete.cases(elss),]
        m <- sample(elss$magnitude[which(elss$focal != elss$neighbour)])
        elss$magnitude[which(elss$focal != elss$neighbour)] <- m
        els <- rbind(els,elss)
      }# for each site
    }else{
      els <- NULL
      sites <- sort(unique(edge_list$site))
      for(i.site in 1:length(sites)){
        elss <- edge_list[edge_list$site == sites[i.site],]
        elss <- elss[complete.cases(elss),]
        m <- sample(elss$magnitude)
        elss$magnitude <- m
        els <- rbind(els,elss)
      }# for each site
    }# if-else keep diagonal
  }# if-else scale
  els
}