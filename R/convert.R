#' Convert from a phylo object to a set of innode, flux, and outnode for a Sankey diagram
#'
#' This will create information for a Sankey diagram. The output is a dataframe: one column for the name of the starting node, one column for the flow from that node to the third column, the name of the ending node.
#' If nodelabels exist, they will be used for node names, otherwise numbers will be used.
#'
#' @param phy The tree in ape phylo format
#' @param tip.weights Optionally a vector of the number of taxa each tip represents. If present, the final width matches this. It should have labels matching the taxon names
#' @param erase.labels If TRUE, erase any node labels; otherwise, they will be used for naming nodes if present
#' @return Data.frame with columns for innode, flux, and outnode
#' @export
convert_phylo_to_sankey <- function(phy, tip.weights=NULL, erase.labels=FALSE) {
  if(!inherits(phy, "phylo")) {
    stop("The input must be a phylo class object (see documentation for the package ape)")
  }
  if(erase.labels) {
    phy$node.label <- NULL
  }
  if(!is.null(phy$node.label)) {
    actual.labels <- phy$node.label[which(nchar(phy$node.label)>0)]
    if(any(duplicated(actual.labels))) {
      warning(paste("There are some duplicate node.labels on your tree: ", paste(actual.labels[duplicated(actual.labels)], collapse=", "), " so we are dropping them all"))
      phy$node.label <- NULL
    }
  }
  if(is.null(tip.weights)) {
    tip.weights <- rep(1, ape::Ntip(phy))
    names(tip.weights) <- phy$tip.label
  }
  if(is.null(names(tip.weights))) {
    warning("The tip.weights vector should have taxon names. We're assuming the names match those in the phy object $tip.label, but be careful")
    names(tip.weights) <- phy$tip.label
  }
  result <- data.frame(innode=phy$edge[,1], flux=rep(0, length(phy$edge[,1])), outnode=phy$edge[,2])
  for(node.row in sequence(length(phy$edge[,1]))) {
    desc <- phangorn::Descendants(phy, result$outnode[node.row], type="tips")[[1]]
    result$flux[node.row] <- sum(tip.weights[phy$tip.label[desc]])
  }
  result$innode <- as.character(result$innode)
  result$outnode <- as.character(result$outnode)
  for(node.id in sequence(max(phy$edge))) {
    if(node.id <= ape::Ntip(phy)) { # we're at a tip
      result$outnode[which(result$outnode==node.id)] <- phy$tip.label[node.id]
    } else {
      if(!is.null(phy$node.label)) {
        local.label <- phy$node.label[node.id - ape::Ntip(phy)]
        if(nchar(local.label)>0) {
          result$outnode[which(result$outnode==node.id)] <- local.label
          result$innode[which(result$innode==node.id)] <- local.label
        }
      }
    }
  }
  return(result)
}

#' Convert flux data.frame to input for http://sankeymatic.com
#'
#' Takes a data.frame from convert_phylo_to_sankey() or similar and creates a text block suitable for pasting into http://sankeymatic.com/build
#' You can then cat() the final object to see it.
#'
#' @param sankey The data.frame with innode, flux, and outnode
#' @return A string suitable for putting into sankeymatic to get a diagram
#' @export
convert_to_sankeymatic <- function(sankey) {
  convert_one_row <- function(x) {
    return(paste0(x["innode"], " [", x["flux"], "] ", x["outnode"]))
  }
  return(paste(apply(sankey, 1, convert_one_row), collapse="\n"))
}
