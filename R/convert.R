#' Convert from a phylo object to a set of source, value, and target for a Sankey diagram
#'
#' This will create information for a Sankey diagram. The output is a dataframe: one column for the name of the starting node, one column for the flow from that node to the third column, the name of the ending node.
#' If nodelabels exist, they will be used for node names, otherwise numbers will be used.
#'
#' @param phy The tree in ape phylo format
#' @param tip.weights Optionally a vector of the number of taxa each tip represents. If present, the final width matches this. It should have labels matching the taxon names
#' @param erase.labels If TRUE, erase any node labels; otherwise, they will be used for naming nodes if present
#' @param convert.numbers If TRUE, convert node numbers to labels if possible
#' @return Data.frame with columns for source, value, and target
#' @export
convert_phylo_to_sankey <- function(phy, tip.weights=NULL, erase.labels=FALSE, convert.numbers = TRUE) {
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
  result <- data.frame(source=phy$edge[,1], value=rep(0, length(phy$edge[,1])), target=phy$edge[,2])
  for(node.row in sequence(length(phy$edge[,1]))) {
    desc <- phangorn::Descendants(phy, result$target[node.row], type="tips")[[1]]
    result$value[node.row] <- sum(tip.weights[phy$tip.label[desc]])
  }
  if(convert.numbers) {
    result$source <- as.character(result$source)
    result$target <- as.character(result$target)
    for(node.id in sequence(max(phy$edge))) {
      if(node.id <= ape::Ntip(phy)) { # we're at a tip
        result$target[which(result$target==node.id)] <- phy$tip.label[node.id]
      } else {
        if(!is.null(phy$node.label)) {
          local.label <- phy$node.label[node.id - ape::Ntip(phy)]
          if(nchar(local.label)>0) {
            result$target[which(result$target==node.id)] <- local.label
            result$source[which(result$source==node.id)] <- local.label
          }
        }
      }
    }
  }
  return(result)
}

#' Convert value data.frame to input for http://sankeymatic.com
#'
#' Takes a data.frame from convert_phylo_to_sankey() or similar and creates a text block suitable for pasting into http://sankeymatic.com/build
#' You can then cat() the final object to see it.
#'
#' @param sankey The data.frame with source, value, and target
#' @return A string suitable for putting into sankeymatic to get a diagram
#' @export
convert_to_sankeymatic <- function(sankey) {
  convert_one_row <- function(x) {
    return(paste0(x["source"], " [", x["value"], "] ", x["target"]))
  }
  return(paste(apply(sankey, 1, convert_one_row), collapse="\n"))
}


#' Convert from a phylo object to a set of Source, Target, Value, Node Label for plotly sankey diagram
#' convert data to data frame for use in plotly sankey diagram
#'
#' @param phy The tree in ape phylo format
#' @param tip.weights Optionally a vector of the number of taxa each tip represents. If present, the final width matches this. It should have labels matching the taxon names
#' @param erase.labels If TRUE, erase any node labels; otherwise, they will be used for naming nodes if present
#' @return A list of links and nodes for sankey plot
#' @export
convert_phylo_to_plotly <- function(phy, tip.weights=NULL, erase.labels=FALSE) {
  links.df <- convert_phylo_to_sankey(phy, tip.weights=tip.weights, erase.labels=FALSE, convert.numbers=FALSE)
  links <- list(source=links.df$source, target=links.df$target, value=links.df$value)
  node.labels <- rep("", ape::Ntip(phy) + ape::Nnode(phy))
  node.labels[1:ape::Ntip(phy)] <- phy$tip.label #node.label produces replacement of zero length

  if(!is.null(phy$node.label)) {
    node.labels[(ape::Ntip(phy)+1) : (ape::Ntip(phy) + ape::Nnode(phy)) ] <- phy$node.label
  }

  #internal <- node.labels[(ape::Ntip(phy)+1) : (ape::Ntip(phy) + ape::Nnode(phy)) ] <- phy$node.label # error: replacement has length zero because no internal nodes have labels, need if else?
  #nodes <- tips + internal
  result_pl <- list(links=links,nodes=node.labels) # problem may be here, this is a list of lists but tips is shorter than the links...need all of the nodes but not all are labeled
  return(result_pl)

}


#' Convert plotly object from convert_phylo_to_plotly to a sankey diagram
#'
#' @param phy The tree in ape phylo format
#' @param color A vector of colors for the nodes
#' @return A plotly santree
#' @export
convert_to_plotly_santree <- function(phy, color="") {
  result_pl <- convert_phylo_to_plotly(phy)
  my.plot <- plot_ly(
    type = "sankey",
    orientation = "h",

    node = list(
      label = result_pl$nodes, #maybe the problem is here - there are fewer labels than nodes
      color = color,
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = result_pl$links
    )
  return(my.plot)
}


#' Convert from a phylo object to an object containing Node, xpos, and edges (ID1, ID2, and Value)
#'
#' @param phy The tree in ape phylo format
#' @return A riverplot object
#' @export
convert_phylo_to_river <- function(phy, tip.weights=NULL, erase.labels=FALSE) {
  plotly <- convert_phylo_to_plotly(phy, tip.weights=tip.weights, erase.labels=FALSE)
  Nodes <- as.character(1:(Ntip(phy) + Nnode(phy))) 
  xpos  <- rep(NA, ape::Ntip(phy) + ape::Nnode(phy))
  for (node.index in sequence(ape::Ntip(phy) + ape::Nnode(phy))) {
    Nodes[node.index] <- phytools::nodeheight(phy, node.index)
  }

  ID1 <- plotly$links$source
  ID2 <- plotly$links$target
  Value <- plotly$links$value
  Edges <- data.frame(ID1=ID1,ID2=ID2,Value=Value)

  river_phy <- riverplot::makeRiver(nodes = Nodes, edges = Edges, node_xpos = xpos)

  return(river_phy)


}
