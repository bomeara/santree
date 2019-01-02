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
    xpos[node.index] <- phytools::nodeheight(phy, node.index)
  }
  #node_labels <- as.character(1:(tip.label(phy) + node.label(phy))) # this doesn't work.

  ID1 <- as.character(plotly$links$source)
  ID2 <- as.character(plotly$links$target)
  Value <- plotly$links$value
  Edges <- data.frame(N1=ID1,N2=ID2,Value=Value)
  ape::plot.phylo(phy, plot=FALSE) # plotting to get yy
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv) # trick learned from phytools
  river_phy <- riverplot::makeRiver(nodes = Nodes, edges = Edges, node_xpos = xpos, node_ypos = obj$yy)

  return(river_phy)
}

#' Uses the Paleobiology Database to get diversity of a lineage through time
#' @param lineage Lineage name
#' @param rank What rank to count at; one of species, genera, genera_plus, families, orders
#' @param time_reso What time resolution to use; one of stage, epoch, period, era
#' @return A data.frame of diversity from PBDB, or NULL if the lineage is not found
get_diversity_for_lineage <- function(lineage, rank="genera", time_reso="stage") {
  diversities <- read.csv(url(paste0("https://paleobiodb.org/data1.2/occs/diversity.txt?base_name=",utils::URLencode(lineage) ,"&recent&count=", utils::URLencode(rank), "&time_reso=", utils::URLencode(time_reso))), stringsAsFactors=FALSE)
  if(ncol(diversities) < 5) {
    diversities <- NULL
  }
  return(diversities)
}

#' Combine diversity of lineages
#'
#' When moving down a tree, merge the diversity of lineages (i.e., common ancestor of acrogymnosperms and angiosperms).
#' @param diversity list List of diversity data.frames
#' @return A data.frame of diversity from PBDB, or NULL if no lineage is not found
merge_diversities_for_lineages <- function(diversity_list) {
  diversities <- diversity_list[[1]]
  if(length(diversity_list)>1) {
    for (additional_lineage_index in sequence(length(diversity_list)-1)) {
      new.diversities <- diversity_list[[additional_lineage_index + 1]]
      diversities <- merge(diversities,new.diversities, by=c('interval_name', 'interval_no', "max_ma", "min_ma"))
      columns <- c("X_Ft","X_bL", "X_FL", "X_bt", "sampled_in_bin", "implied_in_bin", "n_occs")
      for(i in seq_along(columns)) {
        diversities[,columns[i]] <- diversities[,paste0(columns[i], ".x")] +  diversities[,paste0(columns[i], ".y")]
      }
      diversities <- diversities[,-grep("\\.x", colnames(diversities))]
      diversities <- diversities[,-grep("\\.y", colnames(diversities))]
    }
  }
  return(diversities)
}

#' Get spindle diagram
#'
#' A spindle diagram shows the diversity at different time intervals.
#'
#' The Paleobiology Database has information about the number of taxa originating in an interval and persisting after it (X_Ft), the number of taxa originating before the interval and going extinct within it (X_bL), the number of taxa found only within an interval (X_FL), and the number of taxa that originate before an interval and persist after but which are not sampled in an interval (X_bt) [the odd capitalization is in PBDB].
#'
#' We want to use this information for a spindle diagram showing diversity through time. X_bL+X_bt are the number of taxa alive at the interval's start; X_Ft+X_bt are the number of taxa alive at the interval's end. The number of taxa alive at any point within the interval could be as low as X_bt or as high as X_bt+X_FL+max(X_Ft,X_bL). As a compromise, at the midpoint of an interval the width is X_bt+X_FL
#'
#' rank is what rank to count: number of species, genera, families, or orders (there is also a genera_plus option that counts subgenera as genera). While in the study of modern diversity species is the rank most commonly used, in paleontology genera or higher are more usual.
#'
#' time_reso is what resolution to use: stage is finest, era is coarsest
#'
#' This function assumes all taxon names have been resolved to match PBDB's names and that the tree is a chronogram
#' @param phy a chronogram of the groups you want
#' @inheritParams get_diversity_for_lineage
#' @param youngest_taxon How many years ago the youngest tip lives
#' @return a riverplot object
#' @export
convert_phylo_to_spindle_river <- function(phy, rank="genera", time_reso="stage", youngest_taxon=0) {
  phy <- ape::reorder.phylo(phy, "postorder")
  ape::plot.phylo(phy, plot=FALSE, use.edge.length=TRUE) # plotting to get yy and xx
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv) # trick learned from phytools
  #lineage_names <- list()
  lineage_diversities <- list()
  lineage_heights <- phytools::nodeHeights(phy)
  lineage_ranges <- youngest_taxon + max(lineage_heights) - lineage_heights # get depths from present
  nodes.df <- data.frame()
  diversities_list <- lapply(phy$tip.label, merge_diversities_for_lineages) # in node number order
  for (edge_index in sequence(nrow(phy$edge))) {
    my_descendants <- phytools::getDescendants(phy, phy$edge[edge_index,2])
    tip_descendants <- my_descendants[which(my_descendants<=ape::Ntip(phy))]
    if(phy$edge[,2]>ape::Ntip(phy)) {
      diversities_list[[phy$edge[,2]]] <- merge_diversities_for_lineages[tip_descendants]
    }
  }

  # Now traverse the tree, storing heights of each segment and diversity at segment end and midpoint

  # Now convert the x dimensions from time to x positions if needed
  xx_range <- range(obj$xx)

}
