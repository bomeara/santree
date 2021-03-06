% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert_phylo_to_sankey}
\alias{convert_phylo_to_sankey}
\title{Convert from a phylo object to a set of source, value, and target for a Sankey diagram}
\usage{
convert_phylo_to_sankey(phy, tip.weights = NULL, erase.labels = FALSE,
  convert.numbers = TRUE)
}
\arguments{
\item{phy}{The tree in ape phylo format}

\item{tip.weights}{Optionally a vector of the number of taxa each tip represents. If present, the final width matches this. It should have labels matching the taxon names}

\item{erase.labels}{If TRUE, erase any node labels; otherwise, they will be used for naming nodes if present}

\item{convert.numbers}{If TRUE, convert node numbers to labels if possible}
}
\value{
Data.frame with columns for source, value, and target
}
\description{
This will create information for a Sankey diagram. The output is a dataframe: one column for the name of the starting node, one column for the flow from that node to the third column, the name of the ending node.
If nodelabels exist, they will be used for node names, otherwise numbers will be used.
}
