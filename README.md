# santree
Package for making Sankey diagrams for phylogenies

```
library(ape)
devtools::install_github("bomeara/santree")
data("bird.orders")
result <- santree::convert_to_sankeymatic(santree::convert_phylo_to_sankey(bird.orders))
cat(result)
```




And then paste the result to http://sankeymatic.com/build/. You may need to play with figure height and such.

```
library(ape)
library(plotly)
data("bird.orders")
result_pl <- (santree::convert_phylo_to_plotly(bird.orders))
santree::convert_to_plotly_santree(result_pl)

```

```
library(ape)
library(riverplot)
data("bird.orders")
result_rv <- (santree::convert_phylo_to_river(bird.orders))

```

```
library(ape)
library(riverplot)
phy <- rcoal(5)
result_rv <- (santree::convert_phylo_to_river(phy))
plot(result_rv, nodewidth=0)
```
