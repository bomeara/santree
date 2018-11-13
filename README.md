# santree
Package for making Sankey diagrams for phylogenies

```
library(ape)
devtools::install_github("bomeara/santree")
data("bird.orders")
result <- santree::convert_to_sankeymatic(santree::convert_phylo_to_sankey(bird.orders))
cat(result)
```

And then paste the result to http://sankeymatic.com/build/.
