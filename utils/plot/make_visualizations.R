# library(phyloseq)
# library(igraph)
library(ggplot2)
library(network)

if (!require("ggnet")) {
  devtools::install_github("briatte/ggnet")
}
library(ggnet)


# create.layout <- function(net) {
#   # adapted from this answer: https://stackoverflow.com/a/52672660/12931685
#   net.grouped <- net
#   E(net.grouped)$weight <- 1
#
#   ## Add edges with high weight between all nodes in the same group
#   for(i in unique(vertex_attr(net, attr_to_group_on, V(net)))) {
#       group.V <- which(
#         vertex_attr(net, attr_to_group_on, V(net)) == i
#       )
#       if ( length(group.V) > 1 ) {
#         net.grouped <- add_edges(
#           net.grouped, combn(group.V, 2), attr=list(weight=attraction_weight))
#       }
#   }
#
#   ## Now create a layout based on G_Grouped
#   set.seed(420)
#   LO <- layout_method(net.grouped)
# }

transform.weight <- function(x) {
  # return(3 * ((x / max(x)) + 0.2))
  return(log(x) / 3)
}

for ( d in snakemake@input ) {
  tax.df <- read.csv(
    file.path(d, "reduced_taxonomy.tsv"),
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
  # tax.table <- tax_table(as.matrix(tax.df))
  edges <- read.csv(file = file.path(d, "reduced_graph.edgelist"), header = F, as.is = T, sep = "\t")
  edges[, 1] <- as.character(edges[, 1])
  edges[, 2] <- as.character(edges[, 2])
  edges$weight <- transform.weight(edges[, 3])
  net <- network(edges, verbices <- cbind(vertex.names = rownames(tax.df), tax.df))

  palette <- RColorBrewer::brewer.pal(9, "Set1")
  names(palette) <- unique(tax.df$supergroup)
  p <- ggnet2(
    net, edge.size = "weight", color = "supergroup",
    edge.alpha=0.5, size=2.5, alpha=0.8, palette = palette, label = "species", label.size = 2
  )
  ggsave(snakemake@output[[d]], p)


  # net <- graph_from_data_frame(d = edges, directed = F, vertices = cbind(tax = rownames(tax.df), tax.df))
  # E(net)$weight <- E(net)$V3
  # V(net)$name <- as.character(V(net)$name)
  #
  # net <- delete_vertices(net, which(tax.df$supergroup == 'Eukaryota'))
  #
  # dummy_otu_table <- cbind(otu_id = as.vector(V(net)), sample1 = 1:length(V(net)))
  #
  # physeq.obj <- phyloseq(dummy_otu_table, tax.table)
  #
  # hjust <- 1.35
  #
  # plot <- plot_network(
  #   net, physeq.obj, type = "taxa", color = attr_to_group_on,
  #   label = NA, point_size = 1.8, line_weight = 0.8,
  #   layout = layout, line_alpha = 0.2,
  # )# + geom_text(
  #   #aes_string(label="division"), size = 2, hjust = hjust, na.rm = TRUE
  # #)
  # ggsave(snakemake@output[[d]])
}