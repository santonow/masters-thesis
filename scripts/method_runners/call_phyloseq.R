log <- file(snakemake@log[[1]], open="wt")
sink(log, append = TRUE)
sink(log, append = TRUE, type = "message")

library(igraph)
library(phyloseq)
library(Matrix)


otu.table <- otu_table(
  read.csv(snakemake@input[["base"]], header = TRUE, sep = "\t", row.names = 1),
  taxa_are_rows = TRUE
)

configs <- snakemake@config[["phyloseq_configs"]]
for (config.name in names(configs)) {
  if (startsWith(config.name, "config")) {
    config <- configs[[config.name]]
    max.dist <- config[["max.dist"]]
    config[["physeq"]] <- otu.table
    config[["type"]] <- "taxa"
    dists <- distance(otu.table, method = config[["method"]], type = "taxa")
    config[["distance"]] <- dists

    # adapted from phyloseq source code to allow for bigger networks to be processed by igraph
    # see this answer: https://stackoverflow.com/a/58203731/12931685
    if( attributes(dists)$Size != ntaxa(otu.table) ){
      stop("ntaxa(physeq) does not match size of dist object in distance")
    }
    if( !setequal(attributes(dists)$Labels, taxa_names(otu.table)) ){
      stop("taxa_names does not exactly match dist-indices")
    }
    # coerce distance-matrix back into vanilla matrix, Taxa Distance Matrix, TaDiMa
    TaDiMa  <- as.matrix(dists)
    # Add Inf to the diagonal to avoid self-connecting edges (inefficient)
    TaDiMa <- TaDiMa + diag(Inf, ntaxa(otu.table), ntaxa(otu.table))
    # Convert distance matrix to coincidence matrix, CoMa, using max.dist
    CoMa <- TaDiMa < max.dist
    # end of code taken or adapted from phyloseq source code

    net <- graph.adjacency(as(CoMa, "sparseMatrix"), mode="lower")

    for (edge in E(net)) {
      net <- set.edge.attribute(
        graph = net, name = "weight", index = edge, value = 1 - TaDiMa[head_of(net, edge), tail_of(net, edge)]
      )
    }
    write_graph(net, file = snakemake@output[[config.name]], format = "ncol", weights = "weight")
  }
}

sink()
sink(type = "message")