log <- file(snakemake@log[[1]], open="wt")
sink(log, append = TRUE)
sink(log, append = TRUE, type = "message")

library(SpiecEasi)
library(phyloseq)
library(igraph)

otu.table <- otu_table(
  read.csv(snakemake@input[["base"]], header = TRUE, sep = "\t", row.names = 1),
  taxa_are_rows = TRUE
)

configs <- snakemake@config[["spieceasi_configs"]]
for (config.name in names(configs)) {
  if (startsWith(config.name, "config")) {
    config <- configs[[config.name]]
    config[["data"]] <- otu.table
    config[["pulsar.params"]] <- list(ncores=snakemake@threads)
    sf <- do.call(spiec.easi, config)
    net <- adj2igraph(getRefit(sf), vertex.attr=list(name=taxa_names(otu.table)))
    write_graph(net, file = snakemake@output[[config.name]], format = "ncol", weights = "weight")
    # df <- as.data.frame(net)
    # write.table(df, file = snakemake@output[[config.name]], sep = "\t")
  }
}

sink()
sink(type = "message")
