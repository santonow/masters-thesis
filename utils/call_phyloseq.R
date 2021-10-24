log <- file(snakemake@log[[1]], open="wt")
sink(log, append = TRUE)
sink(log, append = TRUE, type = "message")

library(igraph)
library(phyloseq)


otu.table <- otu_table(
  read.csv(snakemake@input[["base"]], header = TRUE, sep = "\t", row.names = 1),
  taxa_are_rows = TRUE
)

configs <- snakemake@config[["phyloseq_configs"]]
for (config.name in names(configs)) {
  if (startsWith(config.name, "config")) {
    config <- configs[[config.name]]
    config[["physeq"]] <- otu.table
    net <- do.call(make_network, config, type="taxa")
    write.table(as_data_frame(net), file=snakemake@output[[config.name]], sep = "\t", )
  }
}

sink()
sink(type = "message")