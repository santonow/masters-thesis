library(SpiecEasi)
library(phyloseq)

log <- file(snakemake@log[[1]], open="wt")
sink(log)

otu.table <- otu_table(
  data.frame(
    read.table(
      snakemake@input[[1]], header = FALSE, sep = "\t"
    )
  ),
  taxa_are_rows = FALSE
)

sf <- spiec.easi(otu.table)

net <- getOptNet(sf)

write.table(net, file=smakemake@input[[1]], sep = "\t", )
