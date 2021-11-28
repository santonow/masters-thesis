from occurence_table import OTUTable

table = OTUTable.from_tsv(
    snakemake.input["base"], snakemake.input.get("meta"), sample_rows=False
)

table.to_hdf5(snakemake.output[0])

