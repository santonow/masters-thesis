import os

from occurence_table import OTUTable

if snakemake.config["input"]["format"] == "tsv":
    table = OTUTable.from_tsv(
        snakemake.input[0], sample_rows=snakemake.config["input"]["sample_rows"]
    )

    table.to_tsv(snakemake.output[0])
