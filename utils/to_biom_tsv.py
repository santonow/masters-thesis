from occurence_table import OTUTable


if snakemake.config["input"]["format"] == "tsv":
    if snakemake.config["input"]["metadata"]:
        metadata_path = snakemake.input["meta"]
    else:
        metadata_path = None
    table = OTUTable.from_tsv(
        snakemake.input["base"],
        sample_rows=snakemake.config["input"]["sample_rows"],
        sample_metadata_fpath=metadata_path
    )

    table.to_tsv(snakemake.output[0], min_tot_abundance=snakemake.config["input"]["min_tot_abundance"])

elif snakemake.config["input"]["format"] == "biom":
    table = OTUTable.from_hdf5(snakemake.input[0])

    table.to_tsv(snakemake.output[0], min_tot_abundance=snakemake.config["input"]["min_tot_abundance"])

