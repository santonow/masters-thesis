from parse_blast_results import prepare_tax_table

prepare_tax_table(snakemake.input[0], snakemake.threads, snakemake.output[0])
