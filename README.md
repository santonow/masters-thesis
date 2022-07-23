# Workflow for inferring microbial interaction networks

[![code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)[![code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)


## Installation

### Workflow

- Clone the repo locally.
- Create snakemake virtual environment (using `conda` or `mamba`): `conda create -c conda-forge -c bioconda -n snakemake snakemake=7.8.0`.
- Run `conda activate snakemake`.

### Visualization

If you want to run an interactive visualization of the workflow results, you'll need to install node: `conda install node=17.8.0`. Then, navigate to the `visualization` directory and run `npm install`.

## Running the workflow

### Configuration

The workflow configuration can be adjusted by editing config.yaml file. 
Most important parameters are in the input dict:
```yaml
input:
  filename: "full_amplicon_occurences.tsv"            # OTU table file
  otu_seqs_filename: "db_w5_seqs.fasta"               # fasta file with consensus/ASV sequences 
                                                      # the identifiers should correspond to OTU ids
  format: "tsv"                                       # file format (tsv/csv/biom) 
  sample_rows: true                                   # whether samples are in rows or in columns
  metadata: false                                     # whether to use metadata file 
                                                      # (currently not used due to only one method allowing it)
  metadata_fname: "amplicon_station_metadata.tsv"     # metadata file
  min_tot_abundance: null                             # minimum total abundance of an OTU
  relative_abundance: 0.00000001                      # minimum relative abundance
  zero_threshold: 0.03                                # controls how many OTU counts can be zero (for 0.03 -> 97%)
```

The workflow can also be customized by including a `taxonomy.tsv` file in `data` directory, with a format like this (column names matter!):
| otu_id	| kingdom	| supergroup	| division	| class	| order	| family	| genus	| species |
| ------    | -------   | ----------    | --------  | ----- | ----- | ------    | ----- | ------- |
| 127368	| Eukaryota	| Opisthokonta	| Metazoa	| Arthropoda	| Crustacea	| Maxillopoda	| Maxillopoda_X	| Maxillopoda_X_sp |
| 122743	| Eukaryota	| Rhizaria	| Radiolaria	| Polycystinea	| Nassellaria	| Sphaerozoidae |
| 143719	| Eukaryota	| Opisthokonta	| Metazoa	| Arthropoda	| Crustacea	| Maxillopoda	| Maxillopoda_X	| Maxillopoda_X_sp |
| 95090	| Eukaryota	| Opisthokonta	| Metazoa	| Arthropoda	| Crustacea	| Maxillopoda	| Maxillopoda_X	| Maxillopoda_X_sp |
| 170886	| Eukaryota	| Rhizaria	| Radiolaria	| Polycystinea	| Nassellaria |

This avoids running the steps of the workflow that involve classifying consensus sequences with `blast` and also removes the requirement to provide consensus sequences.

### Execution

Run `snakemake --use-conda -c n all`, where `n` is the number of cores you want to be used. After finishing (which takes around a day for 16 cores and OTU table with dimensions 13927x334), the results will be available in the `data` directory.

### Results

The results for each method are in the `data/<rule-name>_results` directories. There will be raw results that come out of each method as well as standarized results in a `.tsv` format, with top `100000` most strong associations from each network.

An Excel spreadsheet is also produced in the `data` directory, which has basic statistics for each network, a comparison between the networks and a list of edges for each network.

In the `visualization` directory there will be a new file `data.json`, which is an input to the visualization.

### Visualization

Assuming you installed node as in the [Installation](#installation) section, go to `visualization` directory and run `npm start`. 
The visualization will be available at `localhost:8989`. 
You can customize the port by changing it in the `ws --spa index.html --port 8989` command in the `package.json` file.

In the `visualization` directory there is already a sample `data.json` file, so you can see
the visualization without running the whole workflow.
