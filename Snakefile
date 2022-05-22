import json
import os

import requests


RANKS = [
    "kingdom",
    "supergroup",
    "division",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

from scripts.utils import make_graph_name


configfile: "config.yaml"


PR2 = "https://github.com/pr2database/pr2database/releases/latest/download/pr2_version_4.14.0_SSU_taxo_long.fasta.gz"


def get_pr2_url():
    response = requests.get(
        "https://api.github.com/repos/pr2database/pr2database/releases/latest"
    )
    j = json.loads(response.text)
    for asset in j["assets"]:
        if asset["name"].endswith("SSU_taxo_long.fasta.gz"):
            return asset["browser_download_url"]


METHODS_EXTENSIONS = [
    ("fastspar", ""),
    # ("spieceasi", "edgelist"),
    ("flashweave", "edgelist"),
    # ("phyloseq", "edgelist"),
    ("conet", "tsv"),
]


class PathHandler:
    def __init__(self):
        self.raw_input = config["input"]["filename"]
        self.metadata_input = (
            config["input"].get("metadata_fname", None)
            if config["input"]["metadata"]
            else None
        )
        self.tax_file = "data/taxonomy.tsv"
        self.method2ext = {
            method: extension for method, extension in METHODS_EXTENSIONS
        }
        self.consensus_network_path = "data/consensus_network.edgelist"
        self.trophic_groups = "data/trophic_groups.xlsx"
        self.known_relations = "data/known_relations.tsv"
        self.lima_mendez_relations = "data/lima-mendez_relations.tsv"
        self.pr2_path = "data/blast/pr2.fasta"
        self.blast_results = "data/blast/results.out.gz"

    @property
    def standarized_input(self):
        path, filename = os.path.split(self.raw_input)
        return os.path.join(path, "data", "standarized_" + filename)

    @property
    def standarized_metadata_input(self):
        if self.metadata_input is None:
            return None
        _, filename = os.path.split(self.raw_input)
        name, extensions = filename.split(".", 1)
        metadata_name = name + "_meta"
        return os.path.join("data", "standarized_" + metadata_name + "." + extensions)

    @property
    def standarized_biom_input(self):
        return os.path.splitext(self.standarized_input)[0] + ".biom"

    @property
    def raw_graphs(self):
        filepaths = []
        for method in self.method2ext:
            filepaths.extend(self.make_outputs(method, only_files=True))
        return filepaths

    @property
    def standarized_graphs(self):
        filepaths = []
        for filepath in self.raw_graphs:
            filepaths.append(make_graph_name(filepath))
        return filepaths

    def make_outputs(self, method, only_files=False):
        out = dict()
        ext = self.method2ext[method]
        for k in config[f"{method}_configs"]:
            if k.startswith("config"):
                out[k] = os.path.join(
                    config[f"{method}_configs"]["output_dir"],
                    k + f".{ext}" if ext else k,
                )
                if method == "fastspar":
                    out[k] = directory(out[k])
        if only_files:
            return list(out.values())
        return out

    @property
    def raw_files(self):
        d = {"base": self.raw_input, "tax": self.tax_file}
        if self.metadata_input:
            d["meta"] = self.metadata_input
        return d

    @property
    def standarized_files(self):
        d = {"base": self.standarized_input}
        if self.standarized_metadata_input:
            d["meta"] = self.standarized_metadata_input
        return d

    @property
    def flashweave_input(self):
        l = [self.standarized_biom_input]
        if self.standarized_metadata_input:
            l.append(self.standarized_metadata_input)
        return l

    @property
    def files_for_visualization(self):
        d = dict()
        for method in self.method2ext:
            for config_name, filepath in self.make_outputs(method).items():
                config_id = config_name.removeprefix("config")
                d[
                    f"{method}" + (f"_{config_id}" if config_id else "")
                ] = make_graph_name(filepath)
        d["consensus"] = self.consensus_network_path
        d["trophic_groups"] = self.trophic_groups
        d["taxonomy"] = self.tax_file
        return d


path_handler = PathHandler()


onstart:
    import shutil
    from pathlib import Path

    path = Path().absolute()
    occurence_table_path = path / "scripts/data_preparation/occurence_table.py"
    target_path = path / "custom_conet/occurence_table.py"
    shutil.copy(occurence_table_path, target_path)


onsuccess:
    import os
    from pathlib import Path

    path = Path().absolute()
    target_path = path / "custom_conet/occurence_table.py"
    os.remove(target_path)


rule download_pr2:
    output:
        path_handler.pr2_path,
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/download_pr2.benchmark"
    shell:
        f"""wget {get_pr2_url()} -O {path_handler.pr2_path}.gz
        gunzip {path_handler.pr2_path}.gz
        """


rule run_blast:
    input:
        path_handler.pr2_path,
        config["input"]["otu_seqs_filename"],
    output:
        path_handler.blast_results,
    conda:
        "envs/blast.yaml"
    benchmark:
        "benchmarks/run_blast.benchmark"
    threads: 8
    shell:
        """makeblastdb -in {input[0]} -title pr2 -dbtype nucl -out data/blast/pr2
        blastn -db data/blast/pr2 -query {input[1]} -num_threads {threads} -outfmt 6 \
        | gzip --best -c \
        > {output[0]}
        """


rule get_taxonomy:
    input:
        path_handler.blast_results,
    output:
        path_handler.tax_file,
    conda:
        "envs/blast.yaml"
    benchmark:
        "benchmarks/get_taxonomy.benchmark"
    threads: 8
    shell:
        """python scripts/data_preparation/parse_blast_results.py {input} {threads} {output}"""


rule standarize_input:
    input:
        **path_handler.raw_files,
    output:
        *list(path_handler.standarized_files.values()),
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/standarize_input.benchmark"
    script:
        "scripts/data_preparation/to_biom_tsv.py"


rule get_known_relations:
    input:
        path_handler.tax_file,
    output:
        "data/pida_v1.08.zip",
        path_handler.known_relations,
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/get_known_relations.benchmark"
    script:
        "scripts/data_collection/collect_known_interactions.py"


rule get_lima_mendez_relations:
    input:
        path_handler.tax_file,
    output:
        "data/W7.xlsx",
        "data/Database_W5_OTU_occurences.tsv",
        path_handler.lima_mendez_relations,
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/get_lima_mendez_relations.benchmark"
    script:
        "scripts/data_collection/collect_lima_mendez_interactions.py"


rule fastspar_infer:
    priority: 0
    input:
        path_handler.standarized_input,
    output:
        **path_handler.make_outputs("fastspar"),
    log:
        "logs/fastspar.log",
    benchmark:
        "benchmarks/fastspar.benchmark"
    conda:
        "envs/fastspar.yaml"
    threads: 2
    script:
        "scripts/method_runners/call_fastspar.py"


# rule SpiecEasi_infer:
#     input:
#         **path_handler.standarized_files
#     output:
#         **path_handler.make_outputs("spieceasi")
#     threads: 8
#     log:
#         "logs/spieceasi.log",
#     benchmark:
#         "benchmarks/spieceasi.benchmark"
#     conda:
#         "envs/spieceasi.yaml"
#     script:
#         "scripts/method_runners/call_SpiecEasi.R"


rule make_biom:
    input:
        **path_handler.standarized_files,
    output:
        path_handler.standarized_biom_input,
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/make_biom.benchmark"
    script:
        "scripts/data_preparation/make_biom.py"


rule flashweave_infer:
    input:
        *path_handler.flashweave_input,
    output:
        **path_handler.make_outputs("flashweave"),
    threads: 5
    log:
        "logs/flashweave.log",
    benchmark:
        "benchmarks/flashweave.benchmark"
    conda:
        "envs/flashweave.yaml"
    script:
        "scripts/method_runners/call_flashweave.jl"


# rule phyloseq_infer:
#     input:
#         **path_handler.standarized_files
#     output:
#         **path_handler.make_outputs("phyloseq")
#     log:
#         "logs/phyloseq.log",
#     benchmark:
#         "benchmarks/phyloseq.benchmark"
#     conda:
#         "envs/phyloseq.yaml"
#     script:
#         "scripts/method_runners/call_phyloseq.R"


rule conet_infer:
    input:
        path_handler.standarized_input,
    output:
        *path_handler.make_outputs("conet", only_files=True),
    log:
        "logs/conet.log",
    threads: 8
    benchmark:
        "benchmarks/conet.benchmark"
    conda:
        "envs/custom_conet.yaml"
    shell:
        """python -u custom_conet/run_correlation_inference.py {input} {threads} {output} >> {log}"""


rule standarize_networks:
    input:
        **{"networks": path_handler.raw_graphs, "tax_table": path_handler.tax_file},
    output:
        *path_handler.standarized_graphs,
    log:
        "logs/standarize_networks.log",
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/standarize_networks.benchmark"
    script:
        "scripts/data_preparation/standarize_networks.py"


rule make_consensus_network:
    input:
        *path_handler.standarized_graphs,
    output:
        path_handler.consensus_network_path,
    log:
        "logs/make_consensus_network.log",
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/make_consensus_network.benchmark"
    script:
        "scripts/data_preparation/make_consensus_network.py"


rule generate_vis_file:
    input:
        **path_handler.files_for_visualization,
    output:
        "visualization/data.json",
    conda:
        "envs/prepare_visualization.yaml"
    benchmark:
        "benchmark/generate_vis_file.benchmark"
    script:
        "scripts/prepare_visualization_file.py"


rule generate_stats:
    input:
        path_handler.tax_file,
        path_handler.known_relations,
        path_handler.lima_mendez_relations,
        path_handler.trophic_groups,
        *path_handler.standarized_graphs,
        path_handler.consensus_network_path,
    output:
        "data/stats.xlsx",
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmark/generate_stats.benchmark"
    log:
        "logs/generate_stats.log"
    threads: 4
    shell:
        "python -m scripts.generate_statistics " + " ".join(
            [input[0], input[1], input[2], input[3], threads, output[0], *input[4:]]
        ) >> logs/generate_stats.log


rule all:
    input:
        "data/stats.xlsx",
        "visualization/data.json",
