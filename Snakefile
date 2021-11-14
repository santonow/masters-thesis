import os
import platform
import requests
import json


from utils.utils import make_graph_name

configfile: "config.yaml"

PR2 = "https://github.com/pr2database/pr2database/releases/latest/download/pr2_version_4.14.0_SSU_taxo_long.fasta.gz"

def get_pr2_url():
    response = requests.get('https://api.github.com/repos/pr2database/pr2database/releases/latest')
    j = json.loads(response.text)
    for asset in j['assets']:
        if asset['name'].endswith('SSU_taxo_long.fasta.gz'):
            return asset['browser_download_url']

# PLATFORM = platform.system()
# if PLATFORM == "Darwin":
#     qiime_env_yaml_name = "qiime2-2021.8-py38-osx-conda.yml"
# elif PLATFORM == "Linux":
#     qiime_env_yaml_name = "qiime2-2021.8-py38-linux-conda.yml"
# else:
#     raise SystemError(f"Unsupported system: {PLATFORM}")


METHODS_EXTENSIONS = [
    ("fastspar", ""),
    ("spieceasi", "edgelist"),
    ("flashweave", "edgelist"),
    ("phyloseq", "edgelist")
]



def make_outputs(method, ext, wrapper=None, only_files=False):
    out = dict()
    for k in config[f"{method}_configs"]:
        if k.startswith("config"):
            out[k] = os.path.join(
                config[f"{method}_configs"]["output_dir"],
                k + f".{ext}" if ext else k
            )
            if wrapper is not None:
                out[k] = wrapper(out[k])
    if only_files:
        return list(out.values())
    return out


def pack_input():
    inp = {
        "base": config["input"]["filename"]
    }
    if "metadata_fname" in config["input"]:
        inp["meta"] = config["input"]["metadata_fname"]
    inp["tax"] = "data/taxonomy.tsv"
    return inp

def prepare_filenames(base_fname, prefix="sanitized_"):
    _, filename = os.path.split(base_fname)
    name, extensions = filename.split('.', 1)
    metadata_name = name + "_meta"
    metadata_filepath = os.path.join("data", prefix + metadata_name + '.' + extensions)
    return {
        "base": os.path.join("data", prefix + filename),
        "meta": metadata_filepath
    }

def prepare_biom_filename(base_fname, prefix="sanitized_"):
    names = prepare_filenames(base_fname, prefix)
    filename = os.path.splitext(names["base"])[0]
    return filename + ".biom"

# rule download_qiime2_env_yaml:
#     output:
#         f"envs/{qiime_env_yaml_name}"
#     shell:
#         f"wget -P envs https://data.qiime2.org/distro/core/{qiime_env_yaml_name}"

rule download_pr2:
    output:
        f"data/blast/pr2.fasta"
    conda:
        "envs/data_download.yaml"
    shell:
        f"""echo {get_pr2_url()}
        wget {get_pr2_url()} -O data/blast/pr2.fasta.gz
        gunzip data/blast/pr2.fasta.gz
        """

rule run_blast:
    input:
        "data/blast/pr2.fasta",
        config['input']['otu_seqs_filename']
    output:
        "data/blast/results.out.gz"
    conda:
        "envs/blast.yaml"
    threads: 8
    shell:
        """makeblastdb -in {input[0]} -title pr2 -dbtype nucl -out data/blast/pr2
        blastn -db data/blast/pr2 -query {input[1]} -num_threads {threads} -outfmt 6 \
        | gzip --best -c \
        > data/blast/results.out.gz
        """

rule get_taxonomy:
    input:
        "data/blast/results.out.gz"
    output:
        "data/taxonomy.tsv"
    conda:
        "envs/blast.yaml"
    threads: 8
    shell:
        """python utils/parse_blast_results.py {input} {threads} {output}"""

rule standarize_input:
    input:
        **pack_input()
    output:
        **prepare_filenames(config["input"]["filename"])
    conda:
        "envs/file_manipulation.yaml"
    script:
        "utils/to_biom_tsv.py"

rule fastspar_infer:
    priority: 0
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        **make_outputs("fastspar", "", directory)
    log:
        "logs/fastspar.log"
    benchmark:
        "benchmarks/fastspar.benchmark"
    conda:
        "envs/fastspar.yaml"
    threads: 2
    script:
        "utils/call_fastspar.py"

rule SpiecEasi_infer:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        **make_outputs("spieceasi", "edgelist")
    threads: 8
    log:
        "logs/spieceasi.log"
    benchmark:
        "benchmarks/spieceasi.benchmark"
    conda:
        "envs/spieceasi.yaml"
    script:
        "utils/call_SpiecEasi.R"

rule make_biom:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        prepare_biom_filename(config["input"]["filename"])
    conda:
        "envs/file_manipulation.yaml"
    script:
        "utils/make_biom.py"

rule flashweave_infer:
    input:
        prepare_biom_filename(config["input"]["filename"]),
        prepare_filenames(config["input"]["filename"])["meta"]
    output:
        **make_outputs("flashweave", "edgelist")
    threads: 3
    log:
        "logs/flashweave.log"
    benchmark:
        "benchmarks/flashweave.benchmark"
    conda:
        "envs/flashweave.yaml"
    script:
        "utils/call_flashweave.jl"

rule phyloseq_infer:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        **make_outputs("phyloseq", "edgelist")
    log:
        "logs/phyloseq.log"
    benchmark:
        "benchmarks/phyloseq.benchmark"
    conda:
        "envs/phyloseq.yaml"
    script:
        "utils/call_phyloseq.R"


def prepare_network_inputs(wrapper=lambda x: x):
    d = prepare_filenames(config["input"]["filename"])
    d["networks"] = []
    for method, extension in METHODS_EXTENSIONS:
        d["networks"].extend(make_outputs(method, extension, wrapper, only_files=True))
    return d

rule standarize_networks:
    input:
        **prepare_network_inputs()
    output:
        *prepare_network_inputs(make_graph_name)["networks"]
    log:
        "logs/standarize_networks.log"
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/standarize_networks.benchmark"
    script:
        "utils/standarize_networks.py"



