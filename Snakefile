import os
import requests
import json


RANKS = [
    'kingdom',
    'supergroup',
    'division',
    'class',
    'order',
    'family',
    'genus',
    'species'
]

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
    # ("spieceasi", "edgelist"),
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

def prepare_filenames(base_fname, prefix="standarized_"):
    _, filename = os.path.split(base_fname)
    name, extensions = filename.split('.', 1)
    metadata_name = name + "_meta"
    metadata_filepath = os.path.join("data", prefix + metadata_name + '.' + extensions)
    d = {
        "base": os.path.join("data", prefix + filename)
    }
    if config["input"]["metadata"]:
        d["meta"] = metadata_filepath
    return d

def prepare_biom_filename(base_fname, prefix="standarized_"):
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
        "envs/file_manipulation.yaml"
    shell:
        f"""wget {get_pr2_url()} -O data/blast/pr2.fasta.gz
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

rule get_known_relations:
    input:
        "data/taxonomy.tsv",
    output:
        "data/pida_v1.08.zip",
        "data/known_relations.tsv"
    conda:
        "envs/file_manipulation.yaml"
    script:
        "utils/collect_known_interactions.py"

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
        *(
            [prepare_biom_filename(config["input"]["filename"])] + (
                [prepare_filenames(config["input"]["filename"])["meta"]] if config["input"]["metadata"] else []
            )
        )
    output:
        **make_outputs("flashweave", "edgelist")
    threads: 5
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


def prepare_network_files(wrapper=lambda x: x, input=True):
    d = prepare_filenames(config["input"]["filename"])
    d["networks"] = []
    for method, extension in METHODS_EXTENSIONS:
        d["networks"].extend(make_outputs(method, extension, wrapper, only_files=True))
    d["tax_table"] = "data/taxonomy.tsv"
    dirnames = []
    if not input:
        d["configs"] = []
        for network_fname in d["networks"]:
            fname = os.path.split(network_fname)[1].split(".")[0]
            new_dir = os.path.join(os.path.split(network_fname)[0], fname)
            dirnames.append(directory(new_dir))
            d["configs"].append(
                {
                    "tax_table": os.path.join(new_dir, "reduced_taxonomy.tsv"),
                    "network": os.path.join(new_dir, "reduced_graph.edgelist")
                }
            )
            for rank in RANKS:
                new_dir = os.path.join(os.path.split(network_fname)[0], rank + '_' + fname)
                dirnames.append(directory(new_dir))
                d["configs"].append(
                    {
                        "tax_table": os.path.join(new_dir, "reduced_taxonomy.tsv"),
                        "network": os.path.join(new_dir, "reduced_graph.edgelist")
                    }
                )
        for key in list(d):
            if key not in ["networks", "tax_table", "configs"]:
                del d[key]
    all_files = []
    all_files.extend(d["networks"])
    if not input:
        for conf in d["configs"]:
            all_files.append(conf["tax_table"])
            all_files.append(conf["network"])
    for k, v in d.items():
        if k not in {"networks", "configs", "tax_table"}:
            all_files.append(v)

    return d, all_files, dirnames

rule standarize_networks:
    input:
        **prepare_network_files(input=True)[0]
    output:
        *prepare_network_files(make_graph_name, input=False)[1],
        *prepare_network_files(make_graph_name, input=False)[2]
    log:
        "logs/standarize_networks.log"
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/standarize_networks.benchmark"
    script:
        "utils/standarize_networks.py"

rule make_consensus_network:
    input:
        *prepare_network_files(make_graph_name, input=False)[0]["networks"],
        "data/taxonomy.tsv"
    output:
        "data/consensus_network.edgelist",
        directory("data/consensus_network"),
        "data/consensus_network/reduced_taxonomy.tsv",
        "data/consensus_network/reduced_graph.edgelist"
    log:
        "logs/make_consensus_network.log"
    conda:
        "envs/file_manipulation.yaml"
    benchmark:
        "benchmarks/standarize_networks.benchmark"
    script:
        "utils/make_consensus_network.py"

def prepare_inputs_for_vis():
    res = dict()
    res["networks"] = prepare_network_files(make_graph_name)[0]["networks"]
    res["otu_table"] = prepare_filenames(config["input"]["filename"])["base"]
    res["tax_table"] = "data/taxonomy.tsv"
    return res

def prepare_outputs_for_vis():
    outputs = dict()
    for inp in prepare_network_files(make_graph_name, input=False)[2]:
        base, name = os.path.split(inp)
        outputs[inp] = os.path.join(base, 'vis_' + name + '.png')
    return outputs

rule generate_plots:
    input:
        *prepare_network_files(make_graph_name, input=False)[2],
        "data/consensus_network"
    output:
        **{"data/consensus_network": "data/vis_consensus_network.png"},
        **prepare_outputs_for_vis(),

    conda:
        "envs/phyloseq.yaml"
    script:
        "utils/plot/make_visualizations.R"


