import os

from mt.utils import make_graph_name

configfile: "config.yaml"


METHODS_EXTENSIONS = [
    ("fastspar", ""),
    ("spieceasi", "tsv"),
    ("flashweave", "edgelist"),
    ("phyloseq", "tsv")
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

def prepare_filename_for_flashweave(base_fname, prefix="sanitized_"):
    names = prepare_filenames(base_fname, prefix)
    filename = os.path.splitext(names["base"])[0]
    return filename + '.biom'

rule sanitize_input:
    input:
        **pack_input()
    output:
        **prepare_filenames(config["input"]["filename"])
    conda:
        "envs/file_manipulation.yaml"
    script:
        "mt/to_biom_tsv.py"

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
        "mt/call_fastspar.py"

rule SpiecEasi_infer:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        **make_outputs("spieceasi", "tsv")
    threads: 2
    log:
        "logs/spieceasi.log"
    benchmark:
        "benchmarks/spieceasi.benchmark"
    conda:
        "envs/spieceasi.yaml"
    script:
        "mt/call_SpiecEasi.R"

rule make_biom_for_flashweave:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        prepare_filename_for_flashweave(config["input"]["filename"])
    conda:
        "envs/file_manipulation.yaml"
    script:
        "mt/make_biom.py"

rule flashweave_infer:
    input:
        prepare_filename_for_flashweave(config["input"]["filename"]),
        prepare_filenames(config["input"]["filename"])["meta"]
    output:
        **make_outputs("flashweave", "edgelist")
    threads: 2
    log:
        "logs/flashweave.log"
    benchmark:
        "benchmarks/flashweave.benchmark"
    conda:
        "envs/flashweave.yaml"
    script:
        "mt/call_flashweave.jl"

rule phyloseq_infer:
    input:
        **prepare_filenames(config["input"]["filename"])
    output:
        **make_outputs("phyloseq", "tsv")
    log:
        "logs/phyloseq.log"
    benchmark:
        "benchmarks/phyloseq.benchmark"
    conda:
        "envs/phyloseq.yaml"
    script:
        "mt/call_phyloseq.R"


def prepare_network_inputs(wrapper=lambda x: x):
    d = prepare_filenames(config["input"]["filename"])
    d["networks"] = []
    for method, extension in METHODS_EXTENSIONS:
        d["networks"].extend(make_outputs(method, extension, wrapper, only_files=True))
    return d

rule all:
    input:
        **prepare_network_inputs()
    output:
        *prepare_network_inputs(make_graph_name)["networks"]
    log:
        "logs/standarize_networks.log"
    benchmark:
        "benchmarks/standarize_networks.benchmark"
    script:
        "mt/standarize_networks.py"



