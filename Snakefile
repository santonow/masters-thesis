import os

configfile: "config.yaml"

def make_outputs(method, ext, wrapper=None, only_files=False):
    inputs = dict()
    for k in config[f"{method}_configs"]:
        if k.startswith("config"):
            inputs[k] = os.path.join(
                config[f"{method}_configs"]["output_dir"],
                k + f".{ext}" if ext else ""
            )
            if wrapper is not None:
                inputs[k] = wrapper(inputs[k])
    if only_files:
        return list(inputs.values())
    return inputs


rule fastspar_infer:
    priority: 0
    input:
        config["input"]["filename"]
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
        config["input"]["filename"]
    output:
        **make_outputs("spieceasi", "tsv")
    log:
        "logs/spieceasi.log"
    benchmark:
        "benchmarks/spieceasi.benchmark"
    conda:
        "envs/spieceasi.yaml"
    script:
        "mt/call_SpiecEasi.R"

rule flashweave_infer:
    input:
        config["input"]["filename"]
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
        config["input"]["filename"]
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

rule all:
    input:
        *make_outputs("fastspar", "", directory, only_files=True),
        *make_outputs("spieceasi", "tsv", only_files=True),
        *make_outputs("flashweave", "edgelist", only_files=True),
        *make_outputs("phyloseq", "tsv", only_files=True)

