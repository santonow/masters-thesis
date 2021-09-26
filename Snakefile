import os

configfile: "config.yaml"

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
    print(out)
    return out


def prepare_output_fname(base_fname, method):
    dirname, filename = os.path.split(base_fname)
    return os.path.join(dirname, f"{method}_" + filename)

rule prepare_for_fastspar:
    input:
        config["input"]["filename"]
    output:
        prepare_output_fname(config["input"]["filename"], "fastspar")
    script:
        "mt/to_biom_tsv.py"

rule fastspar_infer:
    priority: 0
    input:
        prepare_output_fname(config["input"]["filename"], "fastspar")
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

