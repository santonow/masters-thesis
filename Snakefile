rule fastspar_infer:
    priority: 0
    input:
        "data/fake_data.tsv"
    output:
        directory("data/fastspar_results")
    log:
        "logs/fastspar.log"
    benchmark:
        "benchmarks/fastspar.benchmark"
    conda:
        "envs/fastspar.yaml"
    threads: 2
    shell:
        """mkdir {output}
        fastspar --otu_table {input} --correlation {output}/correlations.tsv \
        --covariance {output}/covariances.tsv --threads {threads} >> {log}
        """

rule SpiecEasi_infer:
    input:
        "data/fake_data.tsv"
    output:
        "data/spieceasi_results/learned_network.tsv"
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
        "data/fake_data.tsv"
    output:
        "data/julia_results/learned_network.edgelist"
    threads: 2
    log:
        "logs/flashweave.log"
    benchmark:
        "benchmarks/flashweave.benchmark"
    conda:
        "envs/flashweave.yaml"
    script:
        "mt/call_flashweave.jl"

rule all:
    input:
        "data/fastspar_results",
        "data/spieceasi_results/learned_network.tsv",
        "data/julia_results/learned_network.edgelist"

