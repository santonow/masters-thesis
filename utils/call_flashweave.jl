using Pkg

if "FlashWeave" ∉ keys(Pkg.installed())
    Pkg.add("FlashWeave")


using Distributed
addprocs(snakemake.threads)

@everywhere using FlashWeave

open(snakemake.log[1], "a") do out
    redirect_stdout(out) do
        redirect_stderr(out) do
            for (configname, config) in snakemake.config["flashweave_configs"]
                if startswith(configname, "config")
                    config_d = (Symbol(key) => val for (key, val) in config)
                    data_path = snakemake.input[1]
                    metadata_path = snakemake.input[2]
                    netw_results = learn_network(data_path, metadata_path; config_d...)
                    save_network(snakemake.output[configname], netw_results)
                end
            end
        end
    end
end
