using FlashWeave
using Pkg

if "FlashWeave" ∉ keys(Pkg.installed())
    Pkg.add("FlashWeave")
end

open(snakemake.log[1], "w") do out
    redirect_stdout(out) do
        redirect_stderr(out) do
            data_path = snakemake.input[1]
            netw_results = learn_network(
                data_path, sensitive=true, heterogeneous=false, transposed=true
            )
            save_network(snakemake.output[1], netw_results)
        end
    end
end
