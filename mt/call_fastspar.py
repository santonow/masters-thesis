import os
import subprocess

os.makedirs(os.path.abspath(snakemake.config['fastspar_configs']['output_dir']), exist_ok=True)
for config_name, config in snakemake.config['fastspar_configs'].items():

    if config_name.startswith('config'):
        output_path = os.path.abspath(snakemake.output[config_name])
        os.makedirs(output_path, exist_ok=True)
        command = [
            'fastspar',
            '--otu_table',
            os.path.abspath(snakemake.input[0]),
            '--correlation',
            os.path.join(output_path, 'correlations.tsv'),
            '--covariance',
            os.path.join(output_path, 'covariances.tsv'),
            '--threads',
            str(snakemake.threads),
        ]
        for option in ['iterations', 'exclude_iterations', 'threshold']:
            if option in config:
                command.extend(
                    [
                        f'--{option}',
                        str(config[option])
                    ]
                )
        completed_process = subprocess.run(command, capture_output=True)
        with open(snakemake.log[0], 'w') as handle:
            handle.write(f'STDOUT ({config_name})\n')
            handle.write(completed_process.stdout.decode() + '\n')
            handle.write(f'STDERR ({config_name})\n')
            handle.write(completed_process.stderr.decode() + '\n')


