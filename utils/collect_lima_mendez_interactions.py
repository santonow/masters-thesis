import csv
import os
import subprocess
import zipfile

import pandas as pd

LIMA_MENDEZ_URL = 'http://www.raeslab.org/companion/ocean_interactome/tables/W7.xlsx'
LIMA_MENDEZ_PATH = snakemake.output[0]
W5_DB_URL = 'http://taraoceans.sb-roscoff.fr/EukDiv/data/Database_W5_OTU_occurences.tsv.zip'
zipped_W5_DB_path = snakemake.output[1] + '.zip'
W5_DB_path = snakemake.output[1]
TAX_PATH = snakemake.input[0]

subprocess.run(['wget', LIMA_MENDEZ_URL, '-O', snakemake.output[0]])
subprocess.run(['wget', W5_DB_URL, '-O', zipped_W5_DB_path])
with zipfile.ZipFile(zipped_W5_DB_path) as zf:
    zf.extractall(os.path.split(zipped_W5_DB_path)[0])
    os.remove(zipped_W5_DB_path)


md5sum2cid = dict()
with open(W5_DB_path) as handle:
    reader = csv.reader(handle, delimiter='\t')
    next(reader)
    for md5sum, cid, *_ in reader:
        md5sum2cid[md5sum] = cid


existing_cids = set()
with open(TAX_PATH) as handle:
    reader = csv.reader(handle, delimiter='\t')
    next(reader)
    existing_cids.update(cid for cid, *_ in reader)


df = pd.read_excel(LIMA_MENDEZ_PATH, sheet_name='Network')

with open(snakemake.output[2], 'w') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['node1', 'node2', 'kl score', 'spearman score', 'sign'])
    for node1, node2, kl_score, spearman_score, inter_type in zip(
        df['Node 1'], df['Node 2'], df['Best Kullback Leibler score'], df['Best Spearman score '], df['Interaction type']
    ):
        if all(node in md5sum2cid for node in [node1, node2]):
            if md5sum2cid[node1] in existing_cids and md5sum2cid[node2] in existing_cids:
                writer.writerow(
                    [
                        md5sum2cid[node1],
                        md5sum2cid[node2],
                        kl_score,
                        spearman_score,
                        '+' if inter_type == 'copresence' else '-'
                    ]
                )
