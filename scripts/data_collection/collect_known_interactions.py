import csv
import json
import subprocess
import zipfile

import pandas as pd


manual_mapping = {
    "Diatomea": "Bacillariophyta",
    "Ancyromonadidae": "Planomonadida",
    "Bolidomonas": "Triparma",
    "Carpediemonas-like": "Carpediemonas",
    "Centrohelida": "Centroheliozoa",
    "Choanomonada": "Choanoflagellida",
    "Incertae_sedis_rhizaria": "Rhizaria",
    "Incertae_sedis_stramenopilees": "Stramenopiles",
    "Incertae_sedis_stramenopiles": "Stramenopiles",
    "unknown_alveolate": "Alveolata",
    "unknown_amoebozoan": "Amoebozoa",
    "unknown_eukaryote": "Eukaryota",
    "Eukaryote": "Eukaryota",
    "unknown_rhizarian": "Rhizaria",
    "Telonema_subtilis": "Telonema_subtile",
}


def group_to_X(s):
    splitted = s.split("_")
    if len(splitted) < 3:
        return s

    if splitted[-2] == "Group":
        return "_".join(splitted[:-2]) + "_" + "X" * len(splitted[-1])


def match_lineage(lineage, lineages):
    options = []
    if str(lineage[-1]) != "nan":
        options.append("_".join([lineage[-2], lineage[-1]]))
        if options[0] in manual_mapping:
            options.append(manual_mapping[options[0]])
        if lineage[-2] in manual_mapping:
            options.append("_".join([manual_mapping[lineage[-2]], lineage[-1]]))
        species_group_t_X = group_to_X(lineage[-1])
        if species_group_t_X != lineage[-1]:
            options.append(species_group_t_X)
    options.append(lineage[-2])
    if lineage[-2] in manual_mapping:
        options.append(manual_mapping[lineage[-2]])
    genera_group_t_X = group_to_X(lineage[-2])
    if genera_group_t_X != lineage[-2]:
        options.append(genera_group_t_X)
    for option in options:
        if option.lower() in lineages:
            return lineages[option.lower()]


PIDA_URL = "https://github.com/ramalok/PIDA/archive/refs/tags/v1.08.zip"

status = subprocess.run(["wget", PIDA_URL, "-O", snakemake.output[0]])

lineages = dict()
lineages_genus = dict()
with open(snakemake.input[0], "r") as handle:
    reader = csv.reader(handle, delimiter="\t")
    next(reader)
    for _, *lineage in reader:
        for i, elem in enumerate(lineage):
            if elem:
                lineages[elem.lower()] = lineage[: i + 1]
        if len(lineage) == 8:
            lineage = lineage[:-1]
        for i, elem in enumerate(lineage):
            if elem:
                lineages_genus[elem.lower()] = lineage[: i + 1]


with zipfile.ZipFile(snakemake.output[0]) as zfile:
    for fname in zfile.filelist:
        if fname.filename.endswith(".xlsx"):
            with zfile.open(fname) as handle:
                df = pd.read_excel(handle, header=2)


relations = set()
relations_genus_level = set()
for record in df.to_dict(orient="records"):
    lineage_left = [record[f"Taxonomic level {i}: org1"] for i in range(1, 4)] + [
        record["Genus org1"],
        record["Species org1"],
    ]
    lineage_right = [record[f"Taxonomic level {i}: org2"] for i in range(1, 4)] + [
        record["Genus org2"],
        record["Species org2"],
    ]
    if any(x[0] == "Prokaryote" for x in [lineage_left, lineage_right]):
        continue
    left_match = match_lineage(lineage_left, lineages)
    right_match = match_lineage(lineage_right, lineages)
    left_match_genus = left_match[:-1] if left_match is not None and len(left_match) == 8 else left_match
    right_match_genus = right_match[:-1] if right_match is not None and len(right_match) == 8 else right_match

    if left_match is not None and right_match is not None:
        relations.add((tuple(left_match), tuple(right_match), record["Interaction"]))

    if left_match_genus is not None and right_match_genus is not None:
        relations_genus_level.add(
            (tuple(left_match_genus), tuple(right_match_genus), record["Interaction"])
        )

with open(snakemake.output[1], "w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["head", "tail", "interaction"])
    for head, tail, interaction in relations:
        writer.writerow(
            [
                json.dumps(head),
                json.dumps(tail),
                interaction,
            ]
        )

with open(snakemake.output[2], "w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["head", "tail", "interaction"])
    for head, tail, interaction in relations_genus_level:
        writer.writerow(
            [
                json.dumps(head),
                json.dumps(tail),
                interaction,
            ]
        )
