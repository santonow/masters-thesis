import csv
import pkg_resources

import pandas as pd



CHOSEN_COLS = [
    'Latitude (degree N)',
    'Longitude (degree E)',
    'Depth (m)',
    'Temperature (°C)',
    'Salinity (mg/m3)',
    'Oxygen (µmol/L)',
    'Nitrates (µmol/L)',
    'Chlorophyll (µmol/L)',
    'Phosphates (µmol/L)',
    'Distance to coast (km)',
    'NPP',
    'NO2NO3 (mmol/m3)',
    'Iron (µmol/m3)',
    'NH4 (mmol/m3)',
]


CHOSEN_SAMPLES = {
    'amplicon': pkg_resources.resource_filename(
        __name__, 'resources/amplicon_chosen_samples.txt'
    ),
    'metatrans': pkg_resources.resource_filename(
        __name__, 'resources/metatrans_chosen_samples.txt'
    )
}


def read_stations(experiment_type: str) -> set[str]:
    return set(
        x.strip() for x in open(CHOSEN_SAMPLES[experiment_type]).readlines()
    )


def create_metadata_table(output_file: str, experiment_type: str) -> None:
    if experiment_type not in {'amplicon', 'metatrans'}:
        raise ValueError('experiment_type has to be one of [amplicon, metatrans]!')

    metadata = pd.read_excel(
        pkg_resources.resource_filename(__name__, 'resources/metadata/41467_2017_2342_MOESM7_ESM.xlsx'),
        engine='openpyxl', header=1
    )

    chosen_sample_labels = read_stations(experiment_type)

    metadata_records = []

    for record in metadata.to_dict(orient='record'):
        if record['Sample label (Tara Oceans ID)'] in chosen_sample_labels:
            new_record = dict()
            new_record['sample'] = record['Sample label (Tara Oceans ID)']
            for key in CHOSEN_COLS:
                new_record[key] = record[key]
            metadata_records.append(new_record)

    with open(output_file, 'w') as handle:
        keys = list(metadata_records[0].keys())[1:]
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(keys)
        for record in sorted(metadata_records, key=lambda x: x['sample']):
            row = []
            for key in keys:
                row.append(record[key])
            writer.writerow(row)
