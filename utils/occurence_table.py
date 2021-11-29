import csv
import gzip
import os
from typing import Set
from functools import partial

import biom
from biom.util import biom_open
import h5py
import numpy as np


class OTUTable:

    def __init__(self, table: biom.Table, otus_with_lineage: Set[str]):
        self.table = table
        self.otus_with_lineage = otus_with_lineage

    def __getattr__(self, item):
        if item not in self.__dict__:
            return getattr(self.table, item)

    @staticmethod
    def create_opener(gzipped=False):
        if gzipped:
            return partial(gzip.open, mode='rt')
        else:
            return partial(open, mode='r')

    @staticmethod
    def create_writer(gzipped=False):
        if gzipped:
            return partial(gzip.open, mode='wt')
        else:
            return partial(open, mode='w')

    @staticmethod
    def read_metadata(fname: str):
        metadata = dict()
        with open(fname) as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader)
            header = [
                val.replace('/', ' per ') for val in header
            ]
            for line in reader:
                metadata[line[0]] = dict(zip(header[1:], line[1:]))
        return metadata

    def to_csv(
        self, filepath, delimiter=',', gzipped=False,
        metadata_suffix='meta', transpose=False,
        min_tot_abundance=None, filter_tax=True,
        relative_abundance=1e-8, remove_most_zeros=True
    ):
        with self.create_writer(gzipped)(filepath) as writer:
            if transpose:
                tbl: biom.Table = self.table.transpose()
            else:
                tbl = self.table
            chosen_ids = set(self.table.ids(axis='observation'))
            if filter_tax:
                chosen_ids &= self.otus_with_lineage
                filtered_by_tax = len(set(self.table.ids(axis='observation')) - chosen_ids)
            if min_tot_abundance is not None:
                chosen_ids &= {
                    _id for _id in tbl.ids(axis='observation')
                    if sum(tbl.data(_id, axis='observation')) > min_tot_abundance
                }
            filtered_by_total_abundance = len(set(self.table.ids(axis='observation')) - chosen_ids) - filtered_by_tax

            if relative_abundance is not None:
                # based on https://doi.org/10.1126/science.1262073
                tbl_sum = tbl.sum(axis='whole')
                chosen_ids &= {
                    _id for _id in tbl.ids(axis='observation')
                    if sum(tbl.data(_id, axis='observation'))/tbl_sum >= relative_abundance
                }
            filtered_by_relative_abundance = len(
                set(self.table.ids(axis='observation')) - chosen_ids
            ) - filtered_by_tax - filtered_by_total_abundance
            if remove_most_zeros:
                # based on https://doi.org/10.1371/journal.pcbi.1002606, altered to filter only OTUs with >99/100 zeros
                n_samples = tbl.length(axis='sample')
                chosen_ids &= {
                    _id for _id in tbl.ids(axis='observation')
                    if list(tbl.data(_id, axis='observation')).count(0)/n_samples < 99/100
                }
            filtered_by_zero_count = len(
                set(self.table.ids(axis='observation')) - chosen_ids
            ) - filtered_by_tax - filtered_by_total_abundance - filtered_by_relative_abundance

            n_obs = len(tbl.ids(axis='observation'))

            tbl = tbl.filter(
                ids_to_keep=chosen_ids,
                axis='observation'
            )
            print(
                f'Filtering {n_obs - len(chosen_ids)} / {n_obs} OTUs.\n'
                f'- {filtered_by_tax} because no assigned eukaryotic taxonomy.\n'
                f'- {filtered_by_total_abundance} because low total abundance (< {min_tot_abundance}).\n'
                f'- {filtered_by_relative_abundance} because of low relative abundance (< {relative_abundance}).\n'
                f'- {filtered_by_zero_count} because proportion of zeros greater than 99/100'
            )
            for line in tbl.delimited_self(delimiter).split('\n')[1:]:
                writer.write(line + '\n')
        if self.table._sample_metadata is not None:
            pref, fname = os.path.split(filepath)
            name, extensions = fname.split('.', 1)
            metadata_name = name + '_' + metadata_suffix
            if extensions.endswith('.gz'):
                extensions = extensions[:-len('.gz')]
            metadata_filepath = os.path.join(pref, metadata_name + '.' + extensions)
            self.table.metadata_to_dataframe(axis='sample').to_csv(
                metadata_filepath, sep=delimiter, index=False)

    def to_tsv(self, *args, **kwargs):
        return partial(self.to_csv, delimiter='\t')(*args, **kwargs)

    @classmethod
    def from_tsv(
        cls, fpath, sample_metadata_fpath=None, obs_metadata_fpath=None,
        taxonomy_fpath=None, gzipped=False, sample_rows=True
    ):
        dummy_fun = lambda x: x
        otus_with_taxonomy = set()
        if taxonomy_fpath is not None:
            with open(taxonomy_fpath, 'rt') as handle:
                reader = csv.reader(handle, delimiter='\t')
                next(reader)
                for otu_id, *_ in reader:
                    otus_with_taxonomy.add(otu_id)
        with cls.create_opener(gzipped)(fpath) as lines:
            table = biom.Table.from_tsv(lines, None, None, dummy_fun)  # noqa
            if sample_rows:
                table = table.transpose()
            if sample_metadata_fpath is not None:
                sample_metadata = cls.read_metadata(sample_metadata_fpath)
                table.add_metadata(sample_metadata, axis='sample')  # noqa
            if obs_metadata_fpath is not None:
                obs_metadata = cls.read_metadata(sample_metadata_fpath)
                table.add_metadata(obs_metadata, axis='observation')  # noqa
            return cls(table, otus_with_taxonomy)

    def to_hdf5(self, filepath):
        with h5py.File(filepath, 'w') as file:
            self.table.to_hdf5(file, generated_by='unknown')

    @classmethod
    def from_hdf5(cls, filepath):
        with biom_open(filepath) as file:
            return cls(biom.Table.from_hdf5(file))


if __name__ == '__main__':
    table = OTUTable.from_tsv(
        '/Users/ant/masters-thesis/amplicon_occurences.tsv',
        '/Users/ant/masters-thesis/amplicon_station_metadata.tsv'
    )
    print(table._sample_metadata)
    table.to_hdf5('/Users/ant/masters-thesis/amplicon_station_test.biom')
    table = OTUTable.from_hdf5('/Users/ant/masters-thesis/amplicon_station_test.biom')
    print(table.table.head())
    print(table._sample_metadata)
