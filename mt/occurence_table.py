import csv
import gzip
import os
from functools import partial

import biom
import h5py


class OTUTable:

    def __init__(self, table: biom.Table):
        self.table = table

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
            for line in reader:
                metadata[line[0]] = dict(zip(header[1:], line[1:]))
        return metadata

    def to_csv(self, filepath, delimiter=',', gzipped=False, metadata_suffix='meta', transpose=False):
        with self.create_writer(gzipped)(filepath) as writer:
            if transpose:
                tbl = self.table.transpose()
            else:
                tbl = self.table
            for line in tbl.delimited_self(delimiter).split('\n')[1:]:
                writer.write(line + '\n')
        if self.table._sample_metadata is not None:
            pref, fname = os.path.split(filepath)
            name, extensions = fname.split('.', 1)
            metadata_name = name + '_' + metadata_suffix
            if extensions.endswith('.gz'):
                extensions = extensions[:-len('.gz')]
            metadata_filepath = os.path.join(pref, metadata_name + '.' + extensions)
            self.table.metadata_to_dataframe(axis='sample').to_csv(metadata_filepath, sep=delimiter)

    def to_tsv(self, *args, **kwargs):
        return partial(self.to_csv, delimiter='\t')(*args, **kwargs)

    @classmethod
    def from_tsv(cls, fpath, sample_metadata_fpath=None, obs_metadata_fpath=None, gzipped=False, sample_rows=True):
        dummy_fun = lambda x: x
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
            return cls(table)

    def to_hdf5(self, filepath):
        with h5py.File(filepath, 'w') as file:
            self.table.to_hdf5(file, generated_by='unknown')

    # TODO: figure this out
    @classmethod
    def from_hdf5(cls, filepath):
        with h5py.File(filepath, 'r') as file:
            return cls(biom.Table.from_hdf5(file))


if __name__ == '__main__':
    table = OTUTable.from_tsv(
        '/Users/ant/masters-thesis/amplicon_occurences.tsv',
        '/Users/ant/masters-thesis/amplicon_station_metadata.tsv'
    )
    print(table.table.head())
    print(table.table.metadata(id='TARA_A100000393', axis='sample'))
    table.to_tsv('/Users/ant/masters-thesis/amplicon_station_test.tsv')
    table.to_csv('/Users/ant/masters-thesis/amplicon_station_test.csv', delimiter=';')
    table.to_tsv('/Users/ant/masters-thesis/amplicon_station_test.tsv.gz')
    table.to_hdf5('/Users/ant/masters-thesis/amplicon_station_test.biom')
    table = OTUTable.from_hdf5('/Users/ant/masters-thesis/amplicon_station_test.biom')
    print(table.table.head())
    print(table.table.metadata(id='TARA_A100000393', axis='sample'))
