from __future__ import annotations

import csv
import gzip
import os
from typing import Set, Protocol, List, Optional
from functools import partial
from abc import abstractmethod

import biom
from biom.util import biom_open
import h5py
import numpy as np


class TableFilter(Protocol):
    @abstractmethod
    def filter(self, table: biom.Table, current_ids: Set[str]) -> Set[str]:
        raise NotImplementedError


class TaxonomyFilter(TableFilter):
    def __init__(self, otus_with_lineage: Set[str]):
        self.otus_with_lineage = otus_with_lineage

    def filter(self, table: biom.Table, current_ids: Set[str]) -> Set[str]:
        return current_ids & self.otus_with_lineage


class MinTotalAbundanceFilter(TableFilter):
    def __init__(self, min_total_abundance: int):
        self.min_total_abundance = min_total_abundance

    def filter(self, table: biom.Table, current_ids: Set[str]) -> Set[str]:
        return current_ids & {
            _id
            for _id in table.ids(axis="observation")
            if sum(table.data(_id, axis="observation")) > self.min_total_abundance
        }


class RelativeAbundanceFilter(TableFilter):
    def __init__(self, relative_abundance: float):
        self.relative_abundance = relative_abundance

    def filter(self, table: biom.Table, current_ids: Set[str]) -> Set[str]:
        tbl_sum = table.sum(axis="whole")
        return current_ids & {
            _id
            for _id in table.ids(axis="observation")
            if (sum(table.data(_id, axis="observation")) / tbl_sum)
            >= self.relative_abundance
        }


class ZeroPropFilter(TableFilter):
    def __init__(self, zero_threshold: float):
        self.zero_threshold = zero_threshold

    def filter(self, table: biom.Table, current_ids: Set[str]) -> Set[str]:
        n_samples = table.length(axis="sample")
        return current_ids & {
            _id
            for _id in table.ids(axis="observation")
            if (list(table.data(_id, axis="observation")).count(0) / n_samples)
            < 1 - self.zero_threshold
        }


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
            return partial(gzip.open, mode="rt")
        else:
            return partial(open, mode="r")

    @staticmethod
    def create_writer(gzipped=False):
        if gzipped:
            return partial(gzip.open, mode="wt")
        else:
            return partial(open, mode="w")

    @staticmethod
    def read_metadata(fname: str):
        metadata = dict()
        with open(fname) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            header = [val.replace("/", " per ") for val in header]
            for line in reader:
                metadata[line[0]] = dict(zip(header[1:], line[1:]))
        return metadata

    @property
    def matrix(self):
        return np.array(self.table.matrix_data.todense()).astype(np.float64)

    def filter(
        self,
        table: biom.Table,
        min_total_abundance: Optional[int] = None,
        relative_abundance: Optional[float] = 1e-8,
        zero_threshold: float = 0.03,
    ):
        filters: List[TableFilter] = [
            TaxonomyFilter(self.otus_with_lineage),
        ]
        if min_total_abundance is not None:
            filters.append(MinTotalAbundanceFilter(min_total_abundance))
        if relative_abundance is not None:
            filters.append(RelativeAbundanceFilter(relative_abundance))
        if zero_threshold is not None:
            filters.append(ZeroPropFilter(zero_threshold))

        filtered_by_each = dict()
        remaining_ids = set(table.ids(axis="observation"))
        for filter_method in filters:
            remaining_ids = filter_method.filter(table, remaining_ids)
            filtered_by_each[filter_method.__class__.__name__] = len(
                set(table.ids(axis="observation")) - remaining_ids
            ) - sum(filtered_by_each.values())

        n_obs = len(table.ids(axis="observation"))
        filtered_table = table.filter(
            ids_to_keep=remaining_ids, invert=True, axis="observation", inplace=False
        )
        leftover_vector = ["leftover_vector"] + [
            x for x in filtered_table.sum(axis="sample")
        ]
        s = f"Filtering {n_obs - len(remaining_ids)} / {n_obs} OTUs.\n"
        for filter_method, count in filtered_by_each.items():
            s += f"Filter {count} using filtering method {filter_method}.\n"
        print(s)
        return (
            table.filter(ids_to_keep=remaining_ids, axis="observation", inplace=False),
            leftover_vector,
        )

    def to_csv(
        self,
        filepath,
        delimiter=",",
        gzipped=False,
        metadata_suffix="meta",
        transpose=False,
        min_tot_abundance=None,
        relative_abundance=1e-8,
        zero_threshold=0.03,
    ):
        if transpose:
            table: biom.Table = self.table.transpose()
        else:
            table = self.table
        filtered_table, leftover_vec = self.filter(
            table,
            min_total_abundance=min_tot_abundance,
            relative_abundance=relative_abundance,
            zero_threshold=zero_threshold,
        )
        with self.create_writer(gzipped)(filepath) as writer:

            for line in filtered_table.delimited_self(delimiter).split("\n")[1:]:
                writer.write(line + "\n")
            writer.write("\t".join(str(x) for x in leftover_vec))

        if self.table._sample_metadata is not None:
            pref, fname = os.path.split(filepath)
            name, extensions = fname.split(".", 1)
            metadata_name = name + "_" + metadata_suffix
            if extensions.endswith(".gz"):
                extensions = extensions[: -len(".gz")]
            metadata_filepath = os.path.join(pref, metadata_name + "." + extensions)
            self.table.metadata_to_dataframe(axis="sample").to_csv(
                metadata_filepath, sep=delimiter, index=False
            )

    def to_tsv(self, *args, **kwargs):
        return partial(self.to_csv, delimiter="\t")(*args, **kwargs)

    @classmethod
    def from_tsv(
        cls,
        fpath,
        sample_metadata_fpath=None,
        obs_metadata_fpath=None,
        taxonomy_fpath=None,
        gzipped=False,
        sample_rows=True,
    ):
        dummy_fun = lambda x: x
        otus_with_taxonomy = set()
        if taxonomy_fpath is not None:
            with open(taxonomy_fpath, "rt") as handle:
                reader = csv.reader(handle, delimiter="\t")
                next(reader)
                for otu_id, *_ in reader:
                    otus_with_taxonomy.add(otu_id)
        with cls.create_opener(gzipped)(fpath) as lines:
            table = biom.Table.from_tsv(lines, None, None, dummy_fun)  # noqa
            if sample_rows:
                table = table.transpose()
            if sample_metadata_fpath is not None:
                sample_metadata = cls.read_metadata(sample_metadata_fpath)
                table.add_metadata(sample_metadata, axis="sample")  # noqa
            if obs_metadata_fpath is not None:
                obs_metadata = cls.read_metadata(sample_metadata_fpath)
                table.add_metadata(obs_metadata, axis="observation")  # noqa
            return cls(table, otus_with_taxonomy)

    def to_hdf5(self, filepath):
        with h5py.File(filepath, "w") as file:
            self.table.to_hdf5(file, generated_by="unknown")

    @classmethod
    def from_hdf5(cls, filepath):
        with biom_open(filepath) as file:
            return cls(biom.Table.from_hdf5(file))


if __name__ == "__main__":
    table = OTUTable.from_tsv(
        "/Users/ant/masters-thesis/amplicon_occurences.tsv",
        "/Users/ant/masters-thesis/amplicon_station_metadata.tsv",
    )
    print(table._sample_metadata)
    table.to_hdf5("/Users/ant/masters-thesis/amplicon_station_test.biom")
    table = OTUTable.from_hdf5("/amplicon_station_test.biom")
    print(table.table.head())
    print(table._sample_metadata)
