import csv
import gzip
import sys
from collections import Counter
from typing import Optional, List
from io import StringIO

from Bio import SearchIO
from tqdm import tqdm
from multiprocessing import Process, Queue

RANKS = [
    'kingdom',
    'supergroup',
    'division',
    'class',
    'order',
    'family',
    'genus',
    'species'
]


def lines_reader(fname):
    with gzip.open(fname, 'rt') as handle:
        curr_id = ''
        current_lines = []
        for line in handle:
            _id = line.split('\t')[0]
            if curr_id == '':
                curr_id = _id
            if _id != curr_id:
                yield '\n'.join(current_lines)
                current_lines = [line.strip()]
                curr_id = _id
            else:
                current_lines.append(line.strip())
        if current_lines:
            yield '\n'.join(current_lines)


def reader(input_queue: Queue, output_queue: Queue):
    while True:
        item = input_queue.get()
        if item is None:
            output_queue.put(None)
            break
        record = SearchIO.read(StringIO(item), 'blast-tab')
        best_lineage = determine_best_lineage(record.hsps)
        if best_lineage is not None and best_lineage[0] == 'Eukaryota':
            output_queue.put([record.id] + best_lineage)


def writer(output_queue: Queue, output_fname: str, threads: int):
    sentinels_encountered = 0
    with open(output_fname, 'w') as handle:
        csv_writer = csv.writer(handle, delimiter='\t')
        csv_writer.writerow(['otu_id'] + RANKS)
        while True:
            item = output_queue.get()
            if item is None:
                sentinels_encountered += 1
                if sentinels_encountered == threads:
                    break
            else:
                csv_writer.writerow(item)


def get_n_otus(fname):
    otus = set()
    with gzip.open(fname, 'rt') as handle:
        for line in handle:
            otus.add(line.split('\t')[0])
    return len(otus)


def get_lineage_from_description(desc: str) -> List[str]:
    return desc.split('|')[-8:]


def determine_best_lineage(hsps: List[SearchIO.HSP], first_n: int = 10) -> Optional[List[str]]:
    if len(hsps) < first_n:
        return None
    lineages = [
        get_lineage_from_description(hsp.hit_id)
        for hsp in sorted(
            hsps, key=lambda x: x.bitscore, reverse=True
        )[:first_n]
    ]
    i = 8
    lineage_counter = Counter(tuple(lineage) for lineage in lineages)
    most_popular_lineage, count = lineage_counter.most_common(1)[0]
    while not count > first_n / 2:
        i -= 1
        lineage_counter = Counter(tuple(lineage[:i]) for lineage in lineages)
        most_popular_lineage, count = lineage_counter.most_common(1)[0]
    if i <= 0:
        return None
    return list(most_popular_lineage)


def prepare_tax_table(in_file, n_threads, out_file):
    threads = n_threads - 1 if n_threads > 1 else 1
    in_queue = Queue(threads)
    out_queue = Queue(threads)

    processes = []
    for _ in range(threads):
        process = Process(
            target=reader, args=(in_queue, out_queue)
        )
        process.start()
        processes.append(process)
    writer_process = Process(target=writer, args=(out_queue, out_file, threads))
    writer_process.start()
    processes.append(writer_process)

    for lines in tqdm(lines_reader(in_file), total=get_n_otus(in_file)):
        in_queue.put(lines)
    for _ in range(threads):
        in_queue.put(None)

    for process in processes:
        process.join()


if __name__ == '__main__':
    prepare_tax_table(sys.argv[1], int(sys.argv[2]), sys.argv[3])
