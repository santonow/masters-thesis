import csv
import gzip
import sys
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
        try:
            record = SearchIO.read(StringIO(item), 'blast-tab')
        except ValueError:
            print(item)
            continue
        best_lineage = determine_best_lineage(record.hsps)
        if best_lineage is not None and best_lineage[0] == 'Eukaryota':
            output_queue.put([record.id] + best_lineage)


def writer(output_queue: Queue, output_fname: str):
    with open(output_fname, 'w') as handle:
        csv_writer = csv.writer(handle, delimiter='\t')
        csv_writer.writerow(['otu_id'] + RANKS)
        while True:
            item = output_queue.get()
            if item is None:
                break
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
    if len(hsps) < 10:
        return None
    lineages = [
        get_lineage_from_description(hsp.hit_id)
        for hsp in sorted(
            hsps, key=lambda x: x.bitscore, reverse=True
        )[:first_n]
    ]
    i = 8
    while not all(lineage[i-1] == lineages[0][i-1] for lineage in lineages) and i > 0:
        i -= 1
    if i <= 0:
        return None
    return lineages[0][:i]


def prepare_tax_table(in_file, n_threads, out_file):
    in_queue = Queue(n_threads)
    out_queue = Queue(n_threads)

    processes = []
    for _ in range(min(1, n_threads - 1)):
        process = Process(
            target=reader, args=(in_queue, out_queue)
        )
        process.start()
        processes.append(process)
    writer_process = Process(target=writer, args=(out_queue, out_file))
    writer_process.start()
    processes.append(writer_process)

    for lines in tqdm(lines_reader(in_file), total=get_n_otus(in_file)):
        in_queue.put(lines)
    for _ in range(min(1, n_threads - 1)):
        in_queue.put(None)

    for process in processes:
        process.join()


if __name__ == '__main__':
    prepare_tax_table(sys.argv[1], int(sys.argv[2]), sys.argv[3])
