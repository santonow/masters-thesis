import os


def make_graph_name(filepath):
    dirname, fname = os.path.split(filepath)
    fname = os.path.splitext(fname)[0] + ".edgelist"
    return os.path.join(dirname, "standardized_" + fname)
