import argparse
import networkx as nx
import sys

TC_COUNT='tc_count'

def __rec_process_root(G, root):
    G.node[root][TC_COUNT] = 0
    for s in G.neighbors_iter(root):
        if G.out_degree(s) == 0:
            G.node[s][TC_COUNT] = 0
        else:
            __rec_process_root(G, s)
        G.node[root][TC_COUNT] += G.node[s][TC_COUNT] + 1

def __find_root(G):
    for n in G:
        if G.in_degree(n) == 0:
            return n
    return None

def count_transitive_closure(G):
    root = __find_root(G)

    __rec_process_root(G, root)

    N=G.number_of_nodes()

    return float(sum(G.node[n][TC_COUNT] for n in G)) / (N)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    args = parser.parse_args()

    one = nx.read_pajek(args.input)

    sys.stdout.write("%.5e" % count_transitive_closure(one))

