import os
from glob import glob
import networkx as nx
from argparse import ArgumentParser

from typing import Optional

TEXT_DIR = 'perm_txts'
EDGE_TO_ID_MAP = dict[tuple[int], int]
ID_TO_EDGE_MAP = dict[int, tuple[int]]

def perm_from_txt(perm_path: str) -> tuple[int]:
    with open(file_path) as f:
        number_line = f.read()
    number_line = number_line[len('The permutation ('):].split(')')[0]
    number_line = tuple(int(num) for num in number_line.split(', '))
    return number_line

def graph_from_txt(perm_path: str) -> nx.classes.graph.Graph:
    gf = os.path.splitext(os.path.split(perm_path)[-1])[0].split('-')[0]
    return nx.read_edgelist(gf, nodetype=int)

def vstar(graph: nx.classes.graph.Graph, v) -> tuple[tuple[int]]:
    return tuple((min(edge), max(edge)) for edge in graph.edges(v))

def dual_edges(
        vstar: tuple[tuple[int]],
        perm: tuple[int],
        edge_ids: EDGE_TO_ID_MAP,
        id_edges: ID_TO_EDGE_MAP,
    ) -> tuple[tuple[int]]:
    return tuple((id_edges[perm[edge_ids[edge] - 1]] for edge in vstar))

def boundry_walk():
    ...

def get_n_umbrealls(boundry_walks: ...) -> tuple[int]:
    ...

def process_perm(
        perm_path: str,
        graph: Optional[nx.classes.graph.Graph] = None,
        edge_ids: Optional[EDGE_TO_ID_MAP] = None,
        id_edges: Optional[ID_TO_EDGE_MAP] = None
    ) -> None:
    perm = perm_from_txt(perm_path)
    if not graph:
        graph = graph_from_txt(perm_path)
        edge_ids = {edge: id for id, edge in enumerate(graph.edges)}
        id_edges = {id: edge for id, edge in enumerate(graph.edges)}

    boundry_walks = get_boundry_walks(...)
    n_umbrellas = get_n_umbrellas(boundry_walks)
    print(n_umbrellas)

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('-g', dest='graph_file')
    parser.add_argument('-m', dest='max_perms', type=int, default=None)
    args = parser.parse_args()

    perm_pattern = f'{args.graph_file}-*' if args.graph_file else '*'
    perm_files = sorted(glob(
        os.path.join(TEXT_DIR, perm_pattern)
    ))[:args.max_perms]

    if args.graph_file:
        graph = nx.read_edgelist(args.graph_file, nodetype=int)
        edge_ids = {edge: id for id, edge in enumerate(graph.edges)}
        id_edges = {id: edge for id, edge in enumerate(graph.edges)}
    else:
        graph = None
        edge_ids = None
        id_edges = None

    for perm_file in perm_files:
        process_perm(perm_file, graph, edge_ids, id_edges)

if __name__ == '__main__':
    main()
