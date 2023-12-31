import os
from glob import glob
import networkx as nx
from argparse import ArgumentParser
from sys import exit as sys_exit

from typing import Optional


TEXT_DIR = 'perm_txts'
PERM = tuple[int]
NODE = int
EDGE = tuple[NODE]
EDGE_TO_ID_MAP = dict[EDGE, NODE]
ID_TO_EDGE_MAP = dict[NODE, EDGE]


def get_perm(perm_path: str) -> PERM:
    with open(file_path) as f:
        number_line = f.read()
    number_line = number_line[len('The permutation ('):].split(')')[0]
    number_line = tuple(int(num) for num in number_line.split(', '))
    return number_line


class Graph():
    def __init__(file_path: str, file_type: str) -> None:
        if file_type == 'graph':
            self.graph = nx.read_edgelist(file_path, nodetype=int)
        elif file_type == 'perm':
            self.graph = self.graph_from_perm(file_path)
        else:
            sys_exit("Bad file_type. Use 'graph' or 'perm'")

        self.degree_six_nodes = get_degree_sixes()
        self.edge_id_map = get_edge_id_map()
        self.id_edge_map = get_id_edge_map()

    def graph_from_perm(perm_path: str) -> nx.classes.graph.Graph:
        gf = os.path.splitext(os.path.split(perm_path)[-1])[0].split('-')[0]
        return nx.read_edgelist(gf, nodetype=int)

    def get_degree_sixes(self) -> tuple[NODE]:
        return tuple(
            node for node in self.graph.nodes if self.graph.degree(node) == 6
        )

    def get_edge_id_map(self) -> EDGE_TO_ID_MAP:
        return {edge: id for id, edge in enumerate(self.graph.edges)}

    def get_id_edge_map(self) -> ID_TO_EDGE_MAP:
        return {id: edge for id, edge in enumerate(self.graph.edges)}

    def get_vstar(self, node: NODE) -> tuple[EDGE]:
        return tuple((min(edge), max(edge)) for edge in self.graph.edges(node))

    def get_dual_edges(
            self,
            vstar: tuple[EDGE],
            perm: PERM
        ) -> tuple[EDGE]:
        return tuple(
            self.id_edge_map[perm[self.edge_id_map[edge] - 1]]
            for edge in vstar
        )

    def get_dual_face(self, node: NODE, perm: PERM) -> nx.classes.graph.Graph
        return nx.Graph(
            self.get_dual_edges(self.get_vstar(self.graph, node), perm)
        )

    def check_self_dual_perm(perm: PERM): -> bool:
        for node in self.graph.nodes:
            if not nx.is_eulerian(self.get_dual_face(node, perm)):
                return False
        return True

    def check_components(self, perm: PERM) -> bool:
        def inverse_perm(self, perm: PERM) -> PERM:
            return tuple()


    def check_cycle(self, node: NODE, perm: PERM) -> bool:
        dual_face = self.get_dual_face(node, perm)
        node_degrees = dict(dual_face.degree(dual_face.nodes))
        for degree in node_degrees.values():
            if degree != 2:
                return False
        return True

    def check_bowtie(self, node: NODE, perm: PERM) -> bool:
        ...

def process_perm(
        perm_path: str,
        graph: Graph,
    ) -> None:
    perm = get_perm(perm_path)

    if not check_self_dual(graph, perm, edge_id_map, id_edge_map):
        return

    if not degree_sixes:
        return

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

    graph = Graph(args.graph_file, 'graph') if args.graph_file else None

    for perm_file in perm_files:
        if not graph:
            graph = Graph(perm_file, 'perm')
        process_perm(perm_file, graph)

if __name__ == '__main__':
    main()
