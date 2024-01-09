import os
from glob import glob
import networkx as nx
from itertools import product
from argparse import ArgumentParser
from copy import deepcopy
from collections import defaultdict
#from sys import exit as sys_exit
#import matplotlib.pyplot as plt

from typing import Optional
from io import TextIOWrapper


STANDARD_OUT_STR = 'standard_out'
RESULTS_FILE_STR = 'results_file'
RESULTS_DIR = 'results'
PERM = tuple[int]
NODE = int
EDGE = tuple[NODE]
EDGE_TO_ID_MAP = dict[EDGE, NODE]
ID_TO_EDGE_MAP = dict[NODE, EDGE]
DIV_LENGTH = 60

def log_info(
    mssg: str,
    res_file: Optional[TextIOWrapper] = None,
    out: Optional[str] = 'all'
) -> None:
    ALL_STR = 'all'
    if out in (STANDARD_OUT_STR, ALL_STR):
        print(mssg)
    if out in (RESULTS_FILE_STR, ALL_STR):
        res_file.write(mssg + '\n')
    return

def get_perm(perm_path: str) -> PERM:
    with open(perm_path) as f:
        number_line = f.read()
    number_line = number_line[len('The permutation ('):].split(')')[0]
    number_line = tuple(int(num) for num in number_line.split(', '))
    return number_line

def rotation_scheme(seq: list[tuple[NODE]]) -> list[tuple[NODE]]:
    rotations = []
    start_tup = seq.pop(0)
    current = [*start_tup]
    for _ in range(len(seq)):
        last = current[-1]
        for x, y in seq:
            other = (x == last) * y + (y == last) * x
            if other:
                seq.remove((x, y))
                if other == current[0]:
                    rotations.append(tuple(current))
                    if seq:
                        start_tup = seq.pop(0)
                        current = [*start_tup]
                else:
                    current.append(other)
                break
    return rotations

class Graph():
    def __init__(
        self,
        file_path: str,
        res_file: Optional[TextIOWrapper] = None
    ) -> None:
        self.graph = nx.read_edgelist(file_path, nodetype=int)
        self.degree_six_nodes = self.get_degree_six_nodes()
        self.edge_id_map = self.get_edge_id_map()
        self.id_edge_map = self.get_id_edge_map()

        self.res_file = res_file

    def graph_from_perm(perm_path: str) -> nx.classes.graph.Graph:
        gf = os.path.splitext(os.path.split(perm_path)[-1])[0].split('-')[0]
        return nx.read_edgelist(gf, nodetype=int)

    def get_degree_six_nodes(self) -> tuple[NODE]:
        return tuple(
            node for node in self.graph.nodes if self.graph.degree(node) == 6
        )

    def get_edge_id_map(self) -> EDGE_TO_ID_MAP:
        return {edge: id for id, edge in enumerate(self.graph.edges, start=1)}

    def get_id_edge_map(self) -> ID_TO_EDGE_MAP:
        return {id: edge for id, edge in enumerate(self.graph.edges, start=1)}

    def get_vstar(self, node: NODE) -> tuple[EDGE]:
        return tuple((min(edge), max(edge)) for edge in self.graph.edges(node))

    def get_dual_edges(self, vstar: tuple[EDGE], perm: PERM) -> tuple[EDGE]:
        return tuple(
            self.id_edge_map[perm[self.edge_id_map[edge] - 1]]
            for edge in vstar
        )

    def get_dual_face(self, node: NODE, perm: PERM) -> nx.classes.graph.Graph:
        return nx.Graph(self.get_dual_edges(self.get_vstar(node), perm))

    def check_self_dual_perm(self, perm: PERM) -> bool:
        for node in self.graph.nodes:
            if not nx.is_eulerian(self.get_dual_face(node, perm)):
                #if perm == (11, 9, 13, 6, 12, 3, 10, 2, 14, 7, 1, 4, 5, 8):
                #    nx.draw(self.get_dual_face(node, perm), with_labels=True)
                #    plt.show()
                log_info(
                    f'perm(star({node})) was not eulerian',
                    self.res_file
                )
                return False
        return True

    def n_inverse_edges_components(self, degree_six: NODE, perm: PERM) -> int:
        inverse_edges = tuple(
            self.id_edge_map[perm.index(self.edge_id_map[edge]) + 1]
            for edge in self.get_vstar(degree_six)
        )
        return nx.number_connected_components(nx.Graph(inverse_edges))

    def check_cycle(self, node: NODE, perm: PERM) -> bool:
        dual_face = self.get_dual_face(node, perm)
        node_degrees = dict(dual_face.degree(dual_face.nodes))
        for degree in node_degrees.values():
            if degree != 2:
                return False
        return True

    def get_dual_bowtie_edge_options(
        self,
        source_node: NODE,
        perm: PERM
    ) -> list[list[EDGE]]:
        without_cross = self.get_dual_face(source_node, perm)
        for node in without_cross.nodes:
            if without_cross.degree(node) == 4:
                cross_node = node
                break
        without_cross.remove_node(cross_node)
        outer_nodes = tuple(without_cross.nodes)

        # alpha and delta are neighbors & beta and gamma are neighbors
        alpha = outer_nodes[0]
        beta = None
        gamma = None
        delta = None
        for node in outer_nodes[1:]:
            if without_cross.has_edge(alpha, node):
                delta = node
            elif not beta:
                beta = node
            else:
                gamma = node

        return [
            [(alpha, beta), (gamma, delta)],
            [(alpha, gamma), (beta, delta)]
        ]

    def get_pass_sequences(
        self,
        degree_six_node: NODE,
        perm: PERM
    ) -> list[list[EDGE]]:
        pass_seqs = [[]]
        for node in self.graph.nodes:
            dual_face = self.get_dual_face(node, perm)
            if degree_six_node in dual_face.nodes:
                if len(tuple(dual_face.neighbors(degree_six_node))) == 4:
                    options = self.get_dual_bowtie_edge_options(node, perm)
                    if pass_seqs == [[]]:
                        pass_seqs = options
                    else:
                        pass_seqs.extend(deepcopy(pass_seqs))
                        for option_ind, option in enumerate(options):
                            pass_seqs[option_ind].extend(option)
                else:
                    for pass_seq in pass_seqs:
                        pass_seq.append(
                            tuple(dual_face.neighbors(degree_six_node))
                        )
        log_info(
            f'Sequences found for node {degree_six_node}: {len(pass_seqs)}',
            self.res_file
        )
        if len(pass_seqs) > 2:
            log_info(
                f'Found n_pass_seqs > 2 for node {degree_six_node}',
                self.res_file
            )
        return pass_seqs

    def check_pinch_point(
        self,
        degree_six_node: NODE,
        perm: PERM
    ) -> bool:
        seqs = self.get_pass_sequences(degree_six_node, perm)
        for passes in seqs:
            log_info(f'Using pass sequence: {passes}', self.res_file)
            rs = rotation_scheme(passes)
            log_info(f'Rotation Scheme is: {rs}', self.res_file)
            if len(rs) > 1:
                return True

        log_info(
            'None of the bowtie uses produced a pinch point',
            self.res_file
        )
        return False

    def get_n_pinch_points(self, perm: PERM) -> int:
        return sum(
            1 for node in self.degree_six_nodes
            if self.check_pinch_point(node, perm)
        )

    def get_dual_bowtie_nodes(self, perm: PERM) -> dict[NODE, NODE]:
        dual_bowtie_nodes = {}
        for degree_six_node in self.degree_six_nodes:
            dual_face = self.get_dual_face(degree_six_node, perm)
            for df_node in dual_face.nodes:
                if len(tuple(dual_face.neighbors(df_node))) == 4:
                    dual_bowtie_nodes[degree_six_node] = df_node
        return dual_bowtie_nodes

    def check_orientable(self, perm: PERM) -> bool:
        dual_bowtie_nodes = self.get_dual_bowtie_nodes(perm)
        n_dual_bowtie_nodes = len(dual_bowtie_nodes)
        if n_dual_bowtie_nodes == 0:
            options = (None,)
        elif n_dual_bowtie_nodes == 1:
            options = (0, 1)
        elif n_dual_bowtie_nodes == 2:
            options = tuple(product(range(2), repeat=2))
        dual_faces = [
            self.get_dual_face(node, perm) for node in self.graph.nodes
        ]

        def get_option_instruction(node, option) -> int:
            return (
                option if n_dual_bowtie_nodes == 1
                else option[list(dual_bowtie_nodes).index(node)]
            )

        def get_bowtie_walk(center_node, bowtie_passes) -> list[EDGE]:
            return [
                (bowtie_passes[0][0], center_node),
                (center_node, bowtie_passes[0][1]),
                (bowtie_passes[0][1], bowtie_passes[1][0]),
                (bowtie_passes[1][0], center_node),
                (center_node, bowtie_passes[1][1]),
                (bowtie_passes[1][1], bowtie_passes[0][0])
            ]

        def get_cycle_walk(edges) -> list[EDGE]:
            rs = list(rotation_scheme(edges)[0])
            return [tup for tup in zip(rs, rs[1:] + [rs[0]])]

        def get_oriented_walk(edge_count, walk) -> tuple[EDGE]:
            for edge in walk:
                if edge in edge_count:
                    walk = tuple(tuple(reversed(edge)) for edge in walk)
                    break
            return walk

        for option in options:
            edge_count = defaultdict(lambda: 0)
            for node, dual_face in zip(self.graph.nodes, dual_faces):
                if node in dual_bowtie_nodes:
                    bowtie_passes = self.get_dual_bowtie_edge_options(
                        node,
                        perm
                    )[get_option_instruction(node, option)]
                    walk = get_bowtie_walk(
                        dual_bowtie_nodes[node],
                        bowtie_passes
                    )
                else:
                    walk = get_cycle_walk(list(dual_face.edges))
                if not edge_count:
                    for edge in walk:
                        edge_count[edge] = 1
                else:
                    walk = get_oriented_walk(edge_count, walk)
                    for n1, n2 in walk:
                        if (n2, n1) in edge_count:
                            edge_count[(n2, n1)] -= 1
                        else:
                            edge_count[(n1, n2)] += 1
            if all(count == 0 for count in edge_count.values()):
                return True
        return False

def process_perm(perm_path: str, graph: Graph) -> tuple[bool]:
    """ Return: (is_ADC_bool, is_solution_bool) """
    perm = get_perm(perm_path)

    log_info('', graph.res_file)
    log_info('-' * DIV_LENGTH, graph.res_file)
    log_info(f'Perm: {perm}', graph.res_file)
    log_info('-' * DIV_LENGTH, graph.res_file)

    if not graph.check_self_dual_perm(perm):
        log_info(f'Perm is not an ADC', graph.res_file)
        return (False, False)

    n_pinches = graph.get_n_pinch_points(perm)
    is_orientable = graph.check_orientable(perm)
    log_info(f'Is orientable: {is_orientable}', graph.res_file)
    if n_pinches == 1:
        log_info('Solution Found: Pinched Projective Plane', graph.res_file)
        return (True, True)

    if n_pinches == 2:
        log_info('Solution Found: 2-pinch-point Sphere', graph.res_file)
        return (True, True)

    if n_pinches == 0:
        if is_orientable:
            log_info('Solution Found: Torus', graph.res_file)
            return (True, True)
        else:
            log_info('Solution Found: Klien Bottle', graph.res_file)
            return (True, True)

    log_info('No Solution Found', graph.res_file)
    return (True, False)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument('perm_pattern')
    parser.add_argument('-m', dest='max_perms', type=int, default=None)
    parser.add_argument('-r', dest='refresh_results', action='store_true')
    args = parser.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    if args.refresh_results:
        for path in glob(os.path.join(RESULTS_DIR, '*')):
            os.remove(path)

    results_file = os.path.join(
        RESULTS_DIR,
        f'{os.path.split(args.graph_file)[-1]}-results.txt'
    )

    perm_files = sorted(glob(args.perm_pattern))[:args.max_perms]

    with open(results_file, 'w') as res_file_wrapper:
        with open(args.graph_file) as graph_file:
            graph = Graph(graph_file, res_file_wrapper)

        graph_edge_str = ''
        for i, edge in enumerate(graph.graph.edges):
            graph_edge_str += f'Edge {i + 1}: {edge}\n'

        log_info('Graph Edges:', res_file_wrapper)
        log_info(graph_edge_str, res_file_wrapper)

        if not graph.degree_six_nodes:
            log_info('No degree six nodes', res_file_wrapper)
            return
        log_info(
            f'Degree six nodes are: {graph.degree_six_nodes}',
            res_file_wrapper
        )

        n_adc = 0
        n_solution = 0
        for perm_file in perm_files:
            perm_solution = process_perm(perm_file, graph)
            if perm_solution[0]:
                n_adc += 1
                if perm_solution[1]:
                    n_solution += 1

        log_info('', res_file_wrapper)
        log_info('-' * DIV_LENGTH, res_file_wrapper)
        log_info(f'Summary', res_file_wrapper)
        log_info('-' * DIV_LENGTH, res_file_wrapper)
        log_info(f'Permutations checked: {len(perm_files)}', res_file_wrapper)
        log_info(f'ADCs found: {n_adc}', res_file_wrapper)
        log_info(f'Solutions found: {n_solution}', res_file_wrapper)

    return

if __name__ == '__main__':
    main()
