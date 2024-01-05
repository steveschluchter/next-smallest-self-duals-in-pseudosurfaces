import os
from glob import glob
import networkx as nx
from argparse import ArgumentParser
#from sys import exit as sys_exit
import matplotlib.pyplot as plt

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

def rotation_scheme(seq: list[tuple[NODE]]) -> list[list[NODE]]:
    rotations = []
    current = [*seq[0]]
    seq.pop(0)
    for _ in range(len(seq)):
        last = current[-1]
        for x, y in seq:
            other = (x == last) * y + (y == last) * x
            if other:
                seq.remove((x, y))
                if other == current[0]:
                    rotations.append(current)
                    if seq:
                        current = [*seq[0]]
                        seq.pop(0)
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

    def check_bowtie_for_pinch_point(
        self,
        degree_six_node: NODE,
        perm: PERM
    ) -> bool:
        log_info(
            f'perm(star(degree-six node: {degree_six_node})) is a bowtie',
            self.res_file
        )
        without_cross = self.get_dual_face(degree_six_node, perm)
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

        passes = [(alpha, beta), (gamma, delta)]
        for node in self.graph.nodes:
            if node == degree_six_node:
                continue
            dual_face = self.get_dual_face(node, perm)
            for dual_node in dual_face.nodes:
                if dual_node == degree_six_node:
                    passes.append(tuple(dual_face.neighbors(dual_node)))

        def check_passes(current_passes):
            log_info(f'Using pass sequence: {current_passes}', self.res_file)
            rs = rotation_scheme(current_passes)
            log_info(f'Rotation Scheme is: {rs}', self.res_file)
            if len(rs) > 1:
                return True

        if check_passes(passes.copy()):
            return True

        # Other use of bowtie edges
        passes = [(alpha, gamma), (beta, delta), *(passes[2:])]
        if check_passes(passes):
            return True

        log_info(
            'Neither use of the bowtie produced a pinch point',
            self.res_file
        )
        return False

    def get_n_pinch_points(self, perm: PERM) -> int:
        npp = 0
        for node in self.degree_six_nodes:
            if self.check_cycle(node, perm):
                if self.n_inverse_edges_components(node, perm) > 1:
                    npp += 1
            elif self.check_bowtie_for_pinch_point(node, perm):
                npp += 1
        return npp

    def check_orientable(self, perm: PERM) -> bool:
        ...

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
    if n_pinches == 1:
        log_info('Solution Found: Pinched Projective Plane', graph.res_file)
        return (True, True)

    if n_pinches == 2:
        log_info('Solution Found: 2-pinch-point Sphere', graph.res_file)
        return (True, True)

    if n_pinches == 0:
        is_orientable = graph.check_orientable(perm)
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
