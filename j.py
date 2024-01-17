import os
from glob import glob
import networkx as nx
from itertools import product
from argparse import ArgumentParser
from copy import deepcopy
from collections import defaultdict
from itertools import permutations
from tqdm import tqdm
from math import factorial
from re import finditer
from functools import partial

from typing import Optional
from typing import Union
from io import TextIOWrapper


STANDARD_OUT_STR = 'standard_out'
RESULTS_FILE_STR = 'results_file'
PERM_LOG_STR = 'perm_log_file'
RESULTS_DIR = 'results'
PERM = tuple[int]
NODE = int
EDGE = tuple[NODE]
PASS = tuple[NODE]
FACE_WALK = list[EDGE]
GRAPH_WALK = tuple[FACE_WALK]
EDGE_TO_ID_MAP = dict[EDGE, NODE]
ID_TO_EDGE_MAP = dict[NODE, EDGE]
PERM_LOG = list[str]
DIV_LENGTH = 60
PINCHED_PORJECTIVE_PLANE = 'Pinched Projective Plane'
TWO_PINCHPOINT_SPHERE = '2-pinchpoint Sphere'
TORUS = 'Torus'
KLIEN_BOTTLE = 'Klien Bottle'

ALL_PERMS_N = factorial(14)

def log_info(
    mssg: str,
    out: Optional[tuple[str]] = (STANDARD_OUT_STR, RESULTS_FILE_STR),
    res_file: Optional[TextIOWrapper] = None,
    perm_log: Optional[list[str]] = None,
) -> None:
    """ Log to screen, results file, and/or individual perm log """
    if STANDARD_OUT_STR in out:
        print(mssg)
    if RESULTS_FILE_STR in out:
        res_file.write(mssg + '\n')
    if PERM_LOG_STR in out:
        perm_log.append(mssg)
    return

log_perm = partial(log_info, out=PERM_LOG_STR)

def get_perm(perm_path: str) -> PERM:
    """ Parse individual perm from ADC-logging script """
    with open(perm_path) as f:
        number_line = f.read()
    number_line = number_line[len('The permutation ('):].split(')')[0]
    number_line = tuple(int(num) for num in number_line.split(', '))
    return number_line

def rotation_scheme(seq: list[PASS]) -> list[tuple[NODE]]:
    """
    Gets list of cycles from sequence of edges
    Example:
    [(1,2), (2,3), (3,1), (4,5), (5,6), (6,4)] -> [(1,2,3), (4,5,6)]
    """
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
    """ Extends nx.Graph """
    def __init__(
        self,
        file_path: str,
    ) -> None:
        """ Store frequently used values """
        self.graph = nx.read_edgelist(file_path, nodetype=int)
        self.degree_six_nodes = self.get_degree_six_nodes()
        self.edge_id_map = self.get_edge_id_map()
        self.id_edge_map = self.get_id_edge_map()

    def get_degree_six_nodes(self) -> tuple[NODE]:
        return tuple(
            node for node in self.graph.nodes if self.graph.degree(node) == 6
        )

    def get_edge_id_map(self) -> EDGE_TO_ID_MAP:
        return {edge: id for id, edge in enumerate(self.graph.edges, start=1)}

    def get_id_edge_map(self) -> ID_TO_EDGE_MAP:
        return {id: edge for id, edge in enumerate(self.graph.edges, start=1)}

    def get_vstar(self, node: NODE) -> tuple[EDGE]:
        """
        Get all edges including node as (vertex_id_1, vertex_id_2)
        where vertex_id_1 < vertex_id_2
        """
        return tuple((min(edge), max(edge)) for edge in self.graph.edges(node))

    def get_dual_edges(self, vstar: tuple[EDGE], perm: PERM) -> tuple[EDGE]:
        """ Get perm (dual) of vstar """
        return tuple(
            self.id_edge_map[perm[self.edge_id_map[edge] - 1]]
            for edge in vstar
        )

    def get_dual_face(self, node: NODE, perm: PERM) -> nx.classes.graph.Graph:
        """ Make nx.Graph object from perm(vstar(node)) """
        return nx.Graph(self.get_dual_edges(self.get_vstar(node), perm))

    def check_self_dual_perm(self, perm: PERM) -> bool:
        """
        Check that dual_face(vstar(node)) for all nodes in self.graph have a
        closed walk that includes each edge of the graph exactly once.
        """
        for node in self.graph.nodes:
            if not nx.is_eulerian(self.get_dual_face(node, perm)):
                return False
        return True

    def get_dual_bowtie_edge_options(
        self,
        source_node: NODE,
        perm: PERM
    ) -> list[list[PASS]]:
        """
        Get both uses of the four edges including the center of a bowtie that
        produce valid eulerian walks of the bowtie. Note that the tuples here,
        (node_1, node_2), represent a pass through two edges
        (node_1, cross_node) and (cross_node, node_2) where cross_node is the
        node with four edges in the bowtie = dual(vstar(source_node)).
        """
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

    def get_dual_bowtie_nodes(self, perm: PERM) -> dict[NODE, NODE]:
        """
        Get dictionary of degree six node(s) (keys), if the dual of their
        vstar is a bowtie, and the center node of that bowtie (values)
        """
        dual_bowtie_nodes = {}
        for degree_six_node in self.degree_six_nodes:
            dual_face = self.get_dual_face(degree_six_node, perm)
            for df_node in dual_face.nodes:
                if len(tuple(dual_face.neighbors(df_node))) == 4:
                    dual_bowtie_nodes[degree_six_node] = df_node
        return dual_bowtie_nodes

    def get_dual_graph_walks(self, perm: PERM) -> list[GRAPH_WALK]:
        """
        Get
        """
        # Determine how any bowties there are and therefore how many options we
        # have for the use of the bowtie edges
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

        def get_option_instruction(
            node: NODE,
            option: Union[int, tuple[int]]
        ) -> int:
            """
            Return index to pass to self.get_dual_bowtie_edge_options to
            instruct which set of passes to use
            """
            return (
                option if n_dual_bowtie_nodes == 1
                else option[list(dual_bowtie_nodes).index(node)]
            )

        def get_bowtie_walk(center_node, bowtie_passes) -> FACE_WALK:
            """ Unpack bowtie passes into walk edges """
            return [
                (bowtie_passes[0][0], center_node),
                (center_node, bowtie_passes[0][1]),
                (bowtie_passes[0][1], bowtie_passes[1][0]),
                (bowtie_passes[1][0], center_node),
                (center_node, bowtie_passes[1][1]),
                (bowtie_passes[1][1], bowtie_passes[0][0])
            ]

        def get_cycle_walk(edges: list[EDGE]) -> FACE_WALK:
            """ Organize cycle edges with rotation_scheme """
            rs = list(rotation_scheme(edges)[0])
            return [tup for tup in zip(rs, rs[1:] + [rs[0]])]

        dual_graph_walk_options = []
        for option in options:
            graph_walk = []
            for node, dual_face in zip(self.graph.nodes, dual_faces):
                if node in dual_bowtie_nodes:
                    bowtie_passes = self.get_dual_bowtie_edge_options(
                        node,
                        perm
                    )[get_option_instruction(node, option)]
                    face_walk = get_bowtie_walk(
                        dual_bowtie_nodes[node],
                        bowtie_passes
                    )
                else:
                    face_walk = get_cycle_walk(list(dual_face.edges))
                graph_walk.append(face_walk)

            dual_graph_walk_options.append(tuple(graph_walk))

        return dual_graph_walk_options

    def check_pinchpoint(
        self,
        degree_six_node: NODE,
        graph_walk: GRAPH_WALK,
        perm_log: PERM_LOG
    ) -> bool:
        pass_seq = []
        for face_walk in graph_walk:
            for edge_1, edge_2 in zip(
                    face_walk,
                    face_walk[1:] + [face_walk[0]]
                ):
                if edge_1[1] == degree_six_node and edge_2[0] == degree_six_node:
                    pass_seq.append((edge_1[0], edge_2[1]))

        log_perm(
            f'Using pass sequence around node {degree_six_node}: {pass_seq}',
            perm_log=perm_log
        )
        rs = rotation_scheme(pass_seq)
        log_perm(f'Rotation Scheme is: {rs}', perm_log=perm_log)
        return len(rs) > 1

    def get_n_pinchpoints(
        self,
        graph_walk: GRAPH_WALK,
        perm_log: PERM_LOG
    ) -> int:
        return sum(
            1 for node in self.degree_six_nodes
            if self.check_pinchpoint(node, graph_walk, perm_log=perm_log)
        )

    def check_orientable(
        self,
        graph_walk: GRAPH_WALK,
        perm_log: PERM_LOG
    ) -> bool:
        """
        """
        def get_oriented_walk(edge_count, walk) -> tuple[EDGE]:
            """
            If any edge in walk present in same direction in edge_count,
            reverse walk, else return unmodified walk.
            """
            for edge in walk:
                if edge in edge_count:
                    walk = tuple(tuple(reversed(edge)) for edge in walk)
                    break
            return walk

        edge_count = defaultdict(lambda: 0)
        for face_walk in graph_walk:
            face_walk = get_oriented_walk(edge_count, face_walk)
            for n1, n2 in face_walk:
                if (n2, n1) in edge_count:
                    edge_count[(n2, n1)] -= 1
                else:
                    edge_count[(n1, n2)] += 1

        return all(count == 0 for count in edge_count.values())

    def check_dual_graph_walks(self, perm: PERM, perm_log: PERM_LOG) -> None:
        """"""
        log_perm_filled = partial(log_perm, perm_log=perm_log)
        dual_graph_walks = self.get_dual_graph_walks(perm)
        for walk in dual_graph_walks:
            n_pinchpoints = self.get_n_pinchpoints(walk, perm_log=perm_log)
            is_orientable = self.check_orientable(walk, perm_log=perm_log)
            log_perm_filled(f'Number of pinchpoints: {n_pinchpoints}')
            log_perm_filled(f'Is orientable: {is_orientable}')

            if n_pinchpoints == 1:
                log_perm_filled(f'Solution Found: {PINCHED_PORJECTIVE_PLANE}')
            if n_pinchpoints == 2:
                log_perm_filled(f'Solution Found: {TWO_PINCHPOINT_SPHERE}')
            if n_pinchpoints == 0:
                if is_orientable:
                    log_perm_filled(f'Solution Found: {TORUS}')
                else:
                    log_perm_filled(f'Solution Found: {KLIEN_BOTTLE}')


def process_perm(
        perm: Union[str, tuple[int]],
        graph: Graph,
        log_level: int
    ) -> tuple[str, PERM_LOG]:
    """
    Return the 2-complex name, if any, in which graph.graph may be embedded
    where perm is assumed to represent an algebraic dual correspondence (ADC).
    """
    perm_log = []
    log_perm_filled = partial(log_perm, perm_log=perm_log)

    if isinstance(perm, str):
        perm = get_perm(perm)

    if not graph.check_self_dual_perm(perm):
        # Log separately if not an ADC to avoid excessive logging
        if log_level > 1:
            log_perm_filled(
                '\n'
                + '-' * DIV_LENGTH + '\n'
                + f'Perm: {perm}\n'
                + '-' * DIV_LENGTH + '\n'
                + f'Perm is not an ADC'
            )
        return perm_log

    log_perm_filled(
        '\n'
        + '-' * DIV_LENGTH + '\n'
        + f'Perm: {perm}\n'
        + '-' * DIV_LENGTH
    )
    graph.check_dual_graph_walks(perm, perm_log)

    return perm_log


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument(
        'perm_pattern',
        help='Use `all` to evaluate all perms instead of using perm files'
    )
    parser.add_argument('-m', dest='max_perms', type=int, default=None)
    parser.add_argument('-r', dest='refresh_results', action='store_true')
    parser.add_argument('-v', dest='log_level', type=int, default=1)
    parser.add_argument('-n', dest='n_processes', type=int, default=1)
    args = parser.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    if args.refresh_results:
        for path in glob(os.path.join(RESULTS_DIR, '*')):
            os.remove(path)

    results_file = os.path.join(
        RESULTS_DIR,
        f'{os.path.split(args.graph_file)[-1]}-results.txt'
    )

    with open(results_file, 'w') as res_file_handle:
        with open(args.graph_file) as graph_file:
            graph = Graph(graph_file)

        log_main = partial(log_info, res_file=res_file_handle)
        log_main('Graph Edges:')
        for i, edge in enumerate(graph.graph.edges):
            log_main(f'Edge {i + 1}: {edge}', res_file=res_file_handle)

        if not graph.degree_six_nodes:
            log_main('No degree six nodes', res_file=res_file_handle)
            return
        log_main(
            f'Degree six nodes are: {graph.degree_six_nodes}',
            res_file=res_file_handle
        )

        if args.perm_pattern == 'all':
            perms = permutations(range(1, len(graph.graph.edges) + 1))
            total_perms = args.max_perms or ALL_PERMS_N
        else:
            perms = sorted(glob(args.perm_pattern))
            total_perms = len(perms)
            if not perms:
                exit(f'Didn\'t find perms with pattern: {args.perm_pattern}')

        n_adc = 0
        solution_count = {
            PINCHED_PORJECTIVE_PLANE: 0,
            TWO_PINCHPOINT_SPHERE: 0,
            TORUS: 0,
            KLIEN_BOTTLE: 0,
        }
        for perm_n, perm in tqdm(enumerate(perms, start=1), total=total_perms):
                if args.max_perms and perm_n > args.max_perms:
                    perm_n -= 1
                    break
                log = process_perm(perm, graph, args.log_level)
                joined_log = '\n'.join(log)
                if 'Solution Found' in joined_log:
                    n_adc += 1
                    for solution in solution_count:
                        if solution in joined_log:
                            solution_count[solution] += sum(
                                1 for _ in finditer(solution, joined_log)
                            )
                    log_main(joined_log)

        log_main(
            '\n'
            + '-' * DIV_LENGTH + '\n'
            + f'Summary\n'
            + '-' * DIV_LENGTH + '\n'
            + f'Permutations checked: {perm_n}\n'
            + f'ADCs found: {n_adc}'
        )
        for solution, count in solution_count.items():
            log_main(f'{solution} found: {count}')

    return

if __name__ == '__main__':
    main()
