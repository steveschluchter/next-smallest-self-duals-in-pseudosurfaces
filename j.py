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

from typing import Optional
from typing import Union
from io import TextIOWrapper


STANDARD_OUT_STR = 'standard_out'
RESULTS_FILE_STR = 'results_file'
RESULTS_DIR = 'results'
PERM = tuple[int]
NODE = int
EDGE = tuple[NODE]
PASS = tuple[NODE]
EDGE_TO_ID_MAP = dict[EDGE, NODE]
ID_TO_EDGE_MAP = dict[NODE, EDGE]
DIV_LENGTH = 60
PINCHED_PORJECTIVE_PLANE = 'Pinched Projective Plane'
TWO_PINCHPOINT_SPHERE = '2-pinchpoint Sphere'
TORUS = 'Torus'
KLIEN_BOTTLE = 'Klien Bottle'

def log_info(
    mssg: str,
    res_file: Optional[TextIOWrapper] = None,
    out: Optional[str] = 'all'
) -> None:
    """ Log to screena and/or results file """
    ALL_STR = 'all'
    if out in (STANDARD_OUT_STR, ALL_STR):
        print(mssg)
    if out in (RESULTS_FILE_STR, ALL_STR):
        res_file.write(mssg + '\n')
    return

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
        res_file: Optional[TextIOWrapper] = None
    ) -> None:
        """ Store frequently used values """
        self.graph = nx.read_edgelist(file_path, nodetype=int)
        self.degree_six_nodes = self.get_degree_six_nodes()
        self.edge_id_map = self.get_edge_id_map()
        self.id_edge_map = self.get_id_edge_map()

        self.res_file = res_file

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

    def get_pass_sequences(
        self,
        degree_six_node: NODE,
        perm: PERM
    ) -> list[list[PASS]]:
        """
        Get list of lists of passes through degree_six_node from
        dual(vstar(node)) for all nodes in self.groph.
        Pass (node_1, node2)
            = (node_1, degree_six_node), (degree_six_node, node_2)
        """
        pass_seqs = [[]]
        for node in self.graph.nodes:
            dual_face = self.get_dual_face(node, perm)
            if degree_six_node in dual_face.nodes:
                """
                If degree_six_node is the center of the bowtie, make a pass_seq
                for each of the options from get_dual_bowtie_edge_options
                """
                if len(tuple(dual_face.neighbors(degree_six_node))) == 4:
                    options = self.get_dual_bowtie_edge_options(node, perm)
                    if pass_seqs == [[]]:
                        pass_seqs = options
                    else:
                        pass_seqs.extend(deepcopy(pass_seqs))
                        n_seqs = len(pass_seqs)
                        for seq_i, seq in enumerate(pass_seqs):
                            if seq_i < n_seqs // 2:
                                pass_seqs[option_ind].extend(options[0])
                            else:
                                pass_seqs[option_ind].extend(options[1])
                # If degree_six_node isn't at the center of a bowtie, we know
                # it only has two neighbors which respresent a pass
                else:
                    for pass_seq in pass_seqs:
                        pass_seq.append(
                            tuple(dual_face.neighbors(degree_six_node))
                        )
        log_info(
            f'Sequences found for node {degree_six_node}: {len(pass_seqs)}',
            self.res_file
        )
        """
        Note interesting casses where two nodes map to bowties with the same
        center node
        """
        if len(pass_seqs) > 2:
            log_info(
                f'Found n_pass_seqs > 2 for node {degree_six_node}',
                self.res_file
            )
        return pass_seqs

    def check_pinchpoint(
        self,
        degree_six_node: NODE,
        perm: PERM
    ) -> bool:
        """
        Perform rotation_scheme on each pass sequence for degree_six_node. If
        any produce more than 1 cycle (umbrella), we have found a pinchpoint.
        """
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

    def get_n_pinchpoints(self, perm: PERM) -> int:
        return sum(
            1 for node in self.degree_six_nodes
            if self.check_pinchpoint(node, perm)
        )

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

    def check_orientable(self, perm: PERM) -> bool:
        """
        Check if it is possible to impose a consistent orientation on
        dual(vstar(node)) for all nodes in self.graph. This means is it
        possible to use each edge exacty as many times in one direction
        (node_1, node_2) as the other (node_2, node_1) over all of their
        boundry walks.
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

        def get_bowtie_walk(center_node, bowtie_passes) -> list[EDGE]:
            """ Unpack bowtie passes into walk edges """
            return [
                (bowtie_passes[0][0], center_node),
                (center_node, bowtie_passes[0][1]),
                (bowtie_passes[0][1], bowtie_passes[1][0]),
                (bowtie_passes[1][0], center_node),
                (center_node, bowtie_passes[1][1]),
                (bowtie_passes[1][1], bowtie_passes[0][0])
            ]

        def get_cycle_walk(edges: list[EDGE]) -> list[EDGE]:
            """ Organize cycle edges with rotation_scheme """
            rs = list(rotation_scheme(edges)[0])
            return [tup for tup in zip(rs, rs[1:] + [rs[0]])]

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

        for option in options:
            edge_count = defaultdict(lambda: 0)
            for node, dual_face in zip(self.graph.nodes, dual_faces):
                if node in dual_bowtie_nodes:
                    # Get relevant bowtie walk for current node and option
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
                    # If not the first walk, need to orient walk based on
                    # previous use of edges
                    walk = get_oriented_walk(edge_count, walk)
                    for n1, n2 in walk:
                        if (n2, n1) in edge_count:
                            edge_count[(n2, n1)] -= 1
                        else:
                            edge_count[(n1, n2)] += 1
            if (
                    all(count == 0 for count in edge_count.values())
                    and len(edge_count) == 14
                ):
                return True
        return False

def process_perm(
        perm: Union[str, tuple[int]],
        graph: Graph,
        log_level: int
    ) -> str:
    """
    Return the 2-complex name, if any, in which graph.graph may be embedded
    where perm is assumed to represent an algebraic dual correspondence (ADC).
    """
    if isinstance(perm, str):
        perm = get_perm(perm)

    if not graph.check_self_dual_perm(perm):
        # Log separately if not an ADC to avoid excessive logging
        if log_level > 1:
            log_info('', graph.res_file)
            log_info('-' * DIV_LENGTH, graph.res_file)
            log_info(f'Perm: {perm}', graph.res_file)
            log_info('-' * DIV_LENGTH, graph.res_file)
            log_info(f'Perm is not an ADC', graph.res_file)
        return ''

    log_info('', graph.res_file)
    log_info('-' * DIV_LENGTH, graph.res_file)
    log_info(f'Perm: {perm}', graph.res_file)
    log_info('-' * DIV_LENGTH, graph.res_file)

    n_pinchpoints = graph.get_n_pinchpoints(perm)
    is_orientable = graph.check_orientable(perm)
    log_info(f'Number of pinchpoints: {n_pinchpoints}', graph.res_file)
    log_info(f'Is orientable: {is_orientable}', graph.res_file)
    if n_pinchpoints == 1:
        log_info(f'Solution Found: {PINCHED_PORJECTIVE_PLANE}', graph.res_file)
        return PINCHED_PORJECTIVE_PLANE

    if n_pinchpoints == 2:
        log_info(f'Solution Found: {TWO_PINCHPOINT_SPHERE}', graph.res_file)
        return TWO_PINCHPOINT_SPHERE

    if n_pinchpoints == 0:
        if is_orientable:
            log_info(f'Solution Found: {TORUS}', graph.res_file)
            return TORUS
        else:
            log_info(f'Solution Found: {KLIEN_BOTTLE}', graph.res_file)
            return KLIEN_BOTTLE


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
    args = parser.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    if args.refresh_results:
        for path in glob(os.path.join(RESULTS_DIR, '*')):
            os.remove(path)

    results_file = os.path.join(
        RESULTS_DIR,
        f'{os.path.split(args.graph_file)[-1]}-results.txt'
    )

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

        if args.perm_pattern == 'all':
            perms = permutations(range(1, len(graph.graph.edges) + 1))
            n_perms = args.max_perms or factorial(len(graph.graph.edges))
        else:
            perms = sorted(glob(args.perm_pattern))
            n_perms = len(perms)
            if not perms:
                exit(f'Didn\'t find perms with pattern: {args.perm_pattern}')

        n_adc = 0
        solution_count = {
            PINCHED_PORJECTIVE_PLANE: 0,
            TWO_PINCHPOINT_SPHERE: 0,
            TORUS: 0,
            KLIEN_BOTTLE: 0,
        }
        for perm_n, perm in enumerate(tqdm(perms, total=n_perms), start=1):
            if args.max_perms and perm_n > args.max_perms:
                perm_n -= 1
                break
            perm_solution = process_perm(perm, graph, args.log_level)
            if perm_solution:
                n_adc += 1
                solution_count[perm_solution] += 1

        log_info('', res_file_wrapper)
        log_info('-' * DIV_LENGTH, res_file_wrapper)
        log_info(f'Summary', res_file_wrapper)
        log_info('-' * DIV_LENGTH, res_file_wrapper)
        log_info(f'Permutations checked: {perm_n}', res_file_wrapper)
        log_info(f'ADCs found: {n_adc}', res_file_wrapper)
        for solution, count in solution_count.items():
            log_info(f'{solution} found: {count}', res_file_wrapper)

    return

if __name__ == '__main__':
    main()
