import pytest
import networkx as nx

import j

@pytest.fixture(scope='session')
def graph1():
    with open('graoh1') as f:
        graph1 = nx.read_edgelist(f, nodetype=int)
    return graph1

@pytest.mark.parametrize(
    'pass_seq, rotations',
    [
        # graph1 perm: (11, 9, 13, 6, 12, 3, 10, 2, 14, 7, 1, 4, 5, 8)
        ([(7, 1), (3, 6), (7, 5), (5, 6), (3, 4), (1, 4)], [[7, 1, 4, 3, 6, 5]]),
        ([(7, 3), (1, 6), (7, 5), (5, 6), (3, 4), (1, 4)], [[7, 3, 4, 1, 6, 5]]),
        # graph1 perm: (10, 3, 5, 12, 1, 2, 11, 14, 7, 9, 13, 6, 8, 4) 2pp Sphere
        ([(2, 3), (4, 6), (4, 7), (6, 7), (2, 5), (3, 5)], [[2, 3, 5], [4, 6, 7]]),
        ([(6, 5), (3, 7), (6, 1), (7, 4), (3, 4), (1, 5)], [[6, 5, 1], [3, 7, 4]]),
    ]
)
def test_rotation_scheme(pass_seq, rotations):
    assert j.rotation_scheme(pass_seq) == rotations

class TestGraph:
    ...
