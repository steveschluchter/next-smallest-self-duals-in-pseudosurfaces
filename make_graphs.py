#!/anaconda/bin/python3

"""
This progam generates graphs that will be suspended from another vertex
not among the list [2,3,4,5,6,7], vertex 1.
These graphs are an exhaustive list of graphs, up to isomorphism, having
the following properties:
1.) Each graph has 6 vertices.
2.) Each graph is connected.
3.) Each graph has 8 edges.
4.) Each graph has a minumum degree of 2.
5.) The edges of the graph must be incident to all vertices in the list 
    [2,3,4,5,6,7].
The list of graphs is printed at the end, and is stored in the list
named nonIsomorphicGraphs.
"""


from itertools import combinations
import networkx as nx

comb = combinations([2,3,4,5,6,7],2)

graphEdges = combinations(list(comb),8)

nonIsomorphicGraphs = []


print("Steve was here")
print("This will take a while.  Give it less than a minute. :)")

for i in graphEdges:
    
    junkGraphFlag = False
    G = nx.Graph()
    edges = list(i)
    G.add_edges_from(edges)

    if(not list(sorted(G.nodes())) == [2,3,4,5,6,7]):
        continue

    if(not nx.is_connected(G)):
        continue

    degrees = list(nx.degree(G))

    for v in degrees:
        if v[1] == 1:
            continue

    for graph in nonIsomorphicGraphs:

        if (nx.is_isomorphic(G,graph)):
            junkGraphFlag = True
            break
    
    if (junkGraphFlag):
        continue

    nonIsomorphicGraphs.append(G)


print(f"I've got {len(nonIsomorphicGraphs)} distinct nonisomorphic graphs for your computing pleasure.")

graphNum = 1
for graph in nonIsomorphicGraphs:
    print(graphNum, ' ',graph.edges())
    graphNum = graphNum + 1

print("End program!")