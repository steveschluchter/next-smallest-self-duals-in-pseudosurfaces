#!/usr/bin/python

from itertools import permutations
from multiprocessing import Pool
import networkx as nx
import pprint
import time
import sys

def vertexStar(G,vertex):
	
	star = []
	
	for edge in G.edges(vertex):
		if(edge[0] > edge[1]):
			edgeTuple = (edge[1], edge[0])
		else:
			edgeTuple = (edge[0], edge[1])
		star.append(edgeTuple)
	
	return star


#Returns a list of edges incident to a given vertex.
#Note: Since networkx treats the undirected edge (0,1) as being distinct from the undirected edge (1,0), we make the choice to put the lower-numbered vertex in the first entry in each ordered pair.
def vertexStar(G,vertex):
	
	star = []
	
	for edge in G.edges(vertex):
		if(edge[0] > edge[1]):
			edgeTuple = (edge[1], edge[0])
		else:
			edgeTuple = (edge[0], edge[1])
		star.append(edgeTuple)
	
	return star

#Returns the result of permuting the edges in vstar using the permutation perm. 
def dualEdges(vstar,perm):
	
	dualEdges = []
	
	for edge in vstar:
		label = perm[edgestonumbers[edge]-1]
		edge = numberstoedges[label]
		dualEdges.append(edge)
	
	return dualEdges

#Determines if a 6-star permutes to edges including a cycle.
#Note: This is only run under the conditions that perm is an algebraic duality correspondence.
def checkCycles(perm):
	
	starSix = vertexStar(G,degreeSix)
	dualSix = dualEdges(starSix,perm)
	
	H = nx.Graph()
	H.add_edges_from(dualSix)

	#Iterate through the edges and ensure that they are all connected with degree 2.
	for node in H.nodes():
		if(H.degree(node) != 2):
			print ("\nThis 6-star permutes to a bowtie. There is a node (%s) without two degrees (has %s)." % (node,H.degree(node)))
			return False
	return True

#Returns true if a vertex star maps via the inverse permutation to edges inducing more than one component of G, thus detecting multiple umbrellas.
#Notes: This is not an exhaustive test since it is used if starSix maps to edges inducing a cycle.
def checkComponents(perm):
	
	starSix = vertexStar(G,degreeSix)
	inverseEdgeMap = []
	
	for edge in starSix:
		inverseEdgeMap.append(numberstoedges[perm.index(edgestonumbers[edge]) + 1])
	
	print ("\nBuilt Graph from inverse edge map: %s" %(inverseEdgeMap))
	F = nx.Graph()
	F.add_edges_from(inverseEdgeMap)
	
	return (nx.number_connected_components(F) > 1)


def vertex_stars(G, vertex):
	
	star = []
	
	for edge in G.edges(vertex):
		if(edge[0] > edge[1]):
			edgeTuple = (edge[1], edge[0])
		else:
			edgeTuple = (edge[0], edge[1])
		star.append(edgeTuple)
	
	return star	

def is_connected_algebraic_dual(G,perm,fileString):

	#TODO GET VERTICES IN A LIST
	#TODO CHECK EACH VERTEX-STAR IN THE LIST
	    #TODO IF ONE VERTEX FLUNKS IS_CONNECTED OR IS_EULERIAN, then return False

	return True 



if __name__ == "__main__":

    pp = pprint.PrettyPrinter(indent=4)

    #Builds the graph from the structure provided in the input file.
    #Notes: See files: F1.txt, F2.txt, etc. in the GitHub repository for proper notation.
    #Pull the file name from the command line args.
    if(len(sys.argv) == 2):
        filename = sys.argv[1]
    elif(len(sys.argv) > 2):
        print ("Additional arguments supplied. Try: python <<program name>> <<filename>>")
        #sys.exit(1) uncomment this if you want the program to exit here
    else:
        print ("Too few arguments supplied. Try: python <<program name>> <<filename>>")
        #sys.exit(1) uncomment this if you want the program to exit here

    #Open the file and construct the graph from the formatted instructions.

    G = nx.read_edgelist(filename)

    print(G.edges())
    print(type(G.edges()))
    print("Printed G.edges()")
    numbersToEdges = {}
    edgesToNumbers = {}
    list_of_edges = list(G.edges())

    for i in range(0,len(list_of_edges)):
        numbersToEdges[i+1] = (int(list_of_edges[i][0]),int(list_of_edges[i][1]))
        edgesToNumbers[(int(list_of_edges[i][0]),int(list_of_edges[i][1]))] = i+1

    print("checking nodes and shit")

    print(G.nodes)

    print(type(G.nodes))

    print(type(G.nodes[0]))

    print(vertex_stars(G))
	
	#print(vertex_star(G, G.nodes()[0]))
    
    #pp.pprint(vertex_star(G,G.nodes[0]))
    
