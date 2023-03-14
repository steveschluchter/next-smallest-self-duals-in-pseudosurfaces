#!/usr/bin/python

from itertools import permutations
from multiprocessing import Pool
import networkx as nx
import pprint
import time
import sys

FACTORIAL = 11*10*9*8*7*6*5*4*3*2
QUOTIENT = 14*13*12


#Returns a list of edges incident to a given vertex.
#Note: Since networkx treats the undirected edge (0,1) as being distinct from the undirected edge (1,0), we make the choice to put the lower-numbered vertex in the first entry in each ordered pair.
def vertex_star_edges(G,vertex):
	
    star = []
	
    for edge in G.edges(vertex):
   
        edgeSub = (int(edge[0]), int(edge[1]))
		
        if(edge[0] > edge[1]):
		
            edgeTuple = (edgeSub[1], edgeSub[0])
		
        else:
			
            edgeTuple = (edgeSub[0], edgeSub[1])
		
        star.append(edgeTuple)
	
    return star

#Returns the result of permuting the edges in vstar using the permutation perm. 
def dualEdges(vstar,perm):
	
	dualEdges = []
	
	for edge in vstar:

		label = perm[edgesToNumbers[edge]-1]
		edge = numbersToEdges[label]
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


def is_connected_algebraic_dual(G, perm):

	#TODO GET VERTICES IN A LIST
	#TODO CHECK EACH VERTEX-STAR IN THE LIST
	#TODO IF ONE VERTEX FLUNKS IS_CONNECTED OR IS_EULERIAN, then return False

    for v in G.nodes():
        vertex_star = vertex_star_edges(G, v)
        edge_duals = dualEdges(vertex_star, perm) 
                    
        #NOTE: H is the induced graph containing only the permuted edges from star.
        H = nx.Graph()
        H.add_edges_from(edge_duals)
        
        if(not nx.is_eulerian(H)):
            
            return False


    
    return True 


def analyze_perms(perm):
   
    writeThisToFile = f""

    if is_connected_algebraic_dual(G, perm):
        writeThisToFile += f"The permutation {perm} is an algebraic dual permutation of {graphName}."  
        print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

        fileName = f"{graphName}-"
        for i in perm:
            fileName += f"{i}"
        fileName += ".txt"

        with open(fileName, 'w') as f:
            f.write(writeThisToFile)
            print(f"Wrote file {filename}.")

if __name__ == "__main__":
    print("begin program")
    graphName = sys.argv[1]
    print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 
    #pp = pprint.PrettyPrinter(indent=4)
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

    print("checking nodes")

    print(G.nodes)

    print(type(G.nodes))

    print("numbersToEdges " + str(numbersToEdges))
    print("edgesToNumbers " + str(edgesToNumbers))

    fourteen = [digit for digit in range(1, 15)]

    perms = permutations(fourteen, 14)
    
    perm_array = []

    for i in range(QUOTIENT):
        

        print("loading up array")
        print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))        

        perm_array = []

        for j in range(FACTORIAL):

            perm_array.append(next(perms))

        print("Started processing a batch of perms")
        print(f"Started processing with {perm_array[0]}")
        with Pool(processes=7) as pool:
           
            pool.map(analyze_perms, perm_array)
    
    """
    for perm in perms:
        analyze_perms(perm, sys.argv[1])
    """

    print("Finished analyzing perms in batch.") 
    print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

    print("Done analyzing permutations")

    sys.exit(1)
