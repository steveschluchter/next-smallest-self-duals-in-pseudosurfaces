#!/usr/bin/python

from itertools import permutations
from multiprocessing import Pool
import networkx as nx
import time, sys, os, re


FACTORIAL = 11*10*9*8*7*6*5*4*3*2
QUOTIENT = 14*13*12


#Returns a list of edges incident to a given vertex.
#Note: Since networkx treats the undirected edge (0,1) as being distinct from the undirected edge (1,0), we make the choice to put the lower-numbered vertex in the first entry in each ordered pair.
def vertex_star_edges(G, vertex):
    print("vertex_star_edges")
    star = []
	

    for edge in G.edges(str(vertex)):
        print("in loop")
        print(G.edges(str(vertex)))
        edgeSub = (int(edge[0]), int(edge[1]))
		
        if(edge[0] > edge[1]):
		
            edgeTuple = (edgeSub[1], edgeSub[0])
		
        else:
			
            edgeTuple = (edgeSub[0], edgeSub[1])
		
        star.append(edgeTuple)
	
    return star

#Returns the result of permuting the edges in vstar using the permutation perm. 
def dual_edges(vstar,perm):
	
	dual_edges_list = []
	
	for edge in vstar:

		label = perm[edgesToNumbers[edge]-1]
		edge = numbersToEdges[label]
		dual_edges_list.append(edge)
	
	return dual_edges_list

#Determines if a 6-star permutes to edges including a cycle.
#Note: This is only run under the conditions that perm is an algebraic duality correspondence.
def check_cycles(perm):
	

    star_six = vertex_star_edges(G, degree_six)
    dual_six = dual_edges(star_six, perm)

    print(perm)
    print(degree_six)
    print(star_six)
    print(dual_six)

    H = nx.Graph()
    H.add_edges_from(dual_six)

	#Iterate through the edges and ensure that they are all connected with degree 2.
    for node in H.nodes():
        if(H.degree(node) != 2):
            print ("\nThis 6-star permutes to a bowtie. There is a node (%s) without two degrees (has %s)." % (node,H.degree(node)))
            return False        
    print("The 6 stars map to cycles!  Winner")
    return True

#Returns true if a vertex star maps via the inverse permutation to edges inducing more than one component of G, thus detecting multiple umbrellas.
#Notes: This is not an exhaustive test since it is used if starSix maps to edges inducing a cycle.
def check_components(perm):
	
	star_six = vertex_star_edges(G,degree_six)
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
        edge_duals = dual_edges(vertex_star, perm) 
                    
        #NOTE: H is the induced graph containing only the permuted edges from star.
        H = nx.Graph()
        H.add_edges_from(edge_duals)
        
        if(not nx.is_eulerian(H)):
            
            return False


    
    return True 


def analyze_perm(perm):
   
    write_this_to_file = f""

    if is_connected_algebraic_dual(G, perm):
       
        write_this_to_file += f"The permutation {perm} is an algebraic dual permutation of {graph_name}."  
        print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

        filename = f"{graph_name}-"
        for i in perm:
            filename += f"{i}"
        filename += ".txt"

        #TODO Check to see if each vertex star maps to cycles.
        #Else, make all possible choices of boundary walks of bowties.

        if(check_cycles(perm)):
        
            write_this_to_file += "This perm is a winner!  It maps each vertex star to a cycle."                

        #TODO Do check to classify pinchpoint as having one or more umbrellas.
        #There will be at least one special case of a graph having more than one possible pinchpoint.

        #TODO Compute rank of H1(P) via computing the dimension over Z of the facial boundary walks 
        
        """
        with open(filename, 'w') as f:
            #f.write(write_this_to_file)
            print(f"Wrote file {filename}.")
        """


if __name__ == "__main__":
    print("begin program")
    graph_name = sys.argv[1]
    filename = ""    
    #TODO reset time index to local (Pacific) time
    print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 
    #Builds the graph from the structure provided in the input file.
    #Pull the file name from the command line args.
    
    if(len(sys.argv) >= 2):
        filename = sys.argv[1]
    """    
    elif(len(sys.argv) > 2):
        print ("Additional arguments supplied. Try: python <<program name>> <<filename>>")
        #sys.exit(1) uncomment this if you want the program to exit here
    else:
        print ("Too few arguments supplied. Try: python <<program name>> <<filename>>")
        #sys.exit(1) uncomment this if you want the program to exit here

    #Open the file and construct the graph from the formatted instructions.
    """

    G = nx.read_edgelist(filename)

    print(G.edges())
    print(type(G.edges()))
    print("Printed G.edges()")
    numbersToEdges = {}
    edgesToNumbers = {}
    list_of_edges = list(G.edges())

    degree_six = []
    for vertex in G.nodes():    
        if (G.degree[vertex] > 5):
            degree_six.append(int(vertex))

    for i in range(0,len(list_of_edges)):
        numbersToEdges[i+1] = (int(list_of_edges[i][0]),int(list_of_edges[i][1]))
        edgesToNumbers[(int(list_of_edges[i][0]),int(list_of_edges[i][1]))] = i+1

    print("numbersToEdges " + str(numbersToEdges))
    print("edgesToNumbers " + str(edgesToNumbers))
    print("degree_six " + str(degree_six))

    path_to_files = sys.argv[2]
    stream = os.popen(f"ls {path_to_files}")
    graph_files = stream.read().split("\n")
    graph_files = [graph_file for graph_file in graph_files if graph_file != '']
    print(graph_files)
    


    for graph_file in graph_files:

        path_to_file = path_to_files + f"/{graph_file}"
        regex_this = open(path_to_file, 'r')
        regex_this = regex_this.read()
        regexed_checker = re.findall(r"\d+", regex_this)
        tuple_perm = tuple(int(i) for i in regexed_checker)
        analyze_perm(tuple_perm)
        check_cycles(tuple_perm)

    """

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
    
    
    for perm in perms:
        analyze_perms(perm, sys.argv[1])
    """

    #TODO Write code to scrape permutation data out of collections of files
    # pertaining to ingested graph.  There will have to be the use of the os and subprocess commands.

    print("Finished analyzing perms in batch.") 
    print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

    print("Done analyzing permutations")

    sys.exit(1)
