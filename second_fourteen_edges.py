#!/usr/bin/python3

from itertools import permutations
from multiprocessing import Pool
import networkx as nx
import time, sys, os, re

FACTORIAL = 11*10*9*8*7*6*5*4*3*2
QUOTIENT = 14*13*12


#Returns a list of edges incident to a given vertex.
#Note: Since networkx treats the undirected edge (0,1) as being distinct from the undirected edge (1,0), we make the choice to put the lower-numbered vertex in the first entry in each ordered pair.
def vertex_star_edges(graph, vertex):
    
    star = []
    
    for edge in graph.edges(vertex):
        
        edge_sub = (int(edge[0]), int(edge[1]))
		
        if(edge[0] > edge[1]):
		
            edge_tuple = (edge_sub[1], edge_sub[0])
		
        else:
			
            edge_tuple = (edge_sub[0], edge_sub[1])
		
        star.append(edge_tuple)

    return star

def dual_edges(vertexstar, permutation, numberstoedges, edgestonumbers):
	
    dual_edges_list = []
	
    for edge in vertexstar:
        dual_label = permutation[ ( permutation.index(edgestonumbers[edge]) + 1) % len(permutation) ]
        dual_edge = numberstoedges[dual_label]
        dual_edges_list.append(dual_edge)
	
    return dual_edges_list

def is_connected_algebraic_dual(graph, perm, edgestonumbers, numberstoedges):

	#TODO GET VERTICES IN A LIST
	#TODO CHECK EACH VERTEX-STAR IN THE LIST
	#TODO IF ONE VERTEX FLUNKS IS_CONNECTED OR IS_EULERIAN, then return False
    
    vertices = [v for v in list(graph.nodes())]

    for k in vertices:
        vertex_star = vertex_star_edges(graph, k)
        edge_duals = dual_edges(vertex_star, perm, numberstoedges, edgestonumbers)

        H = nx.Graph()
        H.add_edges_from(edge_duals)

        if(not nx.is_eulerian(H)):

            return False

    
    #print('True exiting is_connected_algebraic_dual')

    
    return True


#Determines if a 6-star permutes to edges including a cycle.
#Note: This is only run under the conditions that perm is an algebraic duality correspondence.
def check_cycles(graph, permutation, edgestonumbers, numberstoedges, degreesixvertices):
    

    for i in degreesixvertices:

        starsix = vertex_star_edges(graph, i)
        dualsix = dual_edges(starsix, permutation, numberstoedges, edgestonumbers)

        H = nx.Graph()
        H.add_edges_from(dualsix)

	    #Iterate through the edges and ensure that they are all connected with degree 2.
        for node in H.nodes():

            if(H.degree(node) != 2):

                return False        

    return True

def analyze_perm(graph, permutation, edgestonumbers, numberstoedges, degreesixvertices):
    write_this_to_file = f""

    if is_connected_algebraic_dual(graph, permutation, edgestonumbers, numberstoedges):

        write_this_to_file += f"The permutation {permutation} is an algebraic dual permutation of {graph_name}."  
        print("Write this to file ", write_this_to_file)
        #input("x")
        print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

        filename = f"{graph_name}-"
        for i in permutation:
            filename += f"{i}"
        filename += ".txt"

        #TODO Check to see if each vertex star maps to cycles.
        #Else, make all possible choices of boundary walks of bowties.
        
        if(check_cycles(graph, permutation, edgestonumbers, numberstoedges, degreesixvertices)):

            write_this_to_file += "This perm is a winner!  It maps each vertex star to a cycle."


        #TODO Do check to classify pinchpoint as having one or more umbrellas.
        #There will be at least one special case of a graph having more than one possible pinchpoint.

        #TODO Compute rank of H1(P) via computing the dimension over Z of the facial boundary walks 
        
        with open(filename, 'w') as f:
            f.write(write_this_to_file)
            print(f"Wrote file {filename}.")

        
"""
def analyze_perm(G, perm):

    print("in analyze perm")   
    write_this_to_file = f""

    if is_connected_algebraic_dual(G, perm):

        write_this_to_file += f"The permutation {perm} is an algebraic dual permutation of {graph_name}."  
        #print("Write this to file ", write_this_to_file)
        #input("x")
        #print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

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
        
        
        with open(filename, 'w') as f:
            f.write(write_this_to_file)
            print(f"Wrote file {filename}.")
"""



if __name__ == "__main__":

    print("Begin Program")
    graph_name = sys.argv[1]
    print(time.strftime("Program started at %a, %d %b %Y %H:%M:%S Pacific Time.", time.localtime()))
    master_graph = nx.read_edgelist(graph_name)

    print(master_graph.edges())
    print(type(master_graph.edges()))
    print("Printed G.edges()")
    numbers_to_edges = {}
    edges_to_numbers = {}
    list_of_edges = list(master_graph.edges())

    degree_six = []
    for vertex in master_graph.nodes():    
        if (master_graph.degree[vertex] > 5):
            degree_six.append(vertex)

    for i in range(0,len(list_of_edges)):
        numbers_to_edges[i+1] = (int(list_of_edges[i][0]),int(list_of_edges[i][1]))
        edges_to_numbers[(int(list_of_edges[i][0]),int(list_of_edges[i][1]))] = i+1

    print("numbersToEdges " + str(numbers_to_edges))
    print("edgesToNumbers " + str(edges_to_numbers))
    print("degree_six " + str(degree_six))

    fourteen = [digit for digit in range(1, 15)]
    perms = permutations(fourteen, 14)
    
    
    for permutation in perms:
        analyze_perm(master_graph, permutation, edges_to_numbers, numbers_to_edges, degree_six)

    for i in range(QUOTIENT): 

        print("loading up array")
        print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))        

        perm_array = []

        for j in range(FACTORIAL):

            perm_array.append(next(perms))

        print("Started processing a batch of perms")
        print(f"Started processing with {perm_array[0]}")
        with Pool(processes=7) as pool:
           
            pool.starmap(analyze_perm, perm_array)


    print(time.strftime("Progream ended at %a, %d %b %Y %H:%M:%S Pacific Time."))
    
