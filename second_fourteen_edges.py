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
    print("vertex_star_edges")
    star = []
	

    for edge in graph.edges(vertex):
        #print("in loop")
        #print(G.edges(vertex))
        edge_sub = (int(edge[0]), int(edge[1]))
		
        if(edge[0] > edge[1]):
		
            edge_tuple = (edge_sub[1], edge_sub[0])
		
        else:
			
            edge_tuple = (edge_sub[0], edge_sub[1])
		
        star.append(edge_tuple)
	
    return star

def dual_edges(vstar,perm):
	
    dual_edges_list = []
	
    for edge in vstar:
        dual_label = perm[ ( perm.index(edges_to_numbers[edge]) + 1) % len(perm) ]
        dual_edge = numbers_to_edges[dual_label]
        dual_edges_list.append(dual_edge)
        #print(perm)
        #print(f"edge {edge}")
        #print(f"label {edges_to_numbers[edge]}")
        #print(f"dual_label {dual_label}")
        #print(f"dual_edge {dual_edge}")
	
    return dual_edges_list

def is_connected_algebraic_dual(graph, perm, dual_edges):

    print("in is_connected_algebraic_dual")

	#TODO GET VERTICES IN A LIST
	#TODO CHECK EACH VERTEX-STAR IN THE LIST
	#TODO IF ONE VERTEX FLUNKS IS_CONNECTED OR IS_EULERIAN, then return False
    
    vertices = [v for v in list(graph.nodes())]
    print("Here's vertices!")
    print(vertices)

    for k in vertices:
        print(k)
        print(vertices)
        vertex_star = vertex_star_edges(graph, vertices[k])
        edge_duals = dual_edges(vertex_star, perm)
        H = nx.Graph()
        H.add_edges_from(edge_duals)
        print(H.edges())

        if(not nx.is_eulerian(H)):

            #print("False exiting is_connected_algebraic_dual")
            return False

    
    #print('True exiting is_connected_algebraic_dual')

    return True


#Determines if a 6-star permutes to edges including a cycle.
#Note: This is only run under the conditions that perm is an algebraic duality correspondence.
def check_cycles(graph, permutation):
	

    star_six = vertex_star_edges(graph, degree_six)
    dual_six = dual_edges(star_six, perm)
    
    """
    print(perm)
    print(degree_six)
    print(star_six)
    print(dual_six)
    """

    H = nx.Graph()
    H.add_edges_from(dual_six)

	#Iterate through the edges and ensure that they are all connected with degree 2.
    for node in H.nodes():
        if(H.degree(node) != 2):
            #print ("\nThis 6-star permutes to a bowtie. There is a node (%s) without two degrees (has %s)." % (node,H.degree(node)))
            return False        
    #print("The 6 stars map to cycles!  Winner")
    return True

def analyze_perm(graph, permutation, edgestonumbers, numberstoedges, degreesixvertices):
    print("in analyze perm")   
    write_this_to_file = f""

    if is_connected_algebraic_dual(graph, permutation, edgestonumbers, numberstoedges):

        write_this_to_file += f"The permutation {perm} is an algebraic dual permutation of {graph_name}."  
        #print("Write this to file ", write_this_to_file)
        #input("x")
        #print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())) 

        filename = f"{graph_name}-"
        for i in permutation:
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
            #f.write(write_this_to_file)
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
            #f.write(write_this_to_file)
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
    