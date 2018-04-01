#----------------------------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University                                 #
# Created: 2017-01-24, Updated: 2017-02-12                                         # 
# Descrption: Module for removing cryptic relations based on minimum vertex cover  # 
#----------------------------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import networkx as nx
import itertools as it
from subprocess import call
from networkx.drawing.nx_pydot import write_dot
from networkx import connected_component_subgraphs
from networkx.drawing.nx_pylab import draw_graphviz
from networkx.algorithms.approximation.vertex_cover import min_weighted_vertex_cover

def exact_min_vertex_cover(graph):
    """Returns node subset of connected graph for exact min vertex cover problem
        where ties are broken in lexicographic order of node ids in node subsets."""
    for N in range(1,len(graph.nodes())+1):
        for graph_sub in it.combinations(sorted(graph.nodes(), reverse=True), N):
            graph_temp = graph.copy()
            graph_temp.remove_nodes_from(graph_sub)
            if len(graph_temp.edges()) == 0:
                return list(graph_sub)

def remove_cryptic_relations(relations_list, relations_dot=None, relations_fig=None, threshold=None, verbose=False):
    """"Returns minimum set of individuals to remove for set of cryptic relations"""
    # construct a graph of relations from given iterable of tuples of node ids,
    # then remove individuals corresponding to exact minimum vertex cover,
    # which has exponential-time complexity, but most subgraphs are small
    # unless size of subgraph is greater than the argument 'threshold',
    # then remove individuals corresponding to approximate minimum vertex cover,
    # which has linear-time complexity, and guarantees a factor-2 approximation
    # can handle loop edges, for inbred individuals, where haplotypes are related
    G = nx.Graph()
    for i, j in relations_list:
        G.add_edge(i, j)
    # iterate over each connected component
    removed_list = set()
    if verbose:
        print 'Removing cryptic relations...'
    for g in list(connected_component_subgraphs(G)):
        if verbose:
            print 'Size of current connected component is ' + str(len(g))
        # approximate solution
        if threshold is not None and g.number_of_nodes() >= threshold:
            for i in min_weighted_vertex_cover(g):
                removed_list.add(i)
        # exact solution to mvc
        else:
            for i in exact_min_vertex_cover(g):
                removed_list.add(i)
    # color nodes that were removed
    for i in removed_list:
        G.node[i]['color'] = 'red'
    # save dot file and pdf file of graph
    if relations_dot is not None and relations_fig is not None:
        write_dot(G, relations_dot)
        call('dot -Tpdf {0} -o {1}'.format(relations_dot, relations_fig), shell=True)
    return removed_list
