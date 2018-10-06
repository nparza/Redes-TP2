#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:43:24 2018

@author: sofinico
"""

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from datetime import datetime


#%% EJC CENTRALIDADES


def degrees2dict(graph):
    d = dict(graph.degree)
    return d

def cent_cutoff(graph, centrality):
    c = [  ]            #Fracción de nodos removidos en cada iteración
    mc = [  ]           #Nodos de la componente máxima
    
    maxcomp = max(nx.connected_component_subgraphs(graph),key=len).number_of_nodes()
        
    while maxcomp > 2:
        N = graph.number_of_nodes()                     #Nodos iniciales para normalizar
        quant = list(dict(centrality(graph)).values())  #Calculo la centralidad de mis nodos
        val = max(list(quant))                          #Me quedo con la mas grande
        removed = []                                    #Nodos que voy a sacar
        nodos = list(graph.nodes())
        
        #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima
        for i in range(graph.number_of_nodes()):     
            if quant[i] >= val:
                removed.append(nodos[i])                #Me quedo con el que me voy a sacar
                
        graph.remove_nodes_from(removed)
        #Calculo la longitud de la componente máxima y la fracción de nodos que removí
        if graph.number_of_nodes() > 2:
            maxcomp = max(nx.connected_component_subgraphs(graph),key=len).number_of_nodes()
            c.append(len(removed)/N)
            mc.append(maxcomp/N)
        else:
            break
    return mc, c


centralities = [degrees2dict, nx.betweenness_centrality, nx.subgraph_centrality, 
                nx.eigenvector_centrality] 

Removed_nodes = []
Max_comp = []


t0 = datetime.now()

for i in centralities:
    graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
    mc, c = cent_cutoff(graph,i)
    Removed_nodes.append(c)
    Max_comp.append(mc)    

plt.figure(1)
plt.plot(Max_comp[0],Removed_nodes[0],'b-', label='Degrees')
plt.plot(Max_comp[1],Removed_nodes[1],'r-', label='Betweenness')
plt.plot(Max_comp[2],Removed_nodes[2],'y-', label='Subgraph')
#plt.plot(Max_comp[3],Removed_nodes[3],'g.', label='Eigenvector')
plt.legend()
plt.show(1)

print(datetime.now()-t0)    


#%%



def degree_dist(graph,lista):
    











