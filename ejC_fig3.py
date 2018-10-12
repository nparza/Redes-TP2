#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:43:24 2018

@author: sofinico
"""

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from func import *

#%% Importo los datos y creo los grafos

Essential = ldata('TC02_data/Essential_ORFs_paperHe.txt')
Essential = [Essential[m][1] for m in range(2,len(Essential[2:])-2)]

APMS = ldata('TC02_data/yeast_AP-MS.txt')
LIT = ldata('TC02_data/yeast_LIT.txt')
Y2H = ldata('TC02_data/yeast_Y2H.txt')

LITr = ldata('TC02_data/yeast_LIT_Reguly.txt')
LITr = [row[:2] for row in LITr]
LITr.pop(0)

G_Y2H = nx.Graph(); G_Y2H.add_edges_from(Y2H);

G_LIT = nx.Graph(); G_LIT.add_edges_from(LIT); 

G_LITr = nx.Graph(); G_LITr.add_edges_from(LITr)

G_APMS = nx.Graph(); G_APMS.add_edges_from(APMS); 

graphs = [G_Y2H,G_LIT,G_LITr,G_APMS]

#%% Funciones 

def degrees2dict(graph):
    
    d = dict(graph.degree)
    return d


def cent_cutoff(graph, centrality, connected=False):
    
    c = []          # Fracción de nodos removidos en cada iteración
    mc = []         # Nodos de la componente máxima
    
    maxcomp = max(nx.connected_component_subgraphs(graph),
                  key=len).number_of_nodes()
   
    N = graph.number_of_nodes()     # Nodos iniciales para normalizar    
    
    mc.append(maxcomp/N)
    c.append(0)
    
    while maxcomp > 2:
        
        # Calculo centralidad de los nodos
        
        if not connected:
            quant = list(dict(centrality(graph)).values())
            nodos = list(graph.nodes())
        
        if connected:
            quant = list(dict(centrality(max(nx.connected_component_subgraphs(graph), 
                                             key=len))).values())
            nodos = list(max(nx.connected_component_subgraphs(graph),
                  key=len).nodes())

        
        val = max(list(quant))      # Centralidad más grande
        removed = []                # Nodos que voy a sacar
        
        # Recorro la lista de nodos eligiendo los que tienen centralidad máxima
        # y los guardo en una lista
        
        for i in range(len(nodos)):     
            if quant[i] >= val:
                removed.append(nodos[i])
                
        graph.remove_nodes_from(removed)
        
        # Calculo la longitud de la componente máxima y 
        # la fracción de nodos que removí
        
        if graph.number_of_nodes() > 2:
            maxcomp = max(nx.connected_component_subgraphs(graph),
                          key=len).number_of_nodes()
            c.append(c[-1]+len(removed)/N)
            mc.append(maxcomp/N)
        else:
            break
    
    return mc, c


#%% 

''' Cuidado al correr esto porque redefine 
diccionarios que pueden tener cosas '''

removed_nodes = dict()
max_comp = dict()
centrality_type = dict()


#%% DEGREES 

cent = degrees2dict
label = 'degree'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'local'
removed_nodes[label] = c
max_comp[label] = mc


#%% SHORTEST-PATH BETWEENNESS - 5 min 

cent = nx.betweenness_centrality
label = 'shortest-path'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'betweenness'
removed_nodes[label] = c
max_comp[label] = mc


#%% SUBGRAPH - 1.5 min 

cent = nx.subgraph_centrality
label = 'subgraph'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'local'
removed_nodes[label] = c
max_comp[label] = mc

#%% CLOSENESS - 10 min

cent = nx.closeness_centrality
label = 'closeness'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'NA'
removed_nodes[label] = c
max_comp[label] = mc

#%% CURRENT FLOW BETWEENNESS - 9 min

cent = nx.current_flow_betweenness_centrality
label = 'current flow'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent,connected=True)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'betweenness'
removed_nodes[label] = c
max_comp[label] = mc

#%% EIGENVECTOR - 40 s

cent = nx.eigenvector_centrality_numpy
label = 'eigenvector'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'local'
removed_nodes[label] = c
max_comp[label] = mc


#%% Grafico - FIGURA 3 Zotenko

def plot(label,color):
    plt.plot(removed_nodes[label],
             max_comp[label],
             color = color, 
             label = label,
             linewidth = '1')
    
def applyPlotStyle():
    plt.xlabel('fracción de nodos',weight='bold',fontsize=11)
    plt.ylabel('componente más grande',weight='bold',fontsize=11)
    #plt.grid(linestyle=':')
    plt.legend()
    
plt.figure(10)
plot('degree','r')
plot('shortest-path','magenta')
plot('subgraph','lime')
plot('closeness','gold')
plot('current flow','cyan')
plot('eigenvector','darkgoldenrod')
applyPlotStyle()
plt.show(10)
  

