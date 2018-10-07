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
from func import *

#%% Importo los datos y creo los grafos


Essential = ldata('TC02_data/Essential_ORFs_paperHe.txt')
Essential = [Essential[m][1] for m in range(2,len(Essential[2:])-2)]
APMS = ldata('TC02_data/yeast_AP-MS.txt')
LIT = ldata('TC02_data/yeast_LIT.txt')
LIT_reguly = ldata('TC02_data/yeast_LIT_Reguly.txt')
Y2H = ldata('TC02_data/yeast_Y2H.txt')


Redes_list = [Y2H,LIT,APMS]

G_Y2H=nx.Graph(); G_LIT=nx.Graph(); G_APMS=nx.Graph();
G_Y2H.add_edges_from(Y2H); G_LIT.add_edges_from(LIT); G_APMS.add_edges_from(APMS); 


#%% Funciones 

def degrees2dict(graph):
    d = dict(graph.degree)
    return d

def cent_cutoff(graph, centrality):
    c = [  ]            #Fracción de nodos removidos en cada iteración
    mc = [  ]           #Nodos de la componente máxima
    
    maxcomp = max(nx.connected_component_subgraphs(graph),key=len).number_of_nodes()
    N = graph.number_of_nodes()                     #Nodos iniciales para normalizar    
    
    mc.append(maxcomp/N)
    c.append(0)
    
    while maxcomp > 2:
        
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


#%% DEGREES 

cent = degrees2dict
label = 'degree'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

removed_nodes[label] = c
max_comp[label] = mc


#%% BETWEENNESS 

cent = nx.betweenness_centrality
label = 'betweenness'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

removed_nodes[label] = c
max_comp[label] = mc


#%% SUBGRAPH 

cent = nx.subgraph_centrality
label = 'subgraph'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(G_APMS),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

removed_nodes[label] = c
max_comp[label] = mc


#%% Grafico - FIGURA 3 Zotenko

def plot(label,color):
    plt.plot(removed_nodes[label],
             max_comp[label],
             color = color, 
             label = label)
    
def applyPlotStyle():
    plt.xlabel('fracción de nodos',weight='bold',fontsize=11)
    plt.ylabel('componente más grande',weight='bold',fontsize=11)
    plt.grid(linestyle=':')
    plt.legend()
    
plt.figure(10)
plot('degree','b')
plot('betweenness','r')
plot('subgraph','y')
applyPlotStyle()
plt.show(10)
  

