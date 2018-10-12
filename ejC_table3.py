#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 17:46:16 2018

@author: sofinico
"""

import networkx as nx
import numpy as np
import pandas as pd
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

#%% Label de cada grafo

G_Y2H.graph['label'] = 'Y2H'
G_LIT.graph['label'] = 'LIT'
G_LITr.graph['label'] = 'LIT_r'
G_APMS.graph['label'] = 'APMS'

#%% Agrego atributo de esencialidad a los nodos de cada red

for graph in graphs:
    
    for node in list(graph.nodes()):
        
        essential = False
        i = 0
        
        while essential == False and i < len(Essential):
            
            if Essential[i] == node:
                essential = True
            
            i += 1
        
        graph.nodes[node]['Essential'] = essential
        
#%% 

''' Criterio para remover nodos no-esenciales con la misma distribución 
de grado que los esenciales '''

def remove_random(graph,rm_es):


        
        

#%% Construyo los datos para la TABLA 3

# Tamaño de la componente gigante post remover los nodos esenciales

cg_es = dict()

# Tamaño de la componente gigante post remover nodos random con la misma
# distribución de grado

cg_rand = dict()

# Las dos listas van a estar normalizadas por el número de nodos de la cg original


for graph in graphs:
    
    rm_es = []
    cg = max(nx.connected_component_subgraphs(graph), key=len)
    N = cg.number_of_nodes()  
    
    # Remuevo los esenciales y los guardo en una lista
    
    nodos = list(cg.nodes())
    
    for i in range(len(nodos)):      
       
        if cg.nodes[nodos[i]]['Essential'] == True:
            cg.remove_node(nodos[i])
            rm_es.append(nodos[i])
 
    n = max(nx.connected_component_subgraphs(cg), key=len).number_of_nodes()
    
    cg_es[graph.graph['label']] = n/N

    
    
    




















