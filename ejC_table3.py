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

#%% 

# Agrego atributo de esencialidad a los nodos de cada red

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
        
# Agrego atributos al grafo relacionados a sus nodos esenciales
            
for graph in graphs:
    
    G = graph.copy()
    cg = max(nx.connected_component_subgraphs(G), key=len)
    ess = []
    degree_ess = []
    for node in cg:
        if cg.nodes[node]['Essential'] == True:
            ess.append(node)
            degree_ess.append(cg.degree[node])
    
    # Lista de nodos esenciales de su componente gigante
    
    graph.graph['cg_ess'] = ess
    
    # Distrubución de grado de los nodos esenciales de su comp. gigante
    
    graph.graph['cg_ess_deg'] = degree_ess

    
#%% 
    
import random

''' Criterio para remover nodos no-esenciales con la misma distribución 
de grado que los esenciales '''

def nearest(degree,K):
    
    # Devuelve el grado más cercano a degree 
    # que tiene nodos no esenciales para remover
    
    for key,value in list(K.items()):
        if value == []:
            del K[key]
    
    if degree in list(K.keys()):
        return degree
    else:
        listk = list(K.keys())
        listk.append(degree)
        ordenada = sorted(listk)
        i = ordenada.index(degree)
        if i == len(ordenada)-1:
            return ordenada[i-1]
        dist_izq = ordenada[i] - ordenada[i-1]
        dist_der = ordenada[i+1] - ordenada[i]
        if dist_izq < dist_der:
            return ordenada[i-1]
        elif dist_izq == dist_der:
            return ordenada[random.choice([i-1,i+1])]
        else:
            return ordenada[i+1]
        

def remove_random(graph):
    
    ## Devuelve una lista de nodos no-esenciales
    
    G = graph.copy()
    cg = max(nx.connected_component_subgraphs(G), key=len)
    rand = []
    
    K = dict()
    for i in set(dict(cg.degree).values()):
        K[i] = []
    for node in list(cg.nodes()):
        if cg.nodes[node]['Essential'] == False:
            K[graph.degree[node]].append(node)
    
    Kess = dict()
    for i in set(graph.graph['cg_ess_deg']):
        Kess[i] = 0
    for k in graph.graph['cg_ess_deg']:
        Kess[k] += 1
    
    for degree in Kess: 
        for m in range(Kess[degree]):
            near_deg = nearest(degree,K)
            node = random.choice(K[near_deg])
            rand.append(node)
            K[near_deg].remove(node)
    
    return rand
            
#%% Construyo los datos para la TABLA 3

# Tamaño de la componente gigante post remover los nodos esenciales

cg_es = dict()

# Tamaño de la componente gigante post remover nodos random con la misma
# distribución de grado

cg_rand = dict()

# Las dos listas van a estar normalizadas por el número de nodos de la cg original


for graph in graphs:
    
    CG_orishinal = max(nx.connected_component_subgraphs(graph), key=len)
    N_orishinal = CG_orishinal.number_of_nodes()
    
    
    # Remuevo nodos esenciales
    
    GRAPH = graph.copy()
    CG = max(nx.connected_component_subgraphs(GRAPH), key=len)
    CG.remove_nodes_from(graph.graph['cg_ess'])
    N = max(nx.connected_component_subgraphs(CG), key=len).number_of_nodes()
    
    cg_es[graph.graph['label']] = N/N_orishinal
    
    
    # Remuevo nodos no-esenciales
    
    GRAPH = graph.copy()
    CG = max(nx.connected_component_subgraphs(GRAPH), key=len)
    CG.remove_nodes_from(remove_random(graph))
    N = max(nx.connected_component_subgraphs(CG), key=len).number_of_nodes()
    
    cg_rand[graph.graph['label']] = N/N_orishinal
   

    
    
    
    
    
    



    
    




    
    
    
    
    
    
    
    
    
    
    
    