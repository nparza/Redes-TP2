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
    
    # Lista de nodos esenciales de su componente gigante
    
    graph.graph['cg_ess'] = ess


#%% Funciones 

def degrees2dict(graph):
    
    d = dict(graph.degree)
    return d


def cent_cutoff(graph, centrality, connected=False):
    
    ''' 
        Remueve los nodos con centralidad máxima a la componente más grande
        de un grafo hasta que ésta tenga tres nodos.
    
        INPUTS
        graph: es la componente gigante de algún grafo
        centrality: función centralidad que se utiliza
        connected: True si la función centralidad sólo puede
        ser calculada en un grafo conexo
        
        OUTPUTS
        mc (lista): nodos de la componenta máxima en cada iteración, 
        normalizados por la cantidad de nodos originales
        c (lista): total de nodos removidos a medida que pasan las 
        iteraciones, normalizados por la cantidad de nodos originales
    '''
    
    c = []; mc = []         
   
    N = graph.number_of_nodes() 
    G = graph.copy()
    maxcomp = G.number_of_nodes()
    
    mc.append(maxcomp/N)
    c.append(0)
    
    while maxcomp > 2:
        
        # Calculo centralidad de los nodos
        
        if not connected:
            quant = list(dict(centrality(G)).values())
            nodos = list(G.nodes())
        
        if connected:
            quant = list(dict(centrality(max(nx.connected_component_subgraphs(G), 
                                             key=len))).values())
            nodos = list(max(nx.connected_component_subgraphs(G),
                             key=len).nodes())

        
        val = max(list(quant))      # Centralidad más grande
        removed = []                # Nodos que voy a sacar
        
        # Recorro la lista de nodos eligiendo los que tienen centralidad máxima
        # y los guardo en una lista
        
        for i in range(len(nodos)):     
            if quant[i] >= val:
                removed.append(nodos[i])
                
        G.remove_nodes_from(removed)
        
        # Calculo la longitud de la componente máxima y 
        # la fracción de nodos que removí
        
        if G.number_of_nodes() > 2:
            maxcomp = max(nx.connected_component_subgraphs(G),
                          key=len).number_of_nodes()
            c.append(c[-1]+len(removed)/N)
            mc.append(maxcomp/N)
        else:
            break
    
    return mc, c

#%%
    
import random

def rand_cutoff(graph, step):
    
    ''' 
        Remueve una cantidad (step) de nodos random a la componente más grande
        de un grafo hasta que ésta tenga tres nodos.
    
        INPUTS
        graph: es la componente gigante de algún grafo
        step: cantidad de nodos que remueve en una iteración
        
        OUTPUTS
        list_maxcomp (lista): nodos de la componenta máxima en cada iteración, 
        normalizados por la cantidad de nodos originales
        fraction (lista): total de nodos removidos a medida que pasan las 
        iteraciones, normalizados por la cantidad de nodos originales
    '''
    
    list_maxcomp = []; fraction = []              
    
    N = graph.number_of_nodes()
    maxcomp = graph.number_of_nodes()
    G = graph.copy()
    
    list_maxcomp.append(1)
    fraction.append(0)
    
    while maxcomp > 2:
        
        nodos_disp = list(G.nodes())
        removed = []
        
        for i in range(step):
            node = random.choice(nodos_disp)
            removed.append(node)
            nodos_disp.remove(node)
        
        G.remove_nodes_from(removed)
        maxcomp = max(nx.connected_component_subgraphs(G),
                      key=len).number_of_nodes()
        
        list_maxcomp.append(maxcomp/N)
        fraction.append(fraction[-1]+len(removed)/N)
    
    return list_maxcomp, fraction
        

#%% 

''' Cuidado al correr esto porque redefine 
diccionarios que pueden tener cosas '''

removed_nodes = dict()
max_comp = dict()
centrality_type = dict()

RED = G_Y2H

#%% ESSENTIALS - Remuevo los nodos esenciales

label = 'essentials'

G = RED.copy()
CG = max(nx.connected_component_subgraphs(G), key=len)
N = CG.number_of_nodes()
CG.remove_nodes_from(RED.graph['cg_ess'])
n = max(nx.connected_component_subgraphs(CG), key=len).number_of_nodes()


removed_nodes[label] = [len(graph.graph['cg_ess'])/N]
max_comp[label] = [n/N]

#%% RANDOM - 1.5 min x 10

label = 'random'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len)

C = []
MC = []
for i in range(10):
    mc, c = rand_cutoff(graph, 1)
    C.append(c)
    MC.append(mc)
    print('tarda en correr: ', datetime.now()-ti)
    

l = len(C[0])
for i in range(len(C)):
    if len(C[i]) < l:
        l = len(C[i])
        indice = i

c_prom = []
for j in range(l):
    suma = 0
    for i in range(len(C)):
        suma += C[i][j]
    c_prom.append(suma/len(C))
        
    
centrality_type[label] = 'NA'
removed_nodes[label] = c_prom
max_comp[label] = MC[indice]


#%% DEGREES 

cent = degrees2dict
label = 'degree'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'local'
removed_nodes[label] = c
max_comp[label] = mc


#%% SHORTEST-PATH BETWEENNESS - 5 min 

cent = nx.betweenness_centrality
label = 'shortest-path'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'betweenness'
removed_nodes[label] = c
max_comp[label] = mc


#%% SUBGRAPH - 1.5 min 

cent = nx.subgraph_centrality
label = 'subgraph'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'local'
removed_nodes[label] = c
max_comp[label] = mc

#%% CLOSENESS - 10 min

cent = nx.closeness_centrality
label = 'closeness'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
mc, c = cent_cutoff(graph,cent)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'NA'
removed_nodes[label] = c
max_comp[label] = mc

#%% CURRENT FLOW BETWEENNESS - 9 min

cent = nx.current_flow_betweenness_centrality
label = 'current flow'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
mc, c = cent_cutoff(graph,cent,connected=True)
print('tarda en correr: ', datetime.now()-ti)

centrality_type[label] = 'betweenness'
removed_nodes[label] = c
max_comp[label] = mc

#%% EIGENVECTOR - 40 s

cent = nx.eigenvector_centrality_numpy
label = 'eigenvector'

ti = datetime.now()
graph = max(nx.connected_component_subgraphs(RED),key=len) 
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
             linewidth = '2')
    
def applyPlotStyle():
    plt.xlabel('fracción de nodos',weight='bold',fontsize=11)
    plt.ylabel('componente más grande',weight='bold',fontsize=11)
    #plt.grid(linestyle=':')
    plt.legend()
    
plt.figure(10)

plot('random','gray')
plot('degree','dodgerblue')
plot('shortest-path','orangered')
plot('eigenvector','lime')
plot('subgraph','gold')
#plot('closeness','magenta')
#plot('current flow','r')

plt.plot(removed_nodes['essentials'],
         max_comp['essentials'],
         color = 'darkviolet', 
         label = 'essentials',
         marker = 'h')

plt.xlim([0,0.5])
applyPlotStyle()
plt.show(10)
  
    
