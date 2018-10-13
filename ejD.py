# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 17:33:49 2018

@author: noelp
"""

import networkx as nx
import numpy as np
import random as rnd
from matplotlib import pyplot as plt
from func import *
from datetime import datetime
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
        
##Rewiring
        

def rewire(graph):
    edges = list(dict(graph.edges()))
    
    ##Mezclo los edges
    a = [edges[i][0] for i in range(len(edges))]
    b = [edges[i][1] for i in range(len(edges))]
    rnd.shuffle(a)
    rnd.shuffle(b)
    
    #Creo el nuevo grafo
    newedges = [(a[i],b[i]) for i in range(len(edges))]
    newG = nx.Graph()
    newG.add_edges_from(newedges)
    
    ###Agrego la esencialidad
    for node in list(newG.nodes()):    
        essential = False
        i = 0    
        while essential == False and i < len(Essential):        
            if Essential[i] == node:
                essential = True       
            i += 1    
        newG.nodes[node]['Essential'] = essential
    #Listorti
    
    return newG

def ibeps(graph):
    IBEPS = 0
    edges = list(dict(graph.edges()))
    for e in edges:
        e0 = graph.nodes[e[0]]['Essential']
        e1 = graph.nodes[e[1]]['Essential']
        if e0 == True and e0==e1:
            IBEPS += 1
    return IBEPS

#%%
    
### Esto te da la cantidad de IBEPS para redes generadas random


t0 = datetime.now()

    
times = 100

graph = G_Y2H
IBEPS = []
measuredIBEPS = ibeps(graph)

for t in range(times):
    newGraph = rewire(graph)
    IBEPS.append(ibeps(newGraph))
    
    
#freq, binedges = np.histogram(IBEPS, bins=20)
#norm = sum(freq)
#freq_normed = [i/norm for i in freq]
#bincenters = 0.5*(binedges[1:]+binedges[:-1])
#lins = {'linestyle': 'None'}
#
#plt.figure(1)
#plt.rc('lines', **lins)
#plt.axvline(x=measuredIBEPS, color='k', linestyle='dashed', linewidth=1)
#plt.bar(bincenters,freq_normed, color='orange',edgecolor='black', linewidth=1.2, width = np.diff(binedges))
#plt.grid(True)
#plt.ylabel('Frecuencia')
#plt.xlabel('Fracción de parejas heteronormativas' )
#plt.title('IBEPS medidos: %s' %(round(measuredIBEPS,2) ))
#plt.show(1)
    
print(datetime.now() - t0)    


