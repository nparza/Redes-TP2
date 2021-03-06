# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:36:47 2018

@author: noelp
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


#%% TABLA 1 Zotenko - Calcular datos

redes = [G_Y2H,G_LIT,G_APMS]

redes = [  max(nx.connected_component_subgraphs(r),key=len) for r in redes ]


nodes=[]; edges=[]; kmedio=[]; kmax=[]; kmin=[]; 
clustl=[]; clustg=[]; densidad=[]; diam=[];

pos=0
for red in redes:
    nodes.append(red.number_of_nodes())
    edges.append(red.number_of_edges())
    k=degrees(red, node='All')
    kmedio.append(np.mean(k))
    clustl.append(nx.average_clustering(red))
    pos+=1

#%% TABLA 1 Zotenko - Pandas

caract = pd.DataFrame({ 'red':['Y2H','LIT','APMS'], 
                        'N':nodes,
                        'L':edges,
                        'k_medio':kmedio,
                        'c_local':clustl,
                      })

caract = caract[['red','N','L','k_medio','c_local']]

#Uso cols para copiar y pegar y hacer más rapido la redefinición de lugares
cols = list(caract.columns.values)

print(caract)

#%% TABLA 2 Zotenko - Overlap

Overlap = np.empty([len(redes),len(redes)])

def Over(red,red2):
    ##Overlap Red1 vs Red2 (((((((((ES DISTINTO A Red2 vs Red1)))))))))
    
    o = 0
    for r in list(red.nodes()):
        for n in list(red2.nodes()):
            if r == n:
             o += 1
    o = o/(red.number_of_nodes())
    return o

for i in range(np.shape(Overlap)[0]):
    for j in range(np.shape(Overlap)[1]):
        Overlap[i,j] = Over(redes[i],redes[j])
        
caract = pd.DataFrame(Overlap)
print(caract)


#%% FIGURA 1-a Zotenko

def essential(graph, esential):
    es = []
    for n in list(graph.nodes()):
        for m in  esential:
            if n == m:
                es.append(n)
    return es, len(es)
            


def hub_cutoff(graph, nodes):
    N = graph.number_of_nodes()             #Número de nodos del wacho
    frac_specialhubs = []                   #Fracción de hubs especiales
    frac_hubs = []
    deg = 0                                 #Grado mínimo a partir del cual considero a un nodo un hub
    while deg <= max(degrees(graph, node='All')):
        hubs = []                           #Listita de hubs
        for n in graph.nodes():
            #Me fijo si el grado del nodo es igual o mayor a partir del cual puede ser hub
            
            if degrees(graph,node=n) >= deg:
                hubs.append(n)
        #Ahora voy a contar cuantos de esos hubs es esencial.
        if len(hubs) != 0:
            number = 0      
            for n in hubs:
                for m in nodes:
                    if n == m:
                        number += 1
            frac_specialhubs.append(number/len(hubs))
            frac_hubs.append(len(hubs)/N)
            deg +=1
        else:
            deg +=1
    return frac_specialhubs, frac_hubs

f_Y2H, d_Y2H =hub_cutoff(G_Y2H, essential(G_Y2H,Essential)[0])
f_LIT, d_LIT =hub_cutoff(G_LIT, essential(G_LIT,Essential)[0])
f_APMS, d_APMS =hub_cutoff(G_APMS, essential(G_APMS,Essential)[0])

plt.figure(1)
plt.plot(d_Y2H, f_Y2H,'g.', label='Y2H')
plt.plot(d_LIT, f_LIT,'r.', label='LIT')
plt.plot(d_APMS, f_APMS,'b.', label='APMS')
plt.xlabel('Fracción de hubs', fontsize=12)
plt.ylabel('Fracción de hubs especiales',fontsize=12)
plt.legend()
plt.grid(linestyle=':')
plt.show()
   
 
