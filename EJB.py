# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:36:47 2018

@author: noelp
"""

import networkx as nx
import numpy as np
import pandas as pd


def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)	
    return data


def degrees(grafo, node = 'All'):
    if node == 'All':
        lista=list(dict(grafo.degree).values())
    else:
        lista= grafo.degree(node)
    return lista


#%% Importo los datos y creo los grafos

Directorio = 'C:/Users/noelp/Documents/Git/Redes-TP2'

Essential=ldata(Directorio + '/tc02Data/Essential_ORFs_paperHe.txt')
APMS=ldata(Directorio + '/tc02Data/yeast_AP-MS.txt')
LIT=ldata(Directorio + '/tc02Data/yeast_LIT.txt')
LIT_reguly=ldata(Directorio + '/tc02Data/yeast_LIT_Reguly.txt')
Y2H=ldata(Directorio + '/tc02Data/yeast_Y2H.txt')


Redes_list = [Y2H,LIT,APMS]

G_Y2H=nx.Graph(); G_LIT=nx.Graph(); G_APMS=nx.Graph();
G_Y2H.add_edges_from(Y2H); G_LIT.add_edges_from(LIT); G_APMS.add_edges_from(APMS); 


#%% Calculo cosas
redes=[G_Y2H,G_LIT,G_APMS]

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

#%% Grafico tabla

caract = pd.DataFrame({ 'red':['Y2H','LIT','APMS'], 
                        'N':nodes,
                        'L':edges,
                        'k_medio':kmedio,
                        'c_local':clustl,
                      })

caract = caract[['red','N','L','k_medio','c_local']]

#Uso cols para copiar y pegar y hacer más rapido la redefinición de lugares
cols=list(caract.columns.values)

print(caract)

#%%


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






