# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 15:21:28 2018

@author: noelp
"""

import networkx as nx
import numpy as np
import random as rnd
from matplotlib import pyplot as plt
from func import *
from datetime import datetime
from scipy.optimize import curve_fit
import pandas as pd
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

#%%}


for g in graphs:
    g.remove_edges_from(g.selfloop_edges())


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
        
def essential(graph, esential):
    es = []
    for n in list(graph.nodes()):
        for m in  esential:
            if n == m:
                es.append(n)
    return es, len(es)

#%%
    
#graphs = [  max(nx.connected_component_subgraphs(r),key=len) for r in graphs ]
graphs= [G_Y2H,G_LIT,G_LITr]
t0 = datetime.now()
linear_nonPek = dict()
parameters = dict()
cut = 11

for g in graphs:
    nodes = list(g.nodes())
    dist = np.zeros(cut)
    nodesk = np.zeros(cut)
    
    for n in nodes:
        if g.nodes[n]['Essential'] == True and g.degree[n] < cut:
            dist[g.degree[n]] += 1
            nodesk[g.degree[n]] += 1
        elif g.degree[n] < cut:
            nodesk[g.degree[n]] += 1

    linear_nonPek[g.graph['label']] = [np.log(1-dist[i]/nodesk[i]) for i in range(cut)]
    


def lineal(k,a,b):
    return a*k+b

x= list(range(1,cut))

for g in graphs:
    plt.figure(graphs.index(g))
    l = g.graph['label']
    c = ['g.','r.','y.']
    d = ['g-','r-','y-']
    plt.plot(x,linear_nonPek[l][1:],c[graphs.index(g)], label=l)
    plt.legend()
    plt.show(graphs.index(g))

### Voy a sacar puntos flasheros

linear_nonPek['Y2H'].pop(8)
linear_nonPek['LIT'].pop(-1)


Y2Hx = x.copy() ; Y2Hx.pop(8)
LITx = x.copy(); LITx.pop(-1)
LIT_rx = x.copy(); 

ax = [Y2Hx,LITx,LIT_rx]

for g in graphs:
    P0 = [-0.03,-0.13]
    popt, pcov = curve_fit(lineal, ax[graphs.index(g)], linear_nonPek[g.graph['label']][1:], p0=P0)    
    parameters[g.graph['label']] = (popt[0], popt[1])

for g in graphs:
    l = g.graph['label']
    c = ['g.','r.','y.']
    d = ['g-','r-','y-']
    plt.plot(ax[graphs.index(g)],linear_nonPek[l][1:],c[graphs.index(g)], label=l+' alpha=%s beta=%s'% (round(parameters[l][0],3),round(parameters[l][0],3)))
    plt.plot(ax[graphs.index(g)],[lineal(k,parameters[l][0],parameters[l][1]) for k in ax[graphs.index(g)]],d[graphs.index(g)])
    plt.xlabel('Grado',fontsize=12)
    plt.ylabel('ln(1-Pe)',fontsize=12)
    plt.legend()
    plt.grid(linestyle=':')

#%%

redes = [  max(nx.connected_component_subgraphs(r),key=len) for r in graphs ]

Pares = []
Iguales = []
Esperados = []
 
   
for g in redes:
    nodes = list(g.nodes())
    pairs=0
    equalpairs=0
    for n in nodes:
        for k in nodes:
            if n!=k and n not in list(g.neighbors(k)):
                common = 0
                for i in list(g.neighbors(n)):
                    common += list(g.neighbors(k)).count(i)
                if common >= 3:
                    pairs +=1
                    if g.nodes[n]['Essential'] == g.nodes[k]['Essential']:
                        equalpairs +=1
    Pares.append(pairs/2)
    Iguales.append(equalpairs/2)

    
    
    
    
for g in redes:
    P = 0
    a, b = [1 - np.exp(i) for i in parameters[g.graph['label']]]
    degs = dict(g.degree())
    keys = list(dict(g.degree()).keys())
    N = 0
    
    for n in keys:
        if g.degree(n) < cut:        
            Pn = 1 - ((1-a)**degs[n])*(1-b)
            for m in keys:
                if g.degree(m) < cut and m!=n:
                    Pm = 1 - ((1-a)**degs[m])*(1-b)
                    P += Pn*Pm + (1-Pn)*(1-Pm)
                    N += 2
    Esperados.append(P/(2*N))


Esperados = [Pares[i]*Esperados[i] for i in range(len(Esperados))]
    
 #%%               
caract = pd.DataFrame({ 'red':['Y2H','LIT','LIT_r'], 
                        'Total number of pairs':Pares,
                        'Number of pairs of the same type':Iguales,
                        'Expected equal pairs':Esperados,
                         })

caract = caract[['red','Total number of pairs','Number of pairs of the same type','Expected equal pairs']]

#Uso cols para copiar y pegar y hacer más rapido la redefinición de lugares
cols = list(caract.columns.values)

print(caract)
print(caract.to_latex())           
        
        
#%%








