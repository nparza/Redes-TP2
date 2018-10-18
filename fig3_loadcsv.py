#!/usr/bin/env python3
# -*- coding: utf-8 -*-

label = 'LITr'

entrada = open('fig3_'+label+'.csv')

lines = []

for line in entrada:
    lines.append(line.split(','))
    
for line in lines:
    for i in range(1,len(line)):
        line[i] = float(line[i])
i = 0
x = dict(); y = dict()
while i < len(lines):
    x[lines[i][0]] = lines[i][1:]
    y[lines[i+1][0]] = lines[i+1][1:]
    i += 2

#%% FUNCS

from matplotlib import pyplot as plt

def plot(label,color):
    plt.plot(x[label],
             y[label],
             color = color, 
             label = label, 
             linewidth = '2')
    
def applyPlotStyle():
    plt.xlabel('fracción de nodos',weight='bold',fontsize=12)
    plt.ylabel('componente más grande',weight='bold',fontsize=12)
    plt.title(label)
    #plt.grid(linestyle=':')
    plt.legend()
    

#%%  GRAFICA UNA RED

plt.figure(20)

plot('random','gray')
plot('degree','dodgerblue')
plot('shortest-path','orangered')
plot('eigenvector','lime')
plot('subgraph','gold')

plt.plot(x['essentials'],
         y['essentials'],
         color = 'darkviolet', 
         label = 'essentials',
         marker = 'h')

plt.xlim([0,0.5])
applyPlotStyle()
plt.show(20)    


#%% FIGURA INFORME - CUATRO REDES

def plot(label,color):
    plt.plot(x[label],
             y[label],
             color = color, 
             label = label, 
             linewidth = '1.5')
    
def applyPlotStyle():
    plt.xlabel('fracción de nodos removidos',fontsize=10)
    plt.ylabel('componente más grande',fontsize=10)
    plt.title(label,weight='bold',fontsize=12,loc='left')
    #plt.grid(linestyle=':')
    plt.legend(loc='best',fontsize=7)

redes = ['Y2H','APMS','LIT','LITr']
ejex = [2.5,5,4,3]

plt.figure(101)

for X in range(4):
    
    label = redes[X]
    entrada = open('fig3_'+label+'.csv')
    lines = []
    
    for line in entrada:
        lines.append(line.split(','))
        
    for line in lines:
        for i in range(1,len(line)):
            line[i] = float(line[i])
    i = 0
    x = dict(); y = dict()
    while i < len(lines):
        x[lines[i][0]] = lines[i][1:]
        y[lines[i+1][0]] = lines[i+1][1:]
        i += 2
    
    plt.subplot(2,2,X+1)
    plot('random','gray')
    plot('degree','dodgerblue')
    plot('shortest-path','orangered')
    plot('eigenvector','lime')
    plot('subgraph','gold')
    plt.plot(x['essentials'],
             y['essentials'],
             color = 'darkviolet', 
             label = 'essentials',
             marker = 'h')
    #plt.xlim([0,ejex[X]/10])
    plt.xlim([0,0.5])
    applyPlotStyle()
    if label == 'LITr':
        plt.title('LIT Reguly',weight='bold',fontsize=12,loc='left')
        
plt.show()




