#!/usr/bin/env python3
# -*- coding: utf-8 -*-

label = 'Y2H'

centralities = ['random','degree','shortest-path','eigenvector','subgraph','essentials']

## Dejo los datos en el formato que los deseo guardar

lines_out1 = dict()
lines_out2 = dict()

for c in centralities:
    lines_out1[c]=[]
    lines_out2[c]=[]
    for i in range(len(removed_nodes[c])):
        lines_out1[c].append('%.6f'%removed_nodes[c][i])
        lines_out2[c].append('%.6f'%max_comp[c][i])
    lines_out1[c] = ','.join(lines_out1[c])
    lines_out2[c] = ','.join(lines_out2[c])

f = open('fig3_'+label+'.csv','w')
    
for c in centralities:
    print(c+','+'%s' % lines_out1[c], file=f)
    print(c+','+'%s' % lines_out2[c], file=f)

f.close()