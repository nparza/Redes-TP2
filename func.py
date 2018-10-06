#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:45:47 2018

@author: sofinico
"""

import networkx as nx


def ldata(archive):
    f = open(archive)
    data = []
    for line in f:
        line = line.strip()
        col = line.split()
        data.append(col)	
    return data


def degrees(grafo, node = 'All'):
    if node == 'All':
        lista = list(dict(grafo.degree).values())
    else:
        lista = grafo.degree(node)
    return lista
