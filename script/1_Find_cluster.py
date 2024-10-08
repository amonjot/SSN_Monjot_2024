#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/07/2024

@authors: jeremy & arthur
"""

# Import modules
import igraph as ig
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import sys

# Set working directory
os.chdir("SSN_env/results/")

# Set inputs
condition_name = sys.argv[1]
edges_name = "diamond_ssn_"+condition_name+".edges"
nodes_name = "attributes/Metadata_Unicellular.attributes"
#print(edges_name)
#print(condition_name)
#print(nodes_name)

# Import dataframes containing ...
# ... edges
edges = pd.read_csv(edges_name, sep = ";")
# ... nodes
nodes = pd.read_csv(nodes_name, sep = ";")

# Create graph from the two dataframes
g = ig.Graph.DataFrame(edges, directed=False, use_vids=False, vertices = nodes)
del_vert_isol = [v.index for v in g.vs if v.degree() == 0]

# Create dataframe containing information about the nodes retained for centrality
lst_vs_attrib = g.vertex_attributes()
lst_vs_attrib.insert(0, "CC")
# Centrality dataframe
df_central_seq = pd.DataFrame(columns = lst_vs_attrib)
# Dataframe with all connected component
df_cc = pd.DataFrame(columns=lst_vs_attrib)

## Decomposition into sub-graphs with at least n vertices (here, n=3)
g_sub = g.decompose(minelements=3)
#g_sub = g.decompose()
#print(g.es.attributes())
lst_attribut = ['alignment_len', 'pident', 'bitscore']

count=0
for i in g_sub:

##    recovery of the most central node

#    lst_of_lst_centricity = []
#    for attrib in lst_attribut:
#        lst_centricity = i.eigenvector_centrality(weights = attrib, scale = True)
#        lst_of_lst_centricity.append(lst_centricity)
#    arr = np.array(lst_of_lst_centricity)
#    max_lst = arr.max(axis=1).reshape(-1,1)
#    res_lst_max = list(np.where(arr == max_lst,1,0))
#    lst_sum = np.sum(res_lst_max, axis=0)
#    central_node = lst_sum.argmax(axis=0)

## Create dataframe containing information for each node in each component

    df_sub_cc = i.get_vertex_dataframe()
    df_sub_cc.insert(0, "CC", count)
    df_cc = pd.concat([df_cc, df_sub_cc])

## Create dataframe containing sequences representative of cc

#    df_central_seq = pd.concat([df_central_seq, df_sub_cc.iloc[[central_node]]])
    count += 1


## Statistic files

with open(f"ssn_statistique_{condition_name}.txt", "w") as f:
    f.write("Statistics\n")
    f.write(f"Number of initial sequences: {len(nodes)}\n")
    f.write(f"Number of sequences retained after filtration: {len(df_cc)}\n")
    f.write(f"Number of connected components: {len(g_sub)}\n")
#    f.write(f"Number of reference sequences: {len(df_central_seq)}")


## Save network and connected componant informations

ig.Graph.write_pickle(g, f"ssn_graph_{condition_name}.pickle")
df_cc.to_csv(f"ssn_composantes_connexes_graph_{condition_name}.csv", index=False, sep=";")
#df_central_seq.to_csv(f"ssn_representative_sequence_{condition_name}.csv", index=False, sep=";")


