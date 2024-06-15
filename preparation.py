import networkx, random, copy
import pandas as pd
import os, numpy
import types

import networkx as nx
from matplotlib import pyplot as plt
import csv

"""9606.protein.links.v12.0 파일은 string DB의 전체 파일"""

def get_network_from_excel_file(file_name):
    df = pd.read_excel (file_name)
    g = nx.Graph ()
    g.add_nodes_from (df['protein1'])
    g.add_nodes_from (df['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in df.iterrows ()]
    g.add_edges_from (edges)

    return g

def whole_network_to_csv(_file_name):
    # N = 30
    # with open (_file_name, encoding='utf-8') as myfile:
    #     head = [next (myfile) for x in range (N)]
    # print (len (head))
    # for l in head:
    #     print (l)
    file = open(_file_name, "r")
    txt_dataframe = pd.DataFrame()
    table_data = []
    for line in file:
        tmp_array = line.strip ().split (" ")
        while tmp_array.count (""):
            tmp_array.remove ("")

        table_data.append (tmp_array)
    txt_dataframe = pd.DataFrame(table_data)
    txt_dataframe = txt_dataframe.rename(columns=txt_dataframe.iloc[0])
    txt_dataframe =  txt_dataframe.drop(txt_dataframe.index[0])
    # print(txt_dataframe.shape)
    txt_dataframe['combined_score'] = txt_dataframe['combined_score'].astype(int)

    """경량화, combined score > 400 인 edge만 사용"""
    txt_dataframe = txt_dataframe.loc[txt_dataframe.combined_score>400]
    txt_dataframe.reset_index(drop = True, inplace = True)
    """파일읽는 속도가 빠른 json으로 파일 보내기"""
    txt_dataframe.to_csv("Whole_network.csv",sep= ',',na_rep = 'NaN', index = False)
    return

whole_network_to_csv("./Data/Whole STRING Network/9606.protein.links.v12.0.txt")
"""역할 다 함. json 파일로 STRING 경량화 버젼 저장함"""

_file_name = "test_1.sif"

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """
	Parse sif file into node and edge sets and dictionaries
	returns setNode, setEdge, dictNode, dictEdge
	store_edge_type: if True, dictEdge[(u,v)] = edge_value
	delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set ()
    setEdge = set ()
    dictNode = {}
    dictEdge = {}
    flag = False
    f = open (file_name)
    for line in f:
        # 여기서 line은 'str' class이다.
        words = line.rstrip ("\n").split ()
        if len (words) == 3:
            setNode.add (words[0])
            setEdge.add ((words[0], words[2]))
            score = words[1]
            dictNode[words[0]] = score
            dictEdge[(words[0], words[2])] = words[1]
        else:
            setNode.add (words[0])
    f.close ()
    # print(len(setNode),len(setEdge))#node,edge counting
    return setNode, setEdge, dictNode, dictEdge


## 여기까지는 SIF file reading 완료함

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph ()
    else:
        g = networkx.Graph ()
    return g


def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file (network_file_in_sif,
                                                                               store_edge_type=use_edge_data,
                                                                               delim=delim)
    g = create_graph ()
    if include_unconnected:
        g.add_nodes_from (setNode)
    if use_edge_data:
        for e, w in dictEdge.iteritems ():
            u, v = e
            g.add_edge (u, v, w=w)  # ,{'w':w})
    else:
        g.add_edges_from (setEdge)
    return g

# print(f"Number of nodes: {len(G.nodes())}")
# print(f"Number of edges: {len(G.edges())}")
# average_clustering = nx.average_clustering(G)
# print(f"Average Clustering Coefficient: {average_clustering}")
# density = nx.density(G)
# print(f"Graph Density: {density}")


"""
네트워크 단일 분석 

"""


# # Compute degree centrality
# degree_centrality = nx.degree_centrality(G)
# print(degree_centrality)\
# # Compute closeness centrality
# closeness_centrality = nx.closeness_centrality(G)
# sorted_closeness_centrality = sorted(closeness_centrality.items(), key=lambda x:x[1], reverse=True)
# print("Closeness Centrality:", sorted_closeness_centrality)


"""Distance calculation"""

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length (G, source_id, target_id)


def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                val = get_shortest_path_length_between (network, node_from, node_to)
                values.append (val)
                d = min (values)
                # print d,
                values_outer.append (d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append (val)
                d = min (values)
                values_outer.append (d)
            d = numpy.mean (values_outer)
        # print d
    return d


def Herb_node_checking