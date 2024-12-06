import pandas as pd
import networkx as nx
import random, copy
import TM_MC
import disease_DB
import time
import numpy
from itertools import combinations
from varname import nameof
from collections import defaultdict

def changing_Gene_name_to_ENSP_ID(Gene_list):
    chemical_protein_index = pd.read_csv("./Data/Disease/protein_dict.csv")
    ENSP_code_of_protein = []
    Gene_error_list = []
    for Gene_names in Gene_list:
        ENSP_code_of_protein_1 = chemical_protein_index[chemical_protein_index.preferred_name == Gene_names]["#string_protein_id"].tolist()
        # print(ENSP_code_of_protein_1)
        try :
            ENSP_code_of_protein_1 = ENSP_code_of_protein_1[0]
            # print (ENSP_code_of_protein_1)
            ENSP_code_of_protein = [ENSP_code_of_protein_1] + ENSP_code_of_protein
            # print (ENSP_code_of_protein)
        except :
            Gene_error_list = [Gene_names] + Gene_error_list
            df_1 = pd.read_excel("./Data/Disease/error_gene_list.xlsx")
            _error_list = list(set(df_1["Error_List"].tolist()))
            _error_list = list(set(_error_list + Gene_error_list))
            df= pd.DataFrame({
               "Error_List":  _error_list
            })
            df.to_excel("./Data/Disease/error_gene_list.xlsx")




    print("Encoding Error 명단:",Gene_error_list)
    """ENSP code가 없는 gene들이 있음."""
    # print(ENSP_code_of_protein)
    return ENSP_code_of_protein
"""COL1A1 -> ENSP.0010204으로 변경하는 함수"""
def Barabasi_Herb_list(_Herb_name):
    try:
        a = pd.read_excel("./Data/Herb/Barabasi/TCM_Herb_Barabasi.xlsx")
        herb_target_list = a[a["Korean_name"] == _Herb_name]["Gene Symbol"].tolist()
        # print(_Herb_name + "의 타겟 gene은", herb_target_list)
    except:
        print("맞는 이름의 herb가 없는 것 같습니다.")
    return changing_Gene_name_to_ENSP_ID(herb_target_list)
"""9606.ENSP 출력 list 확보"""
def KIOM_Herb_list(_Herb_name):
    a = TM_MC.Herb_Target_Protein(_Herb_name)
    return a
def Disease_protein_list(_c_code_of_disease):
    d = disease_DB.changing_Gene_name_to_ENSP_ID(_c_code_of_disease)
    return d
# print(Disease_protein_list("C0409959"))
"""Disgenet DB 활용, ENSP ID 출력"""


"""
               Network construction part
"""
def Herb_NETWORK(_herb_name):
    Herb_netwrok = Whole_network_construction().subgraph(Barabasi_Herb_list(_herb_name))
    return Herb_netwrok
def Disease_NETWORK (dis_C_code):
    Dis_network = Whole_network_construction().subgraph(Disease_protein_list(dis_C_code))

    return Dis_network
"""Disgenet DB 기반"""
def get_LCC_Network_list(G_list):

    network_LCC = Whole_network_construction().subgraph(G_list)
    largest_cc = max(nx.connected_components(network_LCC), key=len)
    return  list(set(largest_cc))
"""LCC sub 네트워크를 구성하는 set값을 리턴"""
def Whole_network_construction():
    df = pd.read_csv ('./Data/Whole STRING Network/Whole_network.csv')
    G = nx.Graph ()
    G.add_nodes_from (df['protein1'])
    G.add_nodes_from (df['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in df.iterrows ()]
    G.add_edges_from (edges)

    return G

"""
                               Analysis part 
 - Jaccard Index
 - Proximity Score
 - Separation Score
"""

"""JACCARD index, JACCARD max 값 리턴"""
def overlap_significance(geneids1, geneids2, method = "jaccard"):
    """
    method:  jaccard | jaccard_max
    """
    n1, n2 = set(geneids1), set(geneids2)


    if method == "jaccard":
        val = 1.0 * len(n1 & n2) / len(n1 | n2)
    elif method == "jaccard_max":
        val = 1.0 * len(n1 & n2) / max(map(len, [n1, n2]))

    return val
"""Distance, proximity calculation 함수 모음"""
def get_shortest_path_length_between(G, source_id, target_id):
    try:
        d = nx.shortest_path_length (G, source_id, target_id)
    except:
        d = 0
    return d
"""네트워크에서 두 노드간의 거리를 측정하는 함수"""
def calculate_closest_distance(network, nodes_from, nodes_to):
    ds = defaultdict(dict)

    for i in nodes_from:
        for j in nodes_to:
            if i == j:
                ds[i][j] = 0
            else:
                if nx.has_path (network, i, j):
                    ds[i][j] = nx.shortest_path_length (network, i, j)
                else:
                    ds[i][j] = float ('nan')

    ds = pd.DataFrame.from_dict (ds)
    # nodes_to: rows
    # nodes_from: columns

    dic = {}

    dic['shortest'] = ds.mean().mean()
    dic['closest'] = ds.min().mean()
    return (dic)

def calculate_network_distance(network, nodes_from, nodes_to):
    values_outer = []
    try:
        values = []
        for node_from in nodes_from:
            for node_to in nodes_to:
                val = get_shortest_path_length_between (network, node_from, node_to)
                # print(val)
                values.append (val)
        # print("proximity는",d,"입니다.")
    except:
        pass

    d = sum (values) / len (values)
    return d
"""list 2개의 gene들 사이에서 거리를 측정하여 평균 값을 도출함"""
def get_separation_between_sets(network, nodes_from, nodes_to, lengths=None):
    """
    Calculate dAB in separation metric proposed by Menche et al. 2015
    """
    values = []
    target_to_values = {}
    source_to_values = {}
    for source_id in nodes_from:
        for target_id in nodes_to:
            if lengths is not None:
                d = lengths[source_id][target_id]
            else:
                d = get_shortest_path_length_between(network, source_id, target_id)
            source_to_values.setdefault(source_id, []).append(d)
            target_to_values.setdefault(target_id, []).append(d)
    # Distances to closest node in nodes_to (B) from nodes_from (A)
    for source_id in nodes_from:
        inner_values = source_to_values[source_id]
        values.append(numpy.min(inner_values))
    # Distances to closest node in nodes_from (A) from nodes_to (B)
    for target_id in nodes_to:
        inner_values = target_to_values[target_id]
        values.append(numpy.min(inner_values))
    return values
def get_separation_within_set(network, nodes_from, lengths=None):
    """
    Calculate dAA or dBB in separation metric proposed by Menche et al. 2015
    """
    if len(nodes_from) == 1:
        return [ 0 ]
    values = []
    # Distance to closest node within the set (A or B)
    for source_id in nodes_from:
        inner_values = []
        for target_id in nodes_from:
            if source_id == target_id:
                continue
            if lengths is not None:
                d = lengths[source_id][target_id]
            else:
                d = get_shortest_path_length_between(network, source_id, target_id)
            inner_values.append(d)
        values.append(numpy.min(inner_values))
    return values

def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.iteritems(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values = sorted(values)
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        #print i, low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
             bins.append((low, high, val))
    return bins
def get_degree_equivalents(seeds: object, bins: object, g: object) -> object:
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes
def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                #nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20): # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [ random.choice(nodes) ]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values
def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size)
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed)
    return nodes_random



"""최종 응용 함수"""
def get_separation(network, nodes_from, nodes_to, lengths=None):
    dAA = numpy.mean(get_separation_within_set(network, nodes_from, lengths))
    dBB = numpy.mean(get_separation_within_set(network, nodes_to, lengths))
    dAB = numpy.mean(get_separation_between_sets(network, nodes_from, nodes_to, lengths))
    d = dAB - (dAA + dBB) / 2.0
    return d

"""네트워크, 허브이름, 질병 C_code 넣으면 평균 거리를 측정"""
def LCC_network_proximity_calculate(herb, c_code_disease):

    try:
        z = calculate_network_distance (Whole_network_construction(), get_LCC_Network_list (Herb_NETWORK (herb)),get_LCC_Network_list (Disease_NETWORK (c_code_disease)))
    except:
        z = 0
        print("network calculation ERROR")

    return z
def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set (network.nodes ())
    if len (set (nodes_from) & nodes_network) == 0 or len (set (nodes_to) & nodes_network) == 0:
        return None  # At least one of the node group not in network

    d = calculate_closest_distance (network, nodes_from, nodes_to)

    if n_random:

        if bins is None and (nodes_from_random is None or nodes_to_random is None):
            bins = get_degree_binning (network, min_bin_size,
                                       lengths)  # if lengths is given, it will only use those nodes
        if nodes_from_random is None:
            nodes_from_random = get_random_nodes (nodes_from, network, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        if nodes_to_random is None:
            nodes_to_random = get_random_nodes (nodes_to, network, bins=bins, n_random=n_random,
                                                min_bin_size=min_bin_size, seed=seed)
        random_values_list = list (zip (nodes_from_random, nodes_to_random))
        # values = np.empty(len(nodes_from_random)) #n_random
        null = []
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            res = calculate_closest_distance(network, nodes_from, nodes_to)
            null.append(res)

        null_s = []
        null_c = []
        for i in range (len (null)):
            null_s.append (null[i]['shortest'])
            null_c.append (null[i]['closest'])

        with numpy.errstate (divide='ignore', invalid='ignore'):

            d['avg_shortest'], d['std_shortest'] = numpy.mean (null_s), numpy.std (null_s)
            d['z_shortest'] = (d['shortest'] - d['avg_shortest']) / d['std_shortest']

            d['avg_closest'], d['std_closest'] = numpy.mean (null_c), numpy.std (null_c)
            d['z_closest'] = (d['closest'] - d['avg_closest']) / d['std_closest']


    return (d)

"""엑셀 파일에서 질환-LCC target List를 확보"""
def Dis_skin_LCC_List_from_xlsx(_file_path):
    Dis_list = pd.read_excel(_file_path)["all shared genes"].tolist()
    # print(Dis_list)
    Dis_list_ENSP = changing_Gene_name_to_ENSP_ID(Dis_list)
    dis_LCC = get_LCC_Network_list(Dis_list_ENSP)
    return dis_LCC
