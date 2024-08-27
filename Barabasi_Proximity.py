import pandas as pd
import networkx as nx
import random, copy

import disease_DB
import time
import numpy
from itertools import combinations


"""
Barabasi lab의 reference대로 환경 구축해서 작업하는 파일 
"""
def Herb_list(_Herb_name):
    try:
        a = pd.read_excel("./Data/Herb/Barabasi/TCM_Herb_Barabasi.xlsx")
        herb_target_list = a[a["Korean_name"] == _Herb_name]["Gene Symbol"].tolist()
        # print(_Herb_name + "의 타겟 gene은", herb_target_list)
    except:
        print("맞는 이름의 herb가 없는 것 같습니다.")
    return herb_target_list
"""Target Gene의 list 확보"""
def changing_Gene_name_to_ENSP_ID(Gene_list):
    chemical_protein_index = pd.read_csv("./Data/Disease/protein_dict.csv")
    ENSP_code_of_protein = []
    Gene_error_list = []
    for Gene_names in Gene_list:
        ENSP_code_of_protein_1 = chemical_protein_index.loc[chemical_protein_index.preferred_name == Gene_names]["#string_protein_id"].tolist()
        # print(ENSP_code_of_protein_1)
        try :
            ENSP_code_of_protein_1 = ENSP_code_of_protein_1[0]
            # print (ENSP_code_of_protein_1)
            ENSP_code_of_protein = [ENSP_code_of_protein_1] + ENSP_code_of_protein
            # print (ENSP_code_of_protein)
        except :
            Gene_error_list = [Gene_names] + Gene_error_list
            df_1 = pd.read_excel("./Data/Disease/error_gene_list.xlsx")
            _error_list = df_1["Error_List"].tolist()
            _error_list = _error_list + Gene_error_list
            df= pd.DataFrame({
               "Error_List":  _error_list
            })
            df.to_excel("./Data/Disease/error_gene_list.xlsx")




    print("Encoding Error 명단:",Gene_error_list)
    """ENSP code가 없는 gene들이 있음."""
    # print(ENSP_code_of_protein)
    return ENSP_code_of_protein
"""COL1A1 -> ENSP.0010204으로 변경하는 함수"""
def Disease_protein_list(_c_code_of_disease):
    d = disease_DB.changing_Gene_name_to_ENSP_ID(_c_code_of_disease)
    return d
# print(Disease_protein_list("C0409959"))
"""Disgenet DB 활용"""


"""
                                        Network construction part
"""
def Herb_NETWORK(_herb_name):
    Herb_netwrok = Whole_network_construction().subgraph(Herb_list(_herb_name))
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
    n1, n2 = len(geneids1), len(geneids2)
    n = len(geneids1 & geneids2)


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
"""list 2개의 gene들 사이에서 거리를 측정하여 최소 값을 도출함"""
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
def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None, distance="closest"):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network
    if distance != "closest":
        lengths = get_shortest_path_length_between(network, "temp_n%d_e%d.sif.pcl" % (len(nodes_network), network.number_of_edges()))
        d = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    else:
        d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        random_values_list = zip(nodes_from_random, nodes_to_random)
        values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        if distance != "closest":
            values[i] = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        else:
            values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s) #(z, pval)

"""엑셀 파일에서 질환-LCC target List를 확보"""
def Dis_skin_LCC_List_from_xlsx(_file_path):
    Dis_list = pd.read_excel(_file_path)["all shared genes"].tolist()
    # print(Dis_list)
    Dis_list_ENSP = changing_Gene_name_to_ENSP_ID(Dis_list)
    dis_LCC = get_LCC_Network_list(Dis_list_ENSP)
    return dis_LCC




######### 강활_염증 Pathway 세부 선정 ##########
"""
for i in ["TNFRSF17", "CSF3", "CCR7", "IL1B", "CCR2", "CXCL10", "CXCL8", "LEP", "LEPR", "CCL20", "XCL2", "NGF", "TNFRSF1B", "CXCL1", "TNF", "IL10", "IL17A"
          ]:
    a = calculate_proximity(Whole_network_construction(),changing_Gene_name_to_ENSP_ID(Herb_list("강활")),changing_Gene_name_to_ENSP_ID([i]))
    print(i,a)
"""
""" skin aging_240821"""
"""
#############################          탄력약침 - 약물 조합관계 분석          ###########################

file_50_path = "./Project/Skin_Aging/skin_aging(50).xlsx"
file_150_path = "./Project/Skin_Aging/skin_aging(150).xlsx"

# # get_separation(Whole_network_construction(),get_LCC_Network_list(changing_Gene_name_to_ENSP_ID(Herb_list ("인삼"))),get_LCC_Network_list(changing_Gene_name_to_ENSP_ID(Herb_list ("두시"))))
# print(get_separation(Whole_network_construction(), t_list_1,t_list_2))
# print(separation_score_1)


Skin aging 관련 data 생성, 약재는 총 10가지 (담두시(두시, 제외), 천화분(과루), 구기자(구기), 인삼, 육종용, 백지, 금은화(인동등), 천마, 건율, 부자, 국화)
Herb_candidates_list = ["인삼", "과루", "구기","육종용", "백지", "인동등", "천마", "생부자", "국화"]
combination_list = list(combinations(Herb_candidates_list,2))
a=[]
b=[]
c=[]
for i in combination_list:
    # print(type(i[0]),i[1]) #i[0] : str type
    # print(Herb_list(i[0]))
    for d in [file_50_path, file_150_path]:

        H_list_1 = get_LCC_Network_list(changing_Gene_name_to_ENSP_ID(Herb_list (i[0])))
        H_list_2 = get_LCC_Network_list(changing_Gene_name_to_ENSP_ID(Herb_list (i[1])))
        D_List  = Dis_skin_LCC_List_from_xlsx(d)

        # print(H_list_1, H_list_2)
        separation_score_1 = get_separation(Whole_network_construction(),H_list_1,H_list_2)
        separation_score_2 = get_separation(Whole_network_construction(),H_list_1,D_List)
        separation_score_3 = get_separation(Whole_network_construction(),H_list_2,D_List)

        a.append(i[0])
        a.append(i[0])
        a.append(i[1])
        b.append(i[1])
        b.append(d)
        b.append(d)
        c.append(separation_score_1)
        c.append(separation_score_2)
        c.append(separation_score_3)

Dis_df = pd.DataFrame ({'A': a, 'B': b, 'Separation Score' : c})
Dis_df.to_excel("탄력약침_skin_aging.xlsx")
"""

############## 강활_염증 Pathway 세부 선정 (염증인자들과의 거리 비교)##########
z = calculate_proximity(Whole_network_construction(),["9606.ENSP00000000233","9606.ENSP00000363232"], ["9606.ENSP00000000233","9606.ENSP00000300935"])

print(z)