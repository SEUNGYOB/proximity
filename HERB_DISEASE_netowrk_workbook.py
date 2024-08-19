import pandas as pd
import networkx as nx
import disease_DB
import time
import numpy

"""전체 STRING 네트워크 만들기"""
def Herb_NETWORK(_herb_name):
    Herb_netwrok = Whole_network_construction().subgraph(Herb_list(_herb_name))
    return Herb_netwrok
def Disease_NETWORK (dis_C_code):
    Dis_network = Whole_network_construction().subgraph(Disease_protein_list(dis_C_code))

    return Dis_network
def Whole_network_construction():
    df = pd.read_csv ('./Data/Whole STRING Network/Whole_network.csv')
    G = nx.Graph ()
    G.add_nodes_from (df['protein1'])
    G.add_nodes_from (df['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in df.iterrows ()]
    G.add_edges_from (edges)

    return G

"""Herb의 list를 불러오기, KIOM의 데이터를 확보해서 사용"""
def Herb_list(_herb_name):
    herb_df = pd.read_excel ("./Data/Herb/Herb_Protein_DB_v2/" + _herb_name + "의 protein_list.xlsx")
    herb_protein_list = herb_df["Protein name"].tolist()
    """9606.ENSP0000000으로 변환"""
    herb_gene_list_9606 = herb_protein_list
    for i, ss in enumerate (herb_protein_list):
        herb_gene_list_9606[i] = "9606." + herb_protein_list[i]
    # print (herb_protein_list)
    return herb_gene_list_9606
# Herb_list("가자")

def Disease_protein_list(_c_code_of_disease):
    d = disease_DB.changing_Gene_name_to_ENSP_ID(_c_code_of_disease)
    return d
# print(Disease_protein_list("C0409959"))
"""Disgenet DB 활용"""

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
    print("Encoding Error 명단:",Gene_error_list)
    """ENSP code가 없는 gene들이 있음."""
    # print(ENSP_code_of_protein)
    return ENSP_code_of_protein
"""COL1A11 -> ENSP.0010204으로 변경하는 함수"""

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
    return nx.shortest_path_length (G, source_id, target_id)
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
def get_separation(network, nodes_from, nodes_to, lengths=None):
    dAA = numpy.mean(get_separation_within_set(network, nodes_from, lengths))
    dBB = numpy.mean(get_separation_within_set(network, nodes_to, lengths))
    dAB = numpy.mean(get_separation_between_sets(network, nodes_from, nodes_to, lengths))
    d = dAB - (dAA + dBB) / 2.0
    return d
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

def get_LCC_Network_set(G):
    largest_cc = max(nx.connected_components (G), key=len)
    return  list(largest_cc)
"""LCC sub 네트워크를 구성하는 set값을 리턴"""


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
"""최종 응용 함수"""

def LCC_network_proximity_calculate(herb, c_code_disease):

    try:
        z = calculate_network_distance (Whole_network_construction(), get_LCC_Network_set (Herb_NETWORK (herb)),get_LCC_Network_set (Disease_NETWORK (c_code_disease)))
    except:
        z = 0
        print("network calculation ERROR")

    return z
"""네트워크, 허브이름, 질병 C_code 넣으면 평균 거리를 측정"""



#####시간 측정하는 함수 칸######
start_time = time.time()
# print(Herb_list("강활"))
# print(Disease_protein_list("C0872084"))
# print(calculate_network_distance(Whole_network_construction(),get_LCC_Network_set (Herb_NETWORK ("강활")),get_LCC_Network_set (Disease_NETWORK ("C0409959"))))
# print(get_LCC_Network_set(Herb_NETWORK("강활")))
# print(get_LCC_Network_set(Disease_NETWORK("C0872084")))

# print(LCC_network_proximity_calculate("강활", "C0409959"))

"""강활과 골관절염 proximity 비교"""
print("--- %s seconds ---" % (time.time() - start_time))
"""시간측정"""


