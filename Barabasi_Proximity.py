import pandas as pd
import networkx as nx
import requests
import network_construction_from_string
import network_construction_from_string
import disease_DB
import preparation
import time
import numpy
from itertools import permutations, combinations


"""
Barabasi lab의 reference대로 환경 구축해서 작업하는 파일 
"""
def Herb_list(_Herb_name):
    a = pd.read_excel("./Data/Herb/Barabasi/TCM_Herb_Barabasi.xlsx")
    herb_target_list = a[a["Korean_name"] == _Herb_name]["Gene Symbol"].tolist()
    # print(_Herb_name + "의 타겟 gene은", herb_target_list)
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
def get_LCC_Network_list(G):
    largest_cc = max(nx.connected_components (G), key=len)
    return  list(largest_cc)
"""LCC sub 네트워크를 구성하는 set값을 리턴"""
def Whole_network_construction():
    df = pd.read_csv ('./Data/Whole STRING Network/Whole_network.csv')
    G = nx.Graph ()
    G.add_nodes_from (df['protein1'])
    G.add_nodes_from (df['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in df.iterrows ()]
    G.add_edges_from (edges)

    return G

"""Herb의 list를 불러오기, KIOM의 데이터를 확보해서 사용"""


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

"""최종 응용 함수"""
def get_separation(network, nodes_from, nodes_to, lengths=None):
    dAA = numpy.mean(get_separation_within_set(network, nodes_from, lengths))
    dBB = numpy.mean(get_separation_within_set(network, nodes_to, lengths))
    dAB = numpy.mean(get_separation_between_sets(network, nodes_from, nodes_to, lengths))
    d = dAB - (dAA + dBB) / 2.0
    return d
def LCC_network_proximity_calculate(herb, c_code_disease):

    try:
        z = calculate_network_distance (Whole_network_construction(), get_LCC_Network_list (Herb_NETWORK (herb)),get_LCC_Network_list (Disease_NETWORK (c_code_disease)))
    except:
        z = 0
        print("network calculation ERROR")

    return z
"""네트워크, 허브이름, 질병 C_code 넣으면 평균 거리를 측정"""


######### 강활_염증 Pathway 세부 선정 ##########
#for i in ["TNFRSF17", "CSF3", "CCR7", "IL1B", "CCR2", "CXCL10", "CXCL8", "LEP", "LEPR", "CCL20", "XCL2", "NGF", "TNFRSF1B", "CXCL1", "TNF", "IL10", "IL17A"
#           ]:
#     a = calculate_network_distance(Whole_network_construction(),changing_Gene_name_to_ENSP_ID(Herb_list("강활")),changing_Gene_name_to_ENSP_ID([i]))
#     print(i,a)

#####시간 측정하는 함수 칸######
start_time = time.time()
#
end = time.time()
print(f"{end - start_time:.5f} sec")

"""시간측정"""


def Dis_skin_LCC_List_from_xlsx(_file_path):
    Dis_list = pd.read_excel(_file_path)["all shared genes"].tolist()
    dis_LCC = get_LCC_Network_list(Dis_list)
    return dis_LCC


####### 탄력약침 - 약물 조합관계 분석 #######

file_50_path = "./Project/Skin_Aging/skin_aging(50).xlsx"
file_150_path = "./Project/Skin_Aging/skin_aging(150).xlsx"

"""Skin aging 관련 data 생성, 약재는 총 10가지 (담두시(두시), 천화분(괄루근), 구기자, 인삼, 육종용, 백지, 금은화, 천마, 건율, 부자, 국화)"""
Herb_candidates_list = ["인삼","두시", "괄루근", "구기자","육종용", "백지", "금은화", "천마", "건율", "부자", "국화"]
combination_list = list(combinations(Herb_candidates_list,2))

for i in combination_list:
    for d in [file_50_path, file_150_path]:
        get_separation(Whole_network_construction(),Herb_list("강활"),Dis_skin_LCC_List_from_xlsx(d))