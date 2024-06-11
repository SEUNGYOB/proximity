import pandas as pd
import networkx as nx
import requests
import network_construction_from_string
import network_construction_from_string
import disease_DB
import preparation
import time

def Herb_Network_using_STRING(_herb_name):
    df_herb =pd.read_excel("./Data/Herb/Herb_Protein_DB_v2/"+_herb_name+"의 protein_list.xlsx")
    herb_gene_list = df_herb["Protein name"].tolist()
    # print(herb_gene_list)

    herb_gene_list_9606 = herb_gene_list
    for i,ss in enumerate(herb_gene_list):
        herb_gene_list_9606[i] = "9606." + herb_gene_list[i]
    print(herb_gene_list_9606)

    H_netwrok_DATAFRAME= network_construction_from_string.STRING_API_network_construct_INFO_GET(herb_gene_list_9606)
    H_netwrok_DATAFRAME.columns = ['protein1','protein2']
    """그래프 만들기"""
    g = nx.Graph ()
    g.add_nodes_from (H_netwrok_DATAFRAME['protein1'])
    g.add_nodes_from (H_netwrok_DATAFRAME['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in H_netwrok_DATAFRAME.iterrows ()]
    g.add_edges_from (edges)

    return g
"""STRING API는 많은 유전자를 입력할 경우 해결하지 못 함. 쓸 수 없는 함수"""
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

"""Disgenet DB 활용"""
def Disease_protein_list(_c_code_of_disease):
    d = disease_DB.changing_Gene_name_to_ENSP_ID(_c_code_of_disease)
    return d
# print(Disease_protein_list("C0409959"))

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



"""전체 STRING 네트워크 만들기"""

def Herb_NETWORK(_herb_name):
    Herb_netwrok = Whole_network_construction().subgraph(Herb_list(_herb_name))
    return Herb_netwrok

def Disease_NETWORK (dis_C_code):
    Dis_network = Whole_network_construction().subgraph(Disease_protein_list(dis_C_code))

    return Dis_network


def get_shortest_path_length_between(G, source_id, target_id):
    try:

        spl = nx.shortest_path_length (G, source_id, target_id)
    except:
        spl = 0

    return spl
"""네트워크에서 두 노드간의 거리를 측정하는 함수"""

# print(get_shortest_path_length_between(Whole_network_construction(),"9606.ENSP00000000233","9606.ENSP00000371175"))

def get_LCC_Network_set(G):
    largest_cc = max (nx.connected_components (G), key=len)
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
print(Disease_protein_list("C0872084"))
# print(calculate_network_distance(Whole_network_construction(),get_LCC_Network_set (Herb_NETWORK ("강활")),get_LCC_Network_set (Disease_NETWORK ("C0409959"))))
# print(get_LCC_Network_set(Herb_NETWORK("강활")))
print(get_LCC_Network_set(Disease_NETWORK("C0872084")))

# print(LCC_network_proximity_calculate("강활", "C0409959"))

"""강활과 골관절염 proximity 비교"""
print("--- %s seconds ---" % (time.time() - start_time))
"""시간측정"""
