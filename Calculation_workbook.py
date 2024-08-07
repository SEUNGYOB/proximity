import pandas as pd
import HERB_DISEASE_netowrk_workbook as HDN

def Whole_herb_disease_comparison():
    whole_herb_name_df = pd.read_excel("./Data/Herb/KIOM_DB/medicinal_material.xlsx")
    whole_herb_name_list = whole_herb_name_df["KOREAN"].tolist()
    # whole_herb_name_list = whole_herb_name_list[:2]
    """Herb_list"""
    whole_disease_C_code_df = pd.read_excel ("./Data/Disease/disease_index(Id_to_NID).xlsx")
    whole_disease_C_code = whole_disease_C_code_df["diseaseId"].tolist ()
    # whole_disease_C_code = whole_disease_C_code[:2]
    """Disease_code"""
    """점검"""
    # print (whole_disease_C_code)
    # print(whole_herb_name_list)


    total_prox_table = pd.DataFrame()
    sample_list=[]
    for a in whole_herb_name_list:
        for b in ["C0872084"]:
            c = HDN.LCC_network_proximity_calculate(a, b)
            total_prox_list = [a,b,c]
            print(total_prox_list)
            sample_list.append(total_prox_list)
            total_prox_table=pd.DataFrame(sample_list, columns= ["Herb_name", "Disease_code","Proximity_score"])
    # print(total_prox_table.head(2))
    total_prox_table.to_csv("_HERB_Sarcopenia_Proximity_total.csv", index= False, encoding= 'cp949')
    return
"""리소스 투입이 많이 필요함"""
# Whole_herb_disease_comparison()
"""24.06.05. 13:54 시작"""


def Herb_N_Disease(_H_List,_D_List):
    sample_list = []
    for a in _H_List:
        for b in _D_List:
            c = HDN.LCC_network_proximity_calculate(a, b)
            total_prox_list = [a, b, c]
            # print (total_prox_list)
            sample_list.append (total_prox_list)
            total_prox_table = pd.DataFrame (sample_list, columns=["Herb_name", "Disease_code", "Proximity_score"])
            # print(total_prox_table.head(2))
    total_prox_table.to_csv ("HD_Proximity.csv", index=False)
    return
"""본초 리스트랑 질병 리스트를 받아서 거리 계산해서 CSV로 저장함"""

# H_List  = ["구척", "강활"]
# """본초 이름"""
# D_List = ["C0409959"]
# """질병 이름"""

# Herb_N_Disease(H_List,D_List)

def Herb_complex_Disease_distance(_H_List, _D_name):
    sample_list = []
    for a in _H_List:
        x= HDN.Herb_list(a)
        # print(len(x))
        sample_list = list(set(sample_list + x))
        # print(len(sample_list))
    H_Complex_LCC_Nodes = HDN.get_LCC_Network_set(HDN.Whole_network_construction().subgraph(sample_list))
    D_LCC_nodes = HDN.get_LCC_Network_set (HDN.Disease_NETWORK (_D_name))
    z = HDN.calculate_network_distance(HDN.Whole_network_construction(),H_Complex_LCC_Nodes, D_LCC_nodes)
    print (z)
    return z

# H_Comp = ["홍화"]
# D_name = "C0409959"
# Herb_complex_Disease_distance(H_Comp,D_name)

def Separation_Herb_Herb_Disease(_H_Name_1, _H_Name_2, _D_list) :
    Herb_comp_list_1 = HDN.Herb_list(_H_Name_1)
    Herb_comp_list_2 = HDN.Herb_list(_H_Name_2)

    H_Complex_LCC_Nodes_1 = HDN.get_LCC_Network_set(HDN.Whole_network_construction().subgraph(Herb_comp_list_1))
    H_Complex_LCC_Nodes_2 = HDN.get_LCC_Network_set(HDN.Whole_network_construction().subgraph(Herb_comp_list_2))
    Disease_network = HDN.Whole_network_construction().subgraph (_D_list)
    D_LCC_nodes = HDN.get_LCC_Network_set (Disease_network)

    HERB_1_and_2_separation_score = HDN.get_separation(HDN.Whole_network_construction(),H_Complex_LCC_Nodes_1, H_Complex_LCC_Nodes_2)
    HERB_1_and_Disease_separation_score = HDN.get_separation(HDN.Whole_network_construction(), H_Complex_LCC_Nodes_1, D_LCC_nodes)
    HERB_2_and_Disease_separation_score = HDN.get_separation(HDN.Whole_network_construction(), D_LCC_nodes, H_Complex_LCC_Nodes_2)

    print("Herb_1과 Herb_2의 거리는 ", HERB_1_and_2_separation_score, "Herb_1과 Disease의 거리는 ", HERB_1_and_Disease_separation_score, "Herb_1과 Disease의 거리는 ", HERB_2_and_Disease_separation_score )
    return HERB_1_and_2_separation_score, HERB_1_and_Disease_separation_score, HERB_2_and_Disease_separation_score
"""한약재 2종, Disease list를 받아서 separation Score을 계산"""

D_list = pd.read_excel("./Project/Skin_Aging/skin_aging(50).xlsx")
D_list = D_list["all shared genes"].tolist()
D_list = HDN.changing_Gene_name_to_ENSP_ID(D_list)
# print(D_list)
print(Separation_Herb_Herb_Disease("인삼", "육종용", D_list))