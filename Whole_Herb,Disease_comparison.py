import pandas as pd
import HERB_DISEASE_netowrk_workspace
def Whole_herb_disease_comparison():
    whole_herb_name_df = pd.read_excel("./Data/Herb/KIOM_DB/medicinal_material.xlsx")
    whole_herb_name_list = whole_herb_name_df["KOREAN"].tolist()
    # whole_herb_name_list = whole_herb_name_list[:2]
    """Herb_list"""
    whole_disease_C_code_df = pd.read_excel ("./Data/Disease/disease_index(Id_to_NID).xls")
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
            c = HERB_DISEASE_netowrk_workspace.LCC_network_proximity_calculate(a,b)
            total_prox_list = [a,b,c]
            print(total_prox_list)
            sample_list.append(total_prox_list)
            total_prox_table=pd.DataFrame(sample_list, columns= ["Herb_name", "Disease_code","Proximity_score"])
    # print(total_prox_table.head(2))
    total_prox_table.to_csv("_HERB_Sarcopenia_Proximity_total.csv", index= False, encoding= 'cp949')
    return
"""리소스 투입이 많이 필요함"""
Whole_herb_disease_comparison()
"""24.06.05. 13:54 시작"""


def Herb_N_Disease(_H_List,_D_List):
    sample_list = []
    for a in _H_List:
        for b in _D_List:
            c = HERB_DISEASE_netowrk_workspace.LCC_network_proximity_calculate(a,b)
            total_prox_list = [a, b, c]
            # print (total_prox_list)
            sample_list.append (total_prox_list)
            total_prox_table = pd.DataFrame (sample_list, columns=["Herb_name", "Disease_code", "Proximity_score"])
            # print(total_prox_table.head(2))
    total_prox_table.to_csv ("HD_Proximity.csv", index=False)
    return
"""본초 명단 리스트랑 질병 이름 리스트를 받아서 거리 계산해서 CSV로 저장함"""

H_List  = ["구척", "강활"]
"""본초 이름"""
D_List = ["C0409959"]
"""질병 이름"""

# Herb_N_Disease(H_List,D_List)