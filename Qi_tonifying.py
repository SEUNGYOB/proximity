import pandas as pd
import HERB_DISEASE_netowrk_workbookpy




Qi_herb_list = ["인삼","당삼","황기", "백출", "산약", "백편두", "감초", "대조", "봉밀"]

"""수정 필요한 함수"""
def Group_herb_disease_comparison(_Herb_G,_Disease_G = None):

    whole_disease_C_code_df = pd.read_excel ("./Data/Disease/disease_index(Id_to_NID).xls")
    whole_disease_C_code = whole_disease_C_code_df["diseaseId"].tolist ()
    # whole_disease_C_code = whole_disease_C_code[:2]
    """Disease_code"""
    """점검"""
    # print (whole_disease_C_code)
    # print(whole_herb_name_list)


    total_prox_table = pd.DataFrame()
    sample_list=[]
    for a in _Herb_G:
        for b in ["C0872084"]:
            c = HERB_DISEASE_netowrk_workbookpy.LCC_network_proximity_calculate(a, b)
            total_prox_list = [a,b,c]
            print(total_prox_list)
            sample_list.append(total_prox_list)
            total_prox_table=pd.DataFrame(sample_list, columns= ["Herb_name", "Disease_code","Proximity_score"])
    # print(total_prox_table.head(2))
    total_prox_table.to_csv("_HERB_Sarcopenia_Proximity_total.csv", index= False, encoding= 'cp949')
    return
"""리소스 투입이 많이 필요함"""
Group_herb_disease_comparison(Qi_herb_list,???????)