from multiprocessing import (Pool)
import HERB_DISEASE_netowrk_workbookpy
from itertools import permutations
import HERB_DISEASE_netowrk_workbookpy
import pandas as pd

whole_herb_name_df = pd.read_excel("./Data/Herb/KIOM_DB/medicinal_material.xlsx")
whole_herb_name_list = whole_herb_name_df["KOREAN"].tolist()

whole_disease_C_code_df = pd.read_excel ("./Data/Disease/disease_index(Id_to_NID).xls")
whole_disease_C_code = whole_disease_C_code_df["diseaseId"].tolist ()

def multi_herb_disease_proximity_calculation(h, d):


if __name__ == "__main__":
    # data = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    items = [whole_herb_name_list,whole_disease_C_code]
    data = list(permutations(items,2))
    # data = ["C0872084"]
    cpu = 8
    # data = [(),(),()]

    with Pool(processes=cpu) as pool:
        results = pool.map(HERB_DISEASE_netowrk_workbookpy.LCC_network_proximity_calculate, data)

    print("Results:", results)