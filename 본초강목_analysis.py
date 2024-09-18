import pandas as pd
import openpyxl
import os
import Barabasi_Proximity
import time
def Network_total_analysis(H1,H2):
    whole_network = Barabasi_Proximity.Whole_network_construction()
    H1_target = Barabasi_Proximity.Barabsi_Herb_list(H1)
    H2_target = Barabasi_Proximity.Barabsi_Herb_list(H2)

    """Jaccard index"""
    Jaccard_index = Barabasi_Proximity.overlap_significance(H1_target,H2_target)
    """Overlapping gene"""
    Overlapping_gene_number = len(set(H1_target) & set(H2_target))
    """Proximity Score, Z-score
    d,z,(m,s) 값 도출
    """
    Prox_score = Barabasi_Proximity.calculate_proximity(whole_network,H1_target,H2_target)
    """Separation score"""
    Sep_score = Barabasi_Proximity.get_separation(whole_network,H1_target,H2_target)
    return Jaccard_index,Overlapping_gene_number,Prox_score[0],Prox_score[1],Sep_score

file_dir = "./Project/본초강목 분석/본초강목/"
bcgm_interaction_data = pd.read_excel(file_dir +"본초 강목 정리.xlsx", sheet_name="분석 가능 본초 강목")
herb_pair=bcgm_interaction_data.to_dict("index")
result = {}
for i in range(370):
    herb_pair_dict = {i:[herb_pair[i]["A_korean_name"],herb_pair[i]["관계"],herb_pair[i]["B_korean_name"]]}
    # print(herb_pair_dict)
    result = {**result, **herb_pair_dict}

start_time = time.time()

total_df = pd.DataFrame ()
for key, value in result.items():
    i = Network_total_analysis(value[0],value[2])
    
    bcgm_analysis_df = pd.DataFrame({
        "Herb_1": [value[0]],
        "Herb_2" : [value[2]],
        "Relationship" : [value[1]],
        "Jaccard Index" : [i[0]],
        "number of overlapping gene" : [i[1]],
        "Proximity" : [i[2]],
        "Z-score" : [i[3]],
        "Separation Score": [i[4]]
    })

bcgm_analysis_df.to_excel(file_dir + "본초강목 분석 결과.xlsx")
print("--- %s seconds ---" % (time.time() - start_time))