import Barabasi_Proximity as BP
import pandas as pd


# OK_herb_list = BP.Barabasi_Herb_list("강활")

#세부 pathway gene proximity 계산
"""
data = pd.read_excel("./Result/OA_강활/KEGG pathway gene list.xlsx",header=0)

gene_lists = {}

#dict 형식으로 저장
for col in data.columns:
    gene_lists[col] = data[col].dropna().tolist()
    # c = BP.changing_Gene_name_to_ENSP_ID(gene_lists[col])
    # print(col,":", c)

results = []
for col_1 in gene_lists.keys():
    b = BP.calculate_proximity(BP.Whole_network_construction(),OK_herb_list,BP.changing_Gene_name_to_ENSP_ID(gene_lists[col_1]))
    print(b)
    results.append({"Pathway": col_1, "Proximity": b})

results_df = pd.DataFrame(results)
results_df.to_excel("./Result/OA_강활/Pathway_Proximity_Results.xlsx", index=False, sheet_name="Proximity Results")
"""

##
all_cytokine_list = pd.read_excel("./Result/OA_강활/Cytokine_Genes_Sorted.xlsx")["All Cytokines"].tolist()
# print(all_cytokine_list)

def proximity_calculation(_herb,_list_of_genes):
    b = BP.calculate_proximity(BP.Whole_network_construction(),BP.Barabasi_Herb_list(_herb),BP.changing_Gene_name_to_ENSP_ID(_list_of_genes))
    return b

results=[]
for a in all_cytokine_list:
    z = proximity_calculation("강활",[a])
    print(a,z)
    results.append({"Cytokine": a, "Proximity": z})
results_DF = pd.DataFrame(results)
results_DF.to_excel("강활과 cytokines proximity.xlsx")
