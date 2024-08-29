import  pandas as pd
import openpyxl

# dir = "chemical_property.xlsx"
# df = pd.read_excel(dir)
# compound_list = df["INCHIKEY"].tolist()
#
# Cid_translated_df = pd.read_excel("chemicals.inchikeys.v5.0_filtered.xlsx")
# # print(Cid_translated_df.head(10))
# merge_1 = pd.merge(df,Cid_translated_df,on='INCHIKEY',how='outer')
# merge_1.to_excel("chemical_property_revised(CID).xlsx")
# # print(merge_1.head(100))

###Protein-> Chemical file 기준값 설정 만지기 ###

dir_3 = "./Data/Herb/SY_DB/9606.protein_chemical.links.detailed.v5.0.tsv"
df_3 = pd.read_csv(dir_3, delimiter= '\t')
# condition_1 = (df_3['protein']>0) or (df_3['experimental']>0) or (df_3["prediction"]>0) or (df_3["database"]>0)
condition_2 = df_3["combined_score"]>=400
df_3_revised = df_3[condition_2]
df_3_revised.to_excel("./Data/Herb/SY_DB/9606.protein_chemical.links.detailed.v5.0_(score 400).xlsx")


