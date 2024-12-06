import pandas as pd
import openpyxl
import os
barabasi_db_file_dir = "./Project/본초강목 분석/본초강목/TCM_Herb_Barabasi.xlsx"
df_barabasi_db = pd.read_excel(barabasi_db_file_dir)
xyz = df_barabasi_db["Korean_name"].value_counts()
xyz=pd.DataFrame(xyz)
# print(xyz)

xyz_filtered = xyz[xyz['count']>= 100]
xyz_filtered = xyz_filtered.reset_index()

# print(xyz_filtered.head())
"""Protein Target 갯수 >=100 filtering"""

filtered_barabasi_herb_set = set(xyz_filtered["Korean_name"].tolist())
# print(filtered_barabasi_herb_set)



file_name = "bcgm_240930.xlsx"
project_loc = "./Project/본초강목 분석/본초강목/"
file_dir =project_loc + file_name
df = pd.read_excel(file_dir, sheet_name= "Sheet6")
Herb_1_set = set(df["Herb_1"].tolist())
Herb_2_set = set(df["Herb_2"].tolist())
Total_herb_set = Herb_1_set | Herb_2_set
applicable_list = filtered_barabasi_herb_set & Total_herb_set
"""조건에 부합하는 데이터 프레임 만들기"""
df_final = df[df["Herb_1"].isin(applicable_list)]
df_final = df_final[df_final["Herb_2"].isin(applicable_list)]

df_final.to_excel(project_loc+"cutoffed_list.xlsx")