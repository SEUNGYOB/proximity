import sqlite3
import pandas as pd
import network_construction_from_string
import time
def disgenet_SQL_Datafrme_NID(_C_code_of_disease):

    """엑셀로 저장"""
    conn = sqlite3.connect("./Data/Disease/disgenet_2020.db")
    cur = conn.cursor()

    #table 호출
    # cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    # print(cur.fetchall())
    # [('diseaseAttributes',), ('diseaseClass',), ('disease2class',), ('geneAttributes',), ('geneDiseaseNetwork',), ('variantAttributes',), ('variantGene',), ('variantDiseaseNetwork',)]

    result = cur.execute("SELECT * FROM diseaseAttributes;")
    rows = cur.fetchall()
    cols = [column[0] for column in cur.description]

    data_df = pd.DataFrame.from_records(data=rows, columns=cols)
    # print(data_df.head())
    data_df_wanted_list_1=[]
    df_trash = pd.DataFrame()
    for diseaseIds in getting_disease_type_disgenet():

        data_df_wanted_C_to_N = data_df.loc[data_df.diseaseId == diseaseIds]

        df_trash = pd.concat([df_trash, data_df_wanted_C_to_N], ignore_index=True)
    #     data_df_wanted_list = data_df_wanted_C_to_N["diseaseNID"].tolist()
    #     data_df_wanted_list_1.extend(data_df_wanted_list)
    # data_df_wanted_list = list(set(data_df_wanted_list_1))
    df_trash.to_excel("disease_index(Id_to_NID).xls")
    diseaseNID = df_trash.loc[df_trash.diseaseId == _C_code_of_disease]["diseaseNID"]

    conn.close ()
    return
"""역할을 다 한 함수 disease_index(Id_to_NID).xls 생성"""
def geneNID_to_gene_excel():
    conn = sqlite3.connect ("./Data/Disease/disgenet_2020.db")
    cur = conn.cursor ()
    result = cur.execute ("SELECT * FROM geneAttributes;")
    rows = cur.fetchall ()
    cols = [column[0] for column in cur.description]

    data_df = pd.DataFrame.from_records (data=rows, columns=cols)
    data_df.to_excel("gene_attribution.xls")
"""역할을 다 한 함수(gene_attribution.xls) 생성. geneNID -> Gene name"""
def getting_disease_type_disgenet():
    data = pd.read_csv("./Data/Disease/disease_associations.csv")
    condition = (data.diseaseType=="disease") & ( data.NofGenes > 300)
    data_1 = data.loc[condition]["diseaseId"]

    disease_list_disgenet = data_1.tolist()
    return disease_list_disgenet
"""list type으로 disgenet에 존재하는 "disease"로 분류된 질병 코드를 반환(CXXXXXXX 형식), NofGenes 가 300개 초과인 것들 """

getting_disease_type_disgenet()
def SQL_NID_to_gene_name(_C_code_disease):
    """disease NID 호출"""
    df_trash = pd.read_excel ("./Data/Disease/disease_index(Id_to_NID).xls")
    df_NID = df_trash.loc[df_trash.diseaseId == _C_code_disease]
    # print(df_NID)
    diseaseNID_index = df_NID.index
    diseaseNID = df_NID.loc[diseaseNID_index]["diseaseNID"].tolist()[0]
    # print(diseaseNID)
    conn = sqlite3.connect ("./Data/Disease/disgenet_2020.db")
    cur = conn.cursor ()

    # table 호출
    # cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    # print(cur.fetchall())
    # [('diseaseAttributes',), ('diseaseClass',), ('disease2class',), ('geneAttributes',), ('geneDiseaseNetwork',), ('variantAttributes',), ('variantGene',), ('variantDiseaseNetwork',)]

    result = cur.execute ("SELECT * FROM geneDiseaseNetwork;")
    rows = cur.fetchall ()
    cols = [column[0] for column in cur.description]

    data_df = pd.DataFrame.from_records (data=rows, columns=cols)



    # 여기까지는 문제 없음
    data_geneNID = data_df.loc[data_df.diseaseNID == diseaseNID]["geneNID"].tolist()
    # print(data_geneNID)

    df_gene_attribution = pd.read_excel("./Data/Disease/gene_attribution.xls")
    trash_gene_name_list = []
    for geneNID in data_geneNID:
        # print(geneNID)
        gene_name = df_gene_attribution.loc[df_gene_attribution.geneNID == geneNID]["geneName"].tolist()
        trash_gene_name_list = list(set(gene_name + trash_gene_name_list))

    return trash_gene_name_list
"""Gene code를 넣으면 연관되어 있는 gene name을  리턴함"""
# C0001080
# Gene symbol : ['LEMD3', 'CCND1', 'FGF2', 'FGFR3', 'NPR2', 'TP63', 'EVC', 'STAT1', 'STAT5B', 'IL36RN', 'MAPK3', 'SPRED2', 'GH1', 'IFT20', 'SLC20A1', 'NPPC', 'IL17A', 'MAPK1', 'EPHB2', 'FGFR2', 'AGT', 'COL10A1', 'POU1F1', 'CAV1', 'DCN', 'PTRH1', 'PTHLH', 'COL2A1', 'NLRP3', 'FGF1', 'STAT5A', 'IGF1', 'SNAI1', 'IGHD', 'CYP2C19', 'HDAC6', 'EDN3', 'MSX1', 'TNF', 'CNP', 'ARID1B', 'ACAN', 'SHOX', 'GRB10', 'PTH', 'BCL2']
def STRING_index_make_single_file():
    protein_raw_file = pd.read_excel("./Data/Disease/9606.protein.info.v12.0.xls")
    protein_aliases = pd.read_csv('./Data/Disease/9606.protein.aliases.v12.0.txt', sep = "\t", engine='python')
    # print(protein_aliases.head(1))
    # print(protein_raw_file.head(1))
    protein_raw_file = protein_raw_file.drop(['protein_size','annotation'], axis= 1)
    protein_aliases = protein_aliases.drop('source', axis = 1)
    protein_aliases.columns = protein_raw_file.columns
    # print(protein_aliases.head(1))
    # print(len(protein_aliases), len(protein_raw_file),len(protein_aliases) + len(protein_raw_file) )
    whole_protein_list = pd.concat([protein_raw_file, protein_aliases], axis=0)
    return whole_protein_list
"""STRING DB에서는 alias등이 많아서 파일 2개 합쳐서 index dataframe 생성함,3905590개, excel 저장 안됨 """
"""값을 2개 return 함. ENSP code와 error list """
def changing_Gene_name_to_ENSP_ID(_C_code_disease):
    chemical_protein_index = STRING_index_make_single_file()
    ENSP_code_of_protein = []
    Gene_error_list = []
    for Gene_names in SQL_NID_to_gene_name(_C_code_disease):
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
# changing_Gene_name_to_ENSP_ID("C0001080")



def Disease_dataFrame_from_C_code (_C_code_disease):
    DNf = network_construction_from_string.STRING_API_network_construct_INFO_GET (changing_Gene_name_to_ENSP_ID (_C_code_disease))
    return DNf
"""최종적으로 질병의 C code를 넣으면 Network를 구성하는 nodes들의 정보를 담고 있는 dataframe을 return한다.
쓸모 없음(STRING DB에서는 갯수 제한)"""

# print(Disease_dataFrame_from_C_code(C_code_disease))


start_time = time.time()
# print(changing_Gene_name_to_ENSP_ID("C0409959"))

# print(calculate_closest_distance(Whole_network_construction(),Herb_list("강활"), Disease_protein_list("C0409959")))

print("--- %s seconds ---" % (time.time() - start_time))
"""시간측정"""