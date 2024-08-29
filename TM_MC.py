import pandas as pd

### SY DB 기준으로 HERB의 Target Data 불러옴 ###

def change_korean_to_Latin(korean_herb_name):
    dir = './Data/Herb/KIOM_DB/medicinal_material.xlsx'
    f = pd.read_excel(dir)
    """
    c는 본초 한글 이름들의 list 형식으로 가지고 있음 
    d는 라틴명을 lst형식으로 가지고 있음. 
    """
    c = f["KOREAN"].tolist()
    if korean_herb_name in c:
        a = f[f["KOREAN"] == korean_herb_name]["LATIN"]
        a = a.iloc[0]
        return a
    else :
        print("맞는 이름이 없습니다.")
def Latin_to_CompoundID(_latin_name) :
    comp = pd.read_excel('./Data/Herb/KIOM_DB/medicinal_compound.xlsx')
    comp_revi = comp[comp["LATIN"] == _latin_name]
    comp_List = list(set(comp_revi["ID"].tolist()))

    """
    Compound ID 가 '0' 인 경우는 compound 값이 없기 때문에 list에서 0 값을 지워준다
    """
    remove_set = {0}
    comp_List = [i for i in comp_List if i not in remove_set]
    # print(comp_List)
    print("COMP ID 완료")
    return comp_List
def ADME_Filtering(_compound_list):
    dir_2 = './Data/Herb/KIOM_DB/chemical_property.xlsx'
    chem = pd.read_excel(dir_2)

    wanted_chem_property = chem[chem["ID"].isin(_compound_list)]
    # print(wanted_chem_property)

    ########################################ADME 조건을 변경할 수 있습니다##############################################
    global ADME_ADJ_Chem #밖으로 변수 빼내기#
    ADME_ADJ_Chem = wanted_chem_property[(wanted_chem_property["DL"] >= 0.18) & (wanted_chem_property["OB"] == 'Y')]
    print("Filtering 완료")
    print(ADME_ADJ_Chem["ID"])
    return ADME_ADJ_Chem["ID"].tolist()
    ################################################################################################################
def ID_to_Protein(_filtered_Comp_ID_list):
    dir_3 = 'Data/Herb/KIOM_DB/chemical_protein.xlsx'
    protein = pd.read_excel(dir_3)
    protein_mod2 = []
    for i in _filtered_Comp_ID_list:
        protein_mod1 = protein[protein["ID"] == i]["PREFERRED_NAME"]
        protein_mod1 = protein_mod1.tolist()
        protein_mod2 = list(set(protein_mod2 + protein_mod1))

    print("Target Protein 확보 완료")
    print(protein_mod2)
    return protein_mod2
def Herb_Target_Protein(_herb_name):
    a = ID_to_Protein(ADME_Filtering(Latin_to_CompoundID(change_korean_to_Latin(_herb_name))))

    print("Herb List 확보 최종 완료")
    print(len(a), a)
    return a
# Herb_Target_Protein("인삼")