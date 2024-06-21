import networkx, random, copy
import pandas as pd
import os, numpy
import types

import networkx as nx
from matplotlib import pyplot as plt
import csv

"""9606.protein.links.v12.0 파일은 string DB의 전체 파일"""

def get_network_from_excel_file(file_name):
    df = pd.read_excel (file_name)
    g = nx.Graph ()
    g.add_nodes_from (df['protein1'])
    g.add_nodes_from (df['protein2'])
    edges = [(row['protein1'], row['protein2']) for index, row in df.iterrows ()]
    g.add_edges_from (edges)

    return g

def whole_network_to_csv(_file_name):
    # N = 30
    # with open (_file_name, encoding='utf-8') as myfile:
    #     head = [next (myfile) for x in range (N)]
    # print (len (head))
    # for l in head:
    #     print (l)
    file = open(_file_name, "r")
    txt_dataframe = pd.DataFrame()
    table_data = []
    for line in file:
        tmp_array = line.strip ().split (" ")
        while tmp_array.count (""):
            tmp_array.remove ("")

        table_data.append (tmp_array)
    txt_dataframe = pd.DataFrame(table_data)
    txt_dataframe = txt_dataframe.rename(columns=txt_dataframe.iloc[0])
    txt_dataframe =  txt_dataframe.drop(txt_dataframe.index[0])
    # print(txt_dataframe.shape)
    txt_dataframe['combined_score'] = txt_dataframe['combined_score'].astype(int)

    """경량화, combined score > 400 인 edge만 사용"""
    txt_dataframe = txt_dataframe.loc[txt_dataframe.combined_score>400]
    txt_dataframe.reset_index(drop = True, inplace = True)
    """파일읽는 속도가 빠른 json으로 파일 보내기"""
    txt_dataframe.to_csv("Whole_network.csv",sep= ',',na_rep = 'NaN', index = False)
    return

# whole_network_to_csv("./Data/Whole STRING Network/9606.protein.links.v12.0.txt")
"""역할 다 함. json 파일로 STRING 경량화 버젼 저장함"""

_file_name = "test_1.sif"

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """
	Parse sif file into node and edge sets and dictionaries
	returns setNode, setEdge, dictNode, dictEdge
	store_edge_type: if True, dictEdge[(u,v)] = edge_value
	delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set ()
    setEdge = set ()
    dictNode = {}
    dictEdge = {}
    flag = False
    f = open (file_name)
    for line in f:
        # 여기서 line은 'str' class이다.
        words = line.rstrip ("\n").split ()
        if len (words) == 3:
            setNode.add (words[0])
            setEdge.add ((words[0], words[2]))
            score = words[1]
            dictNode[words[0]] = score
            dictEdge[(words[0], words[2])] = words[1]
        else:
            setNode.add (words[0])
    f.close ()
    # print(len(setNode),len(setEdge))#node,edge counting
    return setNode, setEdge, dictNode, dictEdge


## 여기까지는 SIF file reading 완료함

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph ()
    else:
        g = networkx.Graph ()
    return g


def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file (network_file_in_sif,
                                                                               store_edge_type=use_edge_data,
                                                                               delim=delim)
    g = create_graph ()
    if include_unconnected:
        g.add_nodes_from (setNode)
    if use_edge_data:
        for e, w in dictEdge.iteritems ():
            u, v = e
            g.add_edge (u, v, w=w)  # ,{'w':w})
    else:
        g.add_edges_from (setEdge)
    return g

# print(f"Number of nodes: {len(G.nodes())}")
# print(f"Number of edges: {len(G.edges())}")
# average_clustering = nx.average_clustering(G)
# print(f"Average Clustering Coefficient: {average_clustering}")
# density = nx.density(G)
# print(f"Graph Density: {density}")


"""
네트워크 단일 분석 
"""
# # Compute degree centrality
# degree_centrality = nx.degree_centrality(G)
# print(degree_centrality)\
# # Compute closeness centrality
# closeness_centrality = nx.closeness_centrality(G)
# sorted_closeness_centrality = sorted(closeness_centrality.items(), key=lambda x:x[1], reverse=True)
# print("Closeness Centrality:", sorted_closeness_centrality)


"""Distance calculation"""
def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length (G, source_id, target_id)

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                val = get_shortest_path_length_between (network, node_from, node_to)
                values.append (val)
                d = min (values)
                # print d,
                values_outer.append (d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append (val)
                d = min (values)
                values_outer.append (d)
            d = numpy.mean (values_outer)
        # print d
    return d





def Herb_list_filtering():
    total_herb_list = ['황촉규화', '계골초', '경마자', '오가피', '자오가', '시초', '우슬', '초오제', '백부자', '초오', '부자', '제천오', '천오', '장창포',
                       '석창포', '목천료', '사삼', '사라자', '곽향', '용아초', '저백피', '근골초', '목통', '예지자', '합환화', '합환피', '택사', '대산',
                       '총백', '해백', '구자', '노회', '산강', '초두구', '고량강', '익지', '촉규화', '사인', '백두구', '초과', '천심련', '지모', '양두첨',
                       '시라자', '일당귀', '백지', '당귀', '나포마엽', '침향', '화생의', '독활', '우방자', '우방근', '주사근', '왜지차', '대복피', '빈랑자',
                       '초빈랑', '담남성', '천남성', '마두령', '천선등', '행인', '청호', '유기노', '애엽', '인진호', '한인진', '세신', '천문동', '구향충',
                       '자완', '사원자', '황기', '창출', '백출', '목향', '지각', '등피', '남판람근', '백강잠', '사간', '전가초', '전가근', '동과자', '안식향',
                       '삼과침', '화피', '권삼', '백급', '저마근', '토패모', '영와', '잠사', '용뇌', '우황', '인공우황', '체외배육우황', '우담', '운대자',
                       '개자', '소계', '저실자', '아담자', '수우각', '밀몽화', '섬피', '섬수', '시호', '광동자주', '대엽자주', '차엽', '장뇌', '장유',
                       '능소화', '청과', '도두', '마인', '반묘', '고추', '골담초근', '소두구', '남학슬', '학슬', '홍화', '모정향', '계심', '결명자', '건율',
                       '아차', '계관화', '청상자', '적설초', '아불식초', '봉랍', '녹용', '모과', '명당삼', '서청과', '백굴채', '광조', '감국', '구절초',
                       '국화', '구척', '선퇴', '국거', '승마', '육계', '계지', '육계유', '대계', '아호노', '육종용', '귤홍', '향연', '화귤홍', '불수',
                       '귤핵', '진피', '청피', '천목통', '위령선', '단혈류', '정향유', '사상자', '천궁', '당삼', '의이인', '압척초', '콘두란고', '금용담초',
                       '황련', '황련주자', '동충하초', '운지', '산수유', '고지정', '하천무', '현호색', '산사엽', '산사', '산자고', '번홍화', '파두', '필징가',
                       '선모', '강황', '울금', '아출', '토사자', '천우슬', '백미', '서장경', '백전', '백수오', '쇄양', '향부자', '강향', '양금화', '만타라엽',
                       '석곡', '철피석곡', '광금전초', '구맥', '상산', '백선피', '양지황', '황약자', '황산약', '분비해', '천산룡', '산약', '면비해', '속단',
                       '백편두', '혈갈', '골쇄보', '관중', '우주누로', '한련초', '향유', '합등자', '마황', '마황근', '음양곽', '목적', '등잔세신', '비파엽',
                       '곡정초', '정공등', '두충', '두충초탄', '두충염자', '두충엽', '귀전우', '야마추', '낭독', '비양초', '지금초', '감수', '속수자', '대극',
                       '자충', '검인', '오수유', '금교맥', '관동화', '아위', '황등', '회향', '연교', '용골', '진피', '천패모', '호북패모', '이패모', '절패모',
                       '평패모', '홍두구', '오배자', '아선약', '영지', '등황', '치자', '천마', '소박골', '원화', '겐티아나', '진교', '용담', '현초', '은행엽',
                       '백과', '인삼엽', '인삼', '홍삼', '연전초', '조협', '조각자', '해방풍', '대두황권', '흑두', '두시', '감초초', '감초밀자', '감초',
                       '석류피', '석류', '천수근', '홍기', '소통초', '훤초근', '목근피', '해마', '사극', '수질', '제조', '천년건', '맥아', '어성초', '지구자',
                       '율초', '호프', '대풍자', '감차', '경천', '천선자', '관엽금사도', '사계청', '구골엽', '구필응', '팔각회향', '지풍피', '봉선투골초',
                       '급성자', '모근', '청대', '선복화', '토목향', '금비초', '토근', '천사간', '대청엽', '판람근', '호도', '등심초', '전계혈등', '산내',
                       '시체', '해동피', '홍대극', '지부자', '애편', '건칠', '취령단초', '홍련', '곤포', '독일미', '마발', '익모초', '충위자', '정력자',
                       '고본', '여정실', '백합', '오약', '황매목', '아마인', '노로통', '풍향지', '소합향', '산맥동', '맥문동', '여지핵', '자근', '반변련',
                       '용안육', '산은화', '인동', '금은화', '담죽엽', '상기생', '사과락', '지룡', '구기엽', '구기자', '지골피', '택란', '신근초', '금전초',
                       '석적란', '후박', '신이', '공로목', '동규자', '상표초', '진주', '진주모', '통관등', '옥촉서예', '봉밀', '왕불류행', '고련피', '천련자',
                       '첨과자', '북두근', '박하', '포사엽', '서과상', '목별자', '상엽', '상심자', '상백피', '상지', '파극천', '사향', '목단피', '매화',
                       '오매', '구리향', '육두구', '몰약', '빌베리', '감송향', '석결명', '하엽', '연자심', '연방', '우절', '연자육', '연수', '흑종초자',
                       '삼칠', '백화사설초', '낙화생유', '박하유', '유향', '뇌환', '양총', '우지', '열당', '와송', '목호접', '나도근', '강활', '작약',
                       '죽절삼', '주자삼', '서양삼', '앵속각', '중루', '패장', '자소경', '자소엽', '자소자', '임자', '향가피', '도인', '식방풍', '전호',
                       '견우자', '황백', '한속단', '노근', '여감자', '죽여', '금등롱', '화산삼', '상륙', '고목', '고현삼', '호황련', '청반하', '강반하',
                       '반하', '송엽', '해송자', '송화분', '해풍등', '필발', '후추', '차전초', '차전자', '길경', '광곽향', '과자금', '원지', '옥죽', '황정',
                       '편축', '호장근', '수오등', '하수오', '제하수오', '수홍화자', '강판귀', '저령', '지실', '복령', '복신', '복령피', '마치현', '번백초',
                       '위릉채', '봉교', '하고초', '욱리인', '앵피', '금철쇄', '토형피', '태자삼', '보골지', '익수초', '갈화', '갈근', '분갈', '백두옹',
                       '녹제초', '석위', '박속', '사군자', '동릉초', '합마유', '묘조초', '내복자', '인도사목', '지황', '숙지황', '누로', '대황', '대황초탄',
                       '대황주증', '종대황', '홍경천', '만산홍', '요양화', '칠피', '피마자', '월계화', '금앵자', '금앵근', '영실', '매괴화', '송지', '복분자',
                       '천초근', '양제근', '단삼', '접골목', '지유', '백단향', '자단향', '방풍', '소목', '종절풍', '해조', '대혈등', '삼백초', '천산설련',
                       '호이초', '와릉자', '오미자', '남오미자', '형개', '오공', '낭탕근', '전갈', '현삼', '반지련', '황금', '수분초', '권백', '천리광',
                       '세네가', '번사엽', '사담', '흑지마', '희렴', '희렴주증', '수비계', '방기', '소엽련', '북유기노', '나한과', '토복령', '용규', '일지황화',
                       '괴화', '괴각', '고삼', '산두근', '삼릉', '계혈등', '부평', '야목과', '은시호', '백부근', '반대해', '해삼', '보두', '마전자', '저담',
                       '당약', '청엽담', '해룡', '정향', '맹충', '정류', '포공영', '가자', '모가자', '석명', '측백엽', '사향초', '금과람', '비자', '낙석등',
                       '종려피', '트라가칸타', '질려자', '과루', '과루피', '괄루근', '괄루인', '호로파', '산향원엽', '포황', '유백피', '조구등', '웅담', '웅과엽',
                       '지주향', '길초근', '여로', '마편초', '노봉방', '적소두', '자화지정', '곡기생', '만형엽', '만형자', '모형엽', '황형자', '천목향', '창이자',
                       '산초', '양면침', '건강', '건강초탄', '포강', '생강', '대조', '산조인', '해대']
    unknown_list = ['사라자', '근골초', '유기노', '영와', '수우각', '계관화', '봉랍', '운지', '자충', '진주', '진주모', '포사엽', '석결명', '뇌환', '죽여',
                    '해송자', '송화분', '차전초', '과자금', '욱리인', '세네가', '북유기노', '반대해', '맹충', '트라가칸타']
    already_known_list = list(set(total_herb_list) - set(unknown_list))
    b = pd.DataFrame()
    for a in already_known_list:
        herb_df = pd.read_excel("./Data/Herb/Herb_Protein_DB_v2/"+ a + "의 protein_list.xlsx")
        herb_protein_list = herb_df["Protein name"].tolist()
        herb_protein_number = len(herb_protein_list)
        c = pd.DataFrame({"Herb name" : [a],"protein number": [herb_protein_number]})
        b = pd.concat([b,c])
    # print(b)
    b = b.sort_values("protein number", ascending= False)
    b.to_csv("1234.csv",encoding='euc-kr')
    # mean = b.mean()
    # print(mean)
    return

Herb_list_filtering()