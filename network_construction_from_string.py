import pandas as pd
import requests
import networkx as nx

#참조 : https://string-db.org/help/api/
"""준비물 : gene list"""
def checking_STRING_version():
    request_url = "https://string-db.org/api/tsv-no-header/version"
    response = requests.post(request_url)
    for line in response.text.strip().split("\n"):
        current_server = line.split ("\t")[1]
    return current_server
"""현재의 STRING version을 API로 확인해서 알려줌"""

def STRING_API_network_construct_INFO_GET(_genelist):
    string_api_url = str(checking_STRING_version()) +"/api"
    output_format = "tsv-no-header"
    # output_format = "tsv"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    my_genes = _genelist
    params = {

        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : 9606, # species NCBI identifier
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    nodes_from = []
    nodes_to =[]
    response = requests.post(request_url, data=params)

    for line in response.text.strip ().split ("\n"):
        l = line.strip ().split ("\t")
        # print(l)

        """에러가 떴을 경우 확인해봐야 할 l"""
        p1, p2 = l[0],l[1]
        # filter the interaction according to experimental score
        try :

            combined_score = float (l[5])
            if combined_score >0.4:
                nodes_from.append(p1)
                nodes_to.append(p2)

        except :
            print("ERROR")
    # print (len (nodes_to), len (nodes_from))
    # print ("\t".join ([p1, p2, "experimentally confirmed (prob. %.3f)" % combined_score]))
    df_STRING = pd.DataFrame ({'Nodes_from': nodes_from,
                               'Nodes_to': nodes_to
                               })
    # print(df_STRING)

    return df_STRING
"""combined score >0.4  => medium credibility
   nodes들의 조합들을 통해서 네트워크를 구성한 다음에 dataframe으로 만들어서 리턴해줌
"""


#C0001080의 네트워크 예시
# gene_ENSP_ID_List = ['9606.ENSP00000369889', '9606.ENSP00000480791', '9606.ENSP00000264498']
