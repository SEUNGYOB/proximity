import pandas as pd

df = pd.read_excel("Cytokine_Genes_Sorted.xlsx")
lessnames = df["Gene Name"].tolist()
morenames = df ["295 names"].tolist()

dif = set(morenames) - set(lessnames)

print(dif)