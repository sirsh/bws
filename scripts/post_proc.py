import pandas as pd
import numpy as np


import chardet
rawdata = open("./out.data", 'rb').read()
result = chardet.detect(rawdata)
charenc = result['encoding']

q = input("please provide a UNIQUE run qualifier for output files [e.g RUN07_90_X]:> ")

print("processing the out.data file in the local folder")
template = "#(HIST)#"
with open ("./out.data", newline="\n",encoding=charenc) as f:
    with open("./data/{}_histogram.csv".format(q), "w") as g:
        for l in f:
            if l[:len(template)] == template:
                g.write(l[len(template):])  


data = pd.read_csv("./out.data",sep='\t', comment='#',encoding=charenc)
data["chunk"] = 1.0
data = data[[ "chunk", "L","t", "M0","M1","M2","M3","M4","M5","M6","M7","M8"]]
for m in ["M"+str(i+1) for i in range(8)]:   data[m] = data[m].astype(float) / data["M0"].astype(float)

print("aggregating data...")

grp = data.groupby(["t"])
mn = grp.mean().drop("chunk",1).drop("M0",1).drop("L",1)
std = grp.std().drop("chunk",1).drop("M0",1).drop("L",1) / np.sqrt(grp.count())#count should be the chunk size?
alldata = mn.join(std, rsuffix="error")
cols = np.sort(alldata.columns.values)
alldata = alldata[cols]
alldata.to_csv("./data/{}_moments.csv".format(q))

print("files saved to /data/{}_".format(q))
