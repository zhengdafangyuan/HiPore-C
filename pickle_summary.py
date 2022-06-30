# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 09:08:47 2022

@author: Jiayong Zhong

Nanopore Sequencing summary Pickle
"""
import numpy
import pandas as pd
import pickle
import math
import sys

def N50Len(df):
    df = df.sort_values(by="lengths", ascending=False)
    sumbase = df["lengths"].sum()
    calsumbase = 0
    for ind, rowvalue in df.iterrows():
        calsumbase += rowvalue["lengths"]
        if calsumbase >= 0.5*sumbase :
            N50 = rowvalue["lengths"]
            break
    return (N50)

pfile = sys.argv[1]
fileID = open(pfile, "rb")
Data = pickle.load(fileID)
Data["Minus"] = Data["start_time"].astype("int").div(60000000000)
Data = Data.sort_values(by="Minus")

Readcount = len(Data)
Total_bases = Data["lengths"].sum()
#P = Data["Minus"] <= 60
#f1hreads = sum(P)
#f1hbases = Data.loc[P, "lengths"].sum()
#f1hactivePores = len( set( Data.loc[P, "channelIDs"].to_list() ) )
#print(f1hactivePores, f1hreads, f1hbases)
# Pass reads, Pass bases,  Mean length, N50 read length, q10 ratio, q20 ratio, q30 ratio
PassData =  Data.loc[  Data["quals"] >= 7  , : ]
Passbases = PassData["lengths"].sum()
Passreads = len(PassData)
Menlen = Passbases/Passreads
N50 = N50Len(PassData)
#q10 = PassData.loc[PassData["quals"]>=10, "lengths"].sum()/ Passbases
#q20 = PassData.loc[PassData["quals"]>=20, "lengths"].sum()/ Passbases
#q30 = PassData.loc[PassData["quals"]>=30, "lengths"].sum()/ Passbases
#print("%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"%(Passreads, Passbases, Menlen , N50, q10, q20, q30)  )
keys = ["Readcount", "Total_bases", "Pass_Readcount", "Pass_bases", "Mean_len", "N50_len"]
values = [Readcount, Total_bases, Passreads, Passbases, Menlen, N50]

result_df = pd.DataFrame({"Item":keys,
                          "values":values}, index=keys)
resultfile = pfile+".stat"
result_df.to_csv(resultfile, sep="\t", header=True, index=False)
print("Done!")


