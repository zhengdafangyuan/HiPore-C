#!/usr/bin/env python
# coding: utf-8
# v1 generate juicebox Medium format contact matrix from map.paf
# v2 use datatable import and filter
import sys,gc
import pandas as pd
import numpy as np
import itertools
import datetime
import datatable as dt

import sys,gc
import pandas as pd
import numpy as np
import itertools
import datetime
import datatable as dt

# loadtables 
def LoadTables(filename, usecols):
    dt_df = dt.fread(filename, sep=",", header=True, columns=usecols)
    # MapQual  and LRvdF_pfix filter: if mapqual <= 10 and LRvdF_pfix == 0
    del dt_df[(dt.f.MapQual<=10) & (dt.f.LRvdF_pfix),:]
    df  = dt_df.to_pandas()
    del dt_df
    return(df)

def PrintTime():
    now = datetime.datetime.now()
    otherStyleTime = now.strftime("%Y--%m--%d %H:%M:%S")
    print(otherStyleTime)

def ReAlignIndex(gDF): 
    gDF = gDF.sort_values(by=["chrom", "Position"], ascending=True,  ignore_index=True) 
    return (gDF)

def ContactMatrix(gDF):
    gDF = ReAlignIndex(gDF)
    binidx = BinDict[len(gDF)] 
    i1 = [ t[0] for t in binidx ] 
    i2 = [ t[1] for t in binidx ] 
    contactDF = pd.DataFrame({"read_name":gDF.iloc[i1]["read_name"].to_list(), 
                              "str1":gDF.iloc[i1]["strand"].to_list(), 
                              "chr1":gDF.iloc[i1]["chrom"].to_list(), 
                              "pos1":gDF.iloc[i1]["Position"].values,
                              "frag1":0,
                              "str2":gDF.iloc[i2]["strand"].values, 
                              "chr2":gDF.iloc[i2]["chrom"].to_list(), 
                              "pos2":gDF.iloc[i2]["Position"].values,
                              "frag2":1, 
                               "mapq1":gDF.loc[i1]["MapQual"].values,
                              "mapq2":gDF.loc[i2]["MapQual"].values } ) 
    return(contactDF)

def IniFiles(Exportfile):
    try:
        f=open(Exportfile,'w')
        if f:
            f.truncate()
    except Exception as e:
        print(e)

def ExportFun(Exportfile, gcontactDF, Nloop):
    # export
    ContactDF = pd.DataFrame()
    ContactDF = ContactDF.append(gcontactDF, ignore_index=True) 
    ContactDF = ContactDF.astype({"str1":"int", "str2":"int","pos1":"int","pos2":"int",
                              "frag1":"int","frag2":"int","mapq1":"int","mapq2":"int"})
    ContactDF.to_csv(Exportfile, sep=" ", header=False, index=False, mode="a")
    print("Generated %d pairs contacts in %d reads." %(len(ContactDF), Nloop) )
    return 1

# fragments pairs order
def Binpairs(n):
    lis = [ i for i in range(0, n) ]
    pairs = list(itertools.combinations(list(lis),2))
    return(pairs)

# Loading
RvdFfile = sys.argv[1]
Exportfile = sys.argv[2]
IniFiles(Exportfile) # Initiation files
print("Generate Contact Matrix for %s"%(RvdFfile) )
#read_name,read_length,read_start,read_end,strand,chrom,chrom_length,start,end,Matches,AlignBlock_length,MapQual,subread_length,FragMatchratio,align_idx,Note,LvdF_id,LvdF_start,LvdF_end,RvdF_id,RvdF_start,RvdF_end,rF_start,rF_end,LvdF_pdist,RvdF_pdist,LRvdF_pdist,LvdF_pfix,RvdF_pfix,LRvdF_pfix
usecols = {"read_name","read_start","read_end","strand","chrom","start","end","MapQual","LRvdF_pfix"}
RvdF_DF = LoadTables(RvdFfile, usecols)

# Processing
RvdF_DF["Position"] = (RvdF_DF["start"].to_numpy() + RvdF_DF["end"].to_numpy())/2
RvdF_DF["Position"] =  RvdF_DF["Position"].astype("int")
RvdF_DF = RvdF_DF.loc[:, ["read_name","chrom", "Position", "strand", "MapQual"]]
RvdF_DF = RvdF_DF.drop_duplicates(subset=["read_name”，“chrom", "Position"], keep='first'
selectchrs = [ "chr%d"%i for i in range(1,22+1) ]
selectchrs.extend(["chrX", "chrY"])
RvdF_DF = RvdF_DF[RvdF_DF["chrom"].isin(selectchrs)]
## read name list 
readcount = RvdF_DF.groupby("read_name")["chrom"].count()
readlist = readcount[readcount>=2].index.to_list()
## fragment contact pairs order
binlist = list( set( readcount[readcount>=2].to_list() ) )
BinDict = {}
for n in binlist :
    BinDict[n] = Binpairs(n)
## filter fragment < 2
RvdF_DF.set_index("read_name", drop=False, inplace=True)
RvdF_DF = RvdF_DF.loc[readlist, :]
print("Loading %d reads and %d fragments"%(len(readlist), len(RvdF_DF) ) )
# #Change strand
## forward :0;  reverse : -1
RvdF_DF.loc[RvdF_DF.strand=="-", "strand"] = -1
RvdF_DF.loc[RvdF_DF.strand=="+", "strand"] = 0
Nloop, N, readsNum = 20000, 0, 0
TotalReadsNum = len(readlist)
for N in range(0, TotalReadsNum, Nloop ):   
    gc.collect() 
    Ne = N + Nloop
    if Ne >= TotalReadsNum:
        Ne = TotalReadsNum
    PrintTime() # Time
    print("Processing reads: %d - %d"%(N, Ne) )
    ## get subreads 
    subreads = readlist[N:Ne]
    subRvdF_DF = RvdF_DF.loc[subreads, :]
    # ReadGroup Contacts
    gcontactDF = []
    print("Generate Contact")
    gcontactDF = [ ContactMatrix( subRvdF_DF.loc[read_idx, :] ) for read_idx in  subreads ] # 列表推导式
    status = ExportFun(Exportfile, gcontactDF, Nloop)

readcount = RvdF_DF.groupby(RvdF_DF.index)["chrom"].count()
paircount = sum(readcount*(readcount-1)/2)
print("Finished! \nTherotically it can produce %d contact pairs."%paircount)