#generate juicebox Medium format contact matrix from mapping paf
import sys,gc
import pandas as pd
import os
import itertools
import datetime
import datatable as dt
import argparse

# v1 Complete df
# v2 change dataframe to list export

# Load matrix reader
def LoadMatrixReader(filename, Chunksize, colnamelist, usecols,  sepstr="\t"):
    file_reader = pd.read_table(filename, sep=sepstr,
                                chunksize=Chunksize, iterator=True,
                                header=0, index_col=None, names = colnamelist, usecols = usecols)
    return (file_reader)

def CompleteDF(df, dfhold, Chunksize):
    '''
    Complete the isolated reads alignment records
    '''
    lastread = df.iloc[-1]["read_name"]
    df = pd.concat([dfhold, df]) # concat dfhold and df
    # last read df
    P = df["read_name"] == lastread
    dfhold = df.loc[P, :].copy()
    if len(df) >= Chunksize: # not the last iterally loading
        df = df.drop(df.loc[P].index.to_list() , axis=0)
    return(df, dfhold)

def PrintTime():
    now = datetime.datetime.now()
    otherStyleTime = now.strftime("%Y--%m--%d %H:%M:%S")
    print(otherStyleTime)

def ReAlignIndex(gDF): 
    gDF = gDF.sort_values(by=["chrom", "read_Position","chr_Position"], ascending=True,  ignore_index=True) 
    return (gDF)

def ReAlignIndex_2(gDF):
    gDF = gDF.sort_values(by=["chrom", "chr_Position"], ascending=True,  ignore_index=True) 
    return (gDF)


def Binpairs(n):
    '''
    all pairs
    '''
    lis = [ i for i in range(0, n) ]
    pairs = list(itertools.combinations(list(lis),2))
    return(pairs)

def Binpairs_2(n):
    '''
    adjcent pairs
    '''
    lis_1 = [ i for i in range(0,n) ] 
    lis_2 = [ i for i in range(1,n+1) ]
    pairs = [ (lis_1[i],lis_2[i]) for i in range(0,n-1) ]
    return(pairs)

def Binpairs_3(n):
    '''
    Non-adj pairs
    '''
    pairs_1 = Binpairs(n)
    pairs_2 = Binpairs_2(n)
    pairs = sorted(list(set(pairs_1) ^ set(pairs_2)))
    return(pairs)

def IniFiles(Exportfile):
    '''
    Initiation
    '''
    try:
        f=open(Exportfile,'w')
        if f:
            f.truncate()
    except Exception as e:
        print(e)

def Fragpre(RvdF_DF):
    '''
    mapping Fragment preprocessing
    Input : RvdF_DF
    Output : RvdF_DF
    '''
    # fragment file processing
    RvdF_DF["chr_Position"] = (RvdF_DF["start"].to_numpy() + RvdF_DF["end"].to_numpy())/2
    RvdF_DF["chr_Position"] =  RvdF_DF["chr_Position"].astype("int")
    RvdF_DF["read_Position"] = RvdF_DF["read_start"].copy()
    RvdF_DF = RvdF_DF.loc[:, ["read_name","read_Position", "chrom", "chr_Position", "strand", "MapQual"]]
    ## filter fragment < 2
    readcount = RvdF_DF.groupby("read_name")["chrom"].count()
    readlist = readcount[readcount>=2].index.to_list() ## read name list 
    RvdF_DF.set_index("read_name", drop=False, inplace=True)
    RvdF_DF = RvdF_DF.loc[readlist, :]
    print("Loading %d reads and %d fragments"%(len(readlist), len(RvdF_DF) ) )
    # #Change strand
    ## forward :0;  reverse : -1
    RvdF_DF.loc[RvdF_DF.strand=="-", "strand"] = -1
    RvdF_DF.loc[RvdF_DF.strand=="+", "strand"] = 0
    return(RvdF_DF, len(readlist) )

def Contactrow(binidx, gDF):
    '''
    Generate contact matrix row for each binindx
    '''
    i1 = [ t[0] for t in binidx ] 
    i2 = [ t[1] for t in binidx ]
    read_name = gDF.iloc[i1]["read_name"].to_list()
    str1 = gDF.iloc[i1]["strand"].to_list()
    chr1 = gDF.iloc[i1]["chrom"].to_list()
    pos1 = gDF.iloc[i1]["chr_Position"].to_list()
    frag1 = 0
    str2 = gDF.iloc[i2]["strand"].to_list()
    chr2 = gDF.iloc[i2]["chrom"].to_list()
    pos2 = gDF.iloc[i2]["chr_Position"].to_list()
    frag2 = 1
    mapq1 = 60
    mapq2 = 60
    # all contacts
    if len(binidx) >= 1:
        contactlist = [ '{} {} {} {} {} {} {} {} {} {} {}'.format( read_name[idx], str1[idx], chr1[idx], pos1[idx], frag1,
                                                           str2[idx], chr2[idx], pos2[idx], frag2,
                                                           mapq1, mapq2) for idx in range( len(i1) ) ]
        contactlist = "\n".join(contactlist)
    else:
        contactlist = "-"
    return(contactlist)


def ContactMatrix(gDF,SeparateTag=0):
    '''
    Generate contact matrix
    contact matrix type : all, adj, nonadj
    '''
    gDF = ReAlignIndex(gDF)
    Nfrag = len(gDF)
    if SeparateTag == 0:
        all_binidx = Binpairs(Nfrag)
        contactlist = Contactrow(all_binidx, gDF)
        contactnum = len(all_binidx)
        return(contactlist, contactnum)
    else:
        adj_binidx = Binpairs_2(Nfrag)
        adj_contact = Contactrow(adj_binidx, gDF)
        nonadj_binidx = Binpairs_3(Nfrag)
        nonadj_contact = Contactrow(nonadj_binidx, gDF)
        adj_num, nonadj_num = len(adj_binidx), len(nonadj_binidx)
        return(adj_contact, adj_num, nonadj_contact, nonadj_num)


def ExportFun(Exportfile, contactlist):
    # export
    contactlist = [item+"\n" for item in contactlist if item != "-" ] # remove empty
    with open(Exportfile, "a") as fileID:
        fileID.writelines( "".join(contactlist)  ) 
    fileID.close()
    return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--paf", help="the mapping file",
                        type=str, required=True)
    parser.add_argument('-o',"--out_dir", type=str, required=True,
                        help='the output dir')
    parser.add_argument("--prefix", type=str, required=False, default="",
                        help='the output filename prefix')
    parser.add_argument('-c',"--Chunksize", type=int, required=False, default=500000,
                    help='the chunksize of iterater file reader, default is 500000')
    parser.add_argument('-s', "--separate", type=int, required=False,
                        default=0,
                        help='Whether to export contact matrix for adjcent and non-adjcent pair-contact fragments, defalut is false(0)')
    parser.add_argument('-t','--threads', type=int, required=False,
                        default=5,
                        help='Number of threads to generate contact matrix')

    argv = parser.parse_args()
    paffile = argv.paf 
    Exportdir = argv.out_dir
    prefix = argv.prefix
    Chunksize = argv.Chunksize
    if prefix != "":
        prefix = prefix + "_"
    SeparateTag = argv.separate
    threads =  argv.threads
    # Export filename
    if SeparateTag == 0:
        Exportfile = os.path.join( Exportdir, "%scontact_matrix.txt"%prefix )
        IniFiles(Exportfile)
    else:
        Exportfile_1 = os.path.join( Exportdir, "%sAdj_contact_matrix.txt"%prefix)
        Exportfile_2 = os.path.join( Exportdir, "%sNonadj_contact_matrix.txt"%prefix)
        IniFiles(Exportfile_1) # Initiation files
        IniFiles(Exportfile_2) # Initiation files
    statfile = os.path.join( Exportdir, "generate_pair_contact.stat") # stat summary file
    IniFiles(statfile) # Initiation files

    # Loading
    print("Generate Contact Matrix for %s"%(paffile) )
    ## 3 Laoding RvdF reader
    
    colnames = ['read_name', 'read_length', 'read_start', 'read_end', 
                'strand', 'chrom','chrom_length', 'start', 'end', 'MapQual']
    usecols = [0,1,2,3,4,5,6,7,8,9]
    reader  =  LoadMatrixReader(paffile, Chunksize, colnames, usecols,  sepstr=",") ##  map pafreader
    # Parallle generate pairwise contacts
    from pandarallel import pandarallel
    pandarallel.initialize(nb_workers=threads, progress_bar=False)
    ## stat dict
    Statdict = {"readnum": 0,
                "Pc":0,
                "adj_Pc":0, 
                "nonadj_Pc":0}
    N, Nloop, Ne =  0, 0, 0
    pairs_count = 0
    df_hold = pd.DataFrame()
    for subRvdF_DF in reader:
        gc.collect() # 释放内存
        PrintTime() # Time
        subRvdF_DF, df_hold = CompleteDF(subRvdF_DF, df_hold, Chunksize) 
        subRvdF_DF, readcount = Fragpre(subRvdF_DF) 
        ## get subreads 
        #subreads = list( set( subRvdF_DF["read_name"].to_list() ) )
        Nloop = readcount
        Statdict["readnum"] += Nloop # read num
        Ne = N + Nloop
        print("Processing reads: %d - %d"%(N, Ne-1) )
        N = Ne
        ## generate pair contacts
        subGroup = subRvdF_DF.groupby(subRvdF_DF.read_name)
        if SeparateTag == 0:
            ctype = "all"
            contactlist, contactnumlist = zip( *subGroup.parallel_apply(ContactMatrix, SeparateTag=SeparateTag) )
            status = ExportFun(Exportfile, contactlist)
            all_pc = sum(contactnumlist)
            Statdict["Pc"] += all_pc # read num
            contactlist = None
            print("From %d reads generated %d pairwise contacts." %(readcount, all_pc ) )
        else:
            # adj and nonadj
            adjlist, adjnumlist, nonadjlist, nonadjnumlist = zip( *subGroup.parallel_apply(ContactMatrix, SeparateTag=SeparateTag) )
            
            status_1 = ExportFun(Exportfile_1, adjlist)
            status_2 = ExportFun(Exportfile_2, nonadjlist)
            adj_pc, nonadj_pc = sum(adjnumlist), sum(nonadjnumlist)
            Statdict["adj_Pc"] += adj_pc # pairwise contact contact
            Statdict["nonadj_Pc"] += nonadj_pc 
            adjlist, adjnumlist = None, None
            print("From %d reads generated %d adjcent and %d non adjcent pairwise contacts." %(readcount, adj_pc, nonadj_pc) )
    # Export summary stat
    Statdict["Pc"] = Statdict["adj_Pc"] + Statdict["nonadj_Pc"]
    Statdf  =  pd.DataFrame(Statdict, index=[paffile])
    Statdf.to_csv(statfile, header=True, index=False)
    print("Finished!")