#!/usr/bin/env python
# coding: utf-8

import sys, os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

#Loading Mapping paf file 
def AlignmentLoading(paffile):
    ## loading nglmr alignment paf
    ColNames=["read_name","read_length", "read_start", "read_end", "strand" ,"chrom","chrom_length", "start", "end", "Matches", "AlignBlock_length", "MapQual"]
    read_alignment_DF =  pd.read_table(paffile, sep="\t", header=None, index_col=None, usecols=[x for x in range(0,12)],names= ColNames)
    ## mapq filter
    read_alignment_DF["subread_length"] =  read_alignment_DF["read_end"] - read_alignment_DF["read_start"]
    read_alignment_DF["Identity"] = read_alignment_DF["Matches"] / read_alignment_DF["AlignBlock_length"]
    read_alignment_DF["Identity"] = np.around(read_alignment_DF["Identity"], decimals=3)
    return (read_alignment_DF)


def ReAlignmentLoading(paffile, selectCols = ["read_name","subread_length","read_start", "read_end", "strand" ,"chrom", "start", "end"] ):
    # loading minimap2
    ColNames=["subread_name","subread_length", "subread_start", "subread_end", "strand" ,"chrom","chrom_length", "start", "end", "Matches", "AlignBlock_length", "MapQual"]
    subread_alignment_DF = pd.read_table(paffile, sep="\t", header=None, index_col=None, usecols=[x for x in range(0,12)],names= ColNames)
    ## non empty dataframe
    if ( len(subread_alignment_DF) != 0 ):
        ## extract read_name, read_start, read_end
        ### the first base postion of subread is 1 , transform it to 0
        subread_alignment_DF["read_name"]=subread_alignment_DF["subread_name"].str.split(":", expand=True)[0]
        srposition = subread_alignment_DF["subread_name"].str.split(":", expand=True)[1].str.split("-", expand=True).astype("int")
        # the subreads start  has +1 in seqsk software
        srposition[0] = srposition[0] -1
        subread_alignment_DF["read_start"]= srposition[0] + subread_alignment_DF["subread_start"] 
        subread_alignment_DF["read_end"]= srposition [0] + subread_alignment_DF["subread_end"]
        subread_alignment_DF["subread_length"]= subread_alignment_DF["read_end"] - subread_alignment_DF["read_start"] 
        subread_alignment_DF["subread_name"] = subread_alignment_DF["read_name"] + ":" + srposition[0].astype("string") + "-" + srposition[1].astype("string")
        subread_alignment_DF["Identity"] = subread_alignment_DF["Matches"] / subread_alignment_DF["AlignBlock_length"]
    ## EXPORT
    if selectCols == "all":
        reAlign_DF = subread_alignment_DF
    else:
        reAlign_DF = subread_alignment_DF[selectCols]
    return (reAlign_DF)


# Split Data to samller datasets
def MakeSegmentation(chrlen, Nsection, DF):
    # groupby dataframe of : chr,  start_Section / end_Section
    ## chrsection (len) = chrlen/10  echo chromosome was split to 10 section,  start_Section,   end_Section was id of which chr section was posited 
    chrsectionlen = chrlen/Nsection # split N parts of each chrm, the length of each section
    DF["chrom"] = DF["chrom"].astype("string")
    slen_series = DF["chrom"].apply(lambda x: chrsectionlen[x] )
    ## chr index
    DF["start_Section"] = np.ceil ( (DF["start"].astype(float) + 1) / slen_series.astype(float) ) 
    DF["end_Section"] = np.ceil ( (DF["end"].astype(float)+1) / slen_series.astype(float) )
    Start_dict = DF.groupby(by=["chrom", "start_Section"]) 
    End_dict = DF.groupby(by=["chrom", "end_Section"]) 
    return (Start_dict , End_dict)

# Find reads fragment cover which Genome_virtual digestion fragement
def FindReads_GenomeFragments(rp_df, gp_df, tailtag="start", distance_thred=10):
    ## np.array position match
    Mrow, Ncol = len(gp_df), len(rp_df)
    rp_df["Position"] = rp_df[tailtag] # Find read start or end
    rp_array = rp_df["Position"].values
    gstart_array = gp_df["start"].values
    if tailtag == "start":
        dt = distance_thred
    elif tailtag == "end":
        dt = -1*distance_thred
    p_array = rp_array + dt 
    pindex = np.searchsorted(gstart_array, p_array, side="right") - 1
    align_idx = rp_df.align_idx
    genome_fragment_ids = gp_df['fragment_id'].iloc[ pindex ].values
    # fragment's start and end
    fs = gp_df['start'].iloc[ pindex ].values
    fe =  gp_df['end'].iloc[ pindex ].values
    Fragments_DF = pd.DataFrame({"align_idx":align_idx,
                                 "fragment_id": genome_fragment_ids,
                                "fragment_start":fs,
                                "fragment_end":fe})
    return (Fragments_DF)

# Multi paralllel running
def parallel_run(read_dict, Genome_vd_dict, func, tailtag="start",  cpus=1, distance_thred=10):
    #reads_vdF_df = Parallel(n_jobs=cpus)( delayed(func)(rp_df, Genome_vd_dict[gkey[0]], tailtag, distance_thred) for gkey, rp_df in read_dict.items() )
    reads_vdF_df = Parallel(n_jobs=cpus)( delayed(func)(rp_df, Genome_vd_dict[ gkey[0] ], tailtag, distance_thred) for gkey, rp_df in read_dict )
    return (pd.concat(reads_vdF_df) )

# reads virtual digestion Fragment Positioning
def read_vdFragment_Positioning(rS_vdF_DF, rE_vdF_DF, alignment_DF, pfix_thred=10):
# Use mapping strand (+/-) to determine the left and right genome virtual digestrion fragment of an alignment read's fragment; mapping chr start and end (start, end) 
    rS_vdF_DF = rS_vdF_DF.set_index("align_idx")
    rS_vdF_DF.rename(columns={"fragment_id" : "chrs_vdF_id",
                                   "fragment_start" : "chrs_vdF_start",
                                   "fragment_end" : "chrs_vdF_end"}, inplace=True)
    rE_vdF_DF = rE_vdF_DF.set_index("align_idx")
    rE_vdF_DF.rename(columns={"fragment_id" : "chre_vdF_id",
                                   "fragment_start" : "chre_vdF_start",
                                   "fragment_end" : "chre_vdF_end"}, inplace=True)
    alignment_DF = alignment_DF.set_index("align_idx", drop=False)
    vdF_merge_DF = pd.concat([alignment_DF, rS_vdF_DF, rE_vdF_DF], axis=1)
    
    # reads Left/Right vd Fragment : LvdF / RvdF
    vdF_merge_DF.loc[:,"LvdF_id"] =  vdF_merge_DF["chrs_vdF_id"]
    vdF_merge_DF.loc[:,"LvdF_start"] =  vdF_merge_DF["chrs_vdF_start"]
    vdF_merge_DF.loc[:,"LvdF_end"] =  vdF_merge_DF["chrs_vdF_end"]
    vdF_merge_DF.loc[:,"RvdF_id"] =  vdF_merge_DF["chre_vdF_id"] 
    vdF_merge_DF.loc[:,"RvdF_start"] =  vdF_merge_DF["chre_vdF_start"]
    vdF_merge_DF.loc[:,"RvdF_end"] =  vdF_merge_DF["chre_vdF_end"]
    # read vdF_start/end 
    vdF_merge_DF.loc[:,"rF_start"], vdF_merge_DF.loc[:,"rF_end"] =  vdF_merge_DF["start"], vdF_merge_DF["end"]
    
    ## read fragment mapping to the reverse strand
    reverP = vdF_merge_DF["strand"] == "-"
    vdF_merge_DF.loc[reverP,"LvdF_id"] =  vdF_merge_DF.loc[reverP,"chre_vdF_id"]
    vdF_merge_DF.loc[reverP,"LvdF_start"] =  vdF_merge_DF.loc[reverP,"chre_vdF_end"]
    vdF_merge_DF.loc[reverP,"LvdF_end"] =  vdF_merge_DF.loc[reverP,"chre_vdF_start"]
    vdF_merge_DF.loc[reverP,"RvdF_id"] =  vdF_merge_DF.loc[reverP,"chre_vdF_id"] 
    vdF_merge_DF.loc[reverP,"RvdF_start"] =  vdF_merge_DF.loc[reverP, "chrs_vdF_end"]
    vdF_merge_DF.loc[reverP,"RvdF_end"] =  vdF_merge_DF.loc[reverP, "chrs_vdF_start"]
    vdF_merge_DF.loc[reverP,"rF_start"], vdF_merge_DF.loc[reverP,"rF_end"] =  vdF_merge_DF.loc[reverP, "end"], vdF_merge_DF.loc[reverP,"start"]
    
    # read Left/Right vdF position distance and Fix
    vdF_merge_DF["LvdF_pdist"] = abs( vdF_merge_DF["rF_start"] - vdF_merge_DF["LvdF_start"] ) 
    vdF_merge_DF["RvdF_pdist"]  = abs( vdF_merge_DF["rF_end"] - vdF_merge_DF["RvdF_end"] )
    vdF_merge_DF["LRvdF_pdist"]= vdF_merge_DF["LvdF_pdist"] + vdF_merge_DF["RvdF_pdist"]
    vdF_merge_DF["LvdF_pfix"], vdF_merge_DF["RvdF_pfix"], vdF_merge_DF["LRvdF_pfix"]= False, False, False
    sp_bool = vdF_merge_DF["LvdF_pdist"] <= pfix_thred
    vdF_merge_DF.loc[sp_bool, "LvdF_pfix"] = True
    ep_bool = vdF_merge_DF["RvdF_pdist"] <= pfix_thred
    vdF_merge_DF.loc[ep_bool , "RvdF_pfix"] = True
    vdF_merge_DF.loc[(sp_bool & ep_bool), "LRvdF_pfix"] = True
    # remove some columns
    recols = ['start_Section','end_Section', 'chrs_vdF_id', 'chrs_vdF_start', 'chrs_vdF_end','chre_vdF_id', 'chre_vdF_start', 'chre_vdF_end']
    vdF_merge_DF = vdF_merge_DF.drop(recols, axis=1)
    return (vdF_merge_DF)


def rePosition(vdF_merge_DF, pfix_thred=10):
        # read Left/Right vdF position distance and Fix
    vdF_merge_DF.loc[:,"LvdF_pdist"] = abs( vdF_merge_DF["rF_start"] - vdF_merge_DF["LvdF_start"] ) 
    vdF_merge_DF.loc[:,"RvdF_pdist"]  = abs( vdF_merge_DF["rF_end"] - vdF_merge_DF["RvdF_end"] )
    vdF_merge_DF.loc[:,"LRvdF_pdist"]= vdF_merge_DF["LvdF_pdist"] + vdF_merge_DF["RvdF_pdist"]
    vdF_merge_DF.loc[:,"LvdF_pfix"], vdF_merge_DF.loc[:,"RvdF_pfix"], vdF_merge_DF.loc[:,"LRvdF_pfix"]= False, False, False
    sp_bool = vdF_merge_DF["LvdF_pdist"] <= pfix_thred
    vdF_merge_DF.loc[sp_bool, "LvdF_pfix"] = True
    ep_bool = vdF_merge_DF["RvdF_pdist"] <= pfix_thred
    vdF_merge_DF.loc[ep_bool , "RvdF_pfix"] = True
    vdF_merge_DF.loc[(sp_bool & ep_bool), "LRvdF_pfix"] = True
    return (vdF_merge_DF)

def ModifiedTailLRvdF(DF, taildistThred = 100):
    return(DF)

# Sort and Rearrange alignment index 
def ReAlignIndex(DF):
    DF.sort_values(by=["read_name", "read_start"], ascending=True, inplace=True, ignore_index=True)
    DF.loc[:, 'align_idx'] = DF.index
    return (DF)

# Cover position vdF Df to fragments pair
def ConvertPairsDF(DF):
    (nrow, ncol) = DF.shape
    lF_DF = DF.iloc[0:nrow-1].reset_index(drop=True)
    rF_DF = DF.iloc[1:nrow].reset_index(drop=True)
    ## Join adjacent fragments pair dataframe  rF1-rF2, rF2-rF3
    Fragment_Pairs_DF = lF_DF.join(rF_DF, on=None, lsuffix="_lF", rsuffix="_rF")
    Fragment_Pairs_DF = Fragment_Pairs_DF[ Fragment_Pairs_DF.read_name_lF == Fragment_Pairs_DF.read_name_rF ] ### filter frag pairs of the same read_name
    return (Fragment_Pairs_DF)


# Generate RvdF DataFrame
def GenerateRvdFDataFrame(DF, chrlen, Gvd_dict, InputParaDic, MakeSegmentation=MakeSegmentation, parallel_run=parallel_run, FindReads_GenomeFragments=FindReads_GenomeFragments, read_vdFragment_Positioning=read_vdFragment_Positioning, rePosition = rePosition, ModifiedTailLRvdF = ModifiedTailLRvdF, ReAlignIndex=ReAlignIndex):
    DF = DF.copy()
    DF = ReAlignIndex(DF)
    Nsplit = InputParaDic["Nsplit"]
    cpus= InputParaDic["pthreads"]
    distance_thred = InputParaDic["Position_distance_thred"]
    read_Sdict, read_Edict = MakeSegmentation(chrlen, Nsplit , DF)
    # Map read Fragment to genome virtrual digetstion fragments 
    start_vdF_DF = parallel_run(read_Sdict, Gvd_dict, FindReads_GenomeFragments, "start", cpus, distance_thred)
    end_vdF_DF = parallel_run(read_Edict, Gvd_dict, FindReads_GenomeFragments, "end", cpus, distance_thred)
    ## use a distance threshold to evaluated the whether the align fragment was produced by enzyme digestion
    RvdF_DF =  read_vdFragment_Positioning(start_vdF_DF, end_vdF_DF , DF, distance_thred)
    ## reposition the 5' and 3' tail
    RvdF_DF  = ReAlignIndex(RvdF_DF)
    RvdF_DF  = rePosition(RvdF_DF)
    #RvdF_DF  = ModifiedTailLRvdF(RvdF_DF)
    return(RvdF_DF)

# Filter False Fragments
## if reads with no LRvdF true fragments, then keep the best pdist fragments
def rmFalseFragments(DF, ConvertPairsDF=ConvertPairsDF ):
    # get the pfix fragments count of each reads 
    Truecount = DF["LRvdF_pfix"].groupby(by=DF["read_name"]).sum()
    Fragmentcount= DF[ "align_idx"].groupby(by=DF["read_name"]).count()
    minLRvdF_pdist = DF["LRvdF_pdist"].groupby(by=DF["read_name"]).min()
    DF["Fragment_count"] = DF["read_name"].apply(lambda x : Fragmentcount[x] )
    DF["Read_pfix_count"] = DF["read_name"].apply(lambda x :Truecount[x] )
    DF["min_pdist"] = DF["read_name"].apply(lambda x :minLRvdF_pdist[x] ) 
    # Select double true fragments or the best pdist fragment from no ture LRvdF reads
    Pdture = DF.LRvdF_pfix == True 
    P0 =  DF["Read_pfix_count"] == 0 # reads without true LRvdF
    Pmin = DF["LRvdF_pdist"] == DF["min_pdist"] # keep the best pdist
    DF = DF[ Pdture | (P0 & Pmin)]
    DF = DF.drop(['Fragment_count','Read_pfix_count', 'min_pdist'], axis=1)
    return (DF)

# Overlap Fragments
def FilterIncludedFragments(DF, ConvertPairsDF=ConvertPairsDF):
    '''
    Filter obtained overlap Fragments
    Input: RvdF position dataFrame,  ConvertPairsDF funciton
    Exorpt: RvdF position dataFrame,
    '''
    DF = DF.copy()
    # Cover DF to Fragment Pairs
    Fragment_Pairs_DF = ConvertPairsDF(DF)
    # Included overlap definde
    Included_P = ( Fragment_Pairs_DF["read_start_rF"] >= Fragment_Pairs_DF["read_start_lF"] ) & ( Fragment_Pairs_DF["read_end_rF"] <= Fragment_Pairs_DF["read_end_lF"] )
    # remove the smaller matchratio fragments 
    P1 = Fragment_Pairs_DF["Identity_lF"] >= Fragment_Pairs_DF["Identity_rF"]
    P2 = ~ P1
    lridx = Fragment_Pairs_DF.loc[Included_P & P1, "align_idx_lF"].to_list()
    lridx.extend( Fragment_Pairs_DF.loc[Included_P & P2, "align_idx_rF"].to_list() )
    lridx = list ( set(lridx) )
    if len(lridx) != 0:
        DF = DF.drop(lridx, axis=0)
        print ("Filter %d included Overlap Fragments."%len(lridx) )
        return ( FilterIncludedFragments(DF) )
    return (DF)
    

# Merge overlap Fragments
def Merge_Pairs_Fragments(DF, ConvertPairsD = ConvertPairsDF, ReAlignIndex=ReAlignIndex,  Merge_distance_thred = 50):
    '''
    Merge overlap or gap less than the Merge_distance_thred 
    Input: RvdF position dataFrame,  ConvertPairsDF funciton, Merge_distance_thred
    Exorpt: RvdF position dataFrame
    '''
    DF = DF.copy()
    #DF = DF.drop_duplicates(subset="align_idx", keep='first')
    DF = ReAlignIndex(DF)
    Fragment_Pairs_DF = ConvertPairsDF(DF)
    ## read pair distance
    Fragment_Pairs_DF.loc[:, "read_dist"] =  Fragment_Pairs_DF["read_end_lF"] - Fragment_Pairs_DF["read_start_rF"] 
    Fragment_Pairs_DF.loc[:, "align_dist"] =  Fragment_Pairs_DF["end_lF"] - Fragment_Pairs_DF["start_rF"]  
    ## Overlap Merge 
    ### satisfied all these criteria 
    P1 = Fragment_Pairs_DF["chrom_lF"]  == Fragment_Pairs_DF["chrom_rF"]
    P2 = Fragment_Pairs_DF["strand_lF"]  == Fragment_Pairs_DF["strand_rF"]
    P3 = ( Fragment_Pairs_DF["read_dist"] <= Fragment_Pairs_DF["align_dist"] + Merge_distance_thred ) & ( Fragment_Pairs_DF["read_dist"] >= Fragment_Pairs_DF["align_dist"] - Merge_distance_thred) 
    MP = P1&P2&P3
    # Merge Fragments
    lidx, ridx = Fragment_Pairs_DF.loc[MP, "align_idx_lF"].to_list(), Fragment_Pairs_DF.loc[MP, "align_idx_rF"].to_list()
    # no Merge pairs
    if len(lidx) != 0:
        ## replace right columns end and read_end
        ## Estimate the matches for overlap fragment, 90% error
        DF.loc[lidx, "Matches"] = DF.loc[lidx, "Matches"].to_numpy(dtype=int) + DF.loc[ridx, "Matches"].to_numpy(dtype=int) - 0.9*Fragment_Pairs_DF.loc[MP, "read_dist"].to_numpy(dtype=int)
        DF.loc[lidx, "end"] = DF.loc[ridx, "end"].to_numpy()
        DF.loc[lidx, "read_end"] = DF.loc[ridx,"read_end"].to_numpy()
        DF.loc[lidx,"subread_length"] = DF.loc[lidx, "read_end"].to_numpy() - DF.loc[lidx, "read_start"].to_numpy()
        DF.loc[lidx,"Identity"] = DF.loc[lidx, "Matches"].to_numpy() / (  DF.loc[lidx, "end"].to_numpy() - DF.loc[lidx, "start"].to_numpy()  )
        # merge the overlap fragments 
        lridx = lidx
        lridx.extend(ridx)
        DF = DF.drop(ridx, axis=0)
        print ("Merge %d Pairs Fragments."%len(ridx) )
        return( Merge_Pairs_Fragments(DF)  )
    return( DF )


def FilterOverlapFragments(DF, ConvertPairsD = ConvertPairsDF, overlap_distance_thred=50):
    '''
    Filter Overlap Fragments( overlap length >= overlap_distance_thred )
    Input: RvdF position dataFrame,  ConvertPairsDF funciton， overlap_distance_thred
    Exorpt: RvdF position dataFrame
    '''
    Fragment_Pairs_DF = ConvertPairsDF(DF) # Covert to Fragment_Pairs
    Fragment_Pairs_DF["read_dist"] = Fragment_Pairs_DF["read_end_lF"] - Fragment_Pairs_DF["read_start_rF"]
    P_overlap = Fragment_Pairs_DF["read_dist"] >= overlap_distance_thred
    lridx = []
     # remove the smaller matchratio fragment in the overlap pair 
    P1 = Fragment_Pairs_DF["Identity_lF"] >= Fragment_Pairs_DF["Identity_rF"]
    P2 = ~ P1
    lridx = Fragment_Pairs_DF.loc[ (P_overlap & P1), "align_idx_lF"].to_list()
    lridx.extend( Fragment_Pairs_DF.loc[ (P_overlap & P2), "align_idx_rF"].to_list() )
    lridx = list ( set(lridx) )
    DF = DF.drop(lridx, axis=0)
    print ("Filter %d Overlap Fragments."%len(lridx) )
    if len(lridx) == 0:
        return (DF)
    return ( FilterOverlapFragments(DF) )

def ExportGapBed(DF, ConvertPairsDF = ConvertPairsDF, gap_distance_thred = 50):
    '''
    Export Gap
    '''
    Fragment_Pairs_DF = ConvertPairsDF(DF) # Covert to Fragment_Pairs
    Fragment_Pairs_DF["read_dist"] = Fragment_Pairs_DF["read_start_rF"] - Fragment_Pairs_DF["read_end_lF"]
    P_gap = Fragment_Pairs_DF["read_dist"] >= gap_distance_thred
    Fragment_Pairs_DF = Fragment_Pairs_DF.loc[P_gap]
    ## Get internal subreads
    Fragment_Pairs_DF["read_name"] = Fragment_Pairs_DF["read_name_lF"]   # gap subread start
    Fragment_Pairs_DF["read_length"] = Fragment_Pairs_DF["read_length_lF"]   # gap subread start
    Fragment_Pairs_DF["subread_start"] = Fragment_Pairs_DF["read_end_lF"]   # gap subread start
    Fragment_Pairs_DF["subread_end"] = Fragment_Pairs_DF["read_start_rF"]  # gap subread end
    selectCol = ["read_name", "subread_start", "subread_end" ]
    Gap_DF =  Fragment_Pairs_DF[selectCol]
    return (Gap_DF)

# Generate tail gaps Bed
def tailGap(DF, tail_gap_distance_thred = 100):
    DF = DF.copy()
    ## reads gap in two tails
    DF["tail5_distance"] = DF["read_start"] - 0
    DF["tail3_distance"] = DF["read_length"] - DF["read_end"]
    tailgaps_DF = DF[ ["read_start", "tail5_distance", "tail3_distance"] ].groupby( DF["read_name"] ).min()
    tailgaps_DF["read_end"] =  DF[ "read_end" ].groupby( DF["read_name"] ).max()
    tailgaps_DF["read_length"] =  DF[ "read_length" ].groupby( DF["read_name"] ).max()
    ### tail5 and tail3 export
    tail5_DF = tailgaps_DF[["read_start","tail5_distance", "read_length"]][tailgaps_DF.tail5_distance >= tail_gap_distance_thred ]
    tail5_DF["subread_start"] = 0
    tail5_DF["read_name"] = tail5_DF.index
    tail5_DF = tail5_DF.rename(columns={"read_start":"subread_end"})
    tail3_DF = tailgaps_DF[["read_end","read_length","tail3_distance"]][tailgaps_DF.tail3_distance >= tail_gap_distance_thred ]
    tail3_DF["read_name"] = tail3_DF.index
    tail3_DF = tail3_DF.rename(columns={"read_end":"subread_start", "read_length":"subread_end"})
    # concat
    getcols = ["read_name","subread_start","subread_end"]
    tailgapDF = pd.concat( [tail5_DF[getcols], tail3_DF[getcols] ], axis=0 , ignore_index=True )
    return( tailgapDF )


# Multi Chromosome Checkout
def MultiChromosomeFragmentsCheck(DF):
    DF = DF.astype({"read_name":"string",
                    "chrom":"string"})
    ChrCount = DF["chrom"].groupby(by=DF["read_name"]).apply(lambda x: len( list(set(x)) )   )
    DF.loc[:,"ChrCount"] = DF["read_name"].apply(lambda x : ChrCount[x] )
    ## multi Chrom Fragments read
    ### the sum subread length of major chrom fragments is bigger than the secondary chrom 
    Pmulti = DF["ChrCount"] >= 2 
    mChr_DF = DF[ Pmulti ].copy() 
    sumbase_mChrDF = mChr_DF[['read_name', 'chrom', 'Matches']].groupby(by=['read_name', 'chrom'], as_index=False).sum()
    mChrDF_group = sumbase_mChrDF[['read_name', 'Matches' ]].groupby(by=['read_name']).max()
    sumbase_mChrDF.loc[:,"MaxM"] = sumbase_mChrDF["read_name"].apply(lambda x : mChrDF_group.loc[x,"Matches"] )
    major_ChrDF = sumbase_mChrDF[sumbase_mChrDF["MaxM"] == sumbase_mChrDF["Matches"]].copy()
    major_ChrDF = major_ChrDF.set_index("read_name")
    
    mChr_DF.loc[:, "Major_chrom"] = mChr_DF["read_name"].apply(lambda x: major_ChrDF.loc[x,"chrom"] ).astype("string").to_list()
    # get Secondary_chrom Fragments
    P = np.where( mChr_DF["Major_chrom"] != mChr_DF["chrom"] )
    Seconchrom_DF = mChr_DF.iloc[P]
    return( Seconchrom_DF )


# Exoport multi Chromosome Fragment Mapping
def MultiChromMapping(Seconchrom_DF):
    rawfq = InputParaDic["Rawfq"]
    workdir = InputParaDic["workdir"]
    bashfile = InputParaDic["ChromCheckBash"]
    checkDFlist = []
    for chrom, group in Seconchrom_DF.groupby("Major_chrom"):
        if chrom[0:3] == "chr":
            print("Check %d Fragments in %s" %(group.shape[0], chrom) )
            ## export bed file
            bedfilepath = workdir + "/" + chrom + ".subreads.check.bed"
            paffile = workdir + "/" + chrom + ".subreads.check.paf"
            group[["read_name", "read_start", "read_end"]].to_csv(bedfilepath, header=False, index=False, sep="\t" )
            ## mapping
            cmdstr = "bash %s %s %s %s %s %s"%( bashfile, rawfq, chrom, bedfilepath, paffile, workdir)
            result = os.popen(cmdstr)
            print( result.read() )
            ## loading check mapping results
            checkDFlist.append( ReAlignmentLoading(paffile, "all") )
            ## rm files
            os.system("rm %s"%(bedfilepath) )
            os.system("rm %s"%(paffile) )
    # CheckDF
    CheckDF = pd.concat(checkDFlist, axis=0, ignore_index=True)
    return(CheckDF)

def MultiChromosomeFragmentEvaluation(CheckDF, Mergealign_DF, InputParaDic):
    '''
    Evaluation of the Multichromosome Alignments
    '''
    CheckDF = CheckDF.sort_values(by=["subread_name","Matches"], ascending=False)
    #CheckDF.to_csv("Secondary_Chromosome_Mapping_Check_P1.csv", index=False,header=True)
    # Filter
    Check_filter_DF = CheckDF.drop_duplicates(subset=["subread_name"], keep='first', ignore_index=True).copy()
    Check_filter_DF.loc[:, "read_start1"] = Check_filter_DF["subread_name"].str.split(":", expand=True)[1].str.split("-", expand=True)[0].astype("int")
    Check_filter_DF.loc[:, "read_end1"] = Check_filter_DF["subread_name"].str.split(":", expand=True)[1].str.split("-", expand=True)[1].astype("int")
    IndexDF = Mergealign_DF[["read_name", "read_start", "read_end", "Matches", "align_idx", "read_length"]].set_index(["read_name", "read_start", "read_end"])
    Check_filter_DF.loc[:, 'align_idx'] = Check_filter_DF[["read_name", "read_start1", "read_end1"]].apply(lambda x: IndexDF.loc[(x["read_name"], x["read_start1"]+1,x["read_end1"]), "align_idx"], axis=1)
    Check_filter_DF.loc[:, 'read_length'] = Check_filter_DF[["read_name", "read_start1", "read_end1"]].apply(lambda x: IndexDF.loc[(x["read_name"], x["read_start1"]+1,x["read_end1"]), "read_length"], axis=1)
    Check_filter_DF.loc[:, 'OldMatches'] = Check_filter_DF[["read_name", "read_start1", "read_end1"]].apply(lambda x: IndexDF.loc[(x["read_name"], x["read_start1"]+1,x["read_end1"]), "Matches"], axis=1)
    Check_filter_DF.loc[:, 'Matches'] = Check_filter_DF['Matches'].astype("int")
    # replace to Major chromosome  
    ## to major chromosome mappping there is  0.95 more ercentage of matches in other chromosomes
    Alignreplace_thred = InputParaDic["replace_thred"]
    P_replace = Alignreplace_thred * Check_filter_DF['OldMatches'] <= Check_filter_DF['Matches']
    Check_filter_DF = Check_filter_DF[ P_replace ]
    dropCols =["subread_name", "subread_start", "subread_end","OldMatches", "read_start1", "read_end1"]
    Check_filter_DF = Check_filter_DF.set_index("align_idx", drop=False)
    Check_filter_DF.drop(dropCols, axis=1, inplace=True)
    # replace values
    Mergealign_DF = Mergealign_DF.set_index("align_idx", drop=False).copy()
    Mergealign_DF.loc[Check_filter_DF.index] = Check_filter_DF.copy()
    return(Mergealign_DF)

# Export Fragment Annotation
def GetFragmentDFValue(inDF, outDF, outcol, Pindx, S_incol, E_incol, False_incol):
    P_postive = inDF["strand"] == "+"
    P_negative = inDF["strand"] == "-"
    outDF.loc[Pindx&P_postive, outcol] = inDF.loc[Pindx&P_postive, S_incol]
    outDF.loc[Pindx&P_negative, outcol] = inDF.loc[Pindx&P_negative, E_incol]
    # if false
    outDF.loc[ (~ Pindx), outcol] = inDF.loc[(~ Pindx), False_incol]
    return (outDF)

def ReadsvdFAnnotation(RvdF_DF, GetFragmentDFValue = GetFragmentDFValue ):
    Fragref_DF = RvdF_DF[["read_name", "read_start", "read_end",  "chrom", "strand", "start", "end"]].copy()
    Fragref_DF = Fragref_DF.rename({"chrom":"ref_name"}, axis=1)
    Fragref_DF.loc[:, "rF_start"] = Fragref_DF["start"]
    Fragref_DF.loc[:,"rF_end"] = Fragref_DF["end"] 
    Fragref_DF.loc[:,"Type"] = "ivdF" # ivdF incomplete ,  vdf:complete vdF
    ## left_True and R ture
    PLtrue = RvdF_DF["LvdF_pfix"]
    PRtrue = RvdF_DF["RvdF_pfix"]
    Fragref_DF = GetFragmentDFValue(RvdF_DF, Fragref_DF, "rF_start", PLtrue, "LvdF_start", "RvdF_end", "start")
    ## RvdF_pfix == True,  ref_end = LvdF_start or RvdF_start 
    ## RvdF_pfix == False,  ref_end = chrom end 
    Fragref_DF = GetFragmentDFValue(RvdF_DF, Fragref_DF, "rF_end", PRtrue, "RvdF_end", "LvdF_start", "end")
    Fragref_DF = Fragref_DF.astype({"read_start":"int", "read_end":"int"})
    ## Type
    Fragref_DF.loc[ RvdF_DF["LRvdF_pfix"] ,"Type"] = "vdF"
    ## check  error FragrefDF
    Perror = Fragref_DF["rF_end"] - Fragref_DF["rF_start"] <= 1
    errorReads = Fragref_DF.loc[ Perror , "read_name" ].to_list()
    print( "Check %d  Passes reads."%( sum( ~ Perror ) ) )
    print( "Check %d  Errors reads."%(len(errorReads) ) )
    Fragref_DF = Fragref_DF.loc[ ~ Perror,]
    return(Fragref_DF)

# Initiation Parameter
InputParaDic = {"Rawfq" :sys.argv[1],
                "paffile" :sys.argv[2],
                "workdir" : sys.argv[3],
                "reffa" : sys.argv[4],
                "genome_vd_fragments_path": sys.argv[5],
                "ReAlignBash":sys.argv[6],
                "ChromCheckBash":sys.argv[7],
                "readsinfo":sys.argv[8],
                "remap_paf":"minimap2_subreads_realign.paf",
                "subreadbed" : "subreads.bed",
                "Nsplit" : 10,
                "pthreads" : 5,
                "Merge_distance_thred" : 50,
                "gap_distance_thred" : 30,
                "overlap_distance_thred" : 50,
                "strict_distance_thred" : 30,
                "Position_distance_thred" : 30, 
                "tail_distance_thred" : 50,
                "replace_thred" : 0.95,
                "threads": 5}

os.system( "mkdir -p %s"%InputParaDic["workdir"] )
os.chdir(InputParaDic["workdir"])
# 1 loading geneme Virtual digest data
genome_vd_fragments_path = InputParaDic["genome_vd_fragments_path"]
Genomevd_DF = pd.read_csv(genome_vd_fragments_path) ## columns [start, end, chrom, fragment_length, fragment_id]
chrlen = Genomevd_DF["end"].groupby(by=Genomevd_DF.chrom).max() # chromosome length
### geneome virtual digestion dataframe chromosome groups
Gvd_dict =  dict( list(Genomevd_DF.groupby(by=["chrom"]) ) )

# 2 Loading NGLMR Alignment Fragments
## readsInfo
readsInfo_DF = pd.read_csv(InputParaDic["readsinfo"], sep="\t", header=None, usecols=[0, 1], names=["read_name", "read_length"], index_col=0)
## Loading data
paffile= InputParaDic["paffile"]
alignment_DF = AlignmentLoading(paffile)
alignment_DF["read_length"] = alignment_DF["read_name"].apply(lambda x: readsInfo_DF.loc[x, "read_length"] )
alignment_DF = ReAlignIndex(alignment_DF)
#alignment_DF.drop(['subread_name',  'subread_start', 'subread_end'], axis=1, inplace=True)


# 3 Filter read mappping Fragments
## Filter read fragments
alignment_filter_DF = FilterIncludedFragments(alignment_DF) # inclusion fragment filter
alignment_filter_DF = Merge_Pairs_Fragments(alignment_filter_DF, Merge_distance_thred = InputParaDic["Merge_distance_thred"] ) # Merge Overlap Fragments
alignment_filter_DF = FilterOverlapFragments(alignment_filter_DF,overlap_distance_thred = InputParaDic["overlap_distance_thred"] ) # filter 
alignment_filter_DF = ReAlignIndex(alignment_filter_DF)
alignment_filter_DF["Note"] = "FirstFilter"
# 4 ReAlignment
## Export inappropreated subreads for realignment
### Generate tail gap dataframe  if tail gap >= 100, export and realign
tailGap_DF = tailGap(alignment_filter_DF, InputParaDic["tail_distance_thred"] )
### Generate internal read gap and overlap dataframe 
InternalGap_DF = ExportGapBed(alignment_filter_DF, ConvertPairsDF, InputParaDic["gap_distance_thred"])
### export subreads to bed
subreadbed = os.path.join(InputParaDic["workdir"], InputParaDic["subreadbed"]) 
subreadDF = pd.concat( [tailGap_DF, InternalGap_DF], axis=0 , ignore_index=True )
subreadDF.to_csv(subreadbed, sep="\t", header=False, index=False )
## 5 ReAlignment
paffile=os.path.join(InputParaDic["workdir"], InputParaDic["remap_paf"])
remap_cmdstr = "bash %s %s %s %s %s %s"%(InputParaDic["ReAlignBash"], InputParaDic["workdir"], InputParaDic["Rawfq"], subreadbed, paffile, InputParaDic["reffa"])
result = os.popen(remap_cmdstr)
print( result.read() )


# 5 ReAlignment Fragments Annotation
## Loading data
realignment_DF = ReAlignmentLoading(paffile, "all")
## realignmentDF read length
realignment_DF["read_length"] = realignment_DF["read_name"].apply(lambda x: readsInfo_DF.loc[x, "read_length"] )
realignment_DF = ReAlignIndex(realignment_DF)
realignment_DF.drop(['subread_name',  'subread_start', 'subread_end'], axis=1, inplace=True)
realignment_DF["Note"] = "minimap2_reAlign"

# 6 Merge Alignment and reAlign DF and Filter
Mergealign_DF =  pd.concat([alignment_filter_DF, realignment_DF], axis=0, ignore_index=False)
Mergealign_DF = ReAlignIndex(Mergealign_DF)
#  Merge Inclusion Overlap fragments
Mergealign_filter_DF = FilterIncludedFragments(Mergealign_DF) # inclusion fragment filter
Mergealign_filter_DF = ReAlignIndex(Mergealign_filter_DF)
Mergealign_filter_DF = Merge_Pairs_Fragments(Mergealign_filter_DF, ConvertPairsDF, ReAlignIndex, InputParaDic["overlap_distance_thred"])
Mergealign_filter_DF = ReAlignIndex(Mergealign_filter_DF)
Mergealign_filter_DF = FilterOverlapFragments(Mergealign_filter_DF, InputParaDic["overlap_distance_thred"] ) # filter 
Mergealign_filter_DF = ReAlignIndex(Mergealign_filter_DF)
Mergealign_filter_DF["Matches"] = Mergealign_filter_DF["Matches"].astype("int")
Mergealign_filter_DF["Identity"] = np.around(Mergealign_filter_DF["Identity"], decimals=3 )

# Supp
Checkalign_DF = Mergealign_filter_DF
# Positioning to virtual digest Fragemtns
RvdF_DF = GenerateRvdFDataFrame(Checkalign_DF, chrlen, Gvd_dict, InputParaDic)
RvdF_DF = rePosition(RvdF_DF.copy(), pfix_thred=50)
## ExportRvdF
# 8 Export reads Fragments Composition
RvdF_DF.to_csv("Read_Align_Fragment_RvdF.csv", header=True, index=False, sep=",")