#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:33:32 2025

@author: shawt5
"""

import numpy as np
import pandas as pd
import os, sys
import argparse, datetime
import warnings
import MapATACWithPRO
import glob
from scipy import stats
import subprocess
from scipy.stats import norm



LABEL_HiC='HiC'
LABEL_PRO='Promoter'
LABEL_COL='Label'
FILENAME_COL='Filename'
TARGET_COL='Target'
REF_COL='Reference'
## DEseq
CHROM_COL='CHROM'
START_COL='Start'
END_COL='End'
COORDS_COL='Coords'
Log2FC_COL='log2FoldChange'
FC_COL = 'FoldChange'
FCnorm_COL='NormFoldChange'
GENEID_COL = 'GeneID'
GENENAME_COL='GeneName'
Stat_Score='stat'
Z_Score='Z'
LFCSE_COL='lfcSE'
PValue_COL='pvalue'
Padj_COL='padj'



HiC_A_COL='#chr1'
HiC_B_COL='chr2'
HiC_A_START_COL='x1'
HiC_A_END_COL='x2'
HiC_B_START_COL='y1'
HiC_B_END_COL='y2'
CHECK_COL='Check'
OL_COL='Overlap'
LABEL_PASS='PASS'
LABEL_XY='XY'
LABEL_X='X'
LABEL_Y='Y'
GENEID_X_COL='X_GeneID'
GENEID_Y_COL='Y_GeneID'
DACTpre_COL='DACTpre'
DACT_COL='DACT'
INDEX_COL='INDEX'
STRAND_COL='Strand'

Gnocchi_COL='Gnocchi'
Gnocchi_CHROM_COL='chrom'
Gnocchi_START_COL = 'start'
Gnocchi_END_COL = 'end'
PER_COL='Percentile'
BF_RNA_COL='RNA_BF'
BF_P_COL='P_BF'
BF_R_COL='R_BF'
BF_iga_COL = 'BF_iga'
BF_allComs_COL = 'BF_allComs'
MoDIFI_COL = 'MoDIFI'



def Normalization(data, Col_Name):
    data[Col_Name] = (data[Col_Name] -  data[Col_Name].mean())/ data[Col_Name].std()
    return data



def PathReformat(ResourcesDir):
    if ResourcesDir[-1]=='/':
        ResourcesDir = ResourcesDir
    else:
        ResourcesDir = ResourcesDir+'/'  
    return ResourcesDir
        
    

def importData(RNA_DEseq, dSampleInfo, Target, REF_SET, ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC):
    ## Import data and normalize RNA-seq amd ATAC-seq
    RNA={}
    ATAC={}
    HiC={}
    PRO={}
    REF=[] 
    for ref in REF_SET:
        for label in [LABEL_RNA, LABEL_HiC]:
            if label==LABEL_ATAC:
                Filename ='PRO_ATAC_'+ Target+'_'+ref+'.tsv'
                print (Filename)
                filepath = os.path.join(OutputDir, Filename)
                print(f"Reading file: {filepath}")
                if not os.path.exists(filepath):
                    raise FileNotFoundError(f"File not found: {filepath}")
                elif os.path.getsize(filepath) == 0:
                    raise ValueError(f"File is empty: {filepath}")
                else:
                    with open(filepath) as f:
                        print(f"First line: {f.readline()}")
                
                data = pd.read_csv(filepath, sep='\t')
                #data = pd.read_csv(OutputDir+Filename, sep='\t')
                ATAC[Target+'_'+ref]=data.copy()
            elif label==LABEL_RNA:
                Filename = Target+'_vs_'+ref+'_'+LABEL_RNA+'.txt'
                if RNA_DEseq==None:
                    data = pd.read_csv(OutputDir+Filename, sep='\t')
                else:
                    filename= ''
                    for i in RNA_DEseq:
                        if i.find(Filename)!=-1:
                            filename = i
                    if filename!='':
                        data = pd.read_csv(filename, sep='\t')
                    else:
                        print ('Missing %s' %(Filename))
                        sys.exit(2)
                data = data[data[Log2FC_COL].notna()]
                data = data[[ GENEID_COL, Log2FC_COL, LFCSE_COL, Stat_Score, PValue_COL, Padj_COL]]
                data = data.rename(columns={Stat_Score:Z_Score})
                data[FC_COL] = 2**data[Log2FC_COL]
                data[FCnorm_COL] = data[FC_COL]
                Normalization(data, FCnorm_COL)
                RNA[Target+'_'+ref] = data
            elif label==LABEL_HiC:
                sample = dSampleInfo[dSampleInfo[LABEL_COL]==label]
                sample = sample[sample[TARGET_COL]==Target]
                Filename = sample.iloc[0][FILENAME_COL]
                data = pd.read_csv(Filename, sep='\t')
                data = data[data[HiC_A_COL].str.find('chr')!=-1]
                data = data[data[HiC_B_COL].str.find('chr')!=-1]
                HiC[Target] = data[[HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL, 
                                    HiC_B_COL, HiC_B_START_COL, HiC_B_END_COL]]
    return RNA, HiC


    
#### Step
def OverLap(Target, df_chrom, start, end, OL_COLNAME = OL_COL):
    sS = df_chrom[(df_chrom[START_COL] - start)>0]
    Ss = df_chrom[(df_chrom[START_COL] - start)<=0]
    #
    seSE = sS[(sS[START_COL]-end)>0]
    sSe = sS[(sS[START_COL]-end)<=0]
    sSeE = sSe[(sSe[END_COL] - end)>0]
    sSEe = sSe[(sSe[END_COL] - end)<=0]
    #
    SsE = Ss[start<Ss[END_COL]]
    SseE = SsE[(SsE[END_COL]-end)>0]
    SsEe = SsE[(SsE[END_COL]-end)<=0]
    if sSeE.shape[0]>0:
        sSeE[OL_COLNAME] =  (end - sSeE[START_COL])/(sSeE[END_COL] - sSeE[START_COL])
    if sSEe.shape[0]>0:
        sSEe[OL_COLNAME] =  1
    if SseE.shape[0]>0:
        SseE[OL_COLNAME] =  (end - start)/(SseE[END_COL] - SseE[START_COL])
    if SsEe.shape[0]>0:
        SsEe[OL_COLNAME] =  (SsEe[END_COL] - start)/(SsEe[END_COL] - SsEe[START_COL])
    OL = pd.concat([sSeE, sSEe], axis=0)
    OL = pd.concat([OL, SseE], axis=0)
    OL = pd.concat([OL, SsEe], axis=0)
    return OL
            

                


def mappingWithHIC(Target, atac, hic):
    hic[INDEX_COL] = hic.index
    atac[START_COL] = atac[START_COL].astype(float)
    atac[END_COL] = atac[END_COL].astype(float)
    HIC_ATAC=[]
    atac_out = atac.copy()
    for chrom in range(1,25):
        if chrom<=22:
            chrom = 'chr'+str(chrom)
        elif chrom==23:
            chrom = 'chrX'
        elif chrom==24:
            chrom= 'chrY'
        hic_chrom = hic[hic[HiC_A_COL]==chrom]
        atac_chrom = atac[atac[CHROM_COL]==chrom]
        for index in hic_chrom.index:
            x1 = hic_chrom.loc[index, HiC_A_START_COL]
            x2 = hic_chrom.loc[index, HiC_A_END_COL]
            y1 = hic_chrom.loc[index, HiC_B_START_COL]
            y2 = hic_chrom.loc[index, HiC_B_END_COL]
            atac_chrom_X = OverLap(Target, atac_chrom, x1, x2, OL_COL+'_HiC')
            atac_chrom_Y = OverLap(Target,atac_chrom, y1, y2, OL_COL+'_HiC')
            if atac_chrom_X.shape[0]>0:
                atac_out.loc[atac_chrom_X.index, HiC_A_START_COL] = 1
            if atac_chrom_Y.shape[0]>0:
                atac_out.loc[atac_chrom_Y.index, HiC_B_START_COL] = 1
            if atac_chrom_X.shape[0]>0 and atac_chrom_Y.shape[0]>0:
                #print (index)
                nX = atac_chrom_X.shape[0]
                nY = atac_chrom_Y.shape[0]
                atac_chrom_X = atac_chrom_X.set_index(np.ones(nX)*index)
                atac_chrom_Y = atac_chrom_Y.set_index(np.ones(nY)*index)
                repeated_Y = pd.concat([atac_chrom_Y]*nX, ignore_index=False)
                repeated_X = pd.concat([atac_chrom_X]*nY, ignore_index=False)
                repeated_X = repeated_X.sort_values(by=[COORDS_COL])
                hic_piece = hic_chrom.loc[[index]]
                repeated_hic = pd.concat([hic_piece]*nX*nY, ignore_index=False)
                repeated_hic['X_ATAC_'+atac_chrom_X.columns] = repeated_X[repeated_X.columns]
                repeated_hic['Y_ATAC_'+atac_chrom_Y.columns] = repeated_Y[repeated_Y.columns]
            elif atac_chrom_X.shape[0]>0 and atac_chrom_Y.shape[0]==0:
                nX = atac_chrom_X.shape[0]
                atac_chrom_X = atac_chrom_X.set_index(np.ones(nX)*index)
                repeated_X = atac_chrom_X
                hic_piece = hic_chrom.loc[[index]]
                repeated_hic = pd.concat([hic_piece]*nX, ignore_index=False)
                repeated_hic['X_ATAC_'+atac_chrom_X.columns] = repeated_X[repeated_X.columns]
                repeated_hic['Y_ATAC_'+atac_chrom_Y.columns] = np.nan
            elif atac_chrom_X.shape[0]==0 and atac_chrom_Y.shape[0]>0:
                nY = atac_chrom_Y.shape[0]
                atac_chrom_Y = atac_chrom_Y.set_index(np.ones(nY)*index)
                repeated_Y = atac_chrom_Y
                hic_piece = hic_chrom.loc[[index]]
                repeated_hic = pd.concat([hic_piece]*nY, ignore_index=False)
                repeated_hic['X_ATAC_'+atac_chrom_X.columns] = np.nan
                repeated_hic['Y_ATAC_'+atac_chrom_Y.columns] = repeated_Y[repeated_Y.columns]
            elif atac_chrom_X.shape[0]==0 and atac_chrom_Y.shape[0]==0:
                repeated_hic = hic_chrom.loc[[index]]
                repeated_hic['X_ATAC_'+atac_chrom_X.columns] = np.nan
                repeated_hic['Y_ATAC_'+atac_chrom_Y.columns] = np.nan
            HIC_ATAC.append(repeated_hic)
        print ('Mapping HiC with ATAC-seq: %s finished' %(chrom))  
    HIC_ATAC = pd.concat(HIC_ATAC, ignore_index=False)
    HIC_ATAC = HIC_ATAC.set_index(np.arange(HIC_ATAC.shape[0]))
    return HIC_ATAC, atac_out

def MultipleGeneID(HIC_atac_ALL):
    XY_TABLE=[]
    for XY in [LABEL_X, LABEL_Y]:
        if XY==LABEL_X:
            GENEID = GENEID_X_COL
            LABEL_RV = LABEL_Y
        elif XY==LABEL_Y:
            GENEID = GENEID_Y_COL
            LABEL_RV = LABEL_X
        X = HIC_atac_ALL[HIC_atac_ALL[LABEL_PRO]!=LABEL_RV]
        X_single =  X[X[XY+'_ATAC_'+GENEID_COL].str.find('/')==-1]
        X_single[GENEID]= X_single[XY+'_ATAC_'+GENEID_COL]
        X_multiple = X[X[XY+'_ATAC_'+GENEID_COL].str.find('/')!=-1]
        genes = pd.DataFrame(list(X_multiple[XY+'_ATAC_'+GENEID_COL].str.split('/')))
        genes = genes.set_index(X_multiple.index)
        X_multiple_COL=[X_single]
        for col in genes.columns:
            X_multiple_col = X_multiple.copy()
            X_multiple_col[GENEID] = genes[col]
            X_multiple_COL.append(X_multiple_col)
        X_multiple_COL = pd.concat(X_multiple_COL, axis = 0)
        X_multiple_COL = X_multiple_COL[X_multiple_COL[GENEID].notna()]
        XY_TABLE.append(X_multiple_COL)
    XY_TABLE = pd.concat(XY_TABLE, axis = 0)
    return XY_TABLE
        
def LinkToHiC(ATAC, HiC, Target, OutputDir):
    HIC_ATAC_OUT={}
    HIC_noATAC_OUT = {}
    itern = 0
    for name in list(ATAC.keys()):
        filename = 'atacWithHiCasP_'+name.replace('/','_')+'.tsv'
        filenameNO = 'atacWithHiCnotP_'+name.replace('/','_')+'.tsv'
        if os.path.exists(OutputDir+filename):
            data = pd.read_csv(OutputDir+filename, sep='\t')
            HIC_ATAC_OUT[name] = data.copy()
            data = pd.read_csv(OutputDir+filenameNO, sep='\t')
            HIC_noATAC_OUT[name] = data.copy()
        else:
            atac = ATAC[name].copy()
            if itern==0:
                hic = HiC[Target]
                print ('Start mapping HiC with ATAC: %s (HiC), %s (ATAC)' %(Target, name))
                HIC_ATAC, atac = mappingWithHIC(Target, atac, hic)
                HIC_atac = HIC_ATAC[HIC_ATAC[LABEL_X+'_ATAC_'+OL_COL+'_'+LABEL_HiC].notna()]
                HIC_atac = HIC_atac[HIC_atac[LABEL_Y+'_ATAC_'+OL_COL+'_'+LABEL_HiC].notna()]
                HIC_atac_X = HIC_atac[HIC_atac['X_ATAC_'+GENEID_COL].notna()]
                HIC_atac_XY = HIC_atac_X[HIC_atac_X['Y_ATAC_'+GENEID_COL].notna()]
                HIC_atac_X = HIC_atac_X[HIC_atac_X['Y_ATAC_'+GENEID_COL].isna()]
                HIC_atac_Y = HIC_atac[HIC_atac['Y_ATAC_'+GENEID_COL].notna()]
                HIC_atac_Y = HIC_atac_Y.drop(index=HIC_atac_XY.index)    
                HIC_atac_XY[LABEL_PRO] = LABEL_XY
                HIC_atac_X[LABEL_PRO] = LABEL_X
                HIC_atac_Y[LABEL_PRO] = LABEL_Y
                HIC_atac_ALL = pd.concat([HIC_atac_X, HIC_atac_XY, HIC_atac_Y], axis = 0)
                HIC_atac_noATAC= HIC_ATAC.drop(index=HIC_atac_ALL.index)   ##  cannot calculate DACT
                XY_TABLE = MultipleGeneID(HIC_atac_ALL)
                XY_TABLE = XY_TABLE.sort_values(by=[HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL])
                HIC_atac_noATAC = HIC_atac_noATAC.sort_values(by=[HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL])
                HIC_ATAC_OUT[name] = XY_TABLE.copy()
                HIC_noATAC_OUT[name] = HIC_atac_noATAC
                LOOP_ATAC = XY_TABLE.copy()
                LOOP_nATAC = HIC_atac_noATAC.copy()
                print ('')
            else:
                atac = atac.set_index(atac[COORDS_COL])
                i=0
                for data in [LOOP_ATAC, LOOP_nATAC]:
                    DATA = data.copy()
                    for xy in [LABEL_X, LABEL_Y]:
                        Index_Col=xy+'_ATAC_'+COORDS_COL
                        data = data.set_index(data[Index_Col])
                        data_v = data[data[Index_Col].notna()]
                        data_na = data[data[Index_Col].isna()]
                        for col in [ Log2FC_COL, LFCSE_COL, Z_Score, PValue_COL, Padj_COL, FC_COL, FCnorm_COL]:
                            data_v[xy+'_ATAC_'+col] = atac[col]
                        data =pd.concat([data_v, data_na], axis = 0)
                    data = data.sort_values(by=[HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL])
                    data = data.set_index(DATA.index)
                    if i==0:
                        HIC_ATAC_OUT[name] = data.copy()
                    else:
                        HIC_noATAC_OUT[name] = data.copy()
                    i+=1
            itern+=1
    return HIC_ATAC_OUT, HIC_noATAC_OUT  

def checkRNA(HIC_ATAC_OUT, RNA, LABEL_RNA):
    HIC_ATAC_RNA={}
    for hic_atac_name in list(HIC_ATAC_OUT.keys()):
        name = hic_atac_name
        HIC_ATAC = HIC_ATAC_OUT[hic_atac_name]
        rna = RNA[name]
        rna = rna.set_index(rna[GENEID_COL])
        rna = rna.rename(columns={Stat_Score:Z_Score})
        X = HIC_ATAC[HIC_ATAC[GENEID_X_COL].notna()]
        X = X.set_index(X[GENEID_X_COL])
        X_ol = list(X.index.intersection(rna.index))
        for colname in [FC_COL, Z_Score, PValue_COL]:
            X.loc[X_ol, LABEL_RNA+'_'+colname] = rna.loc[X_ol, colname]
        ##
        Y = HIC_ATAC[HIC_ATAC[GENEID_Y_COL].notna()]
        Y = Y.set_index(Y[GENEID_Y_COL])
        Y_ol = list(Y.index.intersection(rna.index))
        for colname in [FC_COL, Z_Score, PValue_COL]:
            Y.loc[Y_ol, LABEL_RNA+'_'+colname] = rna.loc[Y_ol, colname]
        ##
        RNA_HiC_ATAC = pd.concat([X, Y], axis = 0)
        HIC_ATAC_RNA[hic_atac_name] = RNA_HiC_ATAC.copy()
    return HIC_ATAC_RNA


def CalculateBacon(HIC_ATAC_RNA, Target, REF_SET, Reference, LABEL_RNA, LABEL_ATAC, CTITLE=False,weight=1):
    label = Target+'_'+Reference
    cleaned_values = REF_SET.copy()
    if len(cleaned_values)==1:
        hic_atax_name = Target+'_'+cleaned_values[0] 
        FINAL =  HIC_ATAC_RNA[hic_atax_name].copy()
        FINAL = FINAL[FINAL[LABEL_RNA+'_'+FC_COL].notna()]
        FINAL = combineZscorePE(FINAL, cleaned_values, CTITLE)
    else:
        COMBINE = {}
        Index_All = []
        for  ref in cleaned_values:
            hic_atax_name = Target+'_'+ref 
            data = HIC_ATAC_RNA[hic_atax_name].copy()
            data = data.set_index(data[HiC_A_COL]+data[HiC_A_START_COL].astype(int).astype(str)+'-'+
                                  data[HiC_A_END_COL].astype(int).astype(str)+'-'+
                                  data[HiC_B_END_COL].astype(int).astype(str)+'-'+
                                  data['X_ATAC_Coords'].astype(str)+'-'+
                                   data['Y_ATAC_Coords'].astype(str) +'-'+
                                  data['X_GeneID'].astype(str)+'-'+
                                   data['Y_GeneID'].astype(str))
            Index_All.append((data.index))
            COMBINE[hic_atax_name] = data.copy()
        for i in range(len(Index_All)):
            if i==0: 
                OL = set(Index_All[i])
            else:
                OL = OL.intersection(set(Index_All[i]))
        OL = list(set(OL))
        comUNI = {}
        comOL={}
        for hic_atax_name in list(COMBINE.keys()):
            combineHC = COMBINE[hic_atax_name]
            combineHC_Uni = combineHC.drop(OL)
            combineHC_OL = combineHC.loc[OL]
            combineHC_OL = combineHC_OL[combineHC_OL[LABEL_RNA+'_'+FC_COL].notna()]
            for col in combineHC_OL.columns:
                combineHC_OL = combineHC_OL.rename(columns={col:hic_atax_name+'_'+col})
            comUNI[hic_atax_name] = combineHC_Uni
            comOL[hic_atax_name] = combineHC_OL
        ###
        OL = []
        for hic_atax_name in list(comOL.keys()):
            OL = OL + (list(comOL[hic_atax_name].index))
        OL = list(set(OL))
        itern = 0
        for hic_atax_name in list(comOL.keys()): 
            loop = comOL[hic_atax_name].copy()
            loop = loop.drop_duplicates()
            if itern ==0:
                FINAL = loop.copy()
            else:
                for col in loop.columns:
                    FINAL[col] = loop[col]
            itern+=1
        FINAL = FINAL.set_index(FINAL[hic_atax_name+'_'+INDEX_COL])  
        FINAL = combineZscorePE(FINAL, cleaned_values, CTITLE)
    return FINAL


def oneSidedBayesFactorZ_avoidInf(z, g=None):
    """
    Compute the one-sided Bayes factor (BF10) versus the point null using a g-prior
    and Held (2015)-style summary statistic approximation.
    """
    z = np.asarray(z, dtype=np.float64)
    if g is None:
        g = np.maximum((z ** 2) - 1.0, 0.0)
    g = np.asarray(g, dtype=np.float64)

    sqrtG = np.sqrt(g)
    sqrt1PlusG = np.sqrt(1.0 + g)

    # Posterior tail probability under H2: P(θ > 0 | data)
    a = z * sqrtG / sqrt1PlusG
    tailProb = norm.cdf(a)

    # Two-sided to one-sided correction (≤ 2)
    # Use log-space for stability.
    log_bf12 = np.log(2.0) + norm.logcdf(a)

    # Held (2015) BF20: (1/sqrt(1+g)) * exp((z^2 g)/(2(1+g)))
    # Compute in log-space to avoid overflow.
    log_bf20 = -0.5 * np.log1p(g) + 0.5 * (z * z * g) / (1.0 + g)

    # BF10 = BF12 * BF20  -> logs add
    log_bf10 = log_bf12 + log_bf20

    # Safe exponentiation
    LOG_MAX = np.log(np.finfo(np.float64).max) - 1.0
    LOG_MIN = np.log(np.finfo(np.float64).tiny) + 1.0
    def safe_exp(logx):
        return np.exp(np.clip(logx, LOG_MIN, LOG_MAX))

    bf12 = safe_exp(log_bf12)
    bf20 = safe_exp(log_bf20)
    bf10 = safe_exp(log_bf10)

    return {'BF10': bf10, 'BF20': bf20, 'tailProb': tailProb, 'g': g}



def combineZscorePE(FINAL, cleaned_values, CTITLE):
    data =FINAL.copy()
    data = data.set_index(np.arange(data.shape[0]))
    RNA_COL=data.columns[data.columns.str.find('RNAseq_Z')!=-1]
    gIDX=data.columns[data.columns.str.find('X_GeneID')!=-1]
    gIDY=data.columns[data.columns.str.find('Y_GeneID')!=-1]
    X_COL=data.columns[data.columns.str.find('X_ATAC_Z')!=-1]
    Y_COL=data.columns[data.columns.str.find('Y_ATAC_Z')!=-1]
    X_COL_P=data.columns[data.columns.str.find('X_ATAC_padj')!=-1]
    Y_COL_P=data.columns[data.columns.str.find('Y_ATAC_padj')!=-1]
    CoordX=data.columns[data.columns.str.find('X_ATAC_Coords')!=-1]
    CoordY=data.columns[data.columns.str.find('Y_ATAC_Coords')!=-1]
    for i, gIDx in enumerate(gIDX):
        if len(cleaned_values)==1:
            if CTITLE==False:
                INI=''
            else:
                INI=RNA_COL[i].split('_RNAseq_Z')[0]+'_'
        else:
            INI=RNA_COL[i].split('_RNAseq_Z')[0]+'_'
        data_X = data[data[gIDx].notna()]
        data_Y = data[data[gIDY[i]].notna()]
        data.loc[data_X.index, 'gID']=data_X[gIDx]
        data.loc[data_Y.index, 'gID']=data_Y[gIDY[i]]
        data.loc[data_X.index, INI+'P_RNAseq_Z']=data_X[RNA_COL[i]]
        data.loc[data_Y.index, INI+'P_RNAseq_Z']=data_Y[RNA_COL[i]]
        data.loc[data_X.index, INI+'P_ATAC_Z']=data_X[X_COL[i]]
        data.loc[data_Y.index, INI+'P_ATAC_Z']=data_Y[Y_COL[i]]
        data.loc[data_X.index, INI+'R_ATAC_Z']=data_X[Y_COL[i]]
        data.loc[data_Y.index, INI+'R_ATAC_Z']=data_Y[X_COL[i]]
        data.loc[data_X.index, INI+'P_Coords']=data_X[CoordX[i]]
        data.loc[data_Y.index, INI+'P_Coords']=data_Y[CoordY[i]]
        data.loc[data_X.index, INI+'R_Coords']=data_X[CoordY[i]]
        data.loc[data_Y.index, INI+'R_Coords']=data_Y[CoordX[i]]
        
        data.loc[data_X.index, INI+'P_ATAC_P']=data_X[X_COL_P[i]]
        data.loc[data_Y.index, INI+'P_ATAC_P']=data_Y[Y_COL_P[i]]
        data.loc[data_X.index, INI+'R_ATAC_P']=data_X[Y_COL_P[i]]
        data.loc[data_Y.index, INI+'R_ATAC_P']=data_Y[X_COL_P[i]]
    return data

def saveIntermediateFile(HIC_ATAC_OUT, FileINI):
    for name in list(HIC_ATAC_OUT.keys()):
        name_end = name.replace(':','_vs_').replace('/','_')
        HIC_ATAC_OUT[name].to_csv(FileINI+'_'+name_end+'.tsv', sep='\t', index=False)
        
def RUN(SampleInfo, SamplePair, RNA_DEseq, ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC, PRO_Region, PRO_minOL, CTITLE):
    ATAC = MapATACWithPRO.RUN(SampleInfo, SamplePair, ResourcesDir, OutputDir, PRO_Region, PRO_minOL, LABEL_RNA, LABEL_ATAC, True)
    TMP={}
    for cell_pair in list(ATAC.keys()):
        TMP[cell_pair.replace('/','_')] = ATAC[cell_pair]
    ATAC = TMP.copy()
    warnings.filterwarnings("ignore")
    ResourcesDir = PathReformat(ResourcesDir)
    OutputDir = PathReformat(OutputDir)
    dsamplepair = pd.read_csv(SamplePair, sep='\t')
    Target = dsamplepair.loc[0, 'Target']
    Reference = dsamplepair.loc[0, 'Reference']
    File_endstring=Target+'_vs_'+Reference.replace(' ', '').replace(',','_')
    Check = dsamplepair.loc[0, 'Check']
    if RNA_DEseq!=None:
        RNA_DEseq = RNA_DEseq.replace(' ','').split(',')
    if Check=='PASS':
        REF_SET=[]
        for ref in Reference.split(','):
            REF_SET.append(ref.strip())
        dSampleInfo = pd.read_csv(SampleInfo, sep='\t')
        RNA, HiC = importData(RNA_DEseq, dSampleInfo, Target, REF_SET, ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC)
        HIC_ATAC_OUT, HIC_noATAC_OUT = LinkToHiC(ATAC, HiC, Target, OutputDir)
        HIC_ATAC_RNA = checkRNA(HIC_ATAC_OUT, RNA, LABEL_RNA)
        FINAL = CalculateBacon(HIC_ATAC_RNA,  Target, REF_SET, Reference, LABEL_RNA, LABEL_ATAC,CTITLE,weight=1)
        FINAL.to_csv('ToBacon_'+File_endstring + '.tsv', sep='\t', index=False)
        saveIntermediateFile(HIC_ATAC_RNA, 'atacWithHiCasP')
        saveIntermediateFile(HIC_noATAC_OUT, 'atacWithHiCnotP')
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ResourcesDir", type=str,default='./resources/',
                        help="main directory")
    parser.add_argument("-si", "--SampleInfo", type=str,default='./resources/SampleInfo.tsv',
                        help="import RNA-seq (DESeq formate), ATAC-seq, HiC and GeneList")
    parser.add_argument("-sp", "--SamplePair", type=str,default='./resources/SamplePair.tsv',
                        help="import Sample Combinations")
    parser.add_argument("-o", "--OutputDir", type=str,default='./output/',
                        help="default output folder")
    parser.add_argument("--LABEL_RNA", type=str,default='RNAseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument("--LABEL_ATAC", type=str,default='ATACseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument("--PRO_Region", type=int,default=5000,
                        help="upstream or dowmstream of a gene as promoter region")
    parser.add_argument("--PRO_minOL", type=float,default=0.5,
                        help="the minimal overlap between ATAC-seq and promoter region")  
    parser.add_argument("--RNA_DEseq", type=str,default=None,
                        help="All RNA-seq outptuts after running DEseq")  
    parser.add_argument('--CTITLE', action='store_false',
                        help="Turn off cellline initial title for one pair output")
    parser.set_defaults(CTITLE=False)
    
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.SampleInfo, args.SamplePair, args.RNA_DEseq, args.ResourcesDir, args.OutputDir, 
        args.LABEL_RNA, args.LABEL_ATAC, args.PRO_Region, args.PRO_minOL, args.CTITLE)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
