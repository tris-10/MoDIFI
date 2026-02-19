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



def runCalBF_weighted(cellline, data, Gnocchi, w_rna = 0.5, w_atacrp = 0.25):
    Gnocchi[PER_COL] = Gnocchi['z'].rank(method='min', ascending=True)/Gnocchi.shape[0]
    PRO = data.columns[data.columns.str.find(LABEL_PRO)!=-1]
    INDEX = data.columns[data.columns.str.find(INDEX_COL)!=-1]
    X = list(data.columns[data.columns.str.find('P_ATAC_Z_bacon')!=-1])
    Y = list(data.columns[data.columns.str.find('R_ATAC_Z_bacon')!=-1])
    COLS =  X+Y+['P_RNAseq_Z_bacon']
    for col in COLS:
        if col.find('_baconP')!=-1:
            COLS.remove(col)
    for Z_COL in COLS :
        BF = oneSidedBayesFactorZ_avoidInf(data[Z_COL])
        #BF = oneSidedBayesFactorZ(data[Z_COL])
        BF = BF['BF10']
        ini_COL = Z_COL.split('ATAC_Z')[0]
        if Z_COL!='P_RNAseq_Z_bacon':
            data[ini_COL+'BF'] = BF
        else:
            data[BF_RNA_COL] = BF
    data['BF_iga'] = w_atacrp * np.log(data[BF_R_COL]) + w_atacrp * np.log(data[BF_P_COL]) +  w_rna * np.log(data[BF_RNA_COL])
    data_OUT = data.copy()
    data_OUT = data_OUT.set_index(data_OUT[INDEX_COL])
    for chrom in range(1,24):
        if chrom==23:
            CHR='chrX'
        else:
            CHR = 'chr'+str(chrom)
        Gnocchi_chrom = Gnocchi[Gnocchi[Gnocchi_CHROM_COL]==CHR]
        data_chrom = data[data[HiC_A_COL]==CHR]
        data_chrom_stat = data_chrom.groupby(INDEX_COL)[[HiC_A_START_COL, HiC_A_END_COL, HiC_B_START_COL, HiC_B_END_COL]].agg('mean')
        total = len(data_chrom_stat.index)
        itern = 0
        for index in data_chrom_stat.index:
            GNO_All=[]
            for start_Cols in [[HiC_A_START_COL, HiC_A_END_COL], [HiC_B_START_COL, HiC_B_END_COL]]:
                start = data_chrom_stat.loc[index, start_Cols[0]]
                end = data_chrom_stat.loc[index, start_Cols[1]]
                GNO_A = Gnocchi_chrom[Gnocchi_chrom[Gnocchi_START_COL].between(start, end)]
                GNO_B = Gnocchi_chrom[Gnocchi_chrom[Gnocchi_END_COL].between(start, end)]
                GNO= pd.concat([GNO_A, GNO_B], axis = 0).drop_duplicates()
                GNO_All.append(GNO)
            GNO_All = pd.concat(GNO_All, axis = 0).drop_duplicates()
            GNP_value = GNO_All[PER_COL].mean()
            data_OUT.loc[index, Gnocchi_COL] = GNP_value
            if itern%100==0:
                N = np.round((itern/total)*100, 2)
                print ('%s percent on chrom %s for %s' %(N, chrom, cellline))
            itern+=1
        print ('100 percent on chrom %s for %s' %(chrom, cellline))
        print ('')
    data_OUT = data_OUT.set_index(np.arange(data_OUT.shape[0]))
    return data_OUT



def getMoDIFI_weighted(FINAL_bacon, Target, REF_SET,  Gnocchi):
    cleaned_values = REF_SET.copy()
    data = FINAL_bacon.copy()
    if len(cleaned_values)>1:
        BFiga_COL=[]
        for cellline in REF_SET:
            initial=Target+'_'+cellline
            COLs_BASE = data.columns[data.columns.str.find(initial)!=-1]
            data_BASE = data[COLs_BASE].copy()
            for col in data_BASE.columns:
                new_col = col.split(initial+'_')[1]
                data_BASE = data_BASE.rename(columns={col:new_col})
            data_BASE_OUT = runCalBF_weighted(cellline, data_BASE, Gnocchi, w_rna = 0.5, w_atacrp = 0.25)
            for col in [BF_P_COL, BF_R_COL, BF_RNA_COL, BF_iga_COL, Gnocchi_COL]:
                if col == BF_iga_COL:
                    BFiga_COL.append(initial+'_'+col)
                data[initial+'_'+col] = data_BASE_OUT[col]
        data[BF_iga_COL] = np.exp(data[BFiga_COL].mean(1))
        data[Gnocchi_COL] = data_BASE_OUT[Gnocchi_COL]
    else:
        data = data.drop(columns=["gID"])
        data_BASE_OUT = runCalBF_weighted(REF_SET[0], data, Gnocchi, w_rna = 0.5, w_atacrp = 0.25)
        data[BF_iga_COL] = np.exp(data_BASE_OUT[BF_iga_COL])
        data[Gnocchi_COL] = data_BASE_OUT[Gnocchi_COL]
    data[BF_allComs_COL] = data[BF_iga_COL] * data[Gnocchi_COL]
    data[MoDIFI_COL] = data[BF_allComs_COL]/(data[BF_allComs_COL] + 1 - data[Gnocchi_COL])
    return data

def nonRedundantDACT(MoDIFI):
    FINAL = MoDIFI.copy()
    Index_COLs=FINAL.columns[FINAL.columns.str.find(INDEX_COL)!=-1]
    INDEX = Index_COLs[0]
    FINAL_nonRudant=[]
    for index in FINAL[INDEX].unique():
        loop_index = FINAL[FINAL[INDEX]==index]
        loop_index = loop_index[loop_index[MoDIFI_COL]==loop_index[MoDIFI_COL].max()]
        FINAL_nonRudant.append(loop_index)
    FINAL_nonRudant = pd.concat(FINAL_nonRudant, axis=0)
    return FINAL_nonRudant

        
def RUN(inputFile, OutputDir, PRIOR, CTITLE):
    if os.path.exists(inputFile):
        FINAL_bacon = pd.read_csv(inputFile, sep='\t')
        Gnocchi = pd.read_csv(PRIOR, sep='\t')
        ### Get Filename
        filename = inputFile.split('/')[-1]
        File_endstring = filename.split('ToBacon_')[1].split('_bacon.tsv')[0]
        Target = File_endstring.split('_vs_')[0]
        Reference = File_endstring.split('_vs_')[1]
        REF_SET=[]
        for ref in Reference.split('_'):
            REF_SET.append(ref.strip())
        MoDIFI = getMoDIFI_weighted(FINAL_bacon, Target, REF_SET,  Gnocchi)
        MoDIFI = MoDIFI.drop(columns=[BF_allComs_COL])
        if CTITLE==True:
            initial=Target+'_'+REF_SET[0]
            for col in MoDIFI.columns:
                if col not in [MoDIFI_COL]:
                    MoDIFI = MoDIFI.rename(columns={col:initial+'_'+col})
        MoDIFI_nonRudant = nonRedundantDACT(MoDIFI)        
        MoDIFI_nonRudant.to_csv(OutputDir+'MoDIFI_loop_'+File_endstring+'.tsv', sep='\t', index=False)
        MoDIFI.to_csv(OutputDir+'MoDIFI_all_'+File_endstring+'.tsv', sep='\t', index=False)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", type=str,default='./resources/',
                        help="main directory")
    parser.add_argument("-o", "--OutputDir", type=str,default='./output/',
                        help="default output folder")
    parser.add_argument("-p", "--PRIOR", type=str,default='./resources/Gnocchi.tsv',
                        help="default output folder")
    parser.add_argument('--CTITLE', action='store_false',
                        help="Turn off cellline initial title for one pair output")
    parser.set_defaults(CTITLE=False)
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.inputFile, args.OutputDir, args.PRIOR, args.CTITLE)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
