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
import glob



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


def Normalization(data, Col_Name):
    data[Col_Name] = (data[Col_Name] -  data[Col_Name].mean())/ data[Col_Name].std()
    return data

def generateSampleInfoFromDEseqOutputs(ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC):
    ## RNAseq
    dSampleInfo = pd.DataFrame(columns=[LABEL_COL, FILENAME_COL, TARGET_COL, REF_COL])
    itern=0
    for LABEL in [LABEL_RNA, LABEL_ATAC]:
        files = glob.glob(OutputDir+'*'+LABEL+'.txt')
        if files:
            for f in files:
                dSampleInfo.loc[itern, LABEL_COL] = LABEL
                dSampleInfo.loc[itern, FILENAME_COL] = f
                dSampleInfo.loc[itern, TARGET_COL] = f.split('/')[-1].split('_')[0]
                dSampleInfo.loc[itern, REF_COL] = f.split('/')[-1].split('_')[2]
                itern+=1
    ## HiC
    for cell in dSampleInfo[TARGET_COL].unique():
        pathHiC = ResourcesDir+cell.lower()+'/'+LABEL_HiC.lower()+'/'
        files = glob.glob(pathHiC+'*')
        for f in files:
            dSampleInfo.loc[itern, LABEL_COL] = LABEL_HiC
            dSampleInfo.loc[itern, FILENAME_COL] = f
            dSampleInfo.loc[itern, TARGET_COL] = cell
            itern+=1
    ## promoters
    pathPro = ResourcesDir+'Promoter/'
    files = glob.glob(pathPro+'*')
    if len(files)==1:
        dSampleInfo.loc[itern, LABEL_COL] = LABEL_PRO
        dSampleInfo.loc[itern, FILENAME_COL] = files[0]
        itern+=1
    elif len(files)>1:
        for f in files:
            dSampleInfo.loc[itern, LABEL_COL] = LABEL_PRO
            dSampleInfo.loc[itern, FILENAME_COL] = f
            if f.find('_')!=-1:
                dSampleInfo.loc[itern, TARGET_COL] = f.split('/')[-1].split('_')[1].split('.')[0]
                itern+=1
    return dSampleInfo
        

def PathReformat(ResourcesDir):
    if ResourcesDir[-1]=='/':
        ResourcesDir = ResourcesDir
    else:
        ResourcesDir = ResourcesDir+'/'  
    return ResourcesDir
        
    

def importData(SampleInfo, SamplePair, ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS=False):
    ## Import data and normalize RNA-seq amd ATAC-seq
    STOP=False
    if SampleInfo is not None and os.path.exists(SampleInfo):
        dSampleInfo = pd.read_csv(SampleInfo, sep='\t')
    #if SampleInfo is not None and os.path.exists(SampleInfo):
    #    dSampleInfo = pd.read_csv(SampleInfo), sep='\t')
    else:
        dSampleInfo = generateSampleInfoFromDEseqOutputs(ResourcesDir, OutputDir,  LABEL_RNA, LABEL_ATAC)
    if SamplePair!=None:
        dSamplePair = pd.read_csv(SamplePair, sep='\t')
        #dSamplePair = pd.read_csv(SamplePair, sep='\t')
    else:
        print ('No SamplePair file')
        sys.exit(1)
    RNA={}
    ATAC={}
    HiC={}
    PRO={}
    Target = dSamplePair.iloc[0][TARGET_COL]
    References = dSamplePair.iloc[0][REF_COL]
    REF_SET=[]
    for ref in References.split(','):
        REF_SET.append(ref.strip())
    for label in [LABEL_ATAC, LABEL_PRO]:
        sample = dSampleInfo[dSampleInfo[LABEL_COL]==label]
        if label==LABEL_ATAC:
            sample = sample[sample[TARGET_COL]==Target]
            for index in sample.index:
                ref = sample.loc[index, REF_COL]
                if ref in REF_SET:
                    Filename = sample.loc[index, FILENAME_COL]
                    Target = sample.loc[index, TARGET_COL]
                    REF = sample.loc[index, REF_COL]
                    data = pd.read_csv(Filename, sep='\t')
                    data = data[[COORDS_COL, Log2FC_COL, LFCSE_COL, Stat_Score, PValue_COL, Padj_COL]]
                    data = data.rename(columns={Stat_Score:Z_Score})
                    COORD = pd.DataFrame(list(data[COORDS_COL].str.split(':')))
                    data[CHROM_COL] = COORD[0]
                    COORD = pd.DataFrame(list(COORD[1].str.split('-')))
                    data[START_COL] = COORD[0]
                    data[END_COL] = COORD[1]
                    data[FC_COL] = 2**data[Log2FC_COL]
                    data[FCnorm_COL] = data[FC_COL]
                    Normalization(data, FCnorm_COL)
                    ATAC[Target+'/'+REF]=data
                ## Reverse
        elif label==LABEL_RNA:
            for index in sample.index:
                ref = sample.loc[index, REF_COL]
                if ref in REF_SET:
                    Filename = sample.loc[index, FILENAME_COL]
                    Target = sample.loc[index, TARGET_COL]
                    REF = sample.loc[index, REF_COL]
                    data = pd.read_csv(Filename, sep='\t')
                    data = data[data[Log2FC_COL].notna()]
                    data = data[[ GENEID_COL, Log2FC_COL, LFCSE_COL, Stat_Score, PValue_COL, Padj_COL]]
                    data = data.rename(columns={Stat_Score:Z_Score})
                    data[FC_COL] = 2**data[Log2FC_COL]
                    data[FCnorm_COL] = data[FC_COL]
                    Normalization(data, FCnorm_COL)
                    RNA[Target+'/'+REF] = data
        elif label==LABEL_HiC:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                data = pd.read_csv(Filename, sep='\t')
                data = data[data[HiC_A_COL].str.find('chr')!=-1]
                data = data[data[HiC_B_COL].str.find('chr')!=-1]
                HiC[Target] = data[[HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL, 
                                    HiC_B_COL, HiC_B_START_COL, HiC_B_END_COL]]
        elif label==LABEL_PRO:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                if np.isnan(Target)==True:
                    PRO = pd.read_csv(Filename, sep='\t')
                else:
                    PRO[Target] = pd.read_csv(Filename, sep='\t')
    return ATAC, PRO

    
#### Step
            
def mappingWithPRO(df, PRO, PRO_Region=50000, PRO_minOL=0.5):
    df[START_COL] = df[START_COL].astype(float)
    df[END_COL] = df[END_COL].astype(float)
    df[OL_COL] = np.nan
    for chrom in PRO[CHROM_COL].unique():
        PRO_chrom = PRO[PRO[CHROM_COL]==chrom]
        df_chrom = df[df[CHROM_COL]==chrom]
        total = len(PRO_chrom.index)
        itern = 0
        for index in PRO_chrom.index:
            start = PRO_chrom.loc[index, START_COL]
            end = PRO_chrom.loc[index, END_COL]
            strand = PRO_chrom.loc[index, STRAND_COL]
            gID =  PRO_chrom.loc[index, GENEID_COL]
            if strand=='-':
                start_preion = end
                end_pregion = end + PRO_Region
            elif strand=='+':
                start_preion = start - PRO_Region
                end_pregion = start
            ##
            start = start_preion
            end = end_pregion
            ##
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
                sSeE[OL_COL] =  (end - sSeE[START_COL])/(sSeE[END_COL] - sSeE[START_COL])
            if sSEe.shape[0]>0:
                sSEe[OL_COL] =  1
            if SseE.shape[0]>0:
                SseE[OL_COL] =  (end - start)/(SseE[END_COL] - SseE[START_COL])
            if SsEe.shape[0]>0:
                SsEe[OL_COL] =  (SsEe[END_COL] - start)/(SsEe[END_COL] - SsEe[START_COL])
            OL = pd.concat([sSeE, sSEe], axis=0)
            OL = pd.concat([OL, SseE], axis=0)
            OL = pd.concat([OL, SsEe], axis=0)
            if  OL.shape[0]>0:
                OL = OL[OL[OL_COL]>PRO_minOL]
                ## check whether contain preexist content
                def_temp = df.loc[OL.index]
                def_Nan = def_temp[def_temp[OL_COL].isna()]
                def_NotNan =  def_temp[def_temp[OL_COL].notna()]
                if def_Nan.shape[0]>0:
                    df.loc[def_Nan.index, OL_COL] = OL.loc[def_Nan.index, OL_COL]
                    df.loc[def_Nan.index, STRAND_COL] = strand
                    df.loc[def_Nan.index, GENEID_COL] = gID
                if def_NotNan.shape[0]>0:
                    df.loc[def_NotNan.index, OL_COL] = df.loc[def_NotNan.index, OL_COL].astype(str)+'/'+OL.loc[def_NotNan.index, OL_COL].astype(str)
                    df.loc[def_NotNan.index, STRAND_COL] = df.loc[def_NotNan.index, STRAND_COL].astype(str)+'/'+strand
                    df.loc[def_NotNan.index, GENEID_COL] = df.loc[def_NotNan.index, GENEID_COL].astype(str)+'/'+gID
            if itern%50==0:
                N = np.round((itern/total)*100, 2)
                print ('%s percent for %s to map with promoters' % (N, chrom))
            itern+=1
        print ('100 percent for %s to map with promoters' % (chrom))
        print ('')
    df[OL_COL] = df[OL_COL].astype(str)
    return df

                
def MapATACwithPRO(ATAC, PRO, OutputDir, PRO_Region=5000, PRO_minOL=0.5):
    atac_all=[]
    ATAC_temp={}
    for name in list(ATAC.keys()):
        atac = ATAC[name].copy()
        atac = atac.set_index(atac[COORDS_COL])
        ATAC_temp[name] = atac.copy()
        atac_all.append(atac)
    df = pd.concat(atac_all, ignore_index=False)  
    df = df[[COORDS_COL, CHROM_COL, START_COL, END_COL]]
    df_Coord = df.drop_duplicates()
    #df_Coord.to_csv(OutputDir+'PRO_ATAC_Coord.tsv', sep='\t', index=False)
    if isinstance(PRO, pd.DataFrame):
        df = mappingWithPRO(df_Coord, PRO, PRO_Region, PRO_minOL)
        for name in list(ATAC_temp.keys()):
            atac = ATAC_temp[name].copy()
            for col in [OL_COL, STRAND_COL,GENEID_COL]:
                atac[col] = df.loc[atac.index, col]
            ATAC_temp[name] = atac.copy()
    elif isinstance(PRO, dict):
        for name in list(ATAC.keys()):
            atac = ATAC_temp[name].copy()
            cellline = name.split('/')[0]
            promoter = PRO[cellline]
            atac = mappingWithPRO(atac, promoter, PRO_Region, PRO_minOL)
            ATAC_temp[name] = atac.copy()
    ATAC = ATAC_temp.copy()
    for name in  list(ATAC.keys()):
        ATAC[name].to_csv(OutputDir+'PRO_ATAC_'+name.replace('/','_')+'.tsv', sep='\t', index=False)
    return  ATAC, df_Coord
        
        
def RUN(SampleInfo, SamplePair, ResourcesDir, OutputDir, PRO_Region, PRO_minOL, LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS=False):
    warnings.filterwarnings("ignore")
    ResourcesDir = PathReformat(ResourcesDir)
    OutputDir = PathReformat(OutputDir)
    ATAC, PRO =  importData(SampleInfo, SamplePair, ResourcesDir, OutputDir, LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS)
    ### Annotate ATAC-seq data with promoter regions.
    ATAC, df_Coord = MapATACwithPRO(ATAC, PRO, OutputDir,  PRO_Region, PRO_minOL)
    return ATAC 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ResourcesDir", type=str,default='./resource/',
                        help="main directory")
    parser.add_argument("-si", "--SampleInfo", type=str,default='SampleInfo.tsv',
                        help="import RNA-seq (DESeq formate), ATAC-seq, HiC and GeneList")
    parser.add_argument("-sp", "--SamplePair", type=str,default='SamplePair.tsv',
                        help="import Sample Combinations")
    parser.add_argument("-o", "--OutputDir", type=str,default='./output/',
                        help="default output folder")
    parser.add_argument("-r", "--PRO_Region", type=int,default=5000,
                        help="upstream or dowmstream of a gene as promoter region")
    parser.add_argument("-m", "--PRO_minOL", type=float,default=0.5,
                        help="the minimal overlap between ATAC-seq and promoter region")  
    parser.add_argument("--LABEL_RNA", type=str,default='RNAseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument("--LABEL_ATAC", type=str,default='ATACseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument('--RUN_WITH_WARNINGS', action='store_true',
                        help="Run the script while allowing warnings")
    parser.set_defaults(RUN_WITH_WARNINGS=True)
    
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.SampleInfo, args.SamplePair, args.ResourcesDir, args.OutputDir, args.PRO_Region, args.PRO_minOL, args.LABEL_RNA, args.LABEL_ATAC, args.RUN_WITH_WARNINGS)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
