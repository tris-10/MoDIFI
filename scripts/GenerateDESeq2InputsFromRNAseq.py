#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:26:48 2025

@author: shawt5
"""

import glob
import pandas as pd
import numpy as np
import argparse, datetime


GID_COL = 'GeneID'
geneID_COL='gene_id'
NumRep_COL='Num(rep)'
geneName_COL='GeneName'
GeneType_COL='GeneType'
CHROM_COL='CHROM'
START_COL='Start'
END_COL='End'

def importRNAseq(inputFolder):
    files = glob.glob(inputFolder)
    if files:
        dfs = []
        for f in files:
            df = pd.read_csv(f, sep='\t')  # RSEM files are usually tab-separated
            dfs.append(df)
        dfs = pd.concat(dfs, axis=0, ignore_index=True)
    else:
        print("No matching files found.")
        dfs=[]
    return dfs

def StrinToList(label):
    STRING = label.replace('"', '').replace("'", '').replace("[", '').replace("]", '')
    CELL=[]
    for cell in  STRING.split(','):
        CELL.append(cell.strip())
    return CELL

def RUN(mainDir, label, maxNumRep, filePattern, RNA_COL):
    itern = 0
    RNA_conds = pd.DataFrame(columns=['Sample', 'CellType'])
    CELL=StrinToList(label)
    MaxNum = StrinToList(maxNumRep)
    i=0
    for cellline in CELL:
        for repNum in range(1, int(MaxNum[i])+1):
            inputFolder=mainDir+cellline.lower()+'/rna_seq/rep'+str(repNum)+filePattern
            dfs = importRNAseq(inputFolder)
            if len(dfs)>0:
                dfs = dfs.groupby(geneID_COL)[RNA_COL].agg(['mean']).reset_index()
                dfs = dfs[dfs['mean']>0]
                dfs = dfs.rename(columns={'mean':cellline+'_R'+str(repNum)})
                print ('%s_R%s: %s' %(cellline, repNum, dfs.shape[0]))
                if itern == 0:
                    merged_df = dfs.copy()
                else:
                    merged_df = pd.merge(merged_df, dfs, on=geneID_COL, how='outer')
                RNA_conds.loc[itern, 'Sample'] = cellline+'_R'+str(repNum)
                RNA_conds.loc[itern, 'CellType'] = cellline
                itern+=1
        print ('')
        i+=1
    merged_df = merged_df.fillna(0)
    merged_df[geneID_COL] = pd.DataFrame(list(merged_df[geneID_COL].str.split('.')))[0]
    merged_df = merged_df.rename(columns={geneID_COL:GID_COL})
    
    RNA_ann = merged_df[GID_COL]
    RNA_ann.to_csv('RNA_ann.txt', sep='\t', index=False)
    merged_df.to_csv('RNA_counts.txt', sep='\t', index=False)
    RNA_conds.to_csv('RNA_conds.txt', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--mainDir", type=str,default='./resources/',
                        help="inputFolder")
    parser.add_argument("-l", "--label", type=str, default='["HEK293", "IMR32"]',
                    help="List of labels for cell lines")
    parser.add_argument("-r", "--maxNumRep", type=str,default='[3,3,3,3]',
                        help="maxmal rep #")
    parser.add_argument("-p", "--filePattern", type=str,default='/rsem_quant/*genes.results',
                        help="Labels for cell lines with rep#")
    parser.add_argument("--RNA_COL", type=str,default='TPM',
                        help="RNA-seq quantified column to extract (e.g., TPM, FPKM, expected_count)")
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.mainDir, args.label, args.maxNumRep, args.filePattern, args.RNA_COL)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
    
    
