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
    return dfs, files


def StrinToList(label):
    STRING = label.replace('"', '').replace("'", '').replace("[", '').replace("]", '')
    CELL=[]
    for cell in  STRING.split(','):
        CELL.append(cell.strip())
    return CELL

def RUN(mainDir, label, maxNumRep, filePattern):
    itern = 0
    RNA_conds = pd.DataFrame(columns=['Sample', 'CellType', 'Path'])
    CELL=StrinToList(label)
    MaxNum = StrinToList(maxNumRep)
    i=0
    for cellline in CELL:
        for repNum in range(1, int(MaxNum[i])+1):
            inputFolder=mainDir+cellline.lower()+'/rna_seq/rep'+str(repNum)+filePattern
            dfs, files = importRNAseq(inputFolder)
            if len(dfs)>0:
                RNA_conds.loc[itern, 'Sample'] = cellline+'_R'+str(repNum)
                RNA_conds.loc[itern, 'CellType'] = cellline
                RNA_conds.loc[itern, 'Path'] = files[0]
                itern+=1
        print ('')
        i+=1
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
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.mainDir, args.label, args.maxNumRep, args.filePattern)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
    
    
