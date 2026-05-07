#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 17:26:48 2025

@author: shawt5
"""

import glob, os
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

def importRNAseq(file):
    if os.path.exists(file):
        df = pd.read_csv(file, sep='\t')
    else:
        print("No matching files found.")
        df=[]
    return df


def StrinToList(label):
    STRING = label.replace('"', '').replace("'", '').replace("[", '').replace("]", '')
    CELL=[]
    for cell in  STRING.split(','):
        CELL.append(cell.strip())
    return CELL

def RUN(mainDir, label, FileName):
    itern = 0
    RNA_conds = pd.DataFrame(columns=['Sample', 'CellType', 'Path'])
    CELL=StrinToList(label)
    FileName=StrinToList(FileName)
    for i, cellline in enumerate(CELL):
        files = glob.glob(FileName[i])
        MaxNum = len(files)
        for repNum in range(1, int(MaxNum)+1):
            file =files[repNum-1]
            df = importRNAseq(file)
            if len(df)>0:
                RNA_conds.loc[itern, 'Sample'] = cellline+'_R'+str(repNum)
                RNA_conds.loc[itern, 'CellType'] = cellline
                RNA_conds.loc[itern, 'Path'] = file
                itern+=1
    RNA_conds.to_csv('RNA_conds.txt', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--mainDir", type=str,default='./resources/',
                        help="inputFolder")
    parser.add_argument("-l", "--label", type=str, default='["HEK293", "IMR32"]',
                    help="List of labels for cell lines")
    parser.add_argument("-f", "--FileName", type=str,default='[./resources/hek293/rna_seq/*.genes.results, ./resources/imr32/rna_seq/*.genes.results]',
                        help="maxmal rep #")
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.mainDir, args.label, args.FileName)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
    
    
