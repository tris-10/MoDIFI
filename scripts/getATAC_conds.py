#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 21:48:29 2025

@author: shawt5
"""

import pandas as pd
import argparse, datetime

def getATAC_conds(inputFilePath):
    #ATAC_counts = pd.read_csv(mainDir+"/ATAC_counts.txt", sep='\t')
    ATAC_counts = pd.read_csv(inputFilePath, sep='\t')
    ATAC_conds= pd.DataFrame(columns=['Sample', 'CellType'])
    ATAC_conds['Sample'] = ATAC_counts.columns[1:]
    temp=pd.DataFrame(list(ATAC_conds['Sample'].str.split('_')))
    ATAC_conds['CellType'] = temp[0]
    ATAC_conds.to_csv('ATAC_conds.txt', sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFilePath", type=str,default='./output/',
                        help="main directory")
    
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    getATAC_conds(args.inputFilePath)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
