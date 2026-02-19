library(DESeq2)
library(tidyverse)
library(tximport)

run_deseq <- function(conds, counts, dds, ann, cond1, cond2, out_dir, ann_col, label) {
    res <- results(dds, contrast=c("CellType",cond1, cond2))
    rownames(res) <- str_split_i(rownames(res), "\\.", 1)
    res_tibble <- rownames_to_column(data.frame(res),var=ann_col)
    # prep annotation table
    ann_tibble <- rownames_to_column(ann,var=ann_col)
    #Prep count table
    #counts <- counts(dds,normalized=TRUE)
    rownames(counts) <- str_split_i(rownames(counts), "\\.", 1)
    #colnames(counts) <- conds$Sample
    counts_tibble <- rownames_to_column(data.frame(counts),var=ann_col)
    #Slap everything together
    final_table <-  ann_tibble  %>% left_join(res_tibble) %>% left_join(counts_tibble)
    outfile <- file.path(out_dir, sprintf("%s_vs_%s_%s.txt", cond1, cond2, label))
    write.table(final_table, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
conds_file  <- args[1]
ann_file    <- args[2]
output_dir  <- args[3]
ann_col  <- args[4]
label <- args[5]

# Then use them like:
conds <- read.table(conds_file,header=T,row.names=1)
ann <- read.table(ann_file,header=T,row.names=1, sep="\t")


files <- file.path(paste(conds$Path))
txi.rsem <- tximport(files, type="rsem",txIn=FALSE,txOut=FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, colData=conds, design=~CellType)
dds <- DESeq(ddsTxi)
counts <- counts(dds,normalized=TRUE)

conds$CellType <- as.factor(conds$CellType)
cell_types <- levels(conds$CellType)
for (i in seq_along(cell_types)) {
  for (j in seq_along(cell_types)) {
    if (i != j) {
      run_deseq(conds, counts, dds, ann, cell_types[i], cell_types[j], output_dir, ann_col, label)
    }
  }
}



