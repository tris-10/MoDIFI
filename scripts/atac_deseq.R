library(DESeq2)
library(tidyverse)

run_deseq_atac <- function(dds, ann, counts, cond1, cond2,out_dir,ann_col,label) {
  ann_tibble <- rownames_to_column(data.frame(ann),var=ann_col)
  counts_tibble <- rownames_to_column(data.frame(counts),var=ann_col)
  
  
  res <- results(dds, contrast=c("CellType",cond2,cond1))
  res_tibble <- rownames_to_column(data.frame(res),var=ann_col)
  final_table <-  res_tibble  %>% left_join(counts_tibble) %>% left_join(ann_tibble)
  outfile <- file.path(out_dir, sprintf("%s_vs_%s_%s.txt", cond2, cond1, label))
  write.table(final_table, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
conds_file  <- args[2]
ann_file    <- args[3]
output_dir  <- args[4]
ann_col  <- args[5]
label <- args[6]

# Then use them like:
counts <- read.table(counts_file,header=T,sep="\t",row.names=1)
conds <- read.table(conds_file,header=T,row.names=1)
ann <- read.table(ann_file,header=T,row.names=1)

conds$CellType <- as.factor(conds$CellType)

dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=conds, design=~CellType)
dds <- DESeq(dds)

cell_types <- levels(conds$CellType)
for (i in seq_along(cell_types)) {
  for (j in seq_along(cell_types)) {
    if (i != j) {
      run_deseq_atac(dds, ann, counts, cell_types[i], cell_types[j], output_dir, ann_col, label)
    }
  }
}



