library(data.table)
library(bacon)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]



file_list <- c(input_file)

for (in_path in file_list) {
  
  message("Processing: ", in_path)
  dat <- fread(in_path, sep = "\t", header = TRUE)
  
  # ---- pick Z columns ----
  cols_rna  <- grep("[P]_.*RNAseq_Z$", colnames(dat), value = TRUE)
  cols_atac <- grep("[PR]_.*ATAC_Z$", colnames(dat), value = TRUE)
  cols_coords <- grep("[PR]_.*Coords$$", colnames(dat), value = TRUE)
  cols <- c(cols_rna, cols_atac)
  for (i in seq_along(cols_rna)) {
    cols_coords <- c('gID', cols_coords)
  }
  
  
  for (i in seq_along(cols)) {
    col <- cols[i]
    z   <- dat[[col]]
    gid <- dat[[ cols_coords[i] ]]  # <- change if your gene id column name differs
    
    #idx_keep <- which(is.finite(z) & z != 0 & !is.na(gid))
    idx_keep <- which(is.finite(z)  & z != 0 & !is.na(gid))
    
    if (length(idx_keep) < 10) {
      message("  ", col, ": too few values (", length(idx_keep), "), skip")
      next
    }
    
    z_keep   <- z[idx_keep]
    gid_keep <- gid[idx_keep]
    
    # 1) same gene ID -> same integer group id
    gid_levels <- unique(gid_keep)                # keeps first-seen order
    grp_id     <- match(gid_keep, gid_levels)     # integers 1..K
    K          <- length(gid_levels)
    
    # 2) choose ONE Z per gene ID (unique Z set by gene ID)
    #    Use median if each gene can appear multiple times with slightly different Z
    z_by_gid <- tapply(z_keep, grp_id, median, na.rm = TRUE)
    z_by_gid <- as.numeric(z_by_gid)              # length K
    
    # 3) run bacon on the per-gene Zs
    bc_u   <- bacon::bacon(teststatistics = z_by_gid)
    adjZ_u <- bacon::tstat(bc_u)
    
    # 4) map back: every row gets its gene's bacon-adjusted Z
    adjZ_all <- adjZ_u[grp_id]
    adjP_all <- 2 * pnorm(-abs(adjZ_all))
    
    # 5) write back aligned to original rows
    newZ <- paste0(col, "_bacon")
    newP <- paste0(col, "_baconP")
    
    if (!newZ %in% names(dat)) dat[, (newZ) := NA_real_]
    if (!newP %in% names(dat)) dat[, (newP) := NA_real_]
    
    set(dat, i = idx_keep, j = newZ, value = adjZ_all)
    set(dat, i = idx_keep, j = newP, value = adjP_all)
    
    message("  ", col, ": keep=", length(idx_keep), " unique_gID=", K)
  }
  
  
  # ---- output file ----
  base <- basename(in_path)
  out_base <- if (grepl("\\.txt$", base)) {
    sub("\\.txt$", "_bacon.txt", base)
  } else if (grepl("\\.tsv$", base)) {
    sub("\\.tsv$", "_bacon.tsv", base)
  } else {
    paste0(base, "_bacon")
  }
  
  out_path <- file.path(dirname(in_path), out_base)
  fwrite(dat, file = out_path, sep = "\t", quote = FALSE)
  message("  → Written: ", out_path)
}