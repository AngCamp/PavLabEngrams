## This script saves out a table of 1:1 orthologous genes between mouse and human
## using DIOPT data dump https://www.flyrnai.org/diopt 

library(tidyverse)

pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
homology <- read.delim("/space/grp/DIOPT/DIOPTvs8_export_Sanja Rogic.txt", stringsAsFactors = FALSE)
outfile <- "~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"

# prepare 1:1 high confidence orthologous genes
# Heuristics: keep scores >= 5 (Sanja recommendation) and require mutual best
# score. Remove symbols with more than one match

pc_symbol_mm <- unique(pc_mm$Symbol)
pc_symbol_hg <- unique(pc_hg$Symbol)

homology_filt <- filter(homology,
                        score >= 5 &
                        human_symbol %in% pc_hg$Symbol &
                        symbol2 %in% pc_mm$Symbol &
                        ID_type2 == "MGI" &
                        best_score == "Yes" &
                        best_score_rev == "Yes")

split_hg <- split(homology_filt, homology_filt$human_symbol)
split_mm <- split(homology_filt, homology_filt$symbol2)

which_gt1_hg <- which(unlist(lapply(split_hg, nrow)) > 1)
which_gt1_mm <- which(unlist(lapply(split_mm, nrow)) > 1)

homology_filt <- filter(homology_filt, 
                        !(human_symbol) %in% names(which_gt1_hg) &
                          !(symbol2) %in% names(which_gt1_mm))

stopifnot(all(homology_filt$human_symbol %in% pc_hg$Symbol))
stopifnot(all(homology_filt$symbol2 %in% pc_mm$Symbol))

stopifnot(identical(n_distinct(homology_filt$symbol2), 
                    n_distinct(homology_filt$human_symbol)))

# create a df with a key as not all symbols have an exact 1:1 naming match

symbols <- data.frame(
  Symbol_hg = homology_filt$human_symbol,
  Symbol_mm = homology_filt$symbol2
)
symbols$ID <- paste(symbols$Symbol_hg, symbols$Symbol_mm, sep = "_")

stopifnot(identical(n_distinct(symbols$ID), nrow(symbols)))

write.table(symbols,
            sep = "\t",
            quote = FALSE,
            file = outfile)
