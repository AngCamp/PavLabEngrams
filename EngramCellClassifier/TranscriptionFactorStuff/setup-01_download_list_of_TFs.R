## Download the list of human and mouse TFs from the Animal TF DB
## http://bioinfo.life.hust.edu.cn/AnimalTFDB/
## Note that an older version of this script grabbed these lists from the RIKEN
## website for mouse and the Lambert 2018 publication for Human, but the RIKEN
## link broke so using Animal TFDB for consistency.

library(tidyverse)

outfile_mouse <- "~/Data/Metadata/mouse_tfs.tsv"
outfile_human1 <- "~/Data/Metadata/human_tfs.tsv"
outfile_human2 <- "~/Data/Metadata/human_tfs_lambert2018.tsv"

mouse_url <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF"
human_url1 <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
human_url2 <- "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"

# accessed sept 14 2020

if (!file.exists(outfile_mouse)) {
  download.file(url = mouse_url, destfile = outfile_mouse)
}


if (!file.exists(outfile_human1)) {
  download.file(url = human_url1, destfile = outfile_human1)
}

# Also download the Lambert 2018 list of human TFs, then format table

if (!file.exists(outfile_human2)) {
  download.file(url = human_url2, destfile = outfile_human2)
}


human_tfs <- read.delim(file = outfile_human2,
                        sep = ",",
                        row.names = NULL,
                        col.names = c(
                          "NA",
                          "Ensembl_ID",
                          "Symbol",
                          "DBD",
                          "Is_TF",
                          "Tf_assessment",
                          "Binding_mode",
                          "Motif_status",
                          "Final_Notes",
                          "Final_Comments",
                          "Interpro_IDs.",
                          "EntrezGene_ID",
                          "EntrezGene_Description",
                          "PDB_ID",
                          "TF_tested_by_HT_SELEX",
                          "TF_tested_by_PBM",
                          "Conditional_Binding_Requirements",
                          "Original_Comments",
                          "Vaquerizas_2009_classification",
                          "CisBP_considers_it_a_TF",
                          "TFCat_classification",
                          "Is_a_GO_TF_",
                          "Initial_assessment",
                          "NA",
                          "NA",
                          "TFclass_considers_it_a_TF",
                          "Go_Evidence",
                          "Pfam_Domains_By_ENSP_ID",
                          "Is_C2H2_ZF_KRAB"
                        )
)



# remove irrelevant information, and only take entries with Is.TF == TRUE
human_tfs <- human_tfs[, !str_detect(colnames(human_tfs), "NA")]
human_tfs <- human_tfs[human_tfs$Is_TF == "Yes", ]


write.table(
  human_tfs,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  outfile_human2
)



## deprecated: link no longer works

# scrape the html table of RIKEN's mouse TFs 
# mouse_url <- "http://genome.gsc.riken.jp/TFdb/tf_list.html"
# 
# scrape_table <- mouse_url %>%
#   read_html() %>%
#   html_nodes(xpath = "/html/body/table") %>%
#   html_table(header = TRUE)
# 
# scrape_table <- scrape_table[[1]]
# # remove web specific column and internal Riken ID
# scrape_table <- scrape_table[,3:ncol(scrape_table)]

# write.table(
#   scrape_table,
#   sep = "\t",
#   row.names = FALSE,
#   quote = FALSE,
#   outfile_mouse
# )