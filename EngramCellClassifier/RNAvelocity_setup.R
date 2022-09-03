# ######### RNA velocity draft
#  We are using the velocyto package here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6130801/
#    package is here: http://velocyto.org/
# 
#    
# although later this publication may be of interest as well. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009762
# 
# 
# # first we need fastq fiels downloaded which can then be converted to .bam which can then be converted to .loom for
# # velocyto
# 
# # first we need SRA toolkit to pull from SRA
# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
# - done
# 
# #getting multiple sra runs at once
# https://bioinformatics.stackexchange.com/questions/2644/download-multiple-sra-files
# -done, bash script below is close to what I did double check this
# for acc in $(file and directory here containing .txt with sra accession list);
# do srun fastq-dum $acc done
# 
# here is the SRA run selector for jeager 
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA385870&o=acc_s%3Aa

### ALIGNMENT FOR scRNA-seq

# tutorial here based on illumina pipeline:
# https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/processing-raw-scrna-seq-data.html#trimming-reads
#
# we need to do the alignment, jeager used solexa dynamic trim to trim the files
# this is a potential tutorial on how to do that:
# https://bioinformaticshome.com/tools/rna-seq/descriptions/SolexaQA.html
# - I have emailed the proff asking for their scripts or their .bam files from this study
# 
# # next we can convert fastq to .bam
#  neera shared this guide with me: https://github.com/neerapatadia/BMEG591E_Assignment1/blob/main/BMEG_591E_Assignment1.md
#  
#  there is also this guide here: https://www.youtube.com/watch?v=WLvST3OvAJ0&list=PLvK49Q0ARs91nb8lIpErL2bYU64ADmqc2&index=5
# 
# # finally we can use this guide to get smart seq2 working, jeager was sequenced with smart seq 2
#  
#  http://velocyto.org/velocyto.py/tutorial/cli.html#run-smartseq2-run-on-smartseq2-samples
# 
# Bergen, V., Soldatov, R. A., Kharchenko, P. V., & Theis, F. J. (2021). 
# RNA velocity—current challenges and future perspectives. 
# Molecular systems biology, 17(8), e10282.
# https://www.nature.com/articles/s41587-020-0591-3
# Tool: https://scvelo.readthedocs.io/about/
# 
# Bergen, V., Lange, M., Peidli, S., Wolf, F. A., & Theis, F. J. (2020).
# Generalizing RNA velocity to transient cell states through dynamical modeling. 
# Nature biotechnology, 38(12), 1408-1414.
# https://www.embopress.org/doi/full/10.15252/msb.202110282


# POssible alignment using conda: https://github.com/nikita-telkar/PHB/blob/main/data/BAM-FASTQ_FASTQ-BAM.pdf

