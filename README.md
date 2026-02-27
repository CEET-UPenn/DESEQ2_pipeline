# DESEQ2_pipeline
How we do differential gene expression analysis with DESEQ2



#Methods:
Gene annotations were downloaded from gencode v34 on 01/13/22. 

If files are not aligned, align with run_STAR.py (for this experiment, files provided were already aligned) 

Aligned bam files were sorted by coordinate using samtools 1.201. We then marked duplicate reads using picard 2.23.32 and calculated read counts using RNA-seQC 2.3.63. 
These functions are performed by bam_to_readcount_manager.py, which calls run MarkDuplicates.py and run_rnaseqc.py.

Then the gct files for each experiment and control set were combined using combine_GCTs.py

We identified differentially expressed genes using the DESeq2 1.40.25 R package, excluding genes with <10 reads across all samples and using a paired sample analysis. 
We consider as differentially expressed any gene with Benjamini-Hochberg adjusted p value < 0.05. 
This was done using DESEQ2_paired.R or DESEQ2_paired2.R, depending on experiment naming scheme

#Data Sanity Checks:

PC_expressions_varstable.R and PC_expressions_varstable2.R plot PC 1 vs 2 and 3 vs 4 for experiments 1 and 2, correcting gene expression in the same way the DE analysis did. 
These plots show how the samples cluster, ideally in a way that makes sense across activations, cell types, and controls. 

TPM_plot.R and TPM_plot_OV.R plot gene TPMs across samples to confirm knockdown or overexpression. 
