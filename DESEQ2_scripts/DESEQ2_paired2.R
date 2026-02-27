library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)

#run from within directory containing gct files. will output into same directory 

#use args to get the input files
args = commandArgs(TRUE)

#testing purposes:
#args <- c("OV-Control-LPS,OV-Control-LPS_readcount.gct.gz","OV-Uhrf2-LPS,OV-Uhrf2-LPS_readcount.gct.gz")
#args <- c("OV-Control-PBS,OV-Control-PBS_readcount.gct.gz","OV-Uhrf2-PBS,OV-Uhrf2-PBS_readcount.gct.gz")
#setwd("/project/voightlab_01/lorenzk/CEET/Amatullah/bam")

trait1stuff <- strsplit(args[1], ",")[[1]]
trait2stuff <- strsplit(args[2], ",")[[1]]
#assumes each arguement is traitname,path/to/gctfile

trait1 <- read.delim(trait1stuff[2], header = TRUE, skip=2)
trait2 <- read.delim(trait2stuff[2], header = TRUE, skip=2)

##merge counts files & print
merged <- inner_join(trait1, select(trait2, !Description), "Name")
rownames(merged) <- merged$Name
merged[,2] <- NULL
merged[,1] <- NULL

##make metadata file
meta <- data.frame(cbind(colnames(merged)))
colnames(meta) <- "header"
meta <- separate(meta, header, c("exp", "OV", "type1", "type2", "rep", "repcount"), sep = "[.]", remove = FALSE)

#exp_data <- merged[rowSums(merged[])>=10,]
#exp_data <- rownames_to_column(exp_data)

#run DESeq2

#import data
DE_data <- DESeqDataSetFromMatrix(countData = merged, 
                                  colData = meta,
                                  design = formula(paste("~ repcount + type1")))

#prefilter
keep <- rowSums(counts(DE_data)) >= 10
DE_data <- DE_data[keep,]

#run DESeq
#DE_out <- DESeq(DE_data)
DE_data <- estimateSizeFactors(DE_data)
DE_data <- estimateDispersions(DE_data)
DE_data <- nbinomWaldTest(DE_data, maxit=1000)
res <- results(DE_data, alpha=0.05)
summary(res)

resultsNames(DE_data)

png(paste(trait2stuff[1], "_cooksboxplot.png", sep=''), width=7, height=7, unit="in", res=300);
boxplot(log10(assays(DE_data)[["cooks"]]), range=0, las=2)
dev.off()
#check degrees of freedom
ncol(model.matrix(design(DE_data), colData(DE_data)))

#add gene names
results <- data.frame(res)
results <- rownames_to_column(results)
DE <- separate(results, col = "rowname", into = c("ENSG", "transcript?"), sep = "[.]")
genes <- as.data.frame(genes(EnsDb.Hsapiens.v86))
DE_named <- left_join(DE, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("ENSG" = "gene_id") )

#rearrange things to be pretty
DE_named <- unite(DE_named, ENSG, c("ENSG", "transcript?"), sep = ".")
DE_named <- relocate(DE_named, c("gene_name", "gene_biotype", "seqnames", "start", "end"), .after = "ENSG")
DE_named <- rename(DE_named, "chr" = "seqnames")


#output results files
resultsfile <- str_c(trait2stuff[1], "_DESEQ_results_paired.txt")
write.table(DE_named, file= resultsfile, sep="\t", row.names=FALSE, quote=FALSE)

#pull significant genes
resSig <- subset(DE_named, padj < 0.05)
DE_named$Significant <- ifelse(DE_named$padj < 0.05, TRUE, FALSE)


resultsfile <- str_c(trait2stuff[1], "_DESEQ_results_paired_sig.txt")
write.table(resSig, file= resultsfile, sep="\t", row.names=FALSE, quote=FALSE)

#make volcano plot
forplot <- subset(DE_named, !is.na(Significant))

png(paste(trait2stuff[1], "_volcano.png", sep=''), width=7, height=7, unit="in", res=300);
ggplot(forplot, aes(x=log2FoldChange, y=-log10(padj), color = Significant)) +
  geom_point() +
  ggtitle(paste(trait2stuff[1], " Volcano", sep='')) +
  scale_color_brewer(palette = "Paired")
dev.off()

#png(paste(trait2stuff[1], "_volcano_zoom.png", sep=''), width=7, height=10, unit="in", res=300);
#ggplot(forplot, aes(x=log2FoldChange, y=-log10(padj), color = Significant)) +
#  geom_point() +
#  ggtitle(paste(trait2stuff[1], " Volcano", sep='')) +
#  scale_color_brewer(palette = "Paired")+
#  coord_fixed(ylim = c(0, 12))
#dev.off()
#
