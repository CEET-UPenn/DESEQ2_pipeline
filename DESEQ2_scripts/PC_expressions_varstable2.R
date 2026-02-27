library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)


#run from within directory containing gct files. will output into same directory 


setwd("/project/voightlab_01/lorenzk/CEET/Amatullah/bam")
filename = "EX02_all_readcount.gct.gz"

df <- read.delim(filename, header = TRUE, skip=2)
rownames(df) <- df$Name
df[,2] <- NULL
df[,1] <- NULL

#make meta file
meta <- data.frame(cbind(colnames(df)))
colnames(meta) <- "header"
meta <- separate(meta, header, c("exp", "OV", "type1", "type2", "rep", "repcount"), sep = "[.]", remove = FALSE)

#import data
DE_data <- DESeqDataSetFromMatrix(countData = df, 
                                  colData = meta,
                                  design = formula(paste("~ repcount + type1")))

#prefilter
keep <- rowSums(counts(DE_data)) >= 10
DE_data <- DE_data[keep,]

#run varianceStabilizingTransformation 
var_transformed <- varianceStabilizingTransformation(DE_data, blind = TRUE, fitType = "parametric")
exp_data <- data.frame(assay(var_transformed))


results <- prcomp(exp_data, scale = TRUE)
head(results$rotation)
head(results$x)

plot <- data.frame(results$rotation)
plot <- rownames_to_column(plot)
plot <- separate(plot, rowname, c("exp", "OV", "type1", "type2", "rep", "repcount"), sep = "[.]", remove = FALSE)

png("EX02_PC1PC2_all_samples_varstable.png", width=7, height=7, unit="in", res=300);
ggplot(plot, aes(x=PC1, y=PC2, color = repcount, shape = type1)) +
  geom_point(size = 3, aes(fill=repcount, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

png("EX02_PC3PC4_all_samples_varstable.png", width=7, height=7, unit="in", res=300);
ggplot(plot, aes(x=PC3, y=PC4, color = repcount, shape = type1)) +
  geom_point(size = 3, aes(fill=repcount, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

