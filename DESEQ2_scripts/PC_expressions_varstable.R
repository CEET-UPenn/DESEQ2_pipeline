library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)


#run from within directory containing gct files. will output into same directory 

#use args to get the input files
args = commandArgs(TRUE)

#testing purposes:
#args <- c("Uhrf2_LPS_ctrl,Uhrf2_LPS_ctrl_readcount.gct.gz","Uhrf2_LPS_si,Uhrf2_LPS_si_readcount.gct.gz")
#args <- c("Uhrf2_PBS_ctrl,Uhrf2_PBS_ctrl_readcount.gct.gz","Uhrf2_PBS_si,Uhrf2_PBS_si_readcount.gct.gz")
#args <- c("L3mbtl3_LPS_ctrl,L3mbtl3_LPS_ctrl_readcount.gct.gz","L3mbtl3_LPS_si,L3mbtl3_LPS_si_readcount.gct.gz")
#args <- c("L3mbtl3_PBS_ctrl,L3mbtl3_PBS_ctrl_readcount.gct.gz","L3mbtl3_PBS_si,L3mbtl3_PBS_si_readcount.gct.gz")
setwd("/project/voightlab_01/lorenzk/CEET/Amatullah/bam")
filename = "all_readcount.gct.gz"

df <- read.delim(filename, header = TRUE, skip=2)
rownames(df) <- df$Name
df[,2] <- NULL
df[,1] <- NULL

#make meta file
meta <- data.frame(cbind(colnames(df)))
colnames(meta) <- "header"
meta <- separate(meta, header, c("exp", "type1", "type2", "line"), sep = "[.]", remove = FALSE)

#import data
DE_data <- DESeqDataSetFromMatrix(countData = df, 
                                  colData = meta,
                                  design = formula(paste("~ line + type1")))

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
plot <- separate(plot, rowname, c("exp", "type1", "type2", "line"), sep = "[.]", remove = FALSE)

png("PC1PC2_all_samples_varstable.png", width=7, height=7, unit="in", res=300);
ggplot(plot, aes(x=PC1, y=PC2, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

png("PC3PC4_all_samples_varstable.png", width=7, height=7, unit="in", res=300);
ggplot(plot, aes(x=PC3, y=PC4, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

#separate by gene
ctrls <- filter(plot, type1 == "siControl")
L <- filter(plot, type1 == "siL3M")
U <- filter(plot, type1 == "siUhrf2")

ctrlL <- filter(ctrls, line %in% L$line )
ctrlU <- filter(ctrls, line %in% U$line)

L <- rbind(L, ctrlL)
U <- rbind(U, ctrlU)


png("PC1PC2_all_samples_varstable_L3M.png", width=7, height=7, unit="in", res=300);
ggplot(L, aes(x=PC1, y=PC2, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

png("PC3PC4_all_samples_varstable_L3M.png", width=7, height=7, unit="in", res=300);
ggplot(L, aes(x=PC3, y=PC4, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

png("PC1PC2_all_samples_varstable_Uhrf2.png", width=7, height=7, unit="in", res=300);
ggplot(U, aes(x=PC1, y=PC2, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()

png("PC3PC4_all_samples_varstable_Uhrf2.png", width=7, height=7, unit="in", res=300);
ggplot(U, aes(x=PC3, y=PC4, color = line, shape = type1)) +
  geom_point(size = 3, aes(fill=line, alpha=type2)) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2")
dev.off()


