library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)


#run from within directory containing gct files. will output into same directory 

setwd("/project/voightlab_01/lorenzk/CEET/Amatullah/bam")
filename = "EX02_Uhrf2_tpm.txt"


df <- read.delim(filename, header = FALSE)
df <- separate(df, V1, c("exp", "OV", "type1", "type2", "rep","repcount", "tpm", "gct", "ensg", "transcript"), sep = "[-.]")
df <- rename(df, "TPM" = "V3")
df <- select(df, type1, type2, repcount, TPM)


signals_out <- "Uhrf2_TPM.txt"
df <- rename(df, "Uhrf2_TPM" = "TPM")
write.table(df , file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)

png("EX02_TPM_Uhrf2.png", width=7, height=5, unit="in", res=300);
ggplot(df, aes(x=repcount, y=Uhrf2_TPM, color = type1)) +
  geom_point(size = 3) +
  ylab("Uhrf2 TPM") + xlab("Cell line") +
  labs(title = "Uhrf2 Expression" ) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  facet_wrap(vars(type2))
dev.off()

#boxplot

png("EX02_TPM_Uhrf2_boxplot.png", width=7, height=5, unit="in", res=300);
ggplot(df, aes(x=type1, y=Uhrf2_TPM)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .5)+  
  xlab("Condition") +
  facet_wrap(vars(type2))
dev.off()
