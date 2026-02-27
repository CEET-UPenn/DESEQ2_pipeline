library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)


#run from within directory containing gct files. will output into same directory 

#use args to get the input files
args = commandArgs(TRUE)

setwd("/project/voightlab_01/lorenzk/CEET/Amatullah/bam")
filename = "Uhrf2_tpm.txt"
filename = "L3mbtl3_tpm.txt"

df <- read.delim(filename, header = FALSE)
df <- separate(df, V1, c("exp", "type1", "type2", "line", "tpm", "gct", "ensg", "transcript"), sep = "[-.]")
df <- rename(df, "TPM" = "V3")
df <- select(df, type1, type2, line, TPM)

#separate by gene
ctrls <- filter(df, type1 == "siControl")
L <- filter(df, type1 == "siL3M")
U <- filter(df, type1 == "siUhrf2")

ctrlL <- filter(ctrls, line %in% L$line )
ctrlU <- filter(ctrls, line %in% U$line)

L <- rbind(L, ctrlL)
U <- rbind(U, ctrlU)

signals_out <- "Uhrf2_TPM.txt"
U <- rename(U, "Uhrf2_TPM" = "TPM")
write.table(U , file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)
signals_out <- "L3mbtl3_TPM.txt"
L  <- rename(L , "L3mbtl3_TPM" = "TPM")
write.table(L , file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)



png("TPM_Uhrf2.png", width=7, height=5, unit="in", res=300);
ggplot(U, aes(x=line, y=TPM, color = type1)) +
  geom_point(size = 3) +
  ylab("Uhrf2 TPM") + xlab("Cell line") +
  labs(title = "Uhrf2 Expression" ) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  facet_wrap(vars(type2))
dev.off()

png("TPM_L3mbtl3.png", width=7, height=5, unit="in", res=300);
ggplot(L, aes(x=line, y=TPM, color = type1)) +
  geom_point(size = 3) +
  ylab("L3mbtl3 TPM") + xlab("Cell line") +
  labs(title = "L3mbtl3 Expression" ) +
  scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  facet_wrap(vars(type2))
dev.off()

