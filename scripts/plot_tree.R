#!/usr/bin/env Rscript
library(ggtree)
library(dplyr)
library(tidyverse)
library(treedataverse)
library(ape)
library(treeio)
library(tidytree)
library(ggtreeExtra)
library(ggstar)
library(RColorBrewer)
tree1<-read.tree("strain_tree/RmucY2510_v1.All.SNP.fasttree_reroot.tre") %>% drop.tip(c("Ref_NRRLY2510"))

#tree<-read.tree("strain_tree/Rmuc_v7.All.SNP.poppr.nj.tre")
meta <- read_csv("strain_metdata.csv",col_names=TRUE) %>% 
  select(Strain,SimpleEnv,culture_collection) %>% 
  filter(!is.na(SimpleEnv))

CNV <- read_tsv("plots/strain_median_coverage.txt",col_names=TRUE, col_types = cols(.default = "?", chrlen = "i", CHR = "i"))
strains = colnames(as.data.frame(CNV %>% select(-c("CHR", "chrlen"))))
chroms <- CNV$CHR
tCNV <- t(as.data.frame(CNV %>% select(-c("CHR", "chrlen"))))
colnames(tCNV) <- chroms

tree <- tree1 %>% left_join(meta,by = c("label"="Strain"))
unique(meta$SimpleEnv)
unique(meta$culture_collection)
p <- ggtree(tree, layout="circular", size=0.3, branch.length = "none")
p

colourCount = length(unique(meta$SimpleEnv))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p2<-p + geom_fruit(
  geom=geom_star,
  mapping=aes(y=label, fill=SimpleEnv, starshape=culture_collection),
  starstroke=0.2
) +   scale_fill_manual(values = getPalette(colourCount))

ggsave("strain_tree/cladogram.pdf",p2,width=10,height=10)


p <- ggtree(tree, size=0.5) + geom_tiplab(hjust=-0.5,size=1.5) +
  geom_tippoint(aes(color=SimpleEnv), size=1) +
  scale_fill_manual(values = getPalette(colourCount))
p2 <- gheatmap(p,tCNV,colnames_offset_y=-0.5,font.size = 2,offset=0.5) + scale_fill_continuous(type = "viridis",name="Chrom\nCNV")
ggsave("strain_tree/strain_tree_CNV.pdf",p2,width=10,height=20)

