library(here)
library(vegan)
library(reshape2)
library(ggplot2)
library(gridExtra)

seqtab <- readRDS(here("../../data/seqtab_clean.rds"))
meta <- readRDS(here("../../data/meta_clean_sum.rds"))
tax <- readRDS(here("../../data/tax_clean.rds"))

seqtab <- seqtab[,colSums(seqtab)>10000]

# Rarefy
seqtab_rar <- t(rrarefy(t(seqtab), min(colSums(seqtab))))

#Create relevant subsets
#Capped mound and bare mineral spruce control root samples
cm_spruce <- meta$SampleID_sum[meta$Species=="Spruce"&meta$Treatment=="Control"&meta$Type=="Roots"&meta$Planting_Position=="Capped_Mound"]
em_spruce <- meta$SampleID_sum[meta$Species=="Spruce"&meta$Treatment=="Control"&meta$Type=="Roots"&meta$Planting_Position=="Exposed_Mineral"]

#Capped mound and bare mineral pine control root samples
cm_pine <- meta$SampleID_sum[meta$Species=="Pine"&meta$Treatment=="Control"&meta$Type=="Roots"&meta$Planting_Position=="Capped_Mound"]
em_pine <- meta$SampleID_sum[meta$Species=="Pine"&meta$Treatment=="Control"&meta$Type=="Roots"&meta$Planting_Position=="Exposed_Mineral"]

rownames(meta) <- meta$SampleID_sum
#Run statistical tests (DIFFERENCES BETWEEN 2017 and 2018 for controls?)
#SPRUCE CAPPED MOUND 
meta_cm_spr <- meta[cm_spruce,]
adonis(vegdist(t(seqtab[,cm_spruce])) ~ meta_cm_spr$Year, strata = meta_cm_spr$Block)

#SPRUCE EXPOSED MINERAL
meta_em_spr <- meta[em_spruce,]
adonis2(vegdist(t(seqtab[,em_spruce])) ~ meta_em_spr$Year, strata = meta_em_spr$Block)

#PINE CAPPED MOUND 
meta_cm_pin <- meta[cm_pine,]
adonis2(vegdist(t(seqtab[,cm_pine])) ~ meta_cm_pin$Year, strata = meta_cm_pin$Block)

#PINE EXPOSED MINERAL
meta_em_pin <- meta[em_pine,]
adonis2(vegdist(t(seqtab[,em_pine])) ~ meta_em_pin$Year, strata = meta_em_pin$Block)


### NOW INCLUDE BOTH PLANTING POSITION AND TEST BY YEAR + PLANTING POSITION

#Create relevant subsets
spruce <- meta$SampleID_sum[meta$Species=="Spruce"&meta$Treatment=="Control"&meta$Type=="Roots"]

pine <- meta$SampleID_sum[meta$Species=="Pine"&meta$Treatment=="Control"&meta$Type=="Roots"]

#Run statistical tests
#SPRUCE 
meta_spr <- meta[spruce,]
adonis2(vegdist(t(seqtab[,spruce])) ~ meta_spr$Year*meta_spr$Planting_Position, strata = meta_spr$Block)
adonis2(t(seqtab_rar[,spruce]) ~ meta_spr$Year*meta_spr$Planting_Position, strata = meta_spr$Block)

anova(betadisper(vegdist(t(seqtab_rar[,spruce])), meta_spr$Planting_Position))

#PINE
meta_pin <- meta[pine,]
adonis2(vegdist(t(seqtab_rar[,pine])) ~ meta_pin$Year*meta_pin$Planting_Position, strata = meta_pin$Block)


####MOST ABUNDANT SOTUs ON ROOTS 2017 + 2018
seqtab_rar_p <- prop.table(seqtab_rar, margin = 2)*100

spruce17 <- spruce[grepl("2017", spruce)]
seqtab_rar_spr17 <- seqtab_rar_p[,spruce17]
spruce18 <- spruce[grepl("2018", spruce)]
seqtab_rar_spr18 <- seqtab_rar_p[,spruce18]

pine17 <- pine[grepl("2017", pine)]
seqtab_rar_pin17 <- seqtab_rar_p[,pine17]
pine18 <- pine[grepl("2018", pine)]
seqtab_rar_pin18 <- seqtab_rar_p[,pine18]

###
##Extract the most abundant ones per species and year
top10_f_spr17 <- rownames(seqtab_rar_spr17)[order(rowMeans(seqtab_rar_spr17), decreasing = TRUE)][1:15]
top10_f_spr18 <- rownames(seqtab_rar_spr18)[order(rowMeans(seqtab_rar_spr18), decreasing = TRUE)][1:15]
top10_f_pin17 <- rownames(seqtab_rar_pin17)[order(rowMeans(seqtab_rar_pin17), decreasing = TRUE)][1:15]
top10_f_pin18 <- rownames(seqtab_rar_pin18)[order(rowMeans(seqtab_rar_pin18), decreasing = TRUE)][1:15]

###Let us try to plot this
seqtab_rar_spr <- seqtab_rar_p[,spruce]
seqtab_rar_pin <- seqtab_rar_p[,pine]

#Extract most abundant for spruce and pine, both years
top10_f_spr <- rownames(seqtab_rar_spr)[order(rowMeans(seqtab_rar_spr), decreasing = TRUE)][1:15]
top10_f_pin <- rownames(seqtab_rar_pin)[order(rowMeans(seqtab_rar_pin), decreasing = TRUE)][1:15]

#Create plot for top 15 SOTUs from Norway spruce control roots
seqtab_rar_spr_top10 <- seqtab_rar_spr[top10_f_spr,]
spr_top_plot <- melt(seqtab_rar_spr_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
spr_top_plot$SOTU <- factor(spr_top_plot$SOTU, levels = rev(top10_f_spr))
spr_top_plot$Year <- ifelse(grepl("2017", spr_top_plot$Sample), "2017", "2018")
spr_top_plot$SOTU2 <- paste(spr_top_plot$SOTU,
                                tax$Species[match(as.character(spr_top_plot$SOTU), tax$OTU_ID)],
                                sep = ".")
spr_top_plot$SOTU2 <- factor(spr_top_plot$SOTU2, levels = rev(unique(as.character(spr_top_plot$SOTU2))))
pf1 <- ggplot(spr_top_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Year))+
  ggtitle("Spruce roots control")+
  theme_light()+
  scale_color_manual(values = c("#4886D9", "#50592E"))

#Create plot for top 15 SOTUs from Scots pine control roots
seqtab_rar_pin_top10 <- seqtab_rar_pin[top10_f_pin,]
pin_top_plot <- melt(seqtab_rar_pin_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
pin_top_plot$SOTU <- factor(pin_top_plot$SOTU, levels = rev(top10_f_pin))
pin_top_plot$Year <- ifelse(grepl("2017", pin_top_plot$Sample), "2017", "2018")
pin_top_plot$SOTU2 <- paste(pin_top_plot$SOTU,
                                tax$Species[match(as.character(pin_top_plot$SOTU), tax$OTU_ID)],
                                sep = ".")
pin_top_plot$SOTU2 <- factor(pin_top_plot$SOTU2, levels = rev(unique(as.character(pin_top_plot$SOTU2))))
pf2 <- ggplot(pin_top_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Year))+
  ggtitle("Pine roots control")+
  theme_light()+
  scale_color_manual(values = c("#4886D9", "#50592E"))  

grid.arrange(pf1, pf2, nrow = 1)
ggsave(here("Figure_4.pdf"), grid.arrange(pf1, pf2, nrow = 1), width = 10, height = 5)

#Legend ECM
funguild <- readRDS(here("../../FUNGuild/FUNGuild_UPDATED.rds"))

spr_top_fg <- funguild$ECM[match(as.character(unique(spr_top_plot$SOTU)), funguild$OTU_ID)]

pin_top_fg 

#############################
### Overlap with nursery fungi
nursery <- meta$SampleID_sum[meta$Treatment=="Nursery" & !is.na(meta$Treatment)]

seqtab_n <- seqtab_rar_p[,nursery]
meta_n <- droplevels(meta[nursery,])

seqtab_n_spr <- seqtab_n[,grepl("Spruce", meta_n$Species)]
seqtab_n_pin <- seqtab_n[,grepl("Pine", meta_n$Species)]

##Extract the most abundant ones per species
top10_f_n_spr <- rownames(seqtab_n_spr)[order(rowMeans(seqtab_n_spr), decreasing = TRUE)][1:15]
top10_f_n_pin <- rownames(seqtab_n_pin)[order(rowMeans(seqtab_n_pin), decreasing = TRUE)][1:15]

##Extract all that present with >0.05% average rel ab
nur005_spr <- rownames(seqtab_n_spr)[rowMeans(seqtab_n_spr)>0.01]
nur005_pin <- rownames(seqtab_n_pin)[rowMeans(seqtab_n_pin)>0.01]

##SPRUCE
#Overlap MA nursery - MA control roots 2017+2018
intersect(top10_f_n_spr, top10_f_spr)
#Overlap nursery>0.05% - MA control roots 2017+2018
setdiff(intersect(nur005_spr, top10_f_spr), intersect(top10_f_n_spr, top10_f_spr))

##PINE
#Overlap MA nursery - MA control roots 2017+2018
intersect(top10_f_n_pin, top10_f_pin)
#Overlap nursery>0.05% - MA control roots 2017+2018
setdiff(intersect(nur005_pin, top10_f_pin), intersect(top10_f_n_pin, top10_f_pin))

###PERCENTAGE OF ECM FUNGI ON SEEDLING ROOTS
tax_fg <- readRDS(here("../../FUNGuild/FUNGuild_UPDATED.rds"))
ecm <- tax_fg$OTU_ID[tax_fg$ECM=="ECM"|tax_fg$ECM=="Putative_ECM"]

#Calculate mean percentage of ECM/putative ECM in 2017 control roots
seqtab_2017 <- cbind(seqtab_rar_pin17, seqtab_rar_spr17)
mean(colSums(seqtab_2017[ecm,])/colSums(seqtab_2017))

#Calculate mean percentage of ECM/putative ECM in 2018 control roots
seqtab_2018 <- cbind(seqtab_rar_pin18, seqtab_rar_spr18)
mean(colSums(seqtab_2018[ecm,])/colSums(seqtab_2018))

#Calculate mean percentage of ECM/putative ECM in Norway spruce control roots
mean(colSums(seqtab_rar_spr[ecm,])/colSums(seqtab_rar_spr))

#Calculate mean percentage of ECM/putative ECM in Scots pine control roots
mean(colSums(seqtab_rar_pin[ecm,])/colSums(seqtab_rar_pin))

wilcox.test(colSums(seqtab_rar_spr[ecm,])/colSums(seqtab_rar_spr), colSums(seqtab_rar_pin[ecm,])/colSums(seqtab_rar_pin))
#No signficant difference in % ECM on control roots between spruce and pine

#Number of SOTUs with more than 0.1% rel ab in any sample
seqtab_p <- prop.table(seqtab, margin = 2)*100
seqtab_p[seqtab_p<0.1] <- 0 
seqtab_p2 <- seqtab_p[rowSums(seqtab_p)>0,]
ecm2 <- ecm[ecm%in%rownames(seqtab_p2)]

#### WHICH ONES ARE ALSO PRESENT IN SITE SOILS?
top15_all <- unique(c(top10_f_pin, top10_f_spr))

control_soil_block <- meta$SampleID_sum[grepl("Soil_Block", meta$Type) 
                                        & meta$Treatment=="Control"
                                        & !is.na(meta$Treatment)]
edges <- meta$SampleID_sum[grepl("Edge", meta$Block)][-c(11,12)]

seqtab_csb_e <- seqtab_rar[,c(control_soil_block, edges)]
seqtab_csb_e <- seqtab_csb_e[rowSums(seqtab_csb_e)>10,]

site_fungi <- rownames(seqtab_csb_e)

length(top15_all)
setdiff(top15_all, site_fungi)
length(intersect(top15_all, site_fungi))

