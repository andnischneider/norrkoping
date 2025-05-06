library(here)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(vegan)
library(grid)
library(gridExtra)
library(nlme)
library(phyloseq)
library(FUNGuildR)
library(lefser)
library(rcartocolor)
library(ggsignif)


ggformat <- theme_minimal()+theme(axis.text.x = element_text(angle = 45))
melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("OTU_ID", "SampleID_sum"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}



# Data Import and preparation

seqtab <- readRDS(here("../data/seqtab_clean.rds"))
seqtab <- seqtab[,colSums(seqtab)>1000]

tax <- readRDS(here("../data/tax_clean.rds"))
meta <- readRDS(here("../data/meta_clean_sum.rds"))
meta <- meta[meta$SampleID_sum%in%colnames(seqtab),]
rownames(meta) <- meta$SampleID_sum
rownames(tax) <- tax$OTU_ID
seqtab <- seqtab[,meta$SampleID_sum]
meta$Year[meta$Treatment=="Nursery"] <- "Early_2017"
meta$Year <- factor(meta$Year, levels = c("Early_2017", "2017", "2018"))
meta$Treatment <- factor(meta$Treatment, levels = c("Nursery", "Control", "Osmocote",
                                                    "Argrow", "Control_Site"))



#Soil C/N content data
soil_all_r <- readRDS(here("../data/Soil_data_clean.rds"))

#First Analyses
ggqqplot(soil_all_r$omegaN_Perc)
ggqqplot(soil_all_r$omegaC_Perc)
ggqqplot(soil_all_r$C_N_ratio)

shapiro.test(soil_all_r$omegaN_Perc)
shapiro.test(soil_all_r$C_N_ratio)

#Definitely not normally distributed, so we use PERMANOVA to test significant differences.
rownames(meta) <- meta$SampleID_sum
soil_all_r <- as.data.frame(soil_all_r)
meta_soil <-readRDS(here("../data/meta_soil.rds"))
rownames(meta_soil) <- meta_soil$SampleID_sum
rownames(soil_all_r) <- soil_all_r$SampleID_sum
soil_all_r <- soil_all_r[rownames(meta_soil),]

###################### Edge samples
edges <- meta_soil$SampleID_sum[grepl("Edge", meta_soil$Block)]
soil_all_r_e <- soil_all_r[edges,]
meta_soil_e <- meta_soil[edges,]

shapiro.test(soil_all_r_e$omegaN_Perc)
shapiro.test(soil_all_r_e$C_N_ratio)
#still better to use adonis2
adonis2(soil_all_r_e$omegaC_Perc ~ meta_soil_e$Type*meta_soil_e$Planting_Position)
adonis2(soil_all_r_e$omegaN_Perc ~ meta_soil_e$Type*meta_soil_e$Planting_Position)
adonis2(soil_all_r_e$C_N_ratio ~ meta_soil_e$Type*meta_soil_e$Planting_Position)

soil_all_r_e2 <- merge.data.frame(soil_all_r_e, meta_soil_e, by = "SampleID_sum")

p1=ggplot(soil_all_r_e2, aes(x = Planting_Position, y = omegaC_Perc))+
  geom_boxplot(aes(fill=Planting_Position))+
  facet_wrap(~Type)+
  ggformat+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))
p2=ggplot(soil_all_r_e2, aes(x = Planting_Position, y = omegaN_Perc))+
  geom_boxplot(aes(fill=Planting_Position))+
  facet_wrap(~Type)+
  ggformat+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))
p3=ggplot(soil_all_r_e2, aes(x = Planting_Position, y = C_N_ratio))+
  geom_boxplot(aes(fill=Planting_Position))+
  facet_wrap(~Type)+
  ggformat+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))

grid.arrange(p1, p2, p3, nrow = 1)

ggsave(here("FigureS1_CN.png"), plot = grid.arrange(p1, p2, p3, nrow = 1), width = 5, height = 5)

soil_all_r_e2$Planting_Position <- factor(soil_all_r_e2$Planting_Position, levels = c("North", "East", "South", "West"))

(p_p1=ggplot(soil_all_r_e2, aes(y = omegaC_Perc, x = Type))+
  geom_boxplot(aes(fill=Type), outlier.shape = NA)+
    geom_jitter(aes(col = Planting_Position), size = 2, alpha = 0.6, width = 0.15)+
    scale_x_discrete(labels=c("Soil_CC"="Clear-cut", "Soil_Forest"="Forest"), name=element_blank())+
    scale_y_continuous(name="Carbon mass fraction (%)")+
    scale_color_discrete(name="Location")+
    scale_fill_manual(values = alpha(c("sienna4", "springgreen4"), 0.3))+
    coord_cartesian(ylim=c(0,55))+
    theme_light()+
    theme(legend.position = "none"))

(p_p2=ggplot(soil_all_r_e2, aes(y = omegaN_Perc, x = Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA)+
    geom_jitter(aes(col = Planting_Position), size = 2, alpha = 0.6, width = 0.15)+
    scale_x_discrete(labels=c("Soil_CC"="Clear-cut", "Soil_Forest"="Forest"), name=element_blank())+
    scale_y_continuous(name="Nitrogen mass fraction (%)")+
    scale_color_discrete(name="Location")+
    scale_fill_manual(values = alpha(c("sienna4", "springgreen4"), 0.3))+
    coord_cartesian(ylim=c(0,2.5))+
    theme_light()+
    theme(legend.position = "none"))

(p_p3=ggplot(soil_all_r_e2, aes(y = C_N_ratio, x = Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA)+
    geom_jitter(aes(col = Planting_Position), size = 2, alpha = 0.6, width = 0.15)+
    scale_x_discrete(labels=c("Soil_CC"="Clear-cut", "Soil_Forest"="Forest"), name=element_blank())+
    scale_y_continuous(name="Carbon:Nitrogen ratio")+
    scale_color_discrete(name="Location")+
    scale_fill_manual(values = alpha(c("sienna4", "springgreen4"), 0.3))+
    coord_cartesian(ylim=c(15,35))+
    theme_light()+
    theme(legend.position = "none"))

grid.arrange(p_p1, p_p2, p_p3, nrow = 1)
ggsave(here("FigureS1_CN_NEW0322.pdf"), plot = grid.arrange(p_p1, p_p2, p_p3, nrow = 1), 
       width = 130, 
       height = 100, 
       units = "mm", 
       dpi = 300)


#extract legend
g_legend <- function(a.gplot){
 tmp <- ggplot_gtable(ggplot_build(a.gplot))
 leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
 legend <- tmp$grobs[[leg]]
 return(legend)
}
(p_p1_wl=ggplot(soil_all_r_e2, aes(y = omegaC_Perc, x = Type))+
    geom_boxplot(aes(fill=Type))+
    geom_jitter(aes(col = Planting_Position), size = 3, alpha = 0.6, width = 0.15)+
    scale_x_discrete(labels=c("Soil_CC"="Clear-cut", "Soil_Forest"="Surrounding Forest"), name=element_blank())+
    scale_y_continuous(name="Carbon mass fraction (%)")+
    scale_color_discrete(name="Location")+
    scale_fill_manual(values = alpha(c("sienna4", "springgreen4"), 0.3))+
    theme_minimal())

pdf(here("FigureS1C_legend.pdf"), height = 3, width = 1.5)
grid.newpage()
grid.draw(g_legend(p_p1_wl))
dev.off()

################## Control block soil
control_soil_block <- meta_soil$SampleID_sum[grepl("Soil_Block", meta_soil$Type) 
                                             & meta_soil$Treatment=="Control"
                                             & !is.na(meta_soil$Treatment)]

meta_csb <- meta_soil[control_soil_block,]
soil_all_r_csb <- soil_all_r[control_soil_block,]

#Difference between blocks?
###POINTLESS, since in this part there is only one replicate per block
#adonis2(soil_all_r_csb$omegaC_Perc ~ meta_csb$Block)
#adonis2(soil_all_r_csb$omegaN_Perc ~ meta_csb$Block)
#adonis2(soil_all_r_csb$C_N_ratio ~ meta_csb$Block)
#Plot:
cols_pp <- carto_pal(2, "Earth")
soil_all_r_csb2 <- merge.data.frame(soil_all_r_csb, meta_csb)
(pb_p1 <- ggplot(soil_all_r_csb2, aes(x = Planting_Position, y = omegaC_Perc))+
    geom_boxplot(aes(fill=Planting_Position), outlier.shape = NA)+
    geom_jitter(aes(col=Species), width=0.15, size = 2, alpha = 0.6)+
    scale_fill_manual(values = alpha(cols_pp, 0.8))+
    scale_color_manual(values=c("#C3CC33", "#17BA35"))+
    scale_x_discrete(labels=c("Capped_Mound"="Capped\nMound", "Exposed_Mineral"="Mineral\nSoil"), name=element_blank())+
    scale_y_continuous(name="Carbon mass fraction (%)")+
    coord_cartesian(ylim=c(0,55))+
    theme_light()+
    theme(axis.title.x = element_blank(),
          legend.position = "none"))

(pb_p2 <- ggplot(soil_all_r_csb2, aes(x = Planting_Position, y = omegaN_Perc))+
    geom_boxplot(aes(fill=Planting_Position), outlier.shape = NA)+
    geom_jitter(aes(col=Species), width=0.15, size = 2, alpha = 0.6)+
    scale_fill_manual(values = alpha(cols_pp, 0.8))+
    scale_color_manual(values=c("#C3CC33", "#17BA35"))+
    scale_x_discrete(labels=c("Capped_Mound"="Capped\nMound", "Exposed_Mineral"="Mineral\nSoil"), name=element_blank())+
    scale_y_continuous(name="Nitrogen mass fraction (%)")+
    coord_cartesian(ylim=c(0,2.5))+
    theme_light()+
    theme(axis.title.x = element_blank(),
          legend.position = "none"))

(pb_p3 <- ggplot(soil_all_r_csb2, aes(x = Planting_Position, y = C_N_ratio))+
    geom_boxplot(aes(fill=Planting_Position), outlier.shape = NA)+
    geom_jitter(aes(col=Species), width=0.15, size = 2, alpha = 0.6)+
    scale_fill_manual(values = alpha(cols_pp, 0.8))+
    scale_color_manual(values=c("#C3CC33", "#17BA35"))+
    scale_x_discrete(labels=c("Capped_Mound"="Capped\nMound", "Exposed_Mineral"="Mineral\nSoil"), name=element_blank())+
    scale_y_continuous(name="Carbon:Nitrogen ratio")+
    coord_cartesian(ylim=c(15,35))+
    theme_light()+
    theme(axis.title.x = element_blank(),
          legend.position = "none"))

grid.arrange(pb_p1, pb_p2, pb_p3, nrow = 1)
ggsave(here("Figure1_CN_blocks.pdf"), plot = grid.arrange(pb_p1, pb_p2, pb_p3, nrow = 1), 
       width = 130, 
       height = 100, 
       units = "mm", 
       dpi = 300)

(pb_p1_wl <- ggplot(soil_all_r_csb2, aes(x = Planting_Position, y = omegaC_Perc))+
    geom_boxplot(aes(fill=Planting_Position))+
    geom_jitter(aes(col=Species), width=0.15, size = 3, alpha = 0.6)+
    scale_fill_manual(values = alpha(cols_pp, 0.8))+
    scale_color_manual(values=c("#C3CC33", "#17BA35"))+
    scale_x_discrete(labels=c("Capped_Mound"="Capped\nMound", "Exposed_Mineral"="Mineral\nSoil"), name=element_blank())+
    scale_y_continuous(name="Carbon mass fraction (%)")+
    coord_cartesian(ylim=c(0,55))+
    theme_minimal()+
    theme(axis.title.x = element_blank()))

pdf(here("FigureS1B_legend.pdf"), height = 3, width = 1.5)
grid.newpage()
grid.draw(g_legend(pb_p1_wl))
dev.off()

#Use blocks as strata and test C/N content 
set.seed(555)
adonis2(soil_all_r_csb$omegaC_Perc ~ meta_csb$Species*meta_csb$Planting_Position*meta_csb$Year, strata = meta_csb$Block)
adonis2(soil_all_r_csb$omegaN_Perc ~ meta_csb$Species*meta_csb$Planting_Position*meta_csb$Year, strata = meta_csb$Block)
adonis2(soil_all_r_csb$C_N_ratio ~ meta_csb$Species*meta_csb$Planting_Position*meta_csb$Year, strata = meta_csb$Block)

#########################################################################################
##### Fungal community sequencing data
#Rarefy to same order of magnitude
seqtab_rar <- t(rrarefy(t(seqtab), min(colSums(seqtab))))
div_all <- diversity(seqtab, MARGIN = 2)
div_all <- data.frame(meta, Shannon=div_all)

### Are 2017 and 2018 block soils different?
meta_csb <- meta[control_soil_block,]
ps_csb <- phyloseq(sample_data(meta_csb),
                   otu_table(seqtab_rar[,control_soil_block], taxa_are_rows = TRUE))
ord_csb <- ordinate(ps_csb, "MDS", "bray")
plot_ordination(ps_csb, ord_csb, type = "samples", color = "Block", shape = "Year")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(ord_csb$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(ord_csb$values$Relative_eig[2]*100, digits = 1), "%]"))+
  #scale_color_carto_d(palette = "Earth", direction = 1)+
  theme_minimal()

ggsave(here("Figure1_PCoA.jpg"), width = 5, height = 4)

adonis2(vegdist(t(seqtab_rar[,control_soil_block])) ~ meta_csb$Year*meta_csb$Planting_Position*meta_csb$Species, strata = meta_csb$Block)
##NOTHING SIGNIFICANT.

shannon_model_b_c <- lme(Shannon~Year*Planting_Position*Species, random = ~1|Block, data = div_all[control_soil_block,], method = "ML")
anova(shannon_model_b_c)
## SHANNON NO DIFFERENCES

#Try with summarised blocks
meta_csb$Block2 <- mgsub::mgsub(meta_csb$Block,
                                unique(meta_csb$Block),
                                c("Block_3", "Block_3", "Block_3", "Block_3", "Block_1",
                                  "Block_1", "Block_1", "Block_4", "Block_5", "Block_2",
                                  "Block_2", "Block_2", "Block_4", "Block_5", "Block_2"))
adonis2(vegdist(t(seqtab_rar[,control_soil_block])) ~ meta_csb$Year*meta_csb$Planting_Position*meta_csb$Species, strata = meta_csb$Block2)
div_csb <- merge.data.frame(div_all[control_soil_block,], meta_csb)
shannon_model_b_c2 <- lme(Shannon~Year*Planting_Position*Species, random = ~1|Block2, data = div_csb, method = "ML")
anova(shannon_model_b_c2)

##Most abundant genera
seqtab_rar_c <- seqtab_rar[,control_soil_block]
seqtab_rar_c2 <- t(rrarefy(t(seqtab_rar_c), min(colSums(seqtab_rar_c))))
seqtab_melted_c <- melt_filt2(seqtab_rar_c2)
seqtab_melted_c <- merge.data.frame(seqtab_melted_c, meta[control_soil_block,])
seqtab_melted_c <- merge.data.frame(seqtab_melted_c, tax, by = "OTU_ID")
gen_mat_c <- acast(seqtab_melted_c, Genus~SampleID_sum, value.var = "Count", fun.aggregate = sum)

gen_20_c <- names(sort(rowSums(gen_mat_c), decreasing = TRUE)[1:20])
gen_20_c
##Most abundant ECM genera


##################################################################
######### CLEARCUT EDGE SAMPLES
edges <- meta$SampleID_sum[grepl("Edge", meta$Block)]

meta_e <- droplevels(meta[edges,])
#meta_e$Block <- paste(meta_e$Planting_Position, meta_e$Type, meta_e$Block, sep = ".")

ps_edges <- phyloseq(sample_data(meta_e),
                     otu_table(seqtab_rar[,edges], taxa_are_rows = TRUE))

ord_edges <- ordinate(ps_edges, "MDS", "bray")

plot_ordination(ps_edges, ord_edges, type = "samples", color = "Type", shape = "Planting_Position")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(ord_edges$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(ord_edges$values$Relative_eig[2]*100, digits = 1), "%]"))+
  theme(legend.position = "none")+
  theme_light()+
  scale_color_manual(values = c("sienna4", "springgreen4"))+
  scale_shape_discrete(name="Location")+
  ggtitle("Clearcut Edge Samples")
ggsave(here("FigureS2_PCoA_edges.pdf"), 
       width = 130, 
       height = 100, 
       units = "mm", 
       dpi = 300)

adonis2(dist(t(seqtab_rar[,edges])) ~ meta_e$Type*meta_e$Planting_Position)

div_all_e <- div_all[edges,]
#div_all_e$Block2 <- paste(div_all_e$)
adonis2(div_all_e$Shannon~Type*Planting_Position, data = div_all_e)

ggplot(div_all_e, aes(y = Shannon, x = Type))+
  geom_boxplot(aes(fill = Type))+
  geom_jitter(aes(col = Planting_Position), size = 2, alpha = 0.6, width = 0.15)+
  geom_signif(
    comparisons = list(c("Soil_CC", "Soil_Forest")),
    map_signif_level = TRUE
  )+
  scale_fill_manual(values = alpha(c("sienna4", "springgreen4"), 0.3))+
  scale_y_continuous(name="Shannon Diversity Index")+
  scale_x_discrete(labels=c("Soil_CC"="Clear-cut", "Soil_Forest"="Forest"), name=element_blank())+
  scale_color_discrete(name="Location")+
  theme_light()
ggsave(here("FigureS2A_ShannonEdges.pdf"), 
       width = 80, 
       height = 100, 
       units = "mm", 
       dpi = 300)


#
edges_cc <- edges[grepl("CC", edges)]
meta_cc_ctrl <- meta[c(edges_cc, control_soil_block),]
meta_cc_ctrl$Type <- as.factor(meta_cc_ctrl$Type)
meta_cc_ctrl$Planting_Position <- as.factor(meta_cc_ctrl$Planting_Position)
seqtab_rar_cc_ctrl <- seqtab_rar[,c(edges_cc, control_soil_block)]

ps_cc_ctrl <- phyloseq(sample_data(meta_cc_ctrl),
                       otu_table(seqtab_rar_cc_ctrl, taxa_are_rows = TRUE))
ord_cc_ctrl <- ordinate(ps_cc_ctrl, "MDS", "bray")

plot_ordination(ps_cc_ctrl, ord_cc_ctrl, type = "samples", color = "Species", shape = "Planting_Position")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(ord_edges$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(ord_edges$values$Relative_eig[2]*100, digits = 1), "%]"))+
  theme(legend.position = "none")+
  theme_minimal()+
  scale_color_manual(values = c("#C3CC33", "#17BA35"))+
  ggtitle("Clearcut Edge Samples + Ctrl Block samples")

adonis2(dist(t(seqtab_rar_cc_ctrl))~meta_cc_ctrl$Type)

sort(rowSums(seqtab_rar_cc_ctrl), decreasing = TRUE)[1:10]



####FunGuild
# fung <- get_funguild_db()
# saveRDS(fung, here("funguild.rds"))
# fung <- readRDS(here("funguild.rds"))
# tax2 <- tax
# tax2$Species <- gsub("_", " ", tax2$Species)
# tax2$Taxonomy <- paste(tax2$Kingdom, tax2$Phylum, tax2$Class, tax2$Order, tax2$Family, tax2$Genus, tax2$Species, sep = ";")
# tax2 <- tax2[,c(1,9)]
tax2_guilds <- readRDS(here("../FUNGuild/FUNGuild_UPDATED.rds"))


#Check how many ECM fungi in forest vs cc
seqtab_rar_e <- seqtab_rar[,edges]

tax2_guilds_ecto <- tax2_guilds[tax2_guilds$ECM=="ECM"|tax2_guilds$ECM=="Putative_ECM",]
seqtab_rar_e_ecm <- seqtab_rar_e[tax2_guilds_ecto$OTU_ID,]

#Calculate percentages
#Percentage of ECM fungi and putative ECM fungi
# perc_ecm <- rbind(data.frame(Percent_ECM=colSums(seqtab_rar_e_ecm[,grepl("Soil_CC", colnames(seqtab_rar_e_ecm))])/colSums(seqtab_rar_e[,grepl("Soil_CC", colnames(seqtab_rar_e))])),
#                   data.frame(Percent_ECM=colSums(seqtab_rar_e_ecm[,grepl("Soil_Forest", colnames(seqtab_rar_e_ecm))])/colSums(seqtab_rar_e[,grepl("Soil_Forest", colnames(seqtab_rar_e))])))
perc_ecm <- data.frame(Percent_ECM=colSums(seqtab_rar_e_ecm)/colSums(seqtab_rar_e))
perc_ecm$SampleID_sum <- rownames(perc_ecm)
perc_ecm2 <- merge.data.frame(perc_ecm, meta)
ggplot(perc_ecm2, aes(x = Planting_Position, y = Percent_ECM, fill = Planting_Position))+
  geom_boxplot()+
  facet_wrap(~Type)+
  ggtitle("ECM and Putative ECM")+
  theme_minimal()

ggplot(perc_ecm2, aes(y = Percent_ECM, fill = Type))+
  geom_boxplot()+
  ggtitle("ECM and Putative ECM")+
  scale_fill_manual(values = c("sienna4", "springgreen4"))+
  theme_minimal()

aggregate(perc_ecm2$Percent_ECM, list(perc_ecm2$Type), mean)

#Percentage of Ericoid Mycorrhiza
tax2_guilds_eric <- tax2_guilds[grepl("Ericoid", tax2_guilds$guild),]
seqtab_rar_e_eric <- seqtab_rar_e[tax2_guilds_eric$OTU_ID,]

#Calculate percentages
#Percentage of mycorrhizal fungi
perc_eric <- data.frame(Percent_Eric=colSums(seqtab_rar_e_eric)/colSums(seqtab_rar_e))
perc_eric$SampleID_sum <- rownames(perc_eric)
perc_eric2 <- merge.data.frame(perc_eric, meta)
ggplot(perc_eric2, aes(x = Planting_Position, y = Percent_Eric, fill = Planting_Position))+
  geom_boxplot()+
  facet_wrap(~Type)+
  ggtitle("Percentage of all putative ericoid mycorrhiza")+
  theme_minimal()

#Percentage of any Mycorrhiza
tax2_guilds_myco <- tax2_guilds[grepl("ycorrhiz", tax2_guilds$guild),]
seqtab_rar_e_myco <- seqtab_rar_e[tax2_guilds_myco$OTU_ID,]

#Calculate percentages
#Percentage of mycorrhizal fungi
perc_myco <- data.frame(Percent_Myco=colSums(seqtab_rar_e_myco)/colSums(seqtab_rar_e))
perc_myco$SampleID_sum <- rownames(perc_myco)
perc_myco2 <- merge.data.frame(perc_myco, meta)
ggplot(perc_myco2, aes(x = Planting_Position, y = Percent_Myco, fill = Planting_Position))+
  geom_boxplot()+
  facet_wrap(~Type)+
  ggtitle("Percentage of all mycorrhiza and putative mycorrhiza")+
  theme_minimal()

ggplot(perc_myco2, aes(y = Percent_Myco, fill = Type))+
  geom_boxplot()+
  ggtitle("Percentage of all mycorrhiza and putative mycorrhiza")+
  scale_fill_manual(values = c("sienna4", "springgreen4"))+
  theme_minimal()

aggregate(perc_myco2$Percent_Myco, list(perc_myco2$Type), mean)

#Most abundant genera
seqtab_rar_e2 <- t(rrarefy(t(seqtab_rar_e), min(colSums(seqtab_rar_e))))
seqtab_melted_e <- melt_filt2(seqtab_rar_e2)
seqtab_melted_e <- merge.data.frame(seqtab_melted_e, meta_e)
seqtab_melted_e <- merge.data.frame(seqtab_melted_e, tax, by = "OTU_ID")
gen_mat_e <- acast(seqtab_melted_e, Genus~SampleID_sum, value.var = "Count", fun.aggregate = sum)

#Most abundant genera in forest and CC edge samples
sort(rowSums(gen_mat_e[,grepl("Soil_CC", colnames(gen_mat_e))]), decreasing = TRUE)[1:20]
sort(rowSums(gen_mat_e[,grepl("Soil_Forest", colnames(gen_mat_e))]), decreasing = TRUE)[1:20]


meta_e$Type <- as.factor(meta_e$Type)


###TAXONOMIC BARPLOT...
seqtab_melted_e2 <- seqtab_melted_e
seqtab_melted_e2 <- seqtab_melted_e2[order(seqtab_melted_e2$SampleID_sum),]
seqtab_melted_e2$SampleID_sum <- factor(seqtab_melted_e2$SampleID_sum, levels = sort(unique(as.character(seqtab_melted_e2$SampleID_sum))))

#50 most abundant genera
gen_20 <- names(sort(rowSums(gen_mat_e), decreasing = TRUE)[1:20])
gen_50 <- names(sort(rowSums(gen_mat_e), decreasing = TRUE)[1:50])

seqtab_melted_e2_20 <- seqtab_melted_e2
seqtab_melted_e2_20$Genus[!seqtab_melted_e2_20$Genus%in%gen_20] <- "Others"
seqtab_melted_e2_20$Genus <- factor(seqtab_melted_e2_20$Genus, levels = c(gen_20, "Others"))

ggplot(seqtab_melted_e2_20, aes(x = SampleID_sum, y = Count, fill = Genus))+
  geom_bar(stat = "identity", position = "fill")+
  scale_y_continuous(labels = scales::percent)+
  ggtitle("Barplot, Genus level")+
  theme_minimal()+
  scale_fill_manual(values = c("#28ec9a",
                               "#6b2698",
                               "#64a420",
                               "#cc8eff",
                               "#29a43b",
                               "#ff73be",
                               "#00c686",
                               "#600059",
                               "#a8e48a",
                               "#7892ff",
                               "#b58400",
                               "#ff9cca",
                               "#486d00",
                               "#8f002b",
                               "#344a00",
                               "#ff9388",
                               "#725600",
                               "#870007",
                               "#ffb47e",
                               "#a44f00",
                               "#d3d3d3"))+
  theme(axis.text = element_text(angle = 90))#+
  #scale_x_discrete(labels=seqtab_melted_e2_20)






