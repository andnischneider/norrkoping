library(here)
library(ggplot2)
library(vegan)
library(phyloseq)
library(VennDiagram)
library(lefser)
library(gridExtra)
library(reshape2)
library(rcartocolor)
library(Vennerable)
library(Biostrings)
library(patchwork)
library(ggsignif)


source(here("../data/src/utilsDE.r"))

ggformat <- theme_minimal()+theme(axis.text.x = element_text(angle = 45))
melt_filt2 <- function (mat) {
  melt1 <- reshape2::melt(mat, varnames = c("OTU_ID", "SampleID_sum"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}

plot_venn2 <- function (list) {
  n <- length(list)
  if (n==2) {
    venn <- Venn(SetNames = names(list),
                 Weight = c(`10` = length(setdiff(list[[1]], list[[2]])),
                            `01` = length(setdiff(list[[2]], list[[1]])),
                            `11` = length(intersect(list[[1]], list[[2]]))))
    plot(venn, doWeights = TRUE, type = "circles")
  } else if (n==3) {
    venn <- Venn(SetNames = names(list),
                 Weight = c(`100`= length(setdiff(list[[1]], union(list[[2]], list[[3]]))), 
                            `010`= length(setdiff(list[[2]], union(list[[1]], list[[3]]))), 
                            `001`=length(setdiff(list[[3]], union(list[[2]], list[[1]]))), 
                            `110`= length(setdiff(intersect(list[[1]], list[[2]]), list[[3]])), 
                            `011`= length(setdiff(intersect(list[[3]], list[[2]]), list[[1]])), 
                            `101`= length(setdiff(intersect(list[[1]], list[[3]]), list[[2]])),
                            `111`= length(intersect(list[[1]], intersect(list[[2]], list[[3]])))))
    plot(venn, doWeights = TRUE, type = "circles")
  } else if (n>3) {
    print("This function supports a maximum of three sets")
  }
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



nursery <- meta$SampleID_sum[meta$Treatment=="Nursery" & !is.na(meta$Treatment)]

## diversity comparison split by species (boxplot + lm)
seqtab_n <- seqtab[,nursery]
seqtab_n <- seqtab_n[rowSums(seqtab_n)>3,]
meta_n <- droplevels(meta[nursery,])

div_all <- diversity(seqtab_n, MARGIN = 2)
div_all <- data.frame(meta_n, Shannon=div_all)

ggplot(div_all, aes(x = Type, y = Shannon, fill = Type))+
  geom_boxplot()+
  geom_signif(
    comparisons = list(c("Roots", "Soil_Nursery")),
    map_signif_level = TRUE
  )+
  facet_wrap(~Species)+
  scale_fill_manual(values = c("#4886D9", "#50592E"))+
  theme_light()+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  labs(y = "Shannon Diversity Index")+
  scale_x_discrete(labels=c("Roots", "Peat"))

ggsave(here("FigureS4A_NurDiv.pdf"), 
       width = 70, 
       height = 100, 
       units = "mm", 
       dpi = 300)


wilcox.test(Shannon~Type, data = div_all)
wilcox.test(Shannon~Type, data = div_all[div_all$Species=="Spruce",])
wilcox.test(Shannon~Type, data = div_all[div_all$Species=="Pine",])

# Let's do a PcoA and community difference test as well
set.seed(66)
seqtab_rar <- t(rrarefy(t(seqtab), min(colSums(seqtab))))
seqtab_rar_n <- t(rrarefy(t(seqtab_n), min(colSums(seqtab_n))))
#seqtab_rar <- seqtab_n
seqtab_rar_n <- seqtab_rar_n[rowSums(seqtab_rar_n)>3,]

ps_nur_all <- phyloseq(sample_data(meta_n),
                       otu_table(seqtab_rar_n, taxa_are_rows = TRUE))

ord_nur_all <- ordinate(ps_nur_all, "MDS", "bray")

plot_ordination(ps_nur_all, ord_nur_all, type = "samples", color = "Species", shape = "Type")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(ord_nur_all$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(ord_nur_all$values$Relative_eig[2]*100, digits = 1), "%]"))+
  theme(legend.position = "none")+
  theme_light()+
  scale_color_manual(values = c("#A68E46", "#4C5939"))+
  scale_shape_discrete(labels=c("Roots", "Peat"))

ggsave(here("FigureS4B_NurPCoA.pdf"), 
       width = 120, 
       height = 100, 
       units = "mm", 
       dpi = 300)
# this looks very significant, let's test it.
adonis(vegdist(t(seqtab_rar_n)) ~ meta_n$Type*meta_n$Species)
# Significant effect of both species and sample type on community composition 

########
#Number of common taxa and differential abundance
seqtab_rar_s <- seqtab_rar_n[,grepl("Spruce", colnames(seqtab_rar_n))]
seqtab_rar_p <- seqtab_rar_n[,grepl("Pine", colnames(seqtab_rar_n))]
sotu_list <- list()
sotu_list[["Spruce_SOTUs"]] <- rownames(seqtab_rar_s)[rowSums(seqtab_rar_s)>3]
sotu_list[["Pine_SOTUs"]] <- rownames(seqtab_rar_p)[rowSums(seqtab_rar_p)>3]

#How many are common
length(intersect(sotu_list[[1]], sotu_list[[2]]))
length(setdiff(sotu_list[[1]], sotu_list[[2]]))
length(setdiff(sotu_list[[2]], sotu_list[[1]]))

# grid.newpage()
# futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# grid.draw(venn.diagram(list(Spruce=sotu_list[["Spruce_SOTUs"]],
#                             Pine=sotu_list[["Pine_SOTUs"]]),
#                   NULL))

#####
#differential abundance (Lefser)
#Roots
seqtab_rar_rt <- seqtab_rar_n[,grepl("Roots", colnames(seqtab_rar_n))]
seqtab_rar_rt <- seqtab_rar_rt[rowSums(seqtab_rar_rt)>3,]
sumobj_rt <- SummarizedExperiment(assays <- list(as.matrix(seqtab_rar_rt)),
                               colData = meta_n[meta_n$Type=="Roots",])
res_lefs_rt <- lefser(sumobj_rt, groupCol = "Species")

lefserPlot(res_lefs_rt)+
  theme_minimal()+
  ggtitle("SOTU level LDA scores Sucrose samples (litter vs no litter)")

levels(as.factor(meta_n$Species))

tax_res_rt <- tax[res_lefs_rt$Names,]

#Soil
seqtab_rar_so <- seqtab_rar_n[,grepl("Soil", colnames(seqtab_rar_n))]
seqtab_rar_so <- seqtab_rar_so[rowSums(seqtab_rar_so)>3,]
sumobj_so <- SummarizedExperiment(assays <- list(as.matrix(seqtab_rar_so)),
                                  colData = meta_n[meta_n$Type=="Soil_Nursery",])
res_lefs_so <- lefser(sumobj_so, groupCol = "Species")

lefserPlot(res_lefs_so)+
  theme_minimal()+
  ggtitle("SOTU level LDA scores Sucrose samples (litter vs no litter)")

levels(as.factor(meta_n$Species))

tax_res_so <- tax[res_lefs_so$Names,]

#####
## most abundant species in spruce vs pine
##sum to species level
seqtab_n_m <- melt_filt2(seqtab_n)
#counts_m2 <- merge.data.frame(counts_m, meta_r_its, by = "SampleID")
seqtab_n_m2 <- merge.data.frame(seqtab_n_m, tax, by = "OTU_ID")

seqtab_spec <- acast(seqtab_n_m2, Species~SampleID_sum, value.var = "Count", fun.aggregate = sum)

#propmat
seqtab_spec_p <- prop.table(seqtab_spec, margin = 2)*100

#Remove "Unclassified.Fungi"
seqtab_spec_p <- seqtab_spec_p[rownames(seqtab_spec_p)!="Unclassified.Fungi",]

#Divide into species
seqtab_spec_p_spr <- seqtab_spec_p[,grepl("Spruce", meta_n$Species)]
seqtab_spec_p_pin <- seqtab_spec_p[,grepl("Pine", meta_n$Species)]

##Extract the most abundant ones per species
top10_spr <- rownames(seqtab_spec_p_spr)[order(rowMeans(seqtab_spec_p_spr), decreasing = TRUE)][1:15]
top10_pin <- rownames(seqtab_spec_p_pin)[order(rowMeans(seqtab_spec_p_pin), decreasing = TRUE)][1:15]


seqtab_spec_p_spr_top10 <- seqtab_spec_p_spr[top10_spr,]
spr_top10_plot <- reshape2::melt(seqtab_spec_p_spr_top10, varnames = c("Species", "Sample"), value.name = "Percentage")
spr_top10_plot$Species <- factor(spr_top10_plot$Species, levels = rev(top10_spr))
spr_top10_plot$Type <- ifelse(grepl("Roots", spr_top10_plot$Sample), "Roots", "Soil")
p1 <- ggplot(spr_top10_plot, aes(y = Species, x = Percentage))+
  geom_boxplot(aes(col = Type))+
  ggtitle("Spruce")+
  theme_minimal()

seqtab_spec_p_pin_top10 <- seqtab_spec_p_pin[top10_pin,]
pin_top10_plot <- reshape2::melt(seqtab_spec_p_pin_top10, varnames = c("Species", "Sample"), value.name = "Percentage")
pin_top10_plot$Species <- factor(pin_top10_plot$Species, levels = rev(top10_pin))
pin_top10_plot$Type <- ifelse(grepl("Roots", pin_top10_plot$Sample), "Roots", "Soil")
p2 <- ggplot(pin_top10_plot, aes(y = Species, x = Percentage))+
  geom_boxplot(aes(col = Type))+
  ggtitle("Pine")+
  theme_minimal()

grid.arrange(p1, p2, nrow = 1)


######## THE SAME WITH SOTUs instead of summarizing to Species (doesn't make a big difference apparently)

#Convert to rel ab
seqtab_n_p <- prop.table(seqtab_n, margin = 2)*100
seqtab_p_spr <- seqtab_n_p[,grepl("Spruce", meta_n$Species)]
seqtab_p_pin <- seqtab_n_p[,grepl("Pine", meta_n$Species)]

##Extract the most abundant ones per species
top10_f_spr <- rownames(seqtab_p_spr)[order(rowMeans(seqtab_p_spr), decreasing = TRUE)][1:15]
top10_f_pin <- rownames(seqtab_p_pin)[order(rowMeans(seqtab_p_pin), decreasing = TRUE)][1:15]


seqtab_p_spr_top10 <- seqtab_p_spr[top10_f_spr,]
spr_f_top10_plot <- reshape2::melt(seqtab_p_spr_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
spr_f_top10_plot$SOTU <- factor(spr_f_top10_plot$SOTU, levels = rev(top10_f_spr))
spr_f_top10_plot$Type <- ifelse(grepl("Roots", spr_f_top10_plot$Sample), "Roots", "Soil")
spr_f_top10_plot$SOTU2 <- paste(spr_f_top10_plot$SOTU,
                                tax$Species[match(as.character(spr_f_top10_plot$SOTU), tax$OTU_ID)],
                                sep = ".")
spr_f_top10_plot$SOTU2 <- factor(spr_f_top10_plot$SOTU2, levels = rev(unique(as.character(spr_f_top10_plot$SOTU2))))
pf1 <- ggplot(spr_f_top10_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Type))+
  ggtitle("Norway spruce")+
  theme_light()+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values = c("#4886D9", "#50592E"))

seqtab_p_pin_top10 <- seqtab_p_pin[top10_f_pin,]
pin_f_top10_plot <- reshape2::melt(seqtab_p_pin_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
pin_f_top10_plot$SOTU <- factor(pin_f_top10_plot$SOTU, levels = rev(top10_f_pin))
pin_f_top10_plot$Type <- ifelse(grepl("Roots", pin_f_top10_plot$Sample), "Roots", "Soil")
pin_f_top10_plot$SOTU2 <- paste(pin_f_top10_plot$SOTU,
                                tax$Species[match(as.character(pin_f_top10_plot$SOTU), tax$OTU_ID)],
                                sep = ".")
pin_f_top10_plot$SOTU2 <- factor(pin_f_top10_plot$SOTU2, levels = rev(unique(as.character(pin_f_top10_plot$SOTU2))))
pf2 <- ggplot(pin_f_top10_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Type))+
  ggtitle("Scots pine")+
  theme_light()+
  theme(axis.title.y = element_blank())+
  scale_color_manual(values = c("#4886D9", "#50592E"))  

pf_comb <- pf1 + pf2 & theme(legend.position = "bottom")
pf_comb + plot_layout(guides = "collect")
ggsave(here("FigureS5_MA_Nur.pdf"), 
       width = 190, 
       height = 100, 
       units = "mm", 
       dpi = 300)

#clusters with at least 0.05% relative abundance in spruce and pine
pin_clus_min <- rownames(seqtab_p_spr)[which(rowMeans(seqtab_p_spr)>0.1)]
spr_clus_min <- rownames(seqtab_p_pin)[which(rowMeans(seqtab_p_pin)>0.1)]

intersect(pin_clus_min, spr_clus_min)



pin_clus_min2 <- rownames(seqtab_p_spr)[which(rowMeans(seqtab_p_spr)>0.05)]
spr_clus_min2 <- rownames(seqtab_p_pin)[which(rowMeans(seqtab_p_pin)>0.05)]

###DA SOTUs between planting positions 2017
pp_2017 <- c("cluster254", "cluster258", "cluster925", 
             "cluster158", "cluster370", "cluster405",
             "cluster569")

intersect(pin_clus_min2, pp_2017)
intersect(spr_clus_min2, pp_2017)
intersect(rownames(seqtab_p_spr), pp_2017)
intersect(rownames(seqtab_p_pin), pp_2017)

################################################
###ROOTS ONLY
#Alpha div significantly different?
meta_rt <- meta[colnames(seqtab_rar_rt),]
wilcox.test(diversity(seqtab_rar_rt, MARGIN = 2)~meta_rt$Species)

#Beta div?
set.seed(6)
adonis(vegdist(t(seqtab_rar_rt)) ~ meta_rt$Species)

#Total number of SOTUs
nrow(seqtab_rar_rt)

#How many unique to each species?
seqtab_rar_rt_p <- seqtab_rar_rt[,grepl("Pine", colnames(seqtab_rar_rt))]
seqtab_rar_rt_s <- seqtab_rar_rt[,grepl("Spruce", colnames(seqtab_rar_rt))]
pine_rt_sotus <- rownames(seqtab_rar_rt_p[rowSums(seqtab_rar_rt_p)>3,])
spruce_rt_sotus <- rownames(seqtab_rar_rt_s[rowSums(seqtab_rar_rt_s)>3,])
#common
length(intersect(pine_rt_sotus, spruce_rt_sotus))
length(setdiff(pine_rt_sotus, spruce_rt_sotus))
length(setdiff(spruce_rt_sotus, pine_rt_sotus))

#Total relative abundance of unique SOTUs?
uniq_pine <- setdiff(pine_rt_sotus, spruce_rt_sotus)
median(colSums(seqtab_rar_rt_p[uniq_pine,])/colSums(seqtab_rar_rt_p))*100
uniq_spruce <- setdiff(spruce_rt_sotus, pine_rt_sotus)
median(colSums(seqtab_rar_rt_s[uniq_spruce,])/colSums(seqtab_rar_rt_s))*100

#Most abundant genera in spruce and pine, how do they compare
seqtab_gen <- acast(seqtab_n_m2, Genus~SampleID_sum, value.var = "Count", fun.aggregate = sum)
seqtab_gen_r <- seqtab_gen[,grepl("Roots", colnames(seqtab_gen))]
seqtab_gen_r <- prop.table(seqtab_gen_r, margin = 2)
genera_spruce <- rownames(seqtab_gen_r)[order(rowMeans(seqtab_gen_r[,grepl("Spruce", colnames(seqtab_gen_r))]), decreasing = TRUE)][1:15]
genera_pine <- rownames(seqtab_gen_r)[order(rowMeans(seqtab_gen_r[,grepl("Pine", colnames(seqtab_gen_r))]), decreasing = TRUE)][1:15]

#Most abundant SOTUs, and which are common/unique
seqtab_rar_rt_s_p <- prop.table(seqtab_rar_rt_s, margin = 2)*100
seqtab_rar_rt_p_p <- prop.table(seqtab_rar_rt_p, margin = 2)*100
spruce_rt_sotus_15 <- spruce_rt_sotus[order(rowMeans(seqtab_rar_rt_s_p[spruce_rt_sotus,]), decreasing = TRUE)][1:15]
pine_rt_sotus_15 <- pine_rt_sotus[order(rowMeans(seqtab_rar_rt_p_p[pine_rt_sotus,]), decreasing = TRUE)][1:15]

#clusters with at least 0.05% relative abundance in spruce and pine
pin_clus_rt_min <- rownames(seqtab_rar_rt_p_p)[which(rowMeans(seqtab_rar_rt_p_p)>0.05)]
spr_clus_rt_min <- rownames(seqtab_rar_rt_s_p)[which(rowMeans(seqtab_rar_rt_s_p)>0.05)]

intersect(pin_clus_rt_min, spr_clus_rt_min)
setdiff(spruce_rt_sotus_15, intersect(pin_clus_rt_min, spr_clus_rt_min))
setdiff(pine_rt_sotus_15, intersect(pin_clus_rt_min, spr_clus_rt_min))

###################################################
#### Comparison nursery peat - control block soil
#Create top 15 object with nursery peat samples
seqtab_n_r <- seqtab_rar_n[,nursery]
seqtab_n_peat <- seqtab_n[,grepl("Soil", colnames(seqtab_n_r))]
seqtab_n_peat <- seqtab_n_peat[rowSums(seqtab_n_peat)>3,]
seqtab_n_p_peat <- prop.table(seqtab_n_peat, margin = 2)*100
top10_f_peat_all <- rownames(seqtab_n_p_peat)[order(rowMeans(seqtab_n_p_peat), decreasing = TRUE)][1:15]

seqtab_peat_top10 <- seqtab_n_p_peat[top10_f_peat_all,]
peat_top10_plot <- reshape2::melt(seqtab_peat_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
peat_top10_plot$SOTU <- factor(peat_top10_plot$SOTU, levels = rev(top10_f_peat_all))
peat_top10_plot$Species <- ifelse(grepl("Pine", peat_top10_plot$Sample), "Pine", "Spruce")
peat_top10_plot$SOTU2 <- paste(peat_top10_plot$SOTU,
                               tax$Species[match(as.character(peat_top10_plot$SOTU), tax$OTU_ID)],
                               sep = ".")
peat_top10_plot$SOTU2 <- factor(peat_top10_plot$SOTU2, levels = rev(unique(as.character(peat_top10_plot$SOTU2))))

(p_comp_1 <- ggplot(peat_top10_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Species))+
  ggtitle("Nursery peat")+
  theme_light()+
  scale_color_manual(values = c("#4886D9", "#50592E"))+
  theme(axis.title.y = element_blank()))  

##
#Same with control block soil!
control_soil_block <- meta$SampleID_sum[grepl("Soil_Block", meta$Type) 
                                             & meta$Treatment=="Control"
                                             & !is.na(meta$Treatment)]
meta_csb <- meta[control_soil_block,]
seqtab_csb <- seqtab_rar[,control_soil_block]
seqtab_csb <- seqtab_csb[rowSums(seqtab_csb)>3,]


#Convert to rel ab
seqtab_csb_p <- prop.table(seqtab_csb, margin = 2)*100

##Quickly check which ones from the DA analysis are among the SOTUs
##found on the site
intersect(rownames(seqtab_csb_p), pp_2017)
rowMeans(seqtab_csb_p)[intersect(rownames(seqtab_csb_p), pp_2017)]

#Extract 15 most abundant
top10_csb_p <- rownames(seqtab_csb_p)[order(rowMeans(seqtab_csb_p), decreasing = TRUE)][1:15]
seqtab_csb_top10 <- seqtab_csb_p[top10_csb_p,]
csb_top10_plot <- reshape2::melt(seqtab_csb_top10, varnames = c("SOTU", "Sample"), value.name = "Percentage")
csb_top10_plot$SOTU <- factor(csb_top10_plot$SOTU, levels = rev(top10_csb_p))
csb_top10_plot$Planting_Position <- ifelse(grepl("Capped", csb_top10_plot$Sample), "Capped_Mound", "Exposed_Mineral")
csb_top10_plot$SOTU2 <- paste(csb_top10_plot$SOTU,
                              tax$Species[match(as.character(csb_top10_plot$SOTU), tax$OTU_ID)],
                              sep = ".")
csb_top10_plot$SOTU2 <- factor(csb_top10_plot$SOTU2, levels = rev(unique(as.character(csb_top10_plot$SOTU2))))

(p_comp_2 <- ggplot(csb_top10_plot, aes(y = SOTU2, x = Percentage))+
  geom_boxplot(aes(col = Planting_Position))+
  ggtitle("Control block soil")+
  theme_light()+
  scale_color_manual(values = carto_pal(2, "Earth"), 
                     labels=c("Capped Mound", "Exposed Mineral"),
                     name="Planting Position")+
  theme(axis.title.y = element_blank()))


p_comp_comb <- p_comp_1 + p_comp_2 & theme(legend.position = "bottom")
p_comp_comb + plot_layout(guides = "collect")
ggsave(here("FigureS3_NPeat_ControlBulk.pdf"), 
       width = 200, 
       height = 100, 
       units = "mm", 
       dpi = 300)

#################
#Which of the top15 here are present in the other set?
#Common top15 control block soils - nursery peat all
setdiff(top10_csb_p, rownames(seqtab_n_p_peat))

#Common top15 nursery peat- control block soils all
setdiff(top10_f_peat_all, rownames(seqtab_csb_p))


##########################################
##Venns diagram of SOTUs between nursery peat and block/edge soils
comp_list <- list(Nursery_Peat = rownames(seqtab_n_peat),
                  Control_Blocks = rownames(seqtab_csb))

plot_venn2(comp_list)

###SPLIT BY YEAR
##2017
seqtab_csb_2017 <- seqtab_csb[,grepl("2017", colnames(seqtab_csb))]
seqtab_csb_2017 <- seqtab_csb_2017[rowSums(seqtab_csb_2017)>3,]

seqtab_csb_2018 <- seqtab_csb[,grepl("2018", colnames(seqtab_csb))]
seqtab_csb_2018 <- seqtab_csb_2018[rowSums(seqtab_csb_2018)>3,]

comp_list_2017 <- list(Nursery_Peat = rownames(seqtab_n_peat),
                       Control_Blocks_2017 = rownames(seqtab_csb_2017),
                       Control_Blocks_2018 = rownames(seqtab_csb_2018))

plot_venn2(comp_list_2017)

##Overlap nursery peat, control block (2017+2018) edge samples
edges <- meta$SampleID_sum[grepl("Edge", meta$Block)]

seqtab_e <- seqtab_rar[,edges]
seqtab_e <- seqtab_e[rowSums(seqtab_e)>3,]

comp_list_edg <- list(Nursery_Peat = rownames(seqtab_n_peat),
                      Control_Blocks = rownames(seqtab_csb),
                      Edges = rownames(seqtab_e))

plot_venn2(comp_list_edg)

##alpha diversity
seqtab_comp <- seqtab_rar[,c(control_soil_block, nursery[grepl("Soil", nursery)])]
seqtab_comp <- seqtab_comp[rowSums(seqtab_comp)>3,]

meta_comp <- meta[colnames(seqtab_comp),]

comp_div <- data.frame(Shannon=diversity(seqtab_comp, MARGIN = 2),
                       meta_comp)

comp_div$Treatment <- droplevels(comp_div$Treatment)
boxplot(Shannon~Treatment, data = comp_div)
wilcox.test(Shannon~Treatment, data = comp_div)

##NMDS 
ps_comp <- phyloseq(otu_table(seqtab_comp, taxa_are_rows = T),
                    sample_data(meta_comp))

comp_ord <- ordinate(ps_comp, method = "NMDS", "bray")

plot_ordination(ps_comp, comp_ord, color = "Treatment", shape = "Year")


############## BOTRYTIS???

#The closest taxonomic annotation I could find is Unclassified.Helotiales
helo_clus <- tax$OTU_ID[tax$Genus=="Unclassified.Helotiales"]
seqs <- readDNAStringSet(here("../../Analysis/results/constax/seqs_sw_filt.fasta"))
seqs <- seqs[helo_clus]
writeXStringSet(seqs, here("Helotiales.fasta"), width = 1000)

