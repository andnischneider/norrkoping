library(here)
library(dada2)
library(readxl)
library(tidyverse)
library(decontam)
library(phyloseq)
library(ggplot2)
library(textclean)
library(Biostrings)
source(here("src/featureSelection.R"))

#################
# Data Import

seqtab_nc <- mergeSequenceTables(readRDS(here("results/N2017_ITS_1/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2017_ITS_2/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2017_ITS_3/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2017_ITS_R/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2018_ITS_1/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2018_ITS_2_3/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2018_ITS_4/dada2/seqtab_nc.rds")),
                                 readRDS(here("results/N2018_ITS_5/dada2/seqtab_nc.rds")))

seqtab_nc2 <- t(seqtab_nc)

#Remove failed samples
sp_trash <- c("P15001_1623", "P15001_1624", "P15001_1625", "P15001_1626")
seqtab_nc2 <- seqtab_nc2[,!colnames(seqtab_nc2)%in%sp_trash]

meta_2017 <- as.data.frame(read_excel(here("doc/meta_all.xlsx"), na = "NA"))
meta_2018 <- as.data.frame(read_excel(here("doc/meta_all.xlsx"), sheet = 2, na = "NA"))
meta_2018$Block <- mgsub(meta_2018$Block,
                         c("Block_1", "Block_2", "Block_3", "Block_4"),
                         c("Block_5", "Block_6", "Block_7", "Block_8"))


meta_all <- rbind(meta_2017, meta_2018)

meta_all$Block[grepl("Block", meta_all$Block)] <- ifelse(grepl("Pine", meta_all$Species[grepl("Block", meta_all$Block)]),
                                                         paste(meta_all$Block[grepl("Block", meta_all$Block)], "_P"),
                                                         paste(meta_all$Block[grepl("Block", meta_all$Block)], "_S"))

rownames(meta_all) <- meta_all$SciLifeID

meta_all <- meta_all[colnames(seqtab_nc2),]
meta_all <- meta_all[!grepl("Mock", meta_all$Type),]
seqtab_nc2 <- seqtab_nc2[,rownames(meta_all)]



write_csv(meta_all, here("doc/meta_all.csv"))

##################################
# Cleaning of contaminants
Sam_Con <- grepl("Neg_Con", meta_all$Type)
names(Sam_Con) <- rownames(meta_all)
contam_df <- isContaminant(t(seqtab_nc2), method = "prevalence", neg = Sam_Con, threshold = 0.05)
table(contam_df$contaminant)
contams <- rownames(contam_df)[contam_df$contaminant==TRUE]
non_contams <- rownames(contam_df)[contam_df$contaminant==FALSE]
#Check how much it affects the relative abundance to remove them
boxplot(colSums(seqtab_nc2[non_contams,])/colSums(seqtab_nc2))
boxplot(colSums(seqtab_nc2[non_contams,!Sam_Con])/colSums(seqtab_nc2[,!Sam_Con]))

#For the majority of samples, the vast majority of reads belongs to real ASVs, good.
ps <- phyloseq(otu_table(seqtab_nc2, taxa_are_rows = TRUE),
               sample_data(meta_all))
sample_data(ps)$is.neg <- Sam_Con
sample_data(ps)$Sample_Or_Control <- ifelse(Sam_Con, "Control_Sample", "True_Sample")

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Or_Control=="Control_Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Or_Control=="True_Sample", ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam_df$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#This looks good (threshold of 0.05), and we proceed
meta_rel <- meta_all[!Sam_Con,]
seqtab_nc3 <- seqtab_nc2[non_contams,!Sam_Con]
dim(seqtab_nc3)

##########################
# Merging of replicates

meta_rel$SampleID_sum <- paste(meta_rel$Species, 
                               meta_rel$Type, 
                               meta_rel$Block,
                               meta_rel$Treatment,
                               meta_rel$Planting_Position,
                               meta_rel$Year,
                               sep = ".")

seqtab_sum <- t(apply(seqtab_nc3, 1, function(f){tapply(f, meta_rel$SampleID_sum, sum)}))

meta_rel_sum <- unique(meta_rel[,c(10,5,4,8,6,7,9)])
seqtab_sum <- seqtab_sum[,meta_rel_sum$SampleID_sum]
seqtab_sum <- seqtab_sum[rowSums(seqtab_sum)>0,]
dim(seqtab_sum)

#################################
# Swarm clustering
asv_seqs <- DNAStringSet(rownames(seqtab_sum))
names(asv_seqs) <- paste0("ASV", 1:nrow(seqtab_sum), ";size=", rowSums(seqtab_sum))

dir.create(here("results/swarm"), showWarnings = FALSE)
writeXStringSet(asv_seqs, file = here("results/swarm/asvs.fa"))

system(paste("swarm","-d 3 -z --output-file",
             here("results/swarm/newresults.txt"),"--seeds",
             here("results/swarm/newseeds.fa"), here("results/swarm/asvs.fa")))

# 10004 swarms from 24932 ASVs

seqtab_sum2 <- cbind(ASV_ID=paste0("ASV", 1:nrow(seqtab_sum)), (as.data.frame(data.matrix(seqtab_sum))))

newclus_swarm <- read_lines(here("results/swarm/results.txt")) %>% 
  str_split_fixed(pattern=" ",n = max(str_count(., "ASV"))) %>% 
  gsub(";.*","",.) %>% bind_cols(.,Cluster=paste0("cluster", 1:nrow(.))) %>% 
  pivot_longer(starts_with('...'),values_to="ASV_ID") %>% filter(ASV_ID!="") %>% select(Cluster,ASV_ID) %>% 
  left_join(seqtab_sum2,by="ASV_ID",copy=TRUE) 

seqtab_sw <- newclus_swarm %>% 
  group_by(Cluster) %>% select(starts_with(c("Spruce","Pine", "NA")),"Cluster") %>% 
  summarise(across(everything(), sum), .groups = 'drop') %>% 
  column_to_rownames("Cluster")

seqtab_sw <- seqtab_sw[match(unique(newclus_swarm$Cluster), rownames(seqtab_sw)),]

seqs_sw <- readDNAStringSet(here("results/swarm/seeds.fa"))
newclus_swarm2 <- newclus_swarm$ASV_ID[match(unique(newclus_swarm$Cluster), newclus_swarm$Cluster)]
seqs_sw <- seqs_sw[match(newclus_swarm2, gsub(";size=\\d+;", "", names(seqs_sw)))]
all(gsub(";size=\\d+;", "", names(seqs_sw))==newclus_swarm2)

dir.create(here("results/ITSx"), showWarnings = FALSE)
names(seqs_sw) <- rownames(seqtab_sw)
writeXStringSet(seqs_sw, here("results/ITSx/sw_seeds.fa"))

##################################
# ITSx to remove non-fungal swarms
#Here I run ITSx (to determine fungal sequences), and then import the cleaned results

system(paste("ITSx", "-i", here("results/ITSx/sw_seeds.fa"),
             "-o", here("results/ITSx/swarm_"), 
             "--cpu 1 --multi_thread F --preserve T --partial 50 --minlen 50 --detailed_results T"))

#Import cluster names identified as fungi by ITSx
itsx_res <- read_tsv(here("results/ITSx/swarm.extraction.results"), col_names = F)

itsx_fun <- itsx_res$X1[itsx_res$X3=="F"]
itsx_fun <- itsx_fun[!is.na(itsx_fun)]
seqtab_sw_itsx <- seqtab_sw[itsx_fun,]
boxplot(colSums(seqtab_sw_itsx)/colSums(seqtab_sw))



boxplot(colSums(seqtab_sw_itsx))
#also the vast majority is over 100k, so no problems whatsoever!

#####################
# Abundance filtering
meta_rel_sum$Group <- paste(meta_rel_sum$Species, 
                            meta_rel_sum$Type, 
                            meta_rel_sum$Treatment, 
                            meta_rel_sum$Planting_Position,
                            sep = ".")
ccc <- factor(meta_rel_sum$Group, levels = unique(meta_rel_sum$Group))
names(ccc) <- meta_rel_sum$SampleID_sum

dim(seqtab_sw_itsx)
seqtab_sw_itsx2 <- data.matrix(seqtab_sw_itsx)
seqtab_sw_itsx_filt <- seqtab_sw_itsx2[featureSelectProp(seqtab_sw_itsx2, ccc, 0.00005),]
dim(seqtab_sw_itsx_filt)
seqtab_sw_itsx_filt <- seqtab_sw_itsx2[featureSelect(seqtab_sw_itsx2, ccc, 3, 2),]
dim(seqtab_sw_itsx_filt)

seqs_sw_filt <- seqs_sw[rownames(seqtab_sw_itsx_filt)]
dir.create(here("results/constax"), showWarnings = FALSE)
writeXStringSet(seqs_sw_filt, here("results/constax/seqs_sw_filt.fasta"))

#RUN CONSTAX ON SERVER


################################
# Taxonomy import and cleanup
tax_cs <- read.table(here("results/constax/constax_taxonomy.txt"), sep = "\t", header = TRUE)
tax_cs$Species <- gsub("_$", "", gsub(" ", "_", tax_cs$Species))
#tax_cs$Species <- ifelse(tax_cs$Species==tax_cs$Genus, "", tax_cs$Species)
tax_cs <- as.data.frame(sapply(tax_cs, function(x) gsub("_1$", "", x)))

tax_cs2 <- tax_cs
tax_cs2$Order[grepl("^Gs", tax_cs2$Order)] <- ""
tax_cs2$Class[grepl("^Gs", tax_cs2$Class)] <- ""
tax_cs2[tax_cs2==""] <- NA
#tax_cs2$Family <- ifelse(grepl("Incertae", tax_cs2$Family), paste0("Uncertain.", tax_cs2$Genus), tax_cs2$Family)
#tax_cs2$Family <- gsub("Uncertain.NA", NA, tax_cs2$Family)

# tax_cs2$Order <- ifelse(grepl("Incertae", tax_cs2$Order), paste0("Uncertain.", tax_cs2$Family), tax_cs2$Order)
# tax_cs2$Order <- gsub("Uncertain.NA", NA, tax_cs2$Order)
# tax_cs2$Order <- gsub("Uncertain.Uncertain", "Uncertain", tax_cs2$Order) 

tax_cs2[is.na(tax_cs2)] <- "unidentified"

tax_cs2$Phylum <- ifelse(grepl("unidentified", tax_cs2$Phylum), 
                         paste0("Unclassified.", tax_cs2$Kingdom),
                         tax_cs2$Phylum)

#Let's make a function to automate this
fixNAtax <- function(tax, rank) {
  coln <- which(colnames(tax)==rank)
  namen <- colnames(tax)[coln]
  namen1 <- colnames(tax)[coln-1]
  tax[,namen] <- ifelse(grepl("Unclassified", tax[,namen1]),
                        tax[,namen1],
                        tax[,namen])
  tax[,namen] <- ifelse(grepl("unidentified", tax[,namen]),
                        paste0("Unclassified.", tax[,namen1]),
                        tax[,namen])
  return(tax)
}
tax_cs3 <- fixNAtax(tax_cs2, "Class")
tax_cs3 <- fixNAtax(tax_cs3, "Order")
tax_cs3 <- fixNAtax(tax_cs3, "Family")
tax_cs3 <- fixNAtax(tax_cs3, "Genus")
tax_cs3 <- fixNAtax(tax_cs3, "Species")

#Correct the blocks


saveRDS(tax_cs3, here("results/new_clean/tax_clean.rds"))
saveRDS(seqtab_sw_itsx_filt, here("results/new_clean/seqtab_clean.rds"))
saveRDS(meta_rel_sum, here("results/new_clean/meta_clean_sum.rds"))














