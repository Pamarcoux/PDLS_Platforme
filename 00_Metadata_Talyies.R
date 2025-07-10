library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(Seurat)
library(readxl)
library(reshape2)
library(paletteer)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(paletteer)
library(magick)
library(grid)


setwd('~/Postdoc/R/SC_Flex_TALYIES_July2024/')
data_talyies_full <- read_csv('Metadata/Data_TALYIES_grouped.csv')
source('Figures/Replacement_table.R')  

Sample_order <- c("FL1_PB1","FL1_LN1_1","FL1_LN1_2","FL2_LN1", "FL1_LN2", "FL2_PB1", "FL3_PB1", "FL3_LN1", "FL4_LN1","FL4_BM1", "FL5_PB1_1","FL5_PB1_2","FL5_PB1_3", "FL6_PB1_1","FL6_PB1_2","FL7_LN1","FL8_LN1","FL9_PB1","FL10_PB1","FL10_PB2" ,"FL10_LN2","FL10_BM1", "FL11_PB1", "FL11_LN1", "FL12_PB1",
                  "FL13_PB1", "FL14_LN1","FL15_PB1", "FL16_PB1_1","FL16_PB1_2", "FL17_LN1", "FL18_LN1", "FL19_LN1", "tFL1_LN1","tFL2_LN1","tFL3_PB1","tFL3_LN1","tFL4_PB1", "tFL5_PB1","tFL6_PB1", "tFL7_PB1", "DLBCL1_LN1","DLBCL2_LN1","DLBCL3_PB1","DLBCL4_PB1","DLBCL5_PB1", "DLBCL6_LN1",
                  "MZL1_S1","MZL2_LN1", "MZL3_PB1","MZL3_PB2", "MZL4_S1","MZL5_PB1","tMZL2_LN1")  

data_table_SC_population <- read_csv('Metadata/data_table_SC_population.csv')


data_talyies_full <- data_talyies_full %>% 
  mutate(
    Screening = as.logical(Screening),
    ScRNA_flex = as.logical(ScRNA_flex),
    ScRNA_cite = as.logical(ScRNA_cite),
    D_R = as.factor(D_R),
    Origin = factor(Origin, levels = c("PBMC","LN","BM","Spleen")),
    Treatment = as.factor(Treatment),
    Response_type_total = factor(Response_type_total, levels = c("Low", "Medium", "High")),
    Pop= as.factor(Pop),
    Day = as.factor(Day),
    across(20:length(.),as.numeric)) %>% 
  mutate(Response_Origin = paste(Response_type_total, Origin, sep = "_")) %>% 
  mutate(Disease = ifelse(Disease == "FL_DLBCL","tFL",Disease)) %>% 
  mutate(Disease = ifelse(Disease == "FLt","tFL",Disease)) %>% 
  mutate(Disease = as.factor(Disease)) %>% 
  mutate(Sample_code_paper = Sample_code) %>% 
  mutate(Sample_code_paper = recode(Sample_code_paper, !!!Replacement_table)) %>% 
  mutate(Sample_code_paper = factor(Sample_code_paper, levels = Replacement_table)) %>% 
  mutate(Sample_code = factor(Sample_code, levels = Sample_order))  


######data_talyies_CD20_neg######
data_talyies_CD20_neg <- data_talyies_full %>% 
  filter(Day == "D0" & Screening == T & Pop == "CD22") %>% 
  mutate(ratio_CD20 = Per_CD20/Per_CD22_total) %>% 
  relocate(ratio_CD20, .after = Per_CD20) %>%
  filter(ratio_CD20 < 0.09 | Patient_code %in% c("201326988_20221028", "19T005587-1", "199601089","202200610","202100777","201807661"))

## 19T005587-1 DLBCL1-LN 10% CD22 UT D6
## 201326988_20221028 FL10_PB1 Pas de B D6
## 202100777 / FL11_PB1 CD20 neg
## 202200610 / FL13_PB1 CD20neg
## 201807661 / tFL7 CD20 neg ? 
## 025508129 / tMZL1_PB tMZL
## 19T058974 / tMZL2_LN1
CD20_neg_data <- setdiff(c(data_talyies_CD20_neg$Patient_code, "19T005587-1", "199601089","201807661","025508129","19T058974"), "201908255")
### 201908255 FL1_PB très peu de B CD20 mais reponse quand même

# data_talyies_good <- data_talyies %>% 
#   filter(!Patient_code %in% CD20_neg_data & Treatment %in% c("αCD20-TCB 0,1 nM"))

### data_talyies_filtered ########
data_talyies_filtered <- data_talyies_full %>% 
  filter(Screening == T & ScRNA_flex == T & !Patient_code %in% CD20_neg_data & Disease %in% c("FL","tFL","DLBCL") & Day == "D6") 

Sample_list <- data_talyies_filtered %>% 
  distinct(Sample_code) %>% pull(Sample_code)

Sample_list_extended <- data_talyies_full %>% 
  filter(ScRNA_flex == T & !Patient_code %in% CD20_neg_data & Disease %in% c("FL","tFL","DLBCL")) %>% 
  distinct(Sample_code) %>% pull(Sample_code)

Sample_list_CD20_pos <- data_talyies_full %>% 
  filter(Screening == T & Disease %in% c("FL","tFL","DLBCL") & !Patient_code %in% CD20_neg_data) %>% 
  mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>% 
  distinct(Sample_code) %>% 
  pull(Sample_code)

######Annotation Response_type Single Cell Data ########################
data_talyies_full <- data_talyies_full %>% select(-Response_type_total)
data_response_type_CD10 <- data_talyies_full %>% 
  filter(Treatment %in% c("αCD20-TCB 0,1 nM")) %>% 
  mutate(Response_type_CD10 = case_when(
    Patient_code %in% CD20_neg_data ~ NA,
    B_cell_depletion_CD22_CD10_plus < 45 ~ "Low",
    # B_cell_depletion_CD22_CD10_plus >= 40 & B_cell_depletion_CD22_CD10_plus <= 59 ~ "ZB",  
    B_cell_depletion_CD22_CD10_plus > 45  ~ "High")) %>%
  mutate(Response_type_total = case_when(
    Patient_code %in% CD20_neg_data ~ NA,
    Patient_code %in% c("19T005587-1")  ~ "High",
    B_cell_depletion_total < 45 ~ "Low",
    # B_cell_depletion_total >= 40 & B_cell_depletion_total <= 59 ~ "ZB",  
    B_cell_depletion_total > 45  ~ "High"))%>%
    # B_cell_depletion_CD22_CD10_plus > 89 ~ "VeryHigh")) %>% 

  mutate(B_cell_depletion_CD22_CD10_plus_glofi = B_cell_depletion_CD22_CD10_plus,
         B_cell_depletion_total_glofi = B_cell_depletion_total)%>% 
  reframe(Patient_code,Response_type_CD10,B_cell_depletion_CD22_CD10_plus_glofi,Sample_code,Response_type_total,B_cell_depletion_total_glofi,Sample_code_paper) %>% 
  mutate(Response_type_CD10 = factor(Response_type_CD10, levels = c("Low","Medium", "High")),
         Response_type_total = factor(Response_type_total, levels = c("Low","Medium", "High")))

data_talyies_full <- data_talyies_full %>% full_join(data_response_type_CD10)

data_sample_info_complete <- data_talyies_full %>% 
  select(Patient_code ,  Date_prelev ,  Cevi_code ,  Disease ,
     Origin ,  D_R ,  Date_D0 ,  Screening ,
     ScRNA_flex ,  ScRNA_cite ,  Name_cite ,  Sample_code ,
     Response_type_CD10 ,  B_cell_depletion_CD22_CD10_plus_glofi ,
     Response_type_total ,  B_cell_depletion_total_glofi ,
     Sample_code_paper) %>% 
  distinct(Patient_code, .keep_all =TRUE)

#### Data Papier Glofi #####
colors_response_type <- c("lightskyblue2", "lightseagreen")
sample_colors_LN_PBMC <- setNames(c("#a00000","#1a80bb"), c("LN","PBMC"))

data_talyies <- data_talyies_full %>% 
  filter(Screening == T & Disease %in% c("FL","tFL","DLBCL") & !Patient_code %in% CD20_neg_data) %>% 
  mutate(Sample_code = factor(Sample_code, levels = Sample_order))

data_talyies_glofi <- data_talyies_full %>% 
  filter(Screening == T & Disease %in% c("FL","tFL","DLBCL")) %>% 
  filter(!is.na(B_cell_depletion_total_glofi)) %>%  
  mutate(Sample_code = factor(Sample_code, levels = Sample_order))

data_sample_info_filter <- data_talyies %>% 
  mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>% 
  select(Sample_code) %>%
  distinct() %>% 
  left_join(data_sample_info_complete) 

Sample_list_paper <- data_talyies_glofi %>% 
  distinct(Sample_code_paper) %>% pull(Sample_code_paper)

Sample_list_flex <- data_talyies_filtered %>% 
  distinct(Sample_code_paper) %>% pull(Sample_code_paper)

Sample_left <- setdiff(Sample_list_paper,Sample_list_flex)
colors_sample_left_list <- paletteer_d("ggsci::default_igv", length(Sample_left))
sample_color_left <- setNames(colors_sample_left_list, Sample_left)

colors_sample_flex <- rev(paletteer_d("rcartocolor::Bold"))[1:length(Sample_list_flex)]
sample_colors_flex <- setNames(colors_sample_flex,Sample_list_flex[order(factor(Sample_list_flex))])


sample_colors_paper <- c(sample_colors_flex,sample_color_left)



### Data for Paper Plaftform ####
data_talyies_LN_PBMC <- data_talyies_full %>% 
  filter(Disease %in% c("FL","DLBCL","tFL") & Origin %in% c("LN","PBMC"))

##### FL####
data_talyies_FL <- data_talyies_full %>% 
  filter(Disease %in% c("FL") & Origin %in% c("LN","PBMC")) 

Sample_list_FL <- data_talyies_FL %>%
  distinct(Sample_code) %>% arrange(Sample_code) %>% pull(Sample_code)

# Générer une palette de couleurs
# colors_FL <- paletteer_c("grDevices::Temps", length(Sample_list_FL))
colors_FL <- paletteer_d("ggsci::default_igv", length(Sample_list_FL))


# Associer chaque échantillon à une couleur
sample_colors_FL <- setNames(colors_FL, Sample_list_FL)

####DLBCL_tFL####
data_talyies_tFL_DLBCL <- data_talyies_full %>% 
  filter(Disease %in% c("tFL","DLBCL") & Origin %in% c("PBMC","LN"))

Sample_list_tFL_DLBCL <- data_talyies_tFL_DLBCL %>% 
  distinct(Sample_code) %>% pull(Sample_code)

# Générer une palette de couleurs
colors_tFL_DLBCL <- paletteer_d("ggthemes::Tableau_20", length(Sample_list_tFL_DLBCL))
# Associer chaque échantillon à une couleur
sample_colors_tFL_DLBCL <- setNames(colors_tFL_DLBCL, Sample_list_tFL_DLBCL)

##Combine les 2 palettes###
sample_colors_all <- c(sample_colors_FL,sample_colors_tFL_DLBCL)

####Object Seurat T cells ####
# load(file = "Object_R/Seurat_obj_T_annoted.Robj")

###Regating Tfh
# tfh_Tcells_subset_joined <- subset(Tcells, (CD4 > 0.5 & CD8A < 0.5 & CXCR5 >0.5 & ICOS > 0.5 & IL2RA < 0.5))
# Sample_metadata_tfh <- data.frame(cell_id = rownames(tfh_Tcells_subset_joined@meta.data), Sample_code = tfh_Tcells_subset_joined@meta.data$Sample_code) %>% 
#   mutate(True_Tfh = TRUE) %>% 
#   mutate(cell_id = gsub("[.]", "_", cell_id)) 


Tcells_Sample_metadata <- read_csv('Metadata/Tcells_Sample_metadata.csv') %>%   mutate(
  Screening = as.logical(Screening),
  ScRNA_flex = as.logical(ScRNA_flex),
  ScRNA_cite = as.logical(ScRNA_cite),
  D_R = as.factor(D_R),
  Origin = factor(Origin, levels = c("PBMC","LN","BM","Spleen")),
  Disease = as.factor(Disease),
  Response_type_CD10 = factor(Response_type_CD10, levels = c("Low", "Medium", "High")),
  T_clusters = as.factor(T_clusters),
  Response_type_total = factor(Response_type_total, levels = c("Low", "Medium", "High"))) 
  # select(-True_Tfh) %>% 
  # left_join(Sample_metadata_tfh) %>%
  # mutate(True_Tfh = ifelse(is.na(True_Tfh), FALSE, True_Tfh))
# select(-True_Treg) %>% 
# left_join(Sample_metadata_treg) %>%
#    mutate(True_Treg = ifelse(is.na(True_Treg), FALSE, True_Treg))

dataScore_signature <-read_tsv("Score/T_cells_Tcells_Gene_Data_Expression.tsv") %>%  dplyr::rename(cell_id = id) %>% 
  mutate(cell_id = gsub("\\.", "_", cell_id)) %>% 
  left_join(Tcells_Sample_metadata)


##### Creation table pour Correlation ####
# data_plot_population_for_join <- data_plot_population %>%
#   pivot_wider(
#     names_from = Cell_type,
#     values_from = mean_value,
#     values_fill = list(mean_value = 0)
#   ) %>%
#   mutate(across(everything(), ~ ifelse(is.nan(.), 0, .))) %>%
#    dplyr::rename(
#     SC_Total_Cell = Total,
#     SC_Per_Bcells = Bcells,
#     SC_Per_CD4 = CD4,
#     SC_Per_CD8 = CD8,
#     SC_Per_Myeloid = Myeloid,
#     SC_Per_NK = NK,
#     SC_Per_Other = Other
#   ) %>% 
#   select(-Sample_code_paper)
# 
# 
# data_Plot_cell_pop_CD4_compared_for_join <-  data_Plot_cell_pop_CD4_compared  %>% 
#    dplyr::rename(SC_Per = percentage,
#          Sc_Nb_CD4_total = counttotal,
#          Sc_Nb = count) %>% 
#   pivot_wider(
#     names_from = Cluster_name,
#     values_from = c(SC_Per,Sc_Nb),
#     values_fill = 0
#   ) %>%
#   mutate(across(everything(), ~ ifelse(is.nan(.), 0, .)))
# 
# 
# data_Plot_cell_pop_CD8_compared_for_join <- grouped_data_cluster %>% 
#    dplyr::rename(SC_Per = percentage,
#          Sc_Nb_CD8_total = counttotal,
#          Sc_Nb = count) %>% 
#   pivot_wider(
#     names_from = Cluster_name,
#     values_from = c(SC_Per,Sc_Nb),
#     values_fill = 0
#   ) %>%
#   mutate(across(everything(), ~ ifelse(is.nan(.), 0, .)))
# 
# data_table_SC_population <- data_plot_population_for_join %>% 
#   left_join(data_Plot_cell_pop_CD4_compared_for_join) %>% 
#   left_join(data_Plot_cell_pop_CD8_compared_for_join) %>% 
#   left_join(dataScore_FL20_markers) %>% 
#   relocate(c(Response_type_total,B_cell_depletion_total_glofi,B_Origin),.after = Sample_code)
# 
# # write_csv(data_table_SC_population,"data_table_SC_population.csv")


theme_custom <- function() {
  theme_classic() +  # Fond classique et épuré
    theme(
      axis.text.x = element_text(size = 9, face = "bold", color = "black"),  # Labels des axes X
      axis.text.y = element_text(size = 9, face = "bold", color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold", color = "black", margin = margin(r = 10)),  # Labels des axes X
      legend.title = element_text(size = 8, face = "bold"),  # Titre de la légende
      legend.text = element_text(size = 6,face="bold"),  # Texte de la légende
      plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold"),  # Titre du graphique centré
      strip.placement = "outside",  # Strips à l'extérieur
      strip.background = element_blank(),  # Suppression du cadre autour des strips
      strip.text = element_text(size = 10, face = "bold", angle = 0)  # Style des textes des strips
    )
}

theme_blood <- function() {
  theme_classic() +  # Fond classique et épuré
    theme(
      axis.text.x = element_text(size = 6, face = "bold", color = "black"),  # Labels des axes X
      axis.text.y = element_text(size = 6, face = "bold", color = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold", color = "black"),  # Labels des axes X
      legend.title = element_text(size = 7, face = "bold"),  # Titre de la légende
      legend.text = element_text(size = 6,face="bold"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.size = unit(0.4, "cm"),  
      plot.margin = margin(t = 0.6,  # Top margin
                           r = 0.5,  # Right margin
                           b = 0,  # Bottom margin
                           l = 1), # Left margin 
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),  # Titre du graphique centré
      strip.placement = "outside",  # Strips à l'extérieur
      strip.background = element_blank(),  # Suppression du cadre autour des strips
      strip.text = element_text(size = 6, face = "bold", angle = 0)
    )
}
theme_blood_tfh <- function() {
  theme_classic() +  # Fond classique et épuré
    theme(
      axis.text.x = element_text(size = 6, face = "bold", color = "black"),  # Labels des axes X
      axis.text.y = element_text(size = 6, face = "bold", color = "black"),
      axis.title = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold", color = "black"),  # Labels des axes X
      legend.title = element_text(size = 7, face = "bold"),  # Titre de la légende
      legend.text = element_text(size = 6,face="bold"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.size = unit(0.4, "cm"),  
      plot.margin = margin(t = 0,  # Top margin
                           r = 0.5,  # Right margin
                           b = 0,  # Bottom margin
                           l = 1), # Left margin 
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),  # Titre du graphique centré
      strip.placement = "outside",  # Strips à l'extérieur
      strip.background = element_blank(),  # Suppression du cadre autour des strips
      strip.text = element_text(size = 6, face = "bold", angle = 0)
    )
}
