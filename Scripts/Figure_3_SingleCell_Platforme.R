library(paletteer)

##### Figure 2
###### Filter data ########
sample_colors_LN_PBMC <- setNames(c("#a00000","#1a80bb"), c("LN","PBMC"))
sample_colors_all <- c(sample_colors_FL,sample_colors_tFL_DLBCL)

data_talyies_screened_flex <- data_talyies_full %>% 
  filter(Screening == T & ScRNA_flex == T & Disease %in% c("FL","tFL","DLBCL") ) 

data_talyies_flex <- data_talyies_full %>% 
  filter(ScRNA_flex == T & Disease %in% c("FL","tFL","DLBCL") ) 

Sample_list_screened_flex <- data_talyies_screened_flex %>% 
  distinct(Sample_code) %>% pull(Sample_code)

Sample_list_flex <- data_talyies_flex %>% 
  distinct(Sample_code) %>% pull(Sample_code)

setdiff(Sample_list_flex,Sample_list_screened_flex)

# load(file = "AllSample_v3.Robj")
Sample_list_all <- gsub("_",".",Sample_list_flex)
AllSample_subset <- subset(AllSample, orig.ident %in% c(Sample_list_all,"FL5.PB1","FL6.PB1"))

{Metadata_All <- AllSample_subset@meta.data %>% select(orig.ident,integrated_snn_res.0.3,type,Lymphomatype) %>% 
  mutate(Sample_code = gsub("\\.", "_", orig.ident)) %>% 
  mutate(cell_id = gsub("\\.", "_", rownames(.))) %>% 
    mutate(Sample_code = ifelse(Sample_code == "FL5_PB1","FL5_PB1_3",Sample_code),
           Sample_code = ifelse(Sample_code == "FL6_PB1","FL6_PB1_2",Sample_code)) %>%
  mutate(Sample_code_paper = Sample_code,
         Clusters = integrated_snn_res.0.3) %>% 
  mutate(Sample_code_paper = recode(Sample_code_paper, !!!Replacement_table), 
         Clusters = recode(Clusters, !!!Replacement_table_Sc_all)) %>% 
  mutate(Clusters = factor(Clusters, levels = Replacement_table_Sc_all)) %>% 
  mutate(Sample_code_paper = factor(Sample_code_paper, levels = Replacement_table),
         Sample_code = factor(Sample_code, levels = Sample_order)) %>% 
    mutate(type = if_else(type == "PB", "PBMC", type)) %>% 
  mutate(Cell_type = "Other",
    Cell_type = (if_else(grepl("CD4", Clusters) == TRUE, "CD4", "Other")),
         Cell_type = (if_else(grepl("CD8", Clusters) == TRUE, "CD8", Cell_type)),
         Cell_type = (if_else(grepl("NK", Clusters) == TRUE, "NK", Cell_type)),
         Cell_type = (if_else(grepl("Bcells", integrated_snn_res.0.3) == TRUE, "Bcells", Cell_type)),
         # Cell_type = (if_else(grepl("Healthy", Clusters) == TRUE, "Normal B Cells", Cell_type)),
         Cell_type = (if_else(grepl("myeloid", integrated_snn_res.0.3) == TRUE, "Myeloid", Cell_type))) %>% 
  mutate(Cell_type = factor(Cell_type, levels = c("Bcells","Myeloid", "CD4", "CD8", "NK","Other"))) %>% 
    left_join(data_response_type_CD10)
         
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$orig.ident,col.name = "orig.ident")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Cell_type,col.name = "Cell_type")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Sample_code_paper,col.name = "Sample_code_paper")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Sample_code,col.name = "Sample_code")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Clusters,col.name = "Clusters")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Response_type_total,col.name = "Response_type_total")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$type,col.name = "Origin")
AllSample_subset <- AddMetaData(AllSample_subset, metadata = Metadata_All$Lymphomatype,col.name = "Disease")


}

# AllSample_subset <- AllSample


# AllSample_subset <- AddMetaData(AllSample_subset, metadata = orig.ident,col.name = "orig.ident")


####Umap Celltype####
full_umap_celltype <- DimPlot(AllSample_subset, reduction = "umap_all", label = FALSE, pt.size = 0.4, raster = FALSE, group.by = 'Cell_type', alpha =0.6) +
  scale_color_manual(values = paletteer_d("ggthemes::Classic_Green_Orange_6"))+
  theme_void() +  # Completely removes axes, background, and ticks
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
    legend.title = element_text(size = 12, face = "bold"),  # Legend title styling
    legend.text = element_text(size = 10, face= "bold"),  # Legend text size
    panel.spacing = unit(1, "lines"),  # Space between panels when faceting
    strip.text = element_text(size = 14, face = "bold"),  # Styling facet labels
    # plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Adjust margins
    legend.position = "bottom",  # Hide the legend
    legend.key.size = unit(0.5, 'cm')) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 3)))+
  labs(title = "")  # Add a descriptive title

ggsave(plot = full_umap_celltype, filename = "./Paper_platforme/Figure/Figure2/full_umap.png", device = "png",width = 9, height = 8 , units = "in", dpi = 300)

#####Umap by Sample #####
full_umap_Sample <- DimPlot(AllSample_subset, reduction = "umap_all", label = FALSE, pt.size = 0.4, raster = FALSE, group.by = 'Sample_code', alpha =0.6) +
  scale_color_manual(values = sample_colors_all)+
  theme_void() +  # Completely removes axes, background, and ticks
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
    legend.title = element_text(size = 12, face = "bold"),  # Legend title styling
    legend.text = element_text(size = 10, face= "bold"),  # Legend text size
    panel.spacing = unit(1, "lines"),  # Space between panels when faceting
    strip.text = element_text(size = 14, face = "bold"),  # Styling facet labels
    # plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Adjust margins
    legend.position = "bottom",  # Hide the legend
    legend.key.size = unit(0.5, 'cm')
  ) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 3)))+
  labs(title = "")  # Add a descriptive title

#####Repartion by Sample #####
Sample_repartition_plot <- pochi::MetaDataPlot(AllSample_subset, group.by = "Cell_type", split.by = "Sample_code", as.freq=T)+
  labs(
    title = "",
    x = "",  # Éviter de répéter des informations évidentes
    y = "Frequency",
    fill = "Cell Type"
  ) +
  labs(fill = "Cell Type")+
  scale_fill_manual(values = paletteer_d("ggthemes::Classic_Green_Orange_6"))+
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"))  # Inclinaison des labels

#####Umap by Origin #####
full_umap_Origin <- DimPlot(AllSample_subset, reduction = "umap_all", label = FALSE, pt.size = 0.4, raster = FALSE, group.by = 'Origin', alpha =0.6) +
  scale_color_manual(values = sample_colors_LN_PBMC)+
  theme_void() +  # Completely removes axes, background, and ticks
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
    legend.title = element_text(size = 12, face = "bold"),  # Legend title styling
    legend.text = element_text(size = 10, face= "bold"),  # Legend text size
    panel.spacing = unit(1, "lines"),  # Space between panels when faceting
    strip.text = element_text(size = 14, face = "bold"),  # Styling facet labels
    # plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Adjust margins
    legend.position = "bottom",  # Hide the legend
    legend.key.size = unit(0.5, 'cm')
  ) +
  labs(title = "")  # Add a descriptive title

#####Repartition by Orign #####
Origin_repartition_plot <- pochi::MetaDataPlot(AllSample_subset, group.by = "Cell_type", split.by = "Origin", as.freq=T)+
  labs(
    title = "",
    x = "",  # Éviter de répéter des informations évidentes
    y = "Frequency",
    fill = "Cell Type"
  ) +
  labs(fill = "Cell Type")+
  scale_fill_manual(values = paletteer_d("ggthemes::Classic_Green_Orange_6"))+
  theme_custom()

######Correlation Pop FACS vs SingleCell ######
{graph_list <- list()
  data_facs_population <- data_talyies_flex %>% select(Sample_code_paper,Sample_code,Day,Per_CD22_total,Per_CD11b,Per_CD4,Per_CD8,Per_NK,Per_CD3,Per_CD20,Treatment) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Preparation data
data_facs_population$variable <- factor(data_facs_population$variable, levels = c("CD22","CD20", "CD11b","CD3","CD4", "CD8", "NK", "TGD","NKT"))

data_population_single_cell<- read_csv("Metadata//Repartition_cells.csv") %>% 
  select(-Response_type_CD10,-B_cell_depletion_CD22_CD10_plus_glofi,Sample_code) %>% 
  mutate(Sample_code = ifelse(Sample_code == "FL5_PB1","FL5_PB1_3",Sample_code),
         Sample_code = ifelse(Sample_code == "FL6_PB1","FL6_PB1_2",Sample_code))


data_plot_population <- data_population_single_cell %>%
  mutate(Total = rowSums(across(c(CD8, CD4, Other, NK, Bcells, Myeloid)), na.rm = TRUE)) %>%
  mutate(across(c(CD8, CD4, Other, NK, Bcells, Myeloid)) / Total *100) %>%
  pivot_longer(cols = c(CD8, CD4, Other, NK, Bcells, Myeloid),
               values_to = "Percentage",
               names_to = "Cell_type") %>% 
  group_by(Cell_type,Sample_code,Total) %>% 
  summarize(Percentage_Single_cells = mean(Percentage, na.rm = TRUE),, .groups = "drop")

data_plot_population$Cell_type <- factor(data_plot_population$Cell_type, levels = c("Bcells", "Myeloid","CD4", "CD8", "NK", "TGD","Other"))

list_pop <- setdiff(unique(data_plot_population$Cell_type), "Other")


for (i in (1:(length(list_pop)))) {
  Population <- list_pop[i]
  population_name <- paste0(Population, "+ cells")
  
  
  ###Pop Facs
  data_facs_population_filter <-  data_facs_population %>% filter(Day == "D0" & !variable %in% c("CD3","CD20")) %>% 
    mutate(variable = case_when(
      variable == "CD22" ~ "Bcells",
      variable == "CD11b" ~ "Myeloid",
      TRUE ~ variable # Conserve les autres valeurs inchangées
    )) %>% 
    group_by(Sample_code) %>% 
    mutate(Proportion = value / sum(value, na.rm= TRUE) * 100) %>% 
    filter(variable == Population) %>% 
    select(-variable,-Treatment,-Day) %>% 
    rename(Percentages_Facs = Proportion)
  
###Plot 
  grouped_data_cluster <- subset (data_plot_population, Cell_type == Population) %>% 
    left_join(data_facs_population_filter) %>% 
    mutate(Sample_code = factor(Sample_code, levels = Sample_order))
    # filter(Sample_code %in% Sample_list)
    
  # Calculer le modèle linéaire pour obtenir R^2
  cor_test <- cor.test(grouped_data_cluster$Percentages_Facs, grouped_data_cluster$Percentage_Single_cells)
  r_value <- cor_test$estimate  # Coefficient de corrélation
  p_value <- cor_test$p.value   # Valeur p
  
  p <- ggplot(data = grouped_data_cluster, aes(x = Percentages_Facs, y = Percentage_Single_cells, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    # geom_label( aes(label =Sample_code))+
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_all) +
    labs(title = "",
         x = paste0("% Phenotype D0"),
      y = paste0(" % Single Cell"),
      color = "Sample Code"  # Légende des couleurs
    ) +
    theme_custom()+
    theme(legend.position = "none") # Ajuster la taille du texte de l'axe y
  
  
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 12, col = "black")
    )
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

Plots_correlation_FACS_SC <- plot_grid(plotlist = graph_list, ncol = 5, align = "v", axis = "tblr")+ labs(title = "FL")
}

###################################
#####Montage#####
A <- plot_grid(full_umap_celltype,full_umap_Sample, ncol=2,labels = c("A","B"))
B <- Sample_repartition_plot
C <- Plots_correlation_FACS_SC

Figure_3_total <- plot_grid(
  plot_grid(A, ncol = 1, labels = c("A")),  # Deuxième ligne
  plot_grid(B, ncol = 1, labels = c("C")),
  plot_grid(C, ncol = 1, labels = c("D")),
  nrow=4 )

print(Figure_3_total)
ggsave(
  filename = "Paper_platforme/Figure/Figure_3_total.png",
  plot = Figure_3_total,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 42 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)

#####Montage Suppl #####
Abis <- full_umap_Origin
Bbis <- Origin_repartition_plot

Figure_4_Suppl<- plot_grid(Abis,Bbis, labels = c("A","B"),ncol=2)

print(Figure_4_Suppl)
ggsave(
  filename = "Paper_platforme/Figure/Figure_4_Suppl.png",
  plot = Figure_4_Suppl,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 15 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)



########################################################################################

