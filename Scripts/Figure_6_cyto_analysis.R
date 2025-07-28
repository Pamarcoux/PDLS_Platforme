source('~/Postdoc/R/GlofiResistance/00_Metadata_Talyies.R')
library(flowCore)
library(data.table)
library(impute)
library(pheatmap)
library(scales)
library(ggplotify)



data_metacluster_responder <- read.csv('~/Postdoc/R/Cyto_analysis/Data/Talyies/Metacluster_in_Responder_group.csv',
                                       sep=",") 

####Plot tsne Cytometry #####
# Set your folder
fcs_path <- "~/Postdoc/R/Cyto_analysis/Data/Talyies/FCS_files"
fcs_files <- list.files(fcs_path, pattern = ".fcs$", full.names = TRUE)

# Load and merge
{merged_data_facs <- lapply(fcs_files, function(file) {
  fcs <- read.FCS(file, transformation = FALSE)
  exprs_data <- as.data.frame(exprs(fcs))
  
  # Add sample name
  exprs_data$Sample <- tools::file_path_sans_ext(basename(file))
  return(exprs_data)
}) %>% bind_rows()

# Add metadata
merged_data_facs <- merged_data_facs %>%
  rename("CD3" ='Alexa Fluor 700-A',
         "CD4"='APC-Cy7-A',
         "CD8"="BV421-A",
         "LAG3"="PE-A",
         "TIM3"="PECF 594-A",
         "PD1"="BV605-A",
         "TIGIT"="PerCP-Cy5-5-A",
         "41BB"="APC-A",
         "CD28"="FITC-A",
         "Aqua"="BV510-A"
         )

merged_data_facs$FlowSOM_metacluster_id <- as.factor(merged_data_facs$FlowSOM_metacluster_id)

merged_data_facs <- merged_data_facs %>% 
  mutate(FlowSOM_metacluster_id= recode(FlowSOM_metacluster_id,
         `8` = "CD4 non-exhausted",
         `5` = "CD4 exhausted",
         `2` = "CD8 non-exhausted",
         `3` = "CD8 exhausted",
         '6' = "CD4 exhausted"
  )) %>% 
  dplyr::filter(FlowSOM_metacluster_id %in% c("CD4 non-exhausted", "CD4 exhausted", "CD8 non-exhausted", "CD8 exhausted"))


merged_data_facs$FlowSOM_metacluster_id <- factor(
  merged_data_facs$FlowSOM_metacluster_id,
  levels = c("CD4 non-exhausted", "CD4 exhausted", "CD8 non-exhausted", "CD8 exhausted")
)
}
#####Plot Markers#####
merged_data_facs_pivot <- merged_data_facs %>% 
  select(-CD3,-Aqua) %>% 
  pivot_longer(cols = c("CD4", "CD8", "LAG3", "TIM3", "PD1", "TIGIT", "41BB", "CD28"),
               names_to = "Marker",
               values_to = "Expression") %>%
  mutate(Marker = factor(Marker, levels = c("CD4", "CD8","41BB", "CD28", "LAG3", "TIM3", "PD1", "TIGIT")))

# Define the marker levels
markers <- levels(merged_data_facs_pivot$Marker)

# Define limits per marker (match factor names exactly, including space after "CD8 ")
color_limits <- list(
  "CD4" = c(-300, 33000),
  "CD8" = c(-100, 25000),
  "41BB" = c(-250, 5500),
  "CD28" = c(-800, 828000),
  "LAG3" = c(-700, 3100),
  "TIM3" = c(-700, 27000),
  "PD1" = c(-250, NA),
  "TIGIT" = c(-400, 17300)
)

# Generate one plot per marker with independent color scale
{plot_list <- lapply(markers, function(marker) {
  df_marker <- dplyr::filter(merged_data_facs_pivot, Marker == marker)
  
  ggplot(df_marker, aes(x = opt_SNE_1, y = opt_SNE_2)) +
    geom_point(aes(color = Expression, size = Expression), size = 0.3, alpha = 0.8) +
    scale_color_viridis_c(trans = "log10", na.value = "cadetblue2",option= "turbo",
                          limits = color_limits[[marker]]) +
    ggtitle(marker) +
    theme_blood() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      legend.justification = c("right", "bottom"),
      legend.box.just = "right",
      legend.key.height = unit(0.8, 'cm'),
      legend.key.width = unit(0.3, 'cm'),
      legend.position = "right",
    ) +
    labs(color = "", x = "tSNE 1", y = "tSNE 2")
})

# Combine the plots in a grid layout (e.g. 3 columns)
  # On définit les lignes : 1 à 4 (ligne 1), 5 à 8 (ligne 2)
  p1 <- wrap_plots(plot_list[1:4], ncol = 4)
  p2 <- wrap_plots(plot_list[5:8], ncol = 4)
  
  # Combine les deux avec un espace entre eux
  tsne_markers_plot <- p1 / plot_spacer() / p2 +
    plot_layout(heights = c(1, 0.02, 1))  # 0.1 pour l'espace entre les lignes
tsne_markers_plot 
}
#####Plot Clusters####



# Generate one plot per marker with independent color scale
ncell <- length(merged_data_facs$`SSC-A`)

tsnep <- ggplot(merged_data_facs, aes(x = opt_SNE_1, y = opt_SNE_2)) +
  geom_point(aes(color = FlowSOM_metacluster_id), size = 0.3, alpha = 0.8) +
  # scale_color_manual(values = palette_colors) +
  theme_blood() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    legend.justification = c("right"),
    legend.box.just = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.15, "cm"), 
    legend.position = "none"
  ) +
  labs(color = "Metaclusters", x = "optSNE 1", y = "optSNE 2")+
  guides(color = guide_legend(nrow = 9, override.aes = list(size = 3, alpha = 1)))

tsne_clusters_plot <- tsnep + annotation_custom(
  grob = textGrob(
    label = sprintf("n = %d", ncell),
    x = unit(0.85, "npc"), y = unit(0.03, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 8, col = "black", fontface = "italic")))
 


####HeatMap Cluster Talyies####
marker_order = c("CD3","CD4","CD8","PD1","LAG3","TIGIT","TIM3")


data_heatmap_cluster <- read.csv('~/Postdoc/R/Cyto_analysis/Data/Talyies/Metacluster_Box_Plots_Means.csv',sep=";") %>% 
  rename(Markers = Stats.channel,
         Mean = Mean.of.channel) %>% 
  dplyr::filter(Markers %in% c("CD3","CD4","CD8","PD1","LAG3","TIGIT","TIM3")) %>% 
  mutate(Markers = factor(Markers,
                           levels = rev(marker_order)),
         Population = factor(Population,
                               levels = c("CD4 non-exhausted","CD4 exhausted",
                                          "CD8 non-exhausted","CD8 exhausted")))



data_heatmap_cluster_plot <- data_heatmap_cluster %>%
  group_by(Markers, Population) %>%
  summarise(Mean_marker = mean(Mean), .groups = "drop") %>%
  group_by(Markers) %>% 
  mutate(zscore = scale(Mean_marker)) %>%  # z-score par phenotype (option 1)
  ungroup() %>% 
  mutate(
    zscore_clipped = pmin(pmax(zscore, -1), 1)  # Tronquer entre -2 et 2
  ) %>% 
 dplyr::filter(Population %in% c("CD4 exhausted","CD4 non-exhausted","CD8 exhausted","CD8 non-exhausted"))

heatmap_cluster_plot <- ggplot(data_heatmap_cluster_plot, aes(x = Population , y = Markers, fill = zscore_clipped)) +
  geom_tile(color = "white",width = 0.98, height = 0.98) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1,    # "BlYlRd" n’existe pas, on inverse RdYlBu
                       na.value = "grey90", 
                       name = "Z-score") +
  labs(
    title = "",
    x = "Clusters",
    y = "Markers",
    fill = "Z-score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )+theme_blood()


####Heatmap Response Group pHEATMAP####
sample_list <- data_talyies_full %>% 
  dplyr::filter(Day == "D6", Disease %in% c("FL","DLBCL","tFL")) %>% 
  distinct(Sample_code) 

table_treatment <- read_csv('Paper_platforme/liste_combo.csv') %>%
  mutate(Treatment = apply(.[c("A", "B", "C")], 1, function(row) {
    paste(na.omit(row), collapse = " + ")
  })) %>% 
  mutate(Treatment_reorder = apply(.[c("E", "F", "G")], 1, function(row) {
    paste(na.omit(row), collapse = " + ")
  })) %>% 
  select(Treatment_reorder,Treatment,Treatment_type,Treatment_Cat) %>% 
  # dplyr::filter(!grepl("TCB 10 nM", Treatment_reorder)) %>% #dplyr::filter les Treatments
  dplyr::filter(!grepl("GA101", Treatment_reorder) & !grepl("IL2v", Treatment_reorder) & !grepl("TCB 10 nM", Treatment_reorder) 
                &  !grepl("ZB2", Treatment_reorder) ) %>% #dplyr::filter les Treatments
  mutate(Treatment_reorder = factor(Treatment_reorder, levels = sort(unique(Treatment_reorder)))) 

all_combinations <- crossing(sample_list$Sample_code, table_treatment$Treatment_reorder) %>% 
  rename(Sample_code = 'sample_list$Sample_code', Treatment_reorder = 'table_treatment$Treatment_reorder') %>% 
  left_join(table_treatment,relationship = "many-to-many") %>% 
  left_join(data_sample_info_complete) %>% 
  dplyr::filter(Disease %in% c("FL","DLBCL","tFL") & Screening == TRUE) %>% 
  select(Treatment,Sample_code,Treatment_reorder,Treatment_type,Treatment_Cat) %>% 
  mutate(Treatment_Cat = factor(Treatment_Cat, levels = (c("UT","αCD20-TCB","Inhibiteur_CP","Co_Activator","ADC"))))



cal_z_score <- function(x){
  (x)}

{  
  # filename = "Paper_platforme/Figure/Figure_5_Heatmap/Clusters_by_Treatment_type.png"
  data_treatment_cat <- data_talyies_full %>%
    mutate(B_cell_depletion_total = replace(B_cell_depletion_total, B_cell_depletion_total == "", NA)) %>%
    dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
    right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
    dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>% 
    distinct() %>% 
    dplyr::filter(Treatment_type == "Single" & Treatment_Cat != "UT")
  
  #Annotation Samples
  Metadata_treatment_matrix <- data_treatment_cat %>%
    distinct(Sample_code, .keep_all=TRUE) %>%
    select(Sample_code,Origin,Disease)
  
  df_Metadata_samples_annotations <- as.data.frame(Metadata_treatment_matrix[,-1])
  rownames(df_Metadata_samples_annotations) <- Metadata_treatment_matrix$Sample_code
  
  #Transformation en Matrice
  data_treatment_matrix <- data_treatment_cat %>% select(Sample_code,Treatment_Cat,B_cell_depletion_total) %>%
    dplyr::filter(B_cell_depletion_total > -50 & !Sample_code %in% c("FL3_PB1","FL2_PB1")) %>% 
    group_by(Treatment_Cat,Sample_code) %>% 
    summarize(B_cell_depletion_total = mean(B_cell_depletion_total)) %>% 
    mutate(Treatment_Cat = recode(Treatment_Cat,
                                  Inhibiteur_CP = "ICP Inhibitors",
                                  Co_Activator = "Co Activators"))
    
  
  matrix <- pivot_wider(data_treatment_matrix,
                        names_from = Treatment_Cat,
                        values_from = B_cell_depletion_total) %>% select(-Sample_code)
  
  rownames(matrix) <- unique(data_treatment_matrix$Sample_code)
  Mat_dataExpression_heatmap <- as.matrix(matrix)
  
  ##Voir les NA ##
  na_counts_by_treatment <- colSums(is.na(Mat_dataExpression_heatmap))
  print(na_counts_by_treatment)
  
  na_percentage <- colSums(is.na(Mat_dataExpression_heatmap)) / nrow(Mat_dataExpression_heatmap) * 100
  print(na_percentage)
  
  #### Imputation des datas avec methode KNN###
  imputed_result <- impute.knn(Mat_dataExpression_heatmap, k=5)
  Mat_dataExpression_heatmap_imputed <- imputed_result$data
  
  #Scaling the rows
  Mat_dataExpression_heatmap_norm<- (apply(Mat_dataExpression_heatmap_imputed, 1, cal_z_score))
  
  #Heatmap
  my_colour = list(
    Origin = c(LN = "pink", PBMC = "lightblue"))
  
  par(mar = c(6, 6, 4, 6))  # bottom, left, top, right
  
  resultat_heatmap <- pheatmap(Mat_dataExpression_heatmap_norm,
                               annotation_col = df_Metadata_samples_annotations,
                               # annotation_row = data_gene,
                               legend = TRUE,
                               legend_labels = "Z score",
                               clustering_distance_rows = "euclidean",
                               clustering_distance_cols = "euclidean",
                               clustering_method = "complete",
                               cluster_cols = TRUE,
                               cluster_rows = FALSE,
                               annotation_colors = my_colour,
                               cutree_rows = 2,
                               cutree_cols = 3,
                               angle_col = 90,
                               show_rownames = TRUE,
                               show_colnames = TRUE,  
                               main = "",
                               
                               width = 24,
                               height = 7,
                               display_numbers = TRUE,
                               silent = TRUE)
  dev.off()
  pheatmap_gg <- as_ggplot(resultat_heatmap$gtable)
  k <- 3  
  clusters <- cutree(resultat_heatmap$tree_col, k)
  df_Cluster_response_cat <- data.frame(
    Sample_code = names(clusters),
    Cluster_response_cat = unlist(clusters), row.names = NULL) %>% 
    mutate(Cluster_response_cat = recode(Cluster_response_cat, 
                                         `1` = "Medium Responders", 
                                         `2` = "Low Responders", 
                                         `3` = "High Responders")) 
  
  pheatmap_gg_anot <- pheatmap_gg +
      theme(plot.margin = margin(b = 25)) 

}

  line1 <- linesGrob(
  x = unit(c(0.009, 0.161), "npc"), 
  y = unit(c(0.07, 0.07), "npc"),
  gp = gpar(col = "black", lwd = 3.5)
)

  line2 <- linesGrob(
    x = unit(c(0.168, 0.457), "npc"), 
    y = unit(c(0.07, 0.07), "npc"),
    gp = gpar(col = "black", lwd = 3.5)
  )

line3 <- linesGrob(
  x = unit(c(0.463, 0.83), "npc"), 
  y = unit(c(0.07, 0.07), "npc"),
  gp = gpar(col = "black", lwd = 3.5)
)

# Add vertical text labels
text1 <- textGrob(
  label = "High Responders", rot = 0,
  x = unit(0.09, "npc"), y = unit(0.04, "npc"),
  gp = gpar(col = "black", fontface = "bold", fontsize = 10)
)

text2 <- textGrob(
  label = "Medium Responders", rot = 0,
  x = unit(0.32, "npc"), y = unit(0.04, "npc"),
  gp = gpar(col = "black", fontface = "bold", fontsize = 10)
)

text3 <- textGrob(
  label = "Low Responders", rot = 0,
  x = unit(0.65, "npc"), y = unit(0.04, "npc"),
  gp = gpar(col = "black", fontface = "bold", fontsize = 10)
)

# Combine plot and custom grobs using annotation_custom (trick using ggplotify/grid)

pheatmap_gg_anot <- as.ggplot(~{
  grid.draw(ggplotGrob(pheatmap_gg_anot))
  grid.draw(line1)
  grid.draw(line2)
  grid.draw(line3)
  grid.draw(text1)
  grid.draw(text2)
  grid.draw(text3)
})

####Bar Plot repartition responder cluster ####
data_metacluster_responder <- read.csv('~/Postdoc/R/Cyto_analysis/Data/Talyies/Metacluster_in_Responder_group.csv',
                                       sep=",") 
#Rename and prepare 
data_metacluster_responder <- data_metacluster_responder %>% 
  rename(FlowSOM_metacluster_id = Metacluster.ID,
         Responder_Cluster = clustering.group. ,
         Percentage = X..of.cells) %>% 
  mutate(FlowSOM_metacluster_id= recode(FlowSOM_metacluster_id,
                                        `8` = "CD4 non-exhausted",
                                        `5` = "CD4 exhausted",
                                        `2` = "CD8 non-exhausted",
                                        `3` = "CD8 exhausted",
                                        '6' = "CD4 non-exhausted",
                                        '4' = "NK like cells",
  )) %>% 
  dplyr::filter(FlowSOM_metacluster_id %in% c("CD4 non-exhausted", "CD4 exhausted", "CD8 non-exhausted", "CD8 exhausted")) %>%
  mutate(Responder_Cluster= recode(Responder_Cluster,
                                   `ADC-TCB MED` = "Medium\n Responders", 
                                   `GLOBAL LOW` = "Low\n Responders", 
                                   `TCB HIGH` = "High\n Responders")) %>%
  mutate(Responder_Cluster = factor(Responder_Cluster, 
                                levels = c("Low\n Responders", "Medium\n Responders", "High\n Responders")),
         FlowSOM_metacluster_id = factor(FlowSOM_metacluster_id, 
                                         levels = rev(c("CD4 non-exhausted", "CD4 exhausted", "CD8 non-exhausted", "CD8 exhausted","NK like cells","CD4 Monocytes exhausted")))) %>%
  group_by(FlowSOM_metacluster_id,Responder_Cluster) %>% 
  summarize(MeanPer = mean(Percentage),
            SD= sd(Percentage)) %>% 
  group_by(Responder_Cluster) %>%
  mutate(Normalized_Mean = 100 * MeanPer / sum(MeanPer),
         Normalized_SD = 100 * SD / sum(SD)) %>%
  ungroup()


plot_barplot_reparition_responders <- ggplot(data_metacluster_responder, aes(x = Responder_Cluster, y = Normalized, fill = FlowSOM_metacluster_id)) +
  geom_bar(stat = "identity", color = "black", size = 0.25) + 
  scale_fill_manual(values = rev(hue_pal()(4)))+
  theme_blood()+
  labs(
    title = "",
    x = "",  # Éviter de répéter des informations évidentes
    y = "Proportion (%)",
    fill = "Metaclusters"
  ) +
  theme(legend.position = "left",
        legend.box.margin = margin(0, 10, 0, 0), # Remove margin around legend box
        legend.margin = margin(0.2, 10, 0.2, 0.2))


####Montage####
{empty_plot <- ggplot() + theme_void()
A_cyto <- plot_grid(
  pheatmap_gg_anot, 
  ncol = 1,
  labels = c("A"))
  # vjust = ,)


B_cyto <- plot_grid(
  tsne_markers_plot, 
  ncol = 1,
  labels = c("B"))
  # vjust = ,)

C_D_E <- plot_grid(heatmap_cluster_plot,tsne_clusters_plot,plot_barplot_reparition_responders,
                 rel_widths = c(1.1, 1,0.9),
                 ncol = 3, 
                 labels = c("C","D","E"),
                 axis = "tb",
                 align = "h")

Figure_cyto_total <- plot_grid(
  plot_grid(A_cyto, ncol = 1, labels = c("")),
  plot_grid(B_cyto, ncol = 1, labels = c("")),
  plot_grid(C_D_E, ncol = 1, labels = c("")),
  nrow=3 ,
  rel_heights = c(1, 1,0.8),
  axis = "bt",
  align = "hv")

print(Figure_cyto_total)

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure6/Figure_6_cyto.png",
  plot = Figure_cyto_total,
  device = "png",
  width = 35,        # largeur A4 en cm
  height = 38 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
}


