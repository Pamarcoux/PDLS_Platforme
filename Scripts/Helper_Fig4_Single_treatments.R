{library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggnewscale)
  library(impute)
  library(rstatix)
  library(ggpubr)
  library(grid)
  

  sample_colors_df <- tibble(Sample_code = names(sample_colors_all),
                             sample_color = unname(sample_colors_all))
  
  conditions <- c("B_cell_depletion_total", "Per_CD107a_CD4_normalized", "Per_CD107a_CD8_normalized",
                  "TNFa_normalized", "INFY_normalized", "GrzB_normalized", "IL8_normalized", "IL10_normalized")
  
  Condition = "B_cell_depletion_total"
  
Name <- "All_samples"
Disease_list <- c("FL","tFL","DLBCL")
Origin_list <- c("PBMC","LN")

}

#####Histogram Depletion #####
{Condition <- "Per_CD107a_CD8"
  
  data_treatment <- data_treatment <- data_talyies_full %>%
    select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
    mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
    filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
    right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
    filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
  
  # filter(grepl("TCB 0,01 nM",Treatment_reorder))

# valid_samples <- data_treatment %>%
#   group_by(Sample_code) %>%
#   filter(!is.na(!!sym(Condition))) %>%
#   summarize(treatments_count = n_distinct(Treatment_reorder)) %>%
#   filter(treatments_count == 2) %>%
#   pull(Sample_code)

data_treatment_plot <- data_treatment 
# %>% 
#   filter(Sample_code %in% valid_samples)

ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, colour = "black", size = 0.5) +
  # geom_line(aes(group = Sample_code),
  #           position = position_dodge(0),
  #           color = "grey60",
  #           size = 0.2) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.8, alpha = 0.9) +
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
                     label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
  scale_fill_manual(values = sample_colors_all) + 
  theme_custom()+
  labs(
    title = paste0(Condition),  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = paste0(Condition),
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))
}


###HeapMap ####
#### 8 Heatmaps ####
{data_heatmap <- data_talyies_full %>% 
  select(Origin,Treatment, Disease, Sample_code, Day, all_of(conditions)) %>%
  filter(Day == "D6", Disease %in% Disease_list) %>%
  filter(Origin %in% Origin_list) %>% 
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  filter(Treatment_type == "Single") %>%
  filter(B_cell_depletion_total > -50) %>% 
  filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>% 
  mutate(Per_CD107a_CD8_normalized = ifelse((Treatment == "UT" & Day == "D6" & !is.na(Per_CD107a_CD8_normalized)),1,Per_CD107a_CD8_normalized))


# Compter le nombre de NA dans la colonne Treatment pour chaque Sample_code
No_data_sample_list <- data_heatmap %>%
  group_by(Sample_code) %>%
  summarise(na_count_sample= sum(is.na(B_cell_depletion_total)),
            No_data_sample = ifelse(na_count_sample > (length(unique(data_heatmap$Treatment_reorder))-2), TRUE, FALSE)) %>%
  filter(No_data_sample == TRUE) %>%
  pull(Sample_code)

data_heatmap_pivot <- data_heatmap%>% 
  pivot_longer(cols = all_of(conditions), names_to = "Condition", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Treatment_reorder = factor(Treatment_reorder, levels = rev(levels(factor(Treatment_reorder))))) %>% 
  left_join(sample_colors_df, by = "Sample_code") 

complete_combinations <- expand.grid(
  Sample_code = unique(data_heatmap$Sample_code),
  Treatment_reorder = unique(data_heatmap$Treatment_reorder),
  Condition = unique(data_heatmap_pivot$Condition),
  stringsAsFactors = FALSE)

data_heatmap_pivot_na <- data_heatmap_pivot %>% 
  full_join(complete_combinations) 

# Créer les plots #
plots <- vector("list", length(conditions))
ncol_plot <- 4

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- filter(data_heatmap_pivot_na, Condition == cond) %>% 
    mutate(value_mean = value, na.rm = FALSE)
  
  
  # Heatmap principale
  heatmap_layer <- ggplot(plot_data, aes(x = Sample_code, y = Treatment_reorder)) +
    geom_tile(aes(fill = value_mean), color = "white", size = 0.5) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey80") +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = cond, fill = cond) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom"
    )
  
  # Ajouter la barre colorée des samples
  color_bar <- ggplot(distinct(plot_data, Sample_code, sample_color)) +
    geom_tile(aes(x = Sample_code, y = 1, fill = sample_color), height = 1) +
    scale_fill_identity() +
    theme_void() +
    theme(aspect.ratio = 0.05)
  
  # Empiler le color bar sous la heatmap
  p <- (heatmap_layer + theme(plot.margin = margin(b = 0))) / 
    (color_bar + theme(plot.margin = margin(t = -15)))  
  # Gérer l'affichage des axes y
  col_index <- i %% ncol_plot
  if (col_index != 0) {
    p <- p & theme(axis.text.y = element_blank())
  }
  
  plots[[i]] <- p
}

# Assembler tous les plots #
Plot_heatmap8 <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

#### 8 Heatmap Zscore ####

#Creation Heatmaps
i=0
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- filter(data_heatmap_pivot_na, Condition == cond) %>% 
    filter(!Sample_code %in% No_data_sample_list) %>% 
    group_by(Treatment_reorder) %>%
    mutate(z_score = scale(value), na.rm = TRUE) %>%
    mutate(Sample_code = factor(Sample_code, levels = Sample_order)) %>%
    ungroup()
  
  # Heatmap principale
  heatmap_layer <- ggplot(plot_data, aes(x = Sample_code, y = (Treatment_reorder))) +
    geom_tile(aes(fill = z_score), color = "white", size = 0.5) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey80",limits = c(-2, 3)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = cond, fill = "Z-score") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom"
    )
  
  # Ajouter la barre colorée des samples
  color_bar <- ggplot(distinct(plot_data, Sample_code, sample_color)) +
    geom_tile(aes(x = Sample_code, y = 1, fill = sample_color), height = 1) +
    scale_fill_identity() +
    theme_void() +
    theme(aspect.ratio = 0.05)
  
  # Empiler le color bar sous la heatmap
  p <- (heatmap_layer + theme(plot.margin = margin(b = 0))) / 
    (color_bar + theme(plot.margin = margin(t = -15)))  
  # Gérer l'affichage des axes y
  col_index <- i %% ncol_plot
  if (col_index != 0) {
    p <- p & theme(axis.text.y = element_blank())
  }
  
  plots[[i]] <- p
}

#  Assembler tous les plots
final_plot_heatmapZscore <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

# final_plot_heatmapZscore_1 <- final_plot_heatmapZscore
# final_plot_heatmapZscore_2 
# 
# double_heapmap <- plot_grid(final_plot_heatmapZscore_1,final_plot_heatmapZscore_2, nrow = 2, 
#                            labels = c("A","B"),
#                            rel_heights = c(1, 1))

#### Montage ####
  Figure_4_test <- plot_grid(
    plot_grid(Plot_heatmap8, ncol = 1, labels = c("A")),  
    plot_grid(final_plot_heatmapZscore, ncol = 1, labels = c("B")), 
    nrow =2,
    rel_heights = c(1, 1)) +
  plot_annotation(title = paste0(Name), theme = theme_custom())
  

ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/Heatmaps_",Name,".png"),
       plot = Figure_4_test,
       device = "png",
       width = 29.7,        # largeur A4 en cm
       height = 44 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)
}

#### Correlation ####
Name <- "All_samples"
Disease_list <- c("FL","tFL","DLBCL")
Origin_list <- c("PBMC","LN")



data_correlation <- data_talyies_full %>% 
  filter(Day == "D6", Disease %in% Disease_list) %>%
  filter(Origin %in% Origin_list) %>% 
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  filter(Treatment_type == "Single") %>%
  filter(B_cell_depletion_total > -50) %>% 
  filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))

##### Calcul All Correlation ####
# Sélectionner uniquement les colonnes numériques pour le calcul de la corrélation
numeric_columns <- data_correlation %>%
  select_if(is.numeric)

# Initialiser un dataframe pour stocker les résultats
correlation_results <- data.frame(Variable = character(), Correlation = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Calculer la corrélation et la p-valeur entre B_cell_depletion_total et chaque autre colonne numérique
for (column in names(numeric_columns)) {
  if (column != "B_cell_depletion_total") {
    complete_cases <- complete.cases(numeric_columns[[column]], numeric_columns$B_cell_depletion_total)
    
    if (sum(complete_cases) > 2) { # Vérifier qu'il y a suffisamment de paires de valeurs
      correlation_test <- cor.test(numeric_columns[[column]][complete_cases],
                                   numeric_columns$B_cell_depletion_total[complete_cases])
      correlation_results <- rbind(correlation_results, data.frame(Variable = column, Correlation = correlation_test$estimate, P_Value = correlation_test$p.value, stringsAsFactors = FALSE))
    } else {
      correlation_results <- rbind(correlation_results, data.frame(Variable = column, Correlation = NA, P_Value = NA, stringsAsFactors = FALSE))
    }
  }
}

# Afficher les résultats
print(correlation_results)


##### Plot Correlation #####
condition_list <- correlation_results %>% 
  filter(!is.na(Correlation) & !Variable %in% c("...1","Date_D0")) %>%
  pull(Variable)

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(condition_list))) {
  Condition <- condition_list[i]
  
  population_name <- paste0(Condition)
  
  data_correlation_plot <- data_correlation %>%
    mutate(value = (!!sym(Condition))) %>%
    filter(B_cell_depletion_total > -50) 
  
  # Calculer la corrélation de Pearson entre Area et B_Cell_Depletion
  correlation_result <- cor.test(data_correlation_plot$value, data_correlation_plot$B_cell_depletion_total,use = "pairwise.complete.obs") 
  r_value <- correlation_result$estimate  # Coefficient de corrélation
  p_value <- correlation_result$p.value   # Valeur p
  
  # Créer le graphique
  p <- ggplot(data_correlation_plot, aes(x = B_cell_depletion_total, y = !!sym(Condition),fill)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_all) +
    labs(title = "",
         x = "B_cell_depletion_total (D6)",
         y=paste0(Condition))+
    theme_custom() +
    theme(legend.position = "none")  # Ajuster la taille du texte de l'axe y
  
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D6_Single_1 <- plot_grid(plotlist = graph_list[1:25], ncol = 5,nrow=5, align = "v", axis = "tblr")
plot_correlation_D6_Single_2 <- plot_grid(plotlist = graph_list[26:50], ncol = 5,nrow=5, align = "v", axis = "tblr")
plot_correlation_D6_Single_3 <- plot_grid(plotlist = graph_list[50:57], ncol = 5,nrow=5, align = "v", axis = "tblr")

ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/plot_correlation_D6_Single1.png"),
       plot = plot_correlation_D6_Single_1,
       device = "png",
       width = 29.7,        # largeur A4 en cm
       height = 44 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)

ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/plot_correlation_D6_Single2.png"),
       plot = plot_correlation_D6_Single_2,
       device = "png",
       width = 29.7,        # largeur A4 en cm
       height = 44 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)

ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/plot_correlation_D6_Single3.png"),
       plot = plot_correlation_D6_Single_3,
       device = "png",
       width = 29.7,        # largeur A4 en cm
       height = 44 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)

####Correlation D0 ####
{  # Filtrer les données pour obtenir les valeurs des marqueurs à D0
  markers_d0 <- data_talyies_full %>%
    filter(Day == "D3" & Pop == "CD3_CD4") %>%
    # select(Sample_code, Per_CD3, Per_CD4, Per_CD8, Per_NK,Per_TGD,Per_CD11b) %>%
    # select(Sample_code, Per_CD22_total, Per_CD20, Per_CD19, Per_CD22_CD10_plus) %>%
    # select(Sample_code, Area_mean, Roundness_mean, Viability, Cell_count_PDLS) %>%
    select(Sample_code, Tfr, Tfh, Treg, CD3_Naive,CD3_Central_Memory,CD3_Effector_Memory,CD3_TEMRA) %>%
    rename_with(~ str_replace(., "Per_", ""))
  
  # Préparer les données de déplétion à D6
  data_depletion_d6 <- data_correlation %>%
    select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) 
  
  Sample_list <- unique(data_depletion_d6$Sample_code)
  
  
  ##Retirer les CD20neg##
  # data_talyies <- data_talyies_full %>%
  #   filter(Screening == T & Disease %in% c("FL","tFL","DLBCL") & !Patient_code %in% CD20_neg_data) %>%
  #   mutate(Sample_code = factor(Sample_code, levels = Sample_order))
  # Sample_list <- data_talyies %>%
  #   distinct(Sample_code) %>% pull(Sample_code)
  
  # Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
  data_combined <- data_depletion_d6 %>%
    left_join(markers_d0, by = "Sample_code") %>% 
    filter(Sample_code %in% Sample_list)
  
  # Initialiser un dataframe pour stocker les résultats
  correlation_results <- data.frame(Treatment = character(), Marker = character(), Correlation = numeric(), P_Value = numeric(), Sample_Size = numeric(), stringsAsFactors = FALSE)
  
  # Liste des traitements uniques
  treatments <- unique(data_combined$Treatment_reorder)
  
  # Liste des marqueurs
  markers <- setdiff(colnames(markers_d0), "Sample_code")
  
  # Calculer la corrélation pour chaque traitement et chaque marqueur
  for (treatment in treatments) {
    for (marker in markers) {
      # Filtrer les données pour le traitement et le marqueur actuels
      treatment_data <- data_combined %>%
        filter(Treatment_reorder == treatment)
      
      # Filtrer les valeurs non manquantes
      complete_cases <- complete.cases(treatment_data[[marker]], treatment_data$B_cell_depletion_total)
      
      # Compter le nombre d'échantillons
      sample_size <- sum(complete_cases)
      
      if (sample_size > 1) { # Vérifier qu'il y a suffisamment de paires de valeurs
        # Vérifier la variabilité des données
        if (var(treatment_data[[marker]][complete_cases], na.rm = TRUE) > 0 &
            var(treatment_data$B_cell_depletion_total[complete_cases], na.rm = TRUE) > 0) {
          correlation_test <- cor.test(treatment_data[[marker]][complete_cases],
                                       treatment_data$B_cell_depletion_total[complete_cases])
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = correlation_test$estimate, P_Value = correlation_test$p.value, Sample_Size = sample_size, stringsAsFactors = FALSE))
        } else {
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
        }
      } else {
        correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
      }
    }
  }
  # Afficher les résultats
  print(correlation_results)
  
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_Pop_Effector <- ggplot(correlation_results, aes(x = Marker,
                                                                        y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                        fill = Correlation)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Correlation (R)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("")
    
    # Afficher la heatmap
    print(heatmap_correlation_Pop_Effector)
  }
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_Pop_Effector_nb <- ggplot(correlation_results, aes(x = Marker,
                                                                           y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                           fill = Sample_Size)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Sample_Size, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Sample_Size") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("Nombre de Sample pour calculer les correlations")
    
    # Afficher la heatmap
  }
  
  Figure_correlation_Pop_Effector<- plot_grid(
    plot_grid(heatmap_correlation_Pop_Effector, ncol = 1, labels = c("A")),  
    plot_grid(heatmap_correlation_Pop_Effector_nb, ncol = 1, labels = c("B")), 
    nrow =1,
    rel_heights = c(1, 1)) +
    plot_annotation(title = "Correlation Pop D3 and B Cell_Depletion for each treatment", theme = theme_custom())
  
  # plot_annotation(title = "Correlation between Pop D3 and B_Cell_Depletion for each treatment", theme = theme_custom())
  
  print(Figure_correlation_Pop_Effector)
  
}
#####Correlation entre les marqueurs ICP et la déplétion des cellules B#####
{# Filtrer les données pour obtenir les valeurs des marqueurs à D0
  Day_current <- "D3"
  
  markers_d0 <- data_talyies_full %>%
    filter(Day == "D3", Pop == "CD22") %>%
    # select(Sample_code, TIGIT_Per, PD1_Per, TIM3_Per, LAG3_Per,`PD1/TIM3_Per`,`PD1/LAG3_Per`)
    # select(Sample_code, CD79b_Per, PDL1_Per, PDL2_Per)
    # select(Sample_code, CD79b_Per, PDL1_Per, PDL2_Per)
    select(Sample_code, CD79b_RFI, PDL1_RFI, PDL2_RFI)
  
  # marker_list <- c("TIGIT_Per", "PD1_Per", "TIM3_Per", "LAG3_Per","PD1/TIM3_Per","PD1/LAG3_Per")
  # marker_list <- c("TIGIT_RFI", "PD1_RFI", "TIM3_RFI", "LAG3_RFI","41BB_RFI")
  # marker_list <- c("CD79b_Per", "PDL1_Per", "PDL2_Per")
  marker_list <- c("CD79b_RFI", "PDL1_RFI", "PDL2_RFI")
  
  
  # Préparer les données de déplétion à D6
  data_depletion_d6 <- data_heatmap %>%
    select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
    filter(B_cell_depletion_total > -50) %>%
    distinct()
  
  Sample_list <- unique(data_depletion_d6$Sample_code)
  
  
  ##Retirer les CD20neg##
  # data_talyies <- data_talyies_full %>%
  #   filter(Screening == T & Disease %in% c("FL","tFL","DLBCL") & !Patient_code %in% CD20_neg_data) %>%
  #   mutate(Sample_code = factor(Sample_code, levels = Sample_order))
  # Sample_list <- data_talyies %>%
  #   distinct(Sample_code) %>% pull(Sample_code)
  
  # Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
  data_combined <- data_depletion_d6 %>%
    left_join(markers_d0, by = "Sample_code") %>% 
    filter(Sample_code %in% Sample_list)
  
  # Initialiser un dataframe pour stocker les résultats
  correlation_results <- data.frame(Treatment = character(), Marker = character(), Correlation = numeric(), P_Value = numeric(), Sample_Size = numeric(), stringsAsFactors = FALSE)
  
  # Liste des traitements uniques
  treatments <- unique(data_combined$Treatment_reorder)
  
  # Calculer la corrélation pour chaque traitement et chaque marqueur
  for (treatment in treatments) {
    for (marker in marker_list) {
      # Filtrer les données pour le traitement et le marqueur actuels
      treatment_data <- data_combined %>%
        filter(Treatment_reorder == treatment)
      
      # Filtrer les valeurs non manquantes
      complete_cases <- complete.cases(treatment_data[[marker]], treatment_data$B_cell_depletion_total)
      
      # Compter le nombre d'échantillons
      sample_size <- sum(complete_cases)
      
      if (sample_size > 1) { # Vérifier qu'il y a suffisamment de paires de valeurs
        # Vérifier la variabilité des données
        if (var(treatment_data[[marker]][complete_cases], na.rm = TRUE) > 0 &
            var(treatment_data$B_cell_depletion_total[complete_cases], na.rm = TRUE) > 0) {
          correlation_test <- cor.test(treatment_data[[marker]][complete_cases],
                                       treatment_data$B_cell_depletion_total[complete_cases])
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = correlation_test$estimate, P_Value = correlation_test$p.value, Sample_Size = sample_size, stringsAsFactors = FALSE))
        } else {
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
        }
      } else {
        correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
      }
    }
  }
  
  # Afficher les résultats
  print(correlation_results)
  
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_ICP <- ggplot(correlation_results, aes(x = Marker,
                                                               y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                               fill = Correlation)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Correlation (R)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("")
    
    # Afficher la heatmap
    print(heatmap_correlation_ICP)
  }
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_ICP_nb <- ggplot(correlation_results, aes(x = Marker,
                                                                  y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                  fill = Sample_Size)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Sample_Size, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Sample_Size") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("Nombre de Sample pour calculer les correlations")
    
    # Afficher la heatmap
  }
  
  Figure_correlation_Tcells_ICP <- plot_grid(
    plot_grid(heatmap_correlation_ICP, ncol = 1, labels = c("A")),  
    plot_grid(heatmap_correlation_ICP_nb, ncol = 1, labels = c("B")), 
    nrow =1,
    rel_heights = c(1, 1)) +
    plot_annotation(title = paste0("Correlation between percentage of B-cell expressing ICP at ", Day_current," and B_Cell_Depletion for each treatment"), theme = theme_custom())
  
  print(Figure_correlation_Tcells_ICP)
}

####Correlation avec UT D6 ####
{  # Filtrer les données pour obtenir les valeurs des marqueurs à D0
  markers_d6 <- data_correlation %>%
    filter(Day == "D6" & Treatment == "UT") %>%
    # select(Sample_code, , Per_CD4, Per_CD8, Per_NK,Per_CD22_total,Per_CD22_CD10_plus) %>%
    # select(Sample_code,Per_CD107a_CD4, Per_CD107a_CD8,TNFa, INFY, GrzB, IL8,IL6,IL10) %>% 
    select(Sample_code, Area_mean, Roundness_mean, Viability, Cell_count_PDLS) %>%
    rename_with(~ str_replace(., "Per_", ""))
  
  # Préparer les données de déplétion à D6
  data_depletion_d6 <- data_correlation %>%
    select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) 
  
  Sample_list <- unique(data_depletion_d6$Sample_code)
  
  
  ##Retirer les CD20neg##
  # data_talyies <- data_talyies_full %>%
  #   filter(Screening == T & Disease %in% c("FL","tFL","DLBCL") & !Patient_code %in% CD20_neg_data) %>%
  #   mutate(Sample_code = factor(Sample_code, levels = Sample_order))
  # Sample_list <- data_talyies %>%
  #   distinct(Sample_code) %>% pull(Sample_code)
  
  # Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
  data_combined <- data_depletion_d6 %>%
    left_join(markers_d6, by = "Sample_code") %>% 
    filter(Sample_code %in% Sample_list)
  
  # Initialiser un dataframe pour stocker les résultats
  correlation_results <- data.frame(Treatment = character(), Marker = character(), Correlation = numeric(), P_Value = numeric(), Sample_Size = numeric(), stringsAsFactors = FALSE)
  
  # Liste des traitements uniques
  treatments <- unique(data_combined$Treatment_reorder)
  
  # Liste des marqueurs
  markers <- setdiff(colnames(markers_d6), "Sample_code")
  
  # Calculer la corrélation pour chaque traitement et chaque marqueur
  for (treatment in treatments) {
    for (marker in markers) {
      # Filtrer les données pour le traitement et le marqueur actuels
      treatment_data <- data_combined %>%
        filter(Treatment_reorder == treatment)
      
      # Filtrer les valeurs non manquantes
      complete_cases <- complete.cases(treatment_data[[marker]], treatment_data$B_cell_depletion_total)
      
      # Compter le nombre d'échantillons
      sample_size <- sum(complete_cases)
      
      if (sample_size > 1) { # Vérifier qu'il y a suffisamment de paires de valeurs
        # Vérifier la variabilité des données
        if (var(treatment_data[[marker]][complete_cases], na.rm = TRUE) > 0 &
            var(treatment_data$B_cell_depletion_total[complete_cases], na.rm = TRUE) > 0) {
          correlation_test <- cor.test(treatment_data[[marker]][complete_cases],
                                       treatment_data$B_cell_depletion_total[complete_cases])
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = correlation_test$estimate, P_Value = correlation_test$p.value, Sample_Size = sample_size, stringsAsFactors = FALSE))
        } else {
          correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
        }
      } else {
        correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
      }
    }
  }
  # Afficher les résultats
  print(correlation_results)
  
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_Pop_Effector <- ggplot(correlation_results, aes(x = Marker,
                                                                        y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                        fill = Correlation)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Correlation (R)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("")
    
    # Afficher la heatmap
    print(heatmap_correlation_Pop_Effector)
  }
  {# Créer la heatmap avec les valeurs de corrélation sur les cases
    heatmap_correlation_Pop_Effector_nb <- ggplot(correlation_results, aes(x = Marker,
                                                                           y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                           fill = Sample_Size)) +
      geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
      geom_text(aes(label = round(Sample_Size, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
      scale_fill_distiller(palette = "RdYlBu", direction = -1) +
      labs(x = "Marker",
           y = "Treatment",
           fill = "Sample_Size") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
            axis.ticks = element_blank(),  # Retirer les ticks
            panel.grid = element_blank(),  # Retirer la grille
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.5, "cm"),
            aspect.ratio = 1) +  # Faire des cases carrées
      ggtitle("Nombre de Sample pour calculer les correlations")
    
    # Afficher la heatmap
  }
  
  Figure_correlation_Pop_Effector<- plot_grid(
    plot_grid(heatmap_correlation_Pop_Effector, ncol = 1, labels = c("A")),  
    plot_grid(heatmap_correlation_Pop_Effector_nb, ncol = 1, labels = c("B")), 
    nrow =1,
    rel_heights = c(1, 1)) +
    plot_annotation(title = "Correlation UT D6 and B Cell_Depletion", theme = theme_custom())
  
  # plot_annotation(title = "Correlation between Pop D3 and B_Cell_Depletion for each treatment", theme = theme_custom())
  
  print(Figure_correlation_Pop_Effector)
  
}


#####Correlation Area / Via #####

data_area <- data_talyies_full %>%
  filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>%
  filter(Disease %in% c("FL", "DLBCL", "tFL") ) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment"))


area_d3 <- data_area %>% select(Sample_code,Area_mean,Day,Pop) %>% 
  filter(Day == "D3" & Pop == "CD3_CD4") %>% 
  rename(Area_mean_D3 = Area_mean) %>% 
  select(-Day,-Pop)


area_d6 <- data_area %>% select(Sample_code,Area_mean,Day,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
  filter(Day == "D6" & !is.na(Area_mean)) %>% 
  rename(Area_mean_D6 = Area_mean) %>% 
  select(-Day) %>% 
  left_join(area_d3) %>% 
  mutate( Ratio_mean_area_D6_D3 = Area_mean_D6 / Area_mean_D3)
treatments <- unique(data_combined$Treatment_reorder)
markers <- "ratio_area_D6_D3"
correlation_results <- data.frame(Treatment = character(), Marker = character(), Correlation = numeric(), P_Value = numeric(), Sample_Size = numeric(), stringsAsFactors = FALSE)

# Calculer la corrélation pour chaque traitement et chaque marqueur
for (treatment in treatments) {
  marker <- "ratio_cell_count_D6_D3"
  markers <- "ratio_cell_count_D6_D3"
  treatment_data <- area_d6 %>%
    filter(Treatment_reorder == treatment)
  
  # Filtrer les valeurs non manquantes
  complete_cases <- complete.cases(treatment_data[[marker]], treatment_data$B_cell_depletion_total)
  
  # Compter le nombre d'échantillons
  sample_size <- sum(complete_cases)
  
  if (sample_size > 1) { # Vérifier qu'il y a suffisamment de paires de valeurs
    # Vérifier la variabilité des données
    if (var(treatment_data[[marker]][complete_cases], na.rm = TRUE) > 0 &
        var(treatment_data$B_cell_depletion_total[complete_cases], na.rm = TRUE) > 0) {
      correlation_test <- cor.test(treatment_data[[marker]][complete_cases],
                                   treatment_data$B_cell_depletion_total[complete_cases])
      correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = correlation_test$estimate, P_Value = correlation_test$p.value, Sample_Size = sample_size, stringsAsFactors = FALSE))
    } else {
      correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
    }
  } else {
    correlation_results <- rbind(correlation_results, data.frame(Treatment_reorder = treatment, Marker = marker, Correlation = NA, P_Value = NA, Sample_Size = sample_size, stringsAsFactors = FALSE))
  }
}

# Afficher les résultats
print(correlation_results)

{# Créer la heatmap avec les valeurs de corrélation sur les cases
  heatmap_correlation_Pop_Effector <- ggplot(correlation_results, aes(x = Marker,
                                                                      y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                      fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
    geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    labs(x = "Marker",
         y = "Treatment",
         fill = "Correlation (R)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
          axis.ticks = element_blank(),  # Retirer les ticks
          panel.grid = element_blank(),  # Retirer la grille
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"),
          aspect.ratio = 1) +  # Faire des cases carrées
    ggtitle("")
  
  # Afficher la heatmap
  print(heatmap_correlation_Pop_Effector)
}
{# Créer la heatmap avec les valeurs de corrélation sur les cases
  heatmap_correlation_Pop_Effector_nb <- ggplot(correlation_results, aes(x = Marker,
                                                                         y = factor(Treatment_reorder, levels = rev(sort(unique(data_combined$Treatment_reorder)))),
                                                                         fill = Sample_Size)) +
    geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
    geom_text(aes(label = round(Sample_Size, 2)), color = "black", size = 3) +  # Ajouter les valeurs de corrélation
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    labs(x = "Marker",
         y = "Treatment",
         fill = "Sample_Size") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(hjust = 0, size = 8),  # Aligner le texte de l'axe y à gauche
          axis.ticks = element_blank(),  # Retirer les ticks
          panel.grid = element_blank(),  # Retirer la grille
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"),
          aspect.ratio = 1) +  # Faire des cases carrées
    ggtitle("Nombre de Sample pour calculer les correlations")
  
  # Afficher la heatmap
}

Figure_correlation_Pop_Effector<- plot_grid(
  plot_grid(heatmap_correlation_Pop_Effector, ncol = 1, labels = c("A")),  
  plot_grid(heatmap_correlation_Pop_Effector_nb, ncol = 1, labels = c("B")), 
  nrow =1,
  rel_heights = c(1, 1)) +
  plot_annotation(title = "Correlation UT D6 and B Cell_Depletion", theme = theme_custom())

# plot_annotation(title = "Correlation between Pop D3 and B_Cell_Depletion for each treatment", theme = theme_custom())

print(Figure_correlation_Pop_Effector)

