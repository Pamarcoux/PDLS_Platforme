library(dplyr)
{#####Correlation entre les pops à D0 et la déplétion des cellules B#####
# Filtrer les données pour obtenir les valeurs des marqueurs à D0
markers_d0 <- data_talyies_full %>%
  filter(Day == "D3" & Pop == "CD3_CD4") %>%
  # select(Sample_code, Per_CD3, Per_CD4, Per_CD8, Per_NK,Per_TGD,Per_CD11b) %>% 
  # select(Sample_code, Per_CD22_total, Per_CD20, Per_CD19, Per_CD22_CD10_plus) %>% 
    # select(Sample_code, Area_mean, Roundness_mean, Viability, Cell_count_PDLS) %>% 
  select(Sample_code, Tfr, Tfh, Treg, CD3_Naive,CD3_Central_Memory,CD3_Effector_Memory,CD3_TEMRA) %>% 
  rename_with(~ str_replace(., "Per_", ""))

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_treatment_classification %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
  filter(B_cell_depletion_total > -50) %>%
  filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>%
  filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") ) %>%
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
  plot_annotation(title = "Correlation Pop CD4 D3 and B Cell_Depletion for each treatment", theme = theme_custom())
  
  # plot_annotation(title = "Correlation between Pop D3 and B_Cell_Depletion for each treatment", theme = theme_custom())

print(Figure_correlation_Pop_Effector)

}
#####Correlation entre les marqueurs ICP et la déplétion des cellules B#####
{# Filtrer les données pour obtenir les valeurs des marqueurs à D0
  Day_current <- "D3"
  
markers_d0 <- data_talyies_full %>%
  filter(Day == "D3", Pop == "CD3_CD4") %>%
  select(Sample_code, TIGIT_Per, PD1_Per, TIM3_Per, LAG3_Per,`PD1/TIM3_Per`,`PD1/LAG3_Per`)

marker_list <- c("TIGIT_Per", "PD1_Per", "TIM3_Per", "LAG3_Per","PD1/TIM3_Per","PD1/LAG3_Per")
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

# Liste des marqueurs
markers <- c("TIGIT_Per", "PD1_Per", "TIM3_Per", "LAG3_Per","PD1/TIM3_Per","PD1/LAG3_Per")

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
  plot_annotation(title = paste0("Correlation between percentage of CD4+ cell expressing ICP at ", Day_current," and B_Cell_Depletion for each treatment"), theme = theme_custom())

print(Figure_correlation_Tcells_ICP)
}
##### Correlation ######
data_correlation <- data_talyies_full %>% 
  filter(Day == "D6", Disease %in% Disease_list) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  filter(Treatment_type == "Single") %>%
  filter(B_cell_depletion_total > -50) 


####Calcul all correlations D6 #####
