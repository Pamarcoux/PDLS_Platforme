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
library(ggpubr)
library(grid)
  
  


sample_colors_df <- tibble(Sample_code = names(sample_colors_all),
                           sample_color = unname(sample_colors_all))

conditions <- c("B_cell_depletion_total", "Per_CD107a_CD4_normalized", "Per_CD107a_CD8_normalized",
                "TNFa_normalized", "INFY_normalized", "GrzB_normalized", "IL8_normalized", "IL6_normalized")

Condition = "Response_type_total"

#### Choix des Combos ####
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
}
  
####Heatmap####
######HEAT MAP B Cell Depletion#######
{# Préparer les données
  data_treatment <- data_talyies_full %>%
    select(Treatment, Disease, B_cell_depletion_total, Sample_code, Day, !!sym(Condition),Origin,Response_type_total) %>%
    mutate(B_cell_depletion_total = replace(B_cell_depletion_total, B_cell_depletion_total == "", NA)) %>%
    dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
    right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
    distinct() %>% 
    dplyr::filter(Treatment_type == "Single")
  
  # Compter le nombre de NA dans la colonne Treatment pour chaque Sample_code
  na_count_sample <- data_treatment %>%
    group_by(Sample_code) %>%
    summarise(na_count_sample= sum(is.na(B_cell_depletion_total)))
  
  # Compter le nombre de NA pour les traitements 
  na_counts_treatment <- data_treatment %>%
    group_by(Treatment) %>%
    summarise(na_count_treatment = sum(is.na(B_cell_depletion_total)))
  
  data_treatment_filtered <- data_treatment
  # %>% 
  #   left_join(na_counts_treatment) %>% 
  #   left_join(na_count_sample) %>% 
  #   dplyr::filter(na_count_sample < 48  & na_count_treatment < 30  )
  
  ordered_samples <- levels(data_talyies_full$Sample_code)
  # Trier les Sample_code par nombre décroissant de NA
  # ordered_samples <- na_count_sample %>%
  #   arrange((na_count_sample)) %>%
  #   mutate(Sample_code = as.character(Sample_code)) %>%
  #   pull(Sample_code)

  # Réorganiser data_treatment selon l'ordre des Sample_code
  data_treatment_ordered <- data_treatment_filtered %>%
    dplyr::filter(Sample_code %in% ordered_samples) %>%
    mutate(Sample_code = factor(Sample_code, levels = ordered_samples))
  
  # Créer la heatmap avec les Sample_code triés
  Heatmap_B_cell_Depletion <- ggplot(data_treatment_ordered, 
                                     aes(x = Sample_code,
                                         y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                         fill = B_cell_depletion_total)) +
    geom_tile(color = "white", size = 0.5) +  # Bordures blanches
    # Ajout des valeurs dans chaque case
    # geom_text(aes(label = sprintf("%.1f", B_cell_depletion_total),
    #     color = ifelse(B_cell_depletion_total < 0 | B_cell_depletion_total > 50, "white", "black")), 
    #     size = 3, fontface = "bold") +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-40, 100)) +
    scale_color_identity() +  # Utilise directement les couleurs spécifiées
    labs(x = "",
         y = "",
         fill = "B Cell Depletion (%)") +
    theme_minimal() +
    scale_y_discrete(position = "right") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
          axis.text.y = element_text(hjust = 0, size = 8, face = "bold"),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "left",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.5, "cm"),
          aspect.ratio = 1) +  # Cases carrées
    ggtitle("Heatmap of B Cell Depletion by Treatment Combinations")
  
  # Afficher la heatmap
  print(Heatmap_B_cell_Depletion)
}
#### Figure 8 Heatmaps ####
{data_heatmap <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, all_of(conditions)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  dplyr::filter(Treatment_type == "Single") %>%
  dplyr::filter(B_cell_depletion_total > -50) %>%
  pivot_longer(cols = all_of(conditions), names_to = "Condition", values_to = "value") %>%
  mutate(value = as.numeric(value),
         Treatment_reorder = factor(Treatment, levels = rev(sort(unique(Treatment))))) %>%
  left_join(sample_colors_df, by = "Sample_code") 

##### Créer les plots ----
plots <- vector("list", length(conditions))
ncol_plot <- 4

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap, Condition == cond) %>% 
    group_by(Treatment_Cat,Sample_code) %>% 
    mutate(value_mean = mean(value, na.rm = TRUE))
    
  
  # Heatmap principale
  heatmap_layer <- ggplot(plot_data, aes(x = Sample_code, y = Treatment_Cat)) +
    geom_tile(aes(fill = value_mean), color = "white", size = 0.5) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey90") +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = cond, fill = cond) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
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

#### Assembler tous les plots ####
final_plot_heatmap <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
}

print(final_plot_heatmap)
ggsave(
  filename = "Paper_platforme/Figure/Figure_4_Heatmap.png",
  plot = final_plot_heatmap,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 22 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)



####Heatmap à c/c ######
{Condition = "Per_CD107a_CD8_normalized"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct() %>% 
  dplyr::filter(Treatment_type == "Single") %>% 
  mutate(value = (!!sym(Condition))) %>% 
  dplyr::filter(!Sample_code %in% c("tFL1_LN1","DLBCL5_PB1","tFL5_PB1"))

# Créer la heatmap avec les valeurs normalisées
heatmap_Condition <- ggplot(data_treatment, aes(x = Sample_code,
                                                           y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                                           fill = value)) +
  geom_tile(color = "white", size = 0.5) +  # Bordures blanches
  # Ajout des valeurs dans chaque case
  # geom_text(aes(label = sprintf("%.1f", value),
  #     color = ifelse(value < 0 | value > 50, "white", "black")),
  #     size = 3, fontface = "bold") +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(NA, NA)) +
  scale_color_identity() +  # Utilise directement les couleurs spécifiées
  labs(x = "",
       y = "",
       fill = "B Cell Depletion (%)") +
  theme_minimal() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(hjust = 0, size = 8, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "left",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        aspect.ratio = 1) +  # Cases carrées  # Faire des cases carrées
  labs(x = "",
       y = "",
       fill = paste0("",Condition)) +
  ggtitle(paste0("Heatmap of ",Condition, " Treatment Combinations"))

# Afficher la heatmap
print(heatmap_Condition)
}


##### HeatMap Modulable with Normalization######
{Condition = "INFY"
  
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct() %>% 
  dplyr::filter(Treatment_type == "Single")

data_treatment_normalized <- data_treatment %>%
  left_join(
    data_treatment %>%
      dplyr::filter(Treatment_reorder == "UT") %>%
      select(Sample_code, UT_value = !!sym(Condition)),
    by = "Sample_code"
  ) %>%
  mutate(normalized_value = (!!sym(Condition)) / UT_value)

# Créer la heatmap avec les valeurs normalisées
heatmap_Condition <- ggplot(data_treatment_normalized, aes(x = Sample_code,
                                                      y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                                      fill = normalized_value)) +
  geom_tile(color = "white", size = 1) +  # Bordures blanches
  # Ajout des valeurs dans chaque case
  # geom_text(aes(label = sprintf("%.1f", B_cell_depletion_total),
  #     color = ifelse(B_cell_depletion_total < 0 | B_cell_depletion_total > 50, "white", "black")), 
  #     size = 3, fontface = "bold") +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(NA, NA)) +
  scale_color_identity() +  # Utilise directement les couleurs spécifiées
  labs(x = "",
       y = "",
       fill = "B Cell Depletion (%)") +
  theme_minimal() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(hjust = 0, size = 8, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "left",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        aspect.ratio = 1) +  # Cases carrées  # Faire des cases carrées
  labs(x = "Sample Code",
       y = "Treatment Combination",
       fill = paste0("",Condition)) +
  ggtitle(paste0("Heatmap of ",Condition, " Treatment Combinations (Normalized)"))

# Afficher la heatmap
print(heatmap_Condition)
}

##### HeatMap Modulable without Normalization######
{Condition = "Per_CD107a_CD8_normalized"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct() %>% 
  dplyr::filter(Treatment_type == "Single")

# Calculer le z-score pour chaque Sample
data_treatment_z <- data_treatment %>% #Zscore by Sample
  group_by(Sample_code) %>%
  # mutate(z_score = scale(!!sym(Condition))) %>%
  mutate(z_score = (!!sym(Condition))) %>%
  ungroup()

# Heatmap
heatmap_Condition <- ggplot(data_treatment_z, aes(x = Sample_code,
                                                y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                                fill = z_score)) +
  geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches

  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(x = "Sample Code",
       y = "Treatment Combination",
       fill = paste0("",Condition)) +
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
  ggtitle(paste0("Heatmap of ",Condition, " Treatment Combinations"))

# Afficher la heatmap
print(heatmap_Condition)
}

##### HeatMap Modulable Zscore without Normalization######
{Condition = "B_cell_depletion_total"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct() %>% 
  dplyr::filter(grepl("CD20-TCB",Treatment)) 
  # dplyr::filter(Treatment == "Single")

# Calculer le z-score pour chaque Sample
data_treatment_z <- data_treatment %>% #Zscore by Sample
  group_by(Sample_code) %>%
  mutate(z_score = scale(!!sym(Condition))) %>%
  ungroup()

# Heatmap
heatmap_Condition <- ggplot(data_treatment_z, aes(x = Sample_code,
                                                  y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                                  fill = z_score)) +
  geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(x = "Sample Code",
       y = "Treatment Combination",
       fill = paste0("",Condition)) +
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
  ggtitle(paste0("Heatmap of ",Condition, " Treatment Combinations (Z-score)"))

# Afficher la heatmap
print(heatmap_Condition)
}


#### Correlation ration Area and B Cell Depletion ######
# Préparer les données
{Condition = "Area_mean"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct()


data_area_normalized <- data_treatment %>%
  left_join(
    data_treatment %>%
      dplyr::filter(Treatment_reorder == "UT") %>%
      select(Sample_code, UT_Area_mean = Area_mean),
    by = "Sample_code"
  ) %>%
  mutate(normalized_Area_mean = (Area_mean / UT_Area_mean)) %>%
  dplyr::filter(B_cell_depletion_total > -50 | normalized_Area_mean < 2)

# Calculer la corrélation de Pearson entre Area et B_Cell_Depletion
correlation_result <- cor.test(data_area_normalized$normalized_Area_mean, data_area_normalized$B_cell_depletion_total,method = "pearson", use = "pairwise.complete.obs", alternative = "two.sided")
r_value <- correlation_result$estimate  # Coefficient de corrélation
p_value <- correlation_result$p.value   # Valeur p

# fit_cluster <- lm(B_cell_depletion_total_glofi ~ normalized_Area_mean, data = data_area_normalized)
# summary_fit_cluster <- summary(fit_cluster)
# r_squared_cluster <- summary_fit_cluster$r.squared

##PLOT##
p <- ggplot(data_area_normalized, aes(x = B_cell_depletion_total, y = normalized_Area_mean)) +
  geom_point(aes(fill = Sample_code),
             shape = 21,
             color = "black",
             position = position_dodge(0),
             size = 2.5, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all) +  scale_x_continuous(limits = c(-50, 100)) + # Définir les limites de l'axe x
  geom_smooth(method = "lm", color = "black", fill = "gray80", size= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
  labs(title = "Correlation between Normalized Area and B Cell Depletion",
       x = "B Cell Depletion (%)",
       y = "Normalized Area (compared to UT at D6)") +
  theme_custom()

correlation_area_plot <- p + annotation_custom(
  grob = textGrob(
    label = sprintf("\nR = %.4f\np = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.98, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 12, col = "black")
  ))


# Afficher le plot
print(correlation_area_plot)
}

#### Autres Correlation ####
##### Correlation Condition WITH normalization and B Cell Depletion ######
# Préparer les données
{
  Condition = "IL6"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct()

data_normalized <- data_treatment %>%
  left_join(
    data_treatment %>%
      dplyr::filter(Treatment_reorder == "UT") %>%
      select(Sample_code, UT_value = !!sym(Condition)),
    by = "Sample_code"
  ) %>%
  mutate(normalized_value = (!!sym(Condition) / UT_value)) %>%
  dplyr::filter(B_cell_depletion_total > -50 & !!sym(Condition) < 110  )

# Calculer la corrélation de Pearson entre Area et B_Cell_Depletion
correlation_result <- cor.test(data_normalized$normalized_value, data_normalized$B_cell_depletion_total)
r_value <- correlation_result$estimate  # Coefficient de corrélation
p_value <- correlation_result$p.value   # Valeur p

##PLOT##
p <- ggplot(data_normalized, aes(x = B_cell_depletion_total, y = !!sym(Condition))) +
  geom_point(aes(fill = Sample_code),
             shape = 21,
             color = "black",
             position = position_dodge(0),
             size = 2.5, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all) +  
  scale_x_continuous(limits = c(-50, 100)) + # Définir les limites de l'axe x
  geom_smooth(method = "lm", color = "black", fill = "gray80", size= 0.8) +
  labs(title = paste0("Correlation between ",Condition," and B Cell Depletion"),
       x = "B Cell Depletion (%)",
       y = paste0("",Condition) )+
  theme_custom()

correlation_plot <- p + 
  annotation_custom(grob = textGrob(
      label = sprintf("\nR = %.4f\np = %.2e", r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.98, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 12, col = "black")))

# Afficher le plot
print(correlation_plot)
}


##### Correlation Condition WITHOUT normalization and B Cell Depletion ######
# Préparer les données
{Condition = "Per_CD107a_CD4"
filter_value = 10000

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  dplyr::filter(Treatment == "αTIGIT 10 µg/mL") %>% 
  # dplyr::filter(Treatment_type == "Single") %>% 
  distinct()

data_unnormalized <- data_treatment %>%
  mutate(value = (!!sym(Condition))) %>%
  dplyr::filter(B_cell_depletion_total > -50) %>% 
  dplyr::filter(value < filter_value)

# Calculer la corrélation de Pearson entre Area et B_Cell_Depletion

correlation_result <- cor.test(data_unnormalized$value, data_unnormalized$B_cell_depletion_total)
r_value <- correlation_result$estimate  # Coefficient de corrélation
p_value <- correlation_result$p.value   # Valeur p

##PLOT##
p <- ggplot(data_unnormalized, aes(x = B_cell_depletion_total, y = !!sym(Condition),fill)) +
  geom_smooth(method = "lm", color = "black", fill = "gray80", size= 0.8) +
  geom_point(aes(fill = Sample_code),
             shape = 21,
             color = "black",
             position = position_dodge(0),
             size = 2.5, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all) +  
  scale_x_continuous(limits = c(0, NA)) + # Définir les limites de l'axe x
  labs(title = paste0("Correlation between ",Condition," and B Cell Depletion"),
       x = "B Cell Depletion (%)",
       y = paste0("",Condition) )+
  theme_custom()+
  theme(legend.position = "")

correlation_plot <- p + 
  annotation_custom(grob = textGrob(
    label = sprintf("\nR = %.4f\np = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.90, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 12, col = "black")))

# Afficher le plot
print(correlation_plot)
}

####Calcul Dataframe toutes les correlations a D6 ####
library(dplyr)

# Préparer les données
data_treatment <- data_talyies_full %>%
  select(-...1) %>% 
  dplyr::filter(Day %in% c("D6","D0","D3"), Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  dplyr::filter(B_cell_depletion_total > -50) %>%
  distinct()

# Sélectionner uniquement les colonnes numériques pour le calcul de la corrélation
numeric_columns <- data_treatment %>%
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

#### Single HeatMap Zscore  #####
{Condition = "B_cell_depletion_total"

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  distinct()
# %>% 
#   dplyr::filter(grepl("GA101",Treatment_reorder)) 
# # dplyr::filter(Treatment == "Single")

# Compter le nombre de NA dans la colonne Treatment pour chaque Sample_code
No_data_sample_list <- data_treatment %>%
  group_by(Sample_code) %>%
  summarise(na_count_sample= sum(is.na(B_cell_depletion_total)),
            No_data_sample = ifelse(na_count_sample > (length(unique(data_treatment$Treatment_reorder))-2), TRUE, FALSE)) %>%
  dplyr::filter(No_data_sample == TRUE) %>% 
  pull(Sample_code)

# Calculer le z-score pour chaque Sample
data_treatment_z <- data_treatment %>%
  dplyr::filter(!Sample_code %in% No_data_sample_list) %>% 
  group_by(Sample_code) %>%
  mutate(z_score = scale(!!sym(Condition))) %>%
  ungroup()

# Heatmap
heatmap_Condition <- ggplot(data_treatment_z, aes(x = Sample_code,
                                                  y = factor(Treatment_reorder, levels = rev(sort(unique(Treatment_reorder)))),
                                                  fill = z_score)) +
  geom_tile(color = "white", size = 0.5) +  # Ajouter des bordures blanches
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(x = "Sample Code",
       y = "",
       fill = paste0("Z-Score")) +
  theme_minimal() +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(hjust = 0, size = 8, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "left",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        aspect.ratio = 1) +  # Cases carrées
  ggtitle(paste0("Glofitamab 1 nM ",Condition))

# Afficher la heatmap
print(heatmap_Condition)}


#### HUIT Heatmaps Zscore ####
{data_heatmap <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, all_of(conditions)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("PD1",Treatment_reorder)) 

# Compter le nombre de NA dans la colonne Treatment pour chaque Sample_code
No_data_sample_list <- data_heatmap %>%
  group_by(Sample_code) %>%
  summarise(na_count_sample= sum(is.na(B_cell_depletion_total)),
            No_data_sample = ifelse(na_count_sample > (length(unique(data_heatmap$Treatment_reorder))-2), TRUE, FALSE)) %>%
  dplyr::filter(No_data_sample == TRUE) %>% 
  pull(Sample_code)

data_heatmap_pivot <- data_heatmap %>% 
  pivot_longer(cols = all_of(conditions), names_to = "Condition", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Treatment_reorder = factor(Treatment_reorder, levels = rev(levels(factor(Treatment_reorder))))) %>% 
  left_join(sample_colors_df, by = "Sample_code") 

plots <- vector("list", length(conditions))
ncol_plot <- 4

#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Treatment_reorder) %>%
    dplyr::filter(!Sample_code %in% No_data_sample_list) %>% 
    mutate(z_score = scale(value),
           Sample_code = factor(Sample_code, levels = Sample_order)) %>%
    ungroup()
  
  # Heatmap principale
  heatmap_layer <- ggplot(plot_data, aes(x = Sample_code, y = (Treatment_reorder))) +
    geom_tile(aes(fill = z_score), color = "white", size = 0.5) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey90",limits = c(-2, 3)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = cond, fill = "Z-score") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
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
final_plot_heatmap <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(final_plot_heatmap)
}

ggsave(
  filename = "Paper_platforme/Figure/Figure_4_Heatmap.png",
  plot = final_plot_heatmap,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 22 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
#### HUIT Heatmaps Zscore Categorie treatment ####
{data_heatmap <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, all_of(conditions)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  distinct() %>% 
  dplyr::filter(Treatment_type == "Single")



# Compter le nombre de NA dans la colonne Treatment pour chaque Sample_code
No_data_sample_list <- data_heatmap %>%
  group_by(Sample_code) %>%
  summarise(na_count_sample= sum(is.na(B_cell_depletion_total)),
            No_data_sample = ifelse(na_count_sample > (length(unique(data_heatmap$Treatment_reorder))-2), TRUE, FALSE)) %>%
  dplyr::filter(No_data_sample == TRUE) %>% 
  pull(Sample_code)

data_heatmap_pivot <- data_heatmap %>% 
  pivot_longer(cols = all_of(conditions), names_to = "Condition", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Treatment_reorder = factor(Treatment_reorder, levels = rev(levels(factor(Treatment_reorder))))) %>% 
  left_join(sample_colors_df, by = "Sample_code") 

plots <- vector("list", length(conditions))
ncol_plot <- 4

#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Treatment_Cat,Sample_code) %>% 
    mutate(value_mean = mean(value, na.rm = TRUE)) %>% 
    group_by(Sample_code) %>%
    dplyr::filter(!Sample_code %in% No_data_sample_list) %>% 
    mutate(z_score = scale(value_mean),
           Sample_code = factor(Sample_code, levels = Sample_order)) %>%
    ungroup()
  
  # Heatmap principale
  heatmap_layer <- ggplot(plot_data, aes(x = Sample_code, y = (Treatment_reorder))) +
    geom_tile(aes(fill = z_score), color = "white", size = 0.5) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey90",limits = c(-2, 3)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = cond, fill = "Z-score") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
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
final_plot_heatmap <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(final_plot_heatmap)
}

ggsave(
  filename = "Paper_platforme/Figure/Figure_4_Heatmap.png",
  plot = final_plot_heatmap,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 22 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)


#### Heatmap cluster pHEATMAP####
cal_z_score <- function(x){
  (x)}

{  
  filename = "Paper_platforme/Figure/Figure_5_Heatmap/Clusters_by_Treatment_type.png"
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
    summarize(B_cell_depletion_total = mean(B_cell_depletion_total))
  
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
           filename=filename,
           width = 24,
           height = 7,
           display_numbers = TRUE)
  dev.off()
  print(resultat_heatmap)
}
k <- 3  
clusters <- cutree(resultat_heatmap$tree_col, k)
df_Cluster_response_cat <- data.frame(
  Sample_code = names(clusters),
  Cluster_response_cat = unlist(clusters), row.names = NULL) %>% 
  mutate(Cluster_response_cat = recode(Cluster_response_cat, 
                                       `1` = "ADC_TCB_Medium_responders", 
                                       `2` = "Global_low_responders", 
                                       `3` = "TCB_High_Responders")) 


####Comparaison stat treatment #####
# Créer un dataframe pour stocker les résultats statistiques
{Condition <- "B_cell_depletion_total"
# Cluster_response_cat_filter <- "ADC_TCB_Medium_responders"
Control <- "αTIGIT 1 µg/mL"

#With Cluster Classification
data_treatment_classification <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total,Response_type_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") ) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  dplyr::filter(grepl(Control,Treatment_reorder)) %>%
  left_join(df_Cluster_response_cat, by = "Sample_code") %>%
  # dplyr::filter(Cluster_response_cat %in% Cluster_response_cat_filter) %>% 
  dplyr::filter(!is.na(!!sym(Condition)))


model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment_classification)

# Vérifier la convergence du modèle
summary(model)

# Table d'ANOVA pour le facteur fixe (Treatment)
anova(model)

# Extraire les estimations pour chaque traitement
fixef(model)

# Extraire la variance expliquée par l'effet aléatoire (variabilité inter-patients)
VarCorr(model)

# Calculer l'ICC (Coefficient de corrélation intraclasse) 
# qui indique la proportion de variance due aux patients
icc <- as.data.frame(VarCorr(model))$vcov[1] / 
  (as.data.frame(VarCorr(model))$vcov[1] + attr(VarCorr(model), "sc")^2)
print(paste("ICC =", round(icc, 3)))

# Moyennes estimées pour chaque traitement
emm <- emmeans(model, ~ Treatment_reorder)
print(emm)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = Control)

# # Test si chaque traitement est différent de zéro
# # (si UT a une valeur de 0 et n'est pas dans vos données)
# test(emm)

# # Comparaisons post-hoc entre traitements
# pairs(emm, adjust = "bonferroni")

print(contrasts)
plot(emm, comparisons = TRUE) # Affiche les moyennes et IC95%

# Convertir en data frame pour manipulation
contrast_df <- as.data.frame(contrasts)

# Ajouter une colonne pour les symboles de significativité
contrast_df$signif <- ""
contrast_df$signif[contrast_df$p.value < 0.1] <- "."
contrast_df$signif[contrast_df$p.value < 0.05] <- "*"
contrast_df$signif[contrast_df$p.value < 0.01] <- "**"
contrast_df$signif[contrast_df$p.value < 0.001] <- "***"

# Formater le tableau pour l'affichage
contrast_table <- contrast_df %>%
  mutate(
    p.value = format.pval(p.value, digits = 3),
    p.value_signif = paste0(p.value, " ", signif)
  ) %>%
  select(contrast, estimate, SE, df, t.ratio, p.value_signif)

# Afficher le tableau formaté
print(contrast_table)


# Prepare the contrast data
contrast_df <- as.data.frame(contrasts)

# Clean up labels and compute 95% CI
contrast_df <- contrast_df %>%
  mutate(
    Treatment = sub(paste0(" - ", Control, "$"), "", contrast),  # Enlève " - UT" à la fin
    Treatment = sub("^\\((.*)\\)$", "\\1", Treatment),    
    signif = ifelse(p.value < 0.05, "p < 0.05", "ns"),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) %>% 
  mutate(Treatment = factor(Treatment, levels = rev(unique(Treatment))))

# Forest plot
forest_plot <- ggplot(contrast_df, aes(x = Treatment, 
                                       y = estimate, 
                                       ymin = lower, 
                                       ymax = upper, 
                                       color = signif)) +
  geom_pointrange(size = 0.9, fatten = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  scale_color_manual(values = c("p < 0.05" = "#D73027", "ns" = "gray50")) +
  labs(
    title = paste0("Treatment Effects on ",Condition, " (vs ",Control,")"),
    x = "Treatment",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "top"
  )+
  theme_custom()

####BAR Plot comparaison Treatment ####
barplot <- ggplot(data_treatment_classification, aes(x=Treatment_reorder , y=!!sym(Condition))) +
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
  # stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
  #                    label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
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

final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                        nrow =2,
                        rel_widths = c(1, 0.9))
print(final_plot)
}
