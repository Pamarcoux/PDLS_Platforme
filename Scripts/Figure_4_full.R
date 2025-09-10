# Charger les packages nécessaires
  library(lme4)      # Pour les modèles mixtes
  library(lmerTest)  # Pour obtenir les p-values
  library(emmeans)   # Pour les comparaisons post-hoc

source("../GlofiResistance/00_Metadata_Talyies.R")
  
Control <- "UT"

#### A - B-cell depletion ####
Condition = "B_cell_depletion_total"
  data_treatment <- data_talyies_full %>%
    select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
    mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
    dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
    right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
    dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
  # dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder))
  
  data_treatment_plot <- data_treatment 

  barplot <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
    geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, colour = "black", size = 0.5) +
    # geom_line(aes(group = Sample_code),
    #           position = position_dodge(0),
    #           color = "grey60",
    #           size = 0.2) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 1.4, stroke = 0.8, alpha = 0.9) +
    # stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
    #                    label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
    scale_fill_manual(values = sample_colors_all) + 
    theme_custom()+
    labs(
      title = "",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "B-cell Depletion (%)",
      fill = "Sample_code"  # Titre de la légende
    )  +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))
  
  # Construction du modèle
  model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
                data = data_treatment)
  # Calculer l'ICC (Coefficient de corrélation intraclasse) 
  # qui indique la proportion de variance due aux patients
  icc <- as.data.frame(VarCorr(model))$vcov[1] / 
    (as.data.frame(VarCorr(model))$vcov[1] + attr(VarCorr(model), "sc")^2)
  print(paste("ICC =", round(icc, 3)))
  
  # Moyennes estimées pour chaque traitement
  emm <- emmeans(model, ~ Treatment_reorder)
  print(emm)
  contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = Control)
  
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
    title = paste0("Treatment Effects on B-cell Depletion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.13)
  )

A_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                        ncol =2,
                        rel_widths = c(0.9, 0.7))

#### B - CD107a CD4 ####
Condition = "Per_CD107a_CD4"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
# dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, colour = "black", size = 0.5) +
  # geom_line(aes(group = Sample_code),
  #           position = position_dodge(0),
  #           color = "grey60",
  #           size = 0.2) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 1.4, stroke = 0.8, alpha = 0.9) +
  # stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
  #                    label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
  scale_fill_manual(values = sample_colors_all) + 
  theme_custom()+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "% of CD4+ T-cells CD107a+ (%)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(Per_CD107a_CD4 ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)
# Calculer l'ICC (Coefficient de corrélation intraclasse) 
# qui indique la proportion de variance due aux patients
icc <- as.data.frame(VarCorr(model))$vcov[1] / 
  (as.data.frame(VarCorr(model))$vcov[1] + attr(VarCorr(model), "sc")^2)
print(paste("ICC =", round(icc, 3)))

# Moyennes estimées pour chaque traitement
emm <- emmeans(model, ~ Treatment_reorder)
print(emm)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = Control)

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
    title = paste0("Treatment Effects on % of CD4+ CD107a+ (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.13)
  )

B_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(0.9, 0.7))

#### C - CD107a CD8 ####
Condition = "Per_CD107a_CD8"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
# dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, colour = "black", size = 0.5) +
  # geom_line(aes(group = Sample_code),
  #           position = position_dodge(0),
  #           color = "grey60",
  #           size = 0.2) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 1.4, stroke = 0.8, alpha = 0.9) +
  # stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
  #                    label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
  scale_fill_manual(values = sample_colors_all) + 
  theme_custom()+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "% of CD8+ T-cells CD107a+ (%)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(Per_CD107a_CD8 ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)
# Calculer l'ICC (Coefficient de corrélation intraclasse) 
# qui indique la proportion de variance due aux patients
icc <- as.data.frame(VarCorr(model))$vcov[1] / 
  (as.data.frame(VarCorr(model))$vcov[1] + attr(VarCorr(model), "sc")^2)
print(paste("ICC =", round(icc, 3)))

# Moyennes estimées pour chaque traitement
emm <- emmeans(model, ~ Treatment_reorder)
print(emm)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = Control)

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
    title = paste0("Treatment Effects on % of CD8+ CD107a+ (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.13)
  )

C_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(0.9, 0.7))

#### D - Granzyme B Secretion ####
Condition = "GrzB"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))

data_treatment_plot <- data_treatment 

barplot <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, colour = "black", size = 0.5) +
  # geom_line(aes(group = Sample_code),
  #           position = position_dodge(0),
  #           color = "grey60",
  #           size = 0.2) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 1.4, stroke = 0.8, alpha = 0.9) +
  # stat_compare_means(method = "t.test", paired = TRUE, aes(group = Treatment_reorder), hide.ns = FALSE, 
  #                    label = "p.format", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
  scale_fill_manual(values = sample_colors_all) + 
  theme_custom()+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "Granzyme B Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(GrzB ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)
# Calculer l'ICC (Coefficient de corrélation intraclasse) 
# qui indique la proportion de variance due aux patients
icc <- as.data.frame(VarCorr(model))$vcov[1] / 
  (as.data.frame(VarCorr(model))$vcov[1] + attr(VarCorr(model), "sc")^2)
print(paste("ICC =", round(icc, 3)))

# Moyennes estimées pour chaque traitement
emm <- emmeans(model, ~ Treatment_reorder)
print(emm)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = Control)

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
    title = paste0("Treatment Effects on Granzyme B Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.13)
  )

D_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(0.9, 0.7))


#### Correlation ####
Name <- "All_samples"
Disease_list <- c("FL","tFL","DLBCL")
Origin_list <- c("PBMC","LN")

data_correlation <- data_talyies_full %>% 
  dplyr::filter(Day == "D6", Disease %in% Disease_list) %>%
  dplyr::filter(Origin %in% Origin_list) %>% 
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  dplyr::filter(Treatment_type == "Single") %>%
  dplyr::filter(B_cell_depletion_total > -50) %>% 
  dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))

##### Per_CD8_D0 ######
{Current_day = "D0"
Condition = "Per_CD8"
# Current_treatment = "αTIGIT 0.1 µg/mL"

Current_treatment = "αCD20-TCB 0.1 nM"

markers_d0 <- data_talyies_full %>%
  dplyr::filter(Day %in% Current_day & Pop == "CD3_CD8") %>%
  select(Sample_code,!!sym(Condition)) 

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_correlation %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
  dplyr::filter(Treatment_reorder %in% Current_treatment)

# Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
data_combined <- data_depletion_d6 %>%
  left_join(markers_d0, by = "Sample_code") %>% 
  dplyr::filter(Sample_code %in% Sample_list)

data_unnormalized <- data_combined %>%
  dplyr::filter(Treatment_reorder == Current_treatment) %>% 
  mutate(value = (!!sym(Condition))) %>%
  dplyr::filter(B_cell_depletion_total > -50) 


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
  scale_x_continuous(limits = c(NA, NA)) + # Définir les limites de l'axe x
  labs(title ="",
       x = paste0("B-cell Depletion (%) \n for ",Current_treatment),
       y = ("CD8+ T-cells at D0 (%)"))+
  theme_custom()+
  theme(legend.position = "")

correlation_plot_CD8_D0 <- p + 
  annotation_custom(grob = textGrob(
    label = sprintf("\nR = %.4f\n𝑝 = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.95, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 9, col = "black")))

}

##### Per_CD8_D3 ######
{Current_day = "D3"
Condition = "Per_CD8"
# Current_treatment = "αTIGIT 0.1 µg/mL"

Current_treatment = "αCD20-TCB 0.1 nM"

markers_d0 <- data_talyies_full %>%
  dplyr::filter(Day %in% Current_day & Pop == "CD3_CD8") %>%
  select(Sample_code,!!sym(Condition)) 

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_correlation %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
  dplyr::filter(Treatment_reorder %in% Current_treatment)

# Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
data_combined <- data_depletion_d6 %>%
  left_join(markers_d0, by = "Sample_code") %>% 
  dplyr::filter(Sample_code %in% Sample_list)

data_unnormalized <- data_combined %>%
  dplyr::filter(Treatment_reorder == Current_treatment) %>% 
  mutate(value = (!!sym(Condition))) %>%
  dplyr::filter(B_cell_depletion_total > -50) 


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
  scale_x_continuous(limits = c(NA, NA)) + # Définir les limites de l'axe x
  labs(title ="",
       x = paste0("B-cell Depletion (%) \n for ",Current_treatment),
       y = ("CD8+ T-cells at D3 (%)"))+
  theme_custom()+
  theme(legend.position = "")

correlation_plot_CD8_D3 <- p + 
  annotation_custom(grob = textGrob(
    label = sprintf("\nR = %.4f\n𝑝 = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.95, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 9, col = "black")))

}

##### Ratio_cell_count_D6_D3 ######
{Current_day = "D3"
Condition = "Cell_count_PDLS"
# Current_treatment = "αTIGIT 0.1 µg/mL"

Current_treatment = "αCD20-TCB 0.1 nM"

markers_d0 <- data_talyies_full %>%
  dplyr::filter(Day %in% Current_day & Pop == "CD3_CD8") %>%
  select(Sample_code,!!sym(Condition)) %>% 
  rename(Cell_count_PDLS_D3=Cell_count_PDLS )

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_correlation %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder,Cell_count_PDLS) %>% 
  dplyr::filter(Treatment_reorder %in% Current_treatment) %>% 
  rename(Cell_count_PDLS_D6=Cell_count_PDLS )

# Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
data_combined <- data_depletion_d6 %>%
  left_join(markers_d0, by = "Sample_code") %>% 
  dplyr::filter(Sample_code %in% Sample_list) %>% 
  mutate(Ratio_cell_count_D6_D3 = Cell_count_PDLS_D6 / Cell_count_PDLS_D3)

data_unnormalized <- data_combined %>%
  dplyr::filter(Treatment_reorder == Current_treatment) %>% 
  mutate(value = (Ratio_cell_count_D6_D3)) %>%
  dplyr::filter(B_cell_depletion_total > -50) 


correlation_result <- cor.test(data_unnormalized$value, data_unnormalized$B_cell_depletion_total)
r_value <- correlation_result$estimate  # Coefficient de corrélation
p_value <- correlation_result$p.value   # Valeur p

##PLOT##
p <- ggplot(data_unnormalized, aes(x = B_cell_depletion_total, y = Ratio_cell_count_D6_D3,fill)) +
  geom_smooth(method = "lm", color = "black", fill = "gray80", size= 0.8) +
  geom_point(aes(fill = Sample_code),
             shape = 21,
             color = "black",
             position = position_dodge(0),
             size = 2.5, stroke = 0.4, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all) +  
  scale_x_continuous(limits = c(NA, NA)) + # Définir les limites de l'axe x
  labs(title ="",
       x = paste0("B-cell Depletion (%) \n for ",Current_treatment),
       y = ("cell number fold change \n per PDLS (D6/D3)"))+
  theme_custom()+
  theme(legend.position = "")

correlation_plot_FC_cell<- p + 
  annotation_custom(grob = textGrob(
    label = sprintf("\nR = %.4f\n𝑝 = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.95, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 9, col = "black")))

}
##### Per_CD107a_CD4 ######
{
  Current_day = "D6"
Condition = "Per_CD107a_CD4"
Current_treatment = "αTIGIT 10 µg/mL"

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_correlation %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder,Per_CD107a_CD4) %>% 
  dplyr::filter(Treatment_reorder %in% Current_treatment)

data_unnormalized <- data_depletion_d6 %>%
  dplyr::filter(Treatment_reorder == Current_treatment) %>% 
  mutate(value = (!!sym(Condition))) %>%
  dplyr::filter(B_cell_depletion_total > -50) 


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
  scale_x_continuous(limits = c(NA, NA)) + # Définir les limites de l'axe x
  labs(title ="",
       x = paste0("B-cell Depletion (%) \n for ",Current_treatment),
       y = ("% of CD4+ T-cells CD107a+ (%)"))+
  theme_custom()+
  theme(legend.position = "")

correlation_plot_CD4_Cd107a_tigit <- p + 
  annotation_custom(grob = textGrob(
    label = sprintf("\nR = %.4f\n𝑝 = %.2e", r_value, p_value),
    x = unit(0.05, "npc"), y = unit(0.95, "npc"),
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 9, col = "black")))

}
Correlation_fig4 <- plot_grid(correlation_plot_CD8_D0,correlation_plot_CD8_D3,correlation_plot_FC_cell,correlation_plot_CD4_Cd107a_tigit,
                            ncol = 4)
####Montage ####
Figure_4_total <- plot_grid(A_final_plot,B_final_plot,C_final_plot,D_final_plot,Correlation_fig4,
                            nrow = 5, 
                            labels = c("A","B","C","D","E"),
                            rel_heights = c(1,1,1,1,0.7))

if (CRCT_Share == TRUE) {
ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/Figure_4_total.png"),
       plot = Figure_4_total,
       device = "png",
       width = 36,        # largeur A4 en cm
       height = 48 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)
}

ggsave(
  filename = here("Figure","Figure4","Figure_4_total.png"),
  plot = Figure_4_total,
       device = "png",
       width = 36,        # largeur A4 en cm
       height = 48 ,     # hauteur A4 en cm
       units = "cm",
       dpi = 300)

