# Charger les packages nécessaires
library(lme4)      # Pour les modèles mixtes
library(lmerTest)  # Pour obtenir les p-values
library(emmeans)   # Pour les comparaisons post-hoc

source("../GlofiResistance/00_Metadata_Talyies.R")

Control <- "UT"

#### A - Area_mean ####
Condition = "Area_mean"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
# dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition)/1e7)) +
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
    y = "PDLS Area (mm²)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(Area_mean ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on PDLS Area (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.17)
  )

A_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))

#### B - TNFa ####
Condition = "TNFa"
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
    y = "TFNα Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(TNFa ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on TFNα Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.17)
  )

B_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))

#### C - INFY ####
Condition = "INFY"
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
    y = "INFγ Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(INFY ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on INFγ Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.17)
  )

C_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))

#### D - IL10 Secretion ####
Condition = "IL10"
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
    y = "IL10 Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(IL10 ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on IL10 Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.17)
  )

D_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))

#### E - IL6 Secretion ####
Condition = "IL6"
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
    y = "IL6 Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(IL6 ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on IL6 Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.17)
  )

E_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))

#### F - IL8 Secretion ####
Condition = "IL8"
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
    y = "IL8 Secretion (pg/mL)",
    fill = "Sample_code"  # Titre de la légende
  )  +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"))

# Construction du modèle
model <- lmer(IL8 ~ Treatment_reorder + (1|Sample_code), 
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
    title = paste0("Treatment Effects on IL8 Secretion (vs ",Control,")"),
    x = "",
    y = "Estimated Effect ± 95% CI",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.17)
  )

F_final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                          ncol =2,
                          rel_widths = c(01, 0.6))


####Montage ####
Figure_4_Suppl <- plot_grid(A_final_plot,B_final_plot,C_final_plot,D_final_plot,E_final_plot,F_final_plot,
                            nrow = 6, 
                            labels = c("A","B","C","D","E","F"),
                            rel_heights = c(1,1,1,1,1,1))

if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure4/Figure_4_Suppl.png"),
         plot = Figure_4_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure4","Figure_4_Suppl.png"),
  plot = Figure_4_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)

