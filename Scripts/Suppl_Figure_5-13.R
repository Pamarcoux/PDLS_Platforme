library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggnewscale)
library(impute)
library(ggpubr)
library(grid)
library(lme4)      # Pour les modèles mixtes
library(lmerTest)  # Pour obtenir les p-values
library(emmeans)   # Pour les comparaisons post-hoc


source("../GlofiResistance/00_Metadata_Talyies.R")

sample_colors_df <- tibble(Sample_code = names(sample_colors_all),
                           sample_color = unname(sample_colors_all))

conditions <- c('B-cell Depletion', '% of CD4+ T-cells CD107a+', '% of CD8+ T-cells CD107a+',
                'TNFα Secretion', 'INFγ Secretion', 'Granzyme B Secretion', 'IL8 Secretion', 'IL6 Secretion')

Condition = "Response_type_total"

#### Suppl Fig 5 ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(Treatment_type == "Single")
  # dplyr::filter(grepl("PD-1",Treatment_reorder)) 

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig5_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig5_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig5_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig5_Suppl_Patient)
}
###### Montage #####

Figure_5_Suppl <- plot_grid(Fig5_Suppl_Sample,Fig5_Suppl_Patient,
                            nrow = 2, 
                            labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                            hjust = -0.1)
print(Figure_5_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_5_Suppl.png"),
         plot = Figure_5_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_5_Suppl.png"),
  plot = Figure_5_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)

###Suppl Fig6 Forest plot Glofi ####
##### A - αCD20-TCB 0.01 nM ####
{Control <- "αCD20-TCB 0.01 nM"
Condition = "B_cell_depletion_total"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>% 
  dplyr::filter(grepl("TCB 0.01 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot_0.01 <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
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
    Treatment = sub(paste0("^\\((.*)\\) - \\(", Control, "\\)$"), "\\1", contrast),  
    signif = ifelse(p.value < 0.05, "p < 0.05", "ns"),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) %>% 
  mutate(Treatment = factor(Treatment, levels = rev(unique(Treatment))))

# Forest plot
forest_plot_0.01 <- ggplot(contrast_df, aes(x = Treatment, 
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
    legend.position.inside = c(0.95, 0.93)
  )
}

##### B - αCD20-TCB 0.1 nM ####
{Control <- "αCD20-TCB 0.1 nM"
Condition = "B_cell_depletion_total"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>% 
  dplyr::filter(grepl("TCB 0.1 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot_0.1 <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
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
    Treatment = sub(paste0("^\\((.*)\\) - \\(", Control, "\\)$"), "\\1", contrast),  
    signif = ifelse(p.value < 0.05, "p < 0.05", "ns"),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) %>% 
  mutate(Treatment = factor(Treatment, levels = rev(unique(Treatment))))

# Forest plot
forest_plot_0.1 <- ggplot(contrast_df, aes(x = Treatment, 
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
    legend.position.inside = c(0.95, 0.93)
  )
}
##### C - αCD20-TCB 1 nM ####
{Control <- "αCD20-TCB 1 nM"
Condition = "B_cell_depletion_total"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1")) %>% 
  dplyr::filter(grepl("TCB 1 nM",Treatment_reorder))

data_treatment_plot <- data_treatment 

barplot_1 <- ggplot(data_treatment_plot, aes(x=Treatment_reorder , y=!!sym(Condition))) +
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
    Treatment = sub(paste0("^\\((.*)\\) - \\(", Control, "\\)$"), "\\1", contrast),
    signif = ifelse(p.value < 0.05, "p < 0.05", "ns"),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) %>% 
  mutate(Treatment = factor(Treatment, levels = rev(unique(Treatment))))

# Forest plot
forest_plot_1 <- ggplot(contrast_df, aes(x = Treatment, 
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
    legend.position.inside = c(0.95, 0.93)
  )
}

#####Montage####
Figure_6_Suppl <- plot_grid(forest_plot_0.01,forest_plot_0.1,forest_plot_1,
                            nrow = 3, 
                            labels = c("A","B","C"),
                            rel_heights = c(1,1,1))

if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_6_Suppl.png"),
         plot = Figure_6_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_6_Suppl.png"),
  plot = Figure_6_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 22 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
#### Suppl Fig 8 TIGIT ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
dplyr::filter(grepl("TIGIT",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig8_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig8_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig8_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig8_Suppl_Patient)

###### Montage #####

Figure_8_Suppl <- plot_grid(Fig8_Suppl_Sample,Fig8_Suppl_Patient,
                            nrow = 2, 
                            labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                            hjust = -0.1)
print(Figure_8_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_8_Suppl.png"),
         plot = Figure_8_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_8_Suppl.png"),
  plot = Figure_8_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 9 PD-1 ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("PD-1",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig9_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig9_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig9_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig9_Suppl_Patient)

###### Montage #####

Figure_9_Suppl <- plot_grid(Fig9_Suppl_Sample,Fig9_Suppl_Patient,
                            nrow = 2, 
                            labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                            hjust = -0.1)
print(Figure_9_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_9_Suppl.png"),
         plot = Figure_9_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_9_Suppl.png"),
  plot = Figure_9_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 10 CD19-411B ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("CD19-4-1BB",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig10_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig10_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig10_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig8_Suppl_Patient)

###### Montage #####

Figure_10_Suppl <- plot_grid(Fig10_Suppl_Sample,Fig10_Suppl_Patient,
                            nrow = 2, 
                            labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                            hjust = -0.1)
print(Figure_10_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_10_Suppl.png"),
         plot = Figure_10_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_10_Suppl.png"),
  plot = Figure_10_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 11 CD19-CD28 ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("CD19-CD28",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig11_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig11_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig11_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig11_Suppl_Patient)

###### Montage #####

Figure_11_Suppl <- plot_grid(Fig11_Suppl_Sample,Fig11_Suppl_Patient,
                             nrow = 2, 
                             labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                             hjust = -0.1)
print(Figure_11_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_11_Suppl.png"),
         plot = Figure_11_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_11_Suppl.png"),
  plot = Figure_11_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 12 PD-L1 ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("PD-L1",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig12_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig12_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig12_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig12_Suppl_Patient)

###### Montage #####

Figure_12_Suppl <- plot_grid(Fig12_Suppl_Sample,Fig12_Suppl_Patient,
                             nrow = 2, 
                             labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                             hjust = -0.1)
print(Figure_12_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_12_Suppl.png"),
         plot = Figure_12_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_12_Suppl.png"),
  plot = Figure_12_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 13 CD79 ####
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("CD79",Treatment_reorder))

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

##### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig13_Suppl_Sample <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig13_Suppl_Sample)

####Heatmap per Patient #####
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
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig13_Suppl_Patient <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )
print(Fig13_Suppl_Patient)

###### Montage #####

Figure_13_Suppl <- plot_grid(Fig13_Suppl_Sample,Fig13_Suppl_Patient,
                             nrow = 2, 
                             labels = c("A — Z-scores per Sample","B — Z-scores per Treatment"),
                             hjust = -0.1)
print(Figure_13_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_13_Suppl.png"),
         plot = Figure_13_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_13_Suppl.png"),
  plot = Figure_13_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 48 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)
}

#### Suppl Fig 7  ####
##### αCD20-TCB 0.01 nM #### 
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("αCD20-TCB 0.01",Treatment_reorder))

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


###### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig7_Suppl_Sample_0.01 <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

##### Montage #####

Figure_7_Suppl_0.01 <- plot_grid(Fig7_Suppl_Sample_0.01,barplot_0.01,
                             nrow = 1, 
                             labels = c("A","B"),
                             rel_widths = c(1,0.8),
                             hjust = -0.1)
print(Figure_7_Suppl_0.01)
}

##### αCD20-TCB 0.1 nM #### 
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("αCD20-TCB 0.1",Treatment_reorder))

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


###### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig7_Suppl_Sample_0.1 <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

##### Montage #####

Figure_7_Suppl_0.1 <- plot_grid(Fig7_Suppl_Sample_0.1,barplot_0.1,
                                 nrow = 1, 
                                 labels = c("C","D"),
                                 rel_widths = c(1,0.8),
                                 hjust = -0.1)
print(Figure_7_Suppl_0.1)
}

##### αCD20-TCB 1 nM #### 
{data_heatmap <- data_talyies_full %>%
  rename('B-cell Depletion'=B_cell_depletion_total,
         '% of CD4+ T-cells CD107a+' = Per_CD107a_CD4_normalized,
         '% of CD8+ T-cells CD107a+' = Per_CD107a_CD8_normalized,
         'TNFα Secretion' = TNFa_normalized,
         'INFγ Secretion' = INFY_normalized,
         'Granzyme B Secretion'= GrzB_normalized,
         'IL8 Secretion'=IL8_normalized,
         'IL6 Secretion'=IL6_normalized ) %>% 
  mutate(B_cell_depletion_total='B-cell Depletion') %>% 
  select(Treatment, Disease, Sample_code, Day, all_of(conditions),B_cell_depletion_total) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  dplyr::filter(grepl("αCD20-TCB 1",Treatment_reorder))

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


###### Heatmap per Sample #####
#Creation Heatmaps
for (i in seq_along(conditions)) {
  cond <- conditions[i]
  plot_data <- dplyr::filter(data_heatmap_pivot, Condition == cond) %>% 
    group_by(Sample_code) %>%
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
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold"),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      plot.margin = margin(0, 0, 0, 0)
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
Fig7_Suppl_Sample_1 <- wrap_plots(plots, ncol = ncol_plot, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6, face= "bold"),
    legend.title = element_text(size = 7, face="bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

##### Montage #####

Figure_7_Suppl_1 <- plot_grid(Fig7_Suppl_Sample_1,barplot_1,
                                 nrow = 1, 
                                 labels = c("E","F"),
                                 rel_widths = c(1,0.8),
                                 hjust = -0.1)
print(Figure_7_Suppl_1)
}

#####Montage ####
Figure_7_Suppl <- plot_grid(Figure_7_Suppl_0.01,Figure_7_Suppl_0.1,Figure_7_Suppl_1,
                             nrow = 3)
print(Figure_7_Suppl)


if (CRCT_Share == TRUE) {
  ggsave(filename =  paste0("/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure5/Figure_7_Suppl.png"),
         plot = Figure_7_Suppl,
         device = "png",
         width = 36,        # largeur A4 en cm
         height = 48 ,     # hauteur A4 en cm
         units = "cm",
         dpi = 300)
}

ggsave(
  filename = here("Figure","Figure5","Figure_7_Suppl.png"),
  plot = Figure_7_Suppl,
  device = "png",
  width = 36,        # largeur A4 en cm
  height = 32 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300)


