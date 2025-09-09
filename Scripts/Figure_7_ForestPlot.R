library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggnewscale)
  library(ggpubr)
  library(grid)
  library(lme4)      # Pour les modèles mixtes
  library(lmerTest)  # Pour obtenir les p-values
  library(emmeans)   # Pour les comparaisons post-hoc
  
source("~/Postdoc/R/GlofiResistance/00_Metadata_Talyies.R")

####Calcul the forest plot data ####
#####All Samples####
Control <- "UT"
  # Cluster_response_cat_filter <- "ADC_TCB_Medium_responders"
  
  Condition = "B_cell_depletion_total"
  data_treatment <- data_talyies_full %>%
    select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
    mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
    dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
    right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
    # dplyr::filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
    dplyr::filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
  # dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder))
  
  ##With Cluster Classification
  # data_treatment <- data_talyies_full %>%
  #   select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  #   mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  #   dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  #   right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  #   dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder)) %>% 
  #   left_join(df_Cluster_response_cat, by = "Sample_code") %>% 
  #   dplyr::filter(Cluster_response_cat %in% Cluster_response_cat_filter)
  
    # Construction du modèle
  model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
                data = data_treatment)

  emm <- emmeans(model, ~ Treatment_reorder)

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
  
  # # Afficher le tableau formaté
  # print(contrast_table)

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
  mutate(Treatment = factor(Treatment, levels = unique(Treatment)))

order_treatment <- contrast_df %>%
  select(Treatment) %>%
  distinct() %>%
  pull(Treatment)

All_treatment <- as_tibble(unique(contrast_df$Treatment)) %>% 
  rename(Treatment = value) %>% 
  mutate(Treatment = factor(Treatment, levels = order_treatment))


# Forest plot

forest_plot_all <- ggplot(contrast_df, aes(x = Treatment, 
                                       y = estimate, 
                                       ymin = lower, 
                                       ymax = upper, 
                                       color = signif)) +
  geom_pointrange(size = 0.9, fatten = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")+
  scale_color_manual(values = c("p < 0.05" = "#D73027", "ns" = "gray50"),
                     na.translate = FALSE ) +
  labs(
    x = "Treatment",
    color = "Significance",
    y = "All Samples"
  ) +
  theme_custom()+
  theme(
    axis.text.x = element_blank(),
    legend.position="inside",
    legend.position.inside = c(0.92,0.85),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_blank(),
    plot.margin = margin(t = 0.5,  # Top margin
                         r = 2,  # Right margin
                         b = 1,  # Bottom margin
                         l = 1))+ # Left margin 
  ylim(-30,125) +
  annotation_custom(
    grob = textGrob("Your Text Here", gp = gpar(fontsize = 20)),  # size in points
    xmin = -Inf, xmax = Inf,
    ymin = unit(0.7, "npc"), ymax = unit(0.7, "npc")              # 90% from bottom
  )
#####Low Responder Forest Plot####
{Control <- "UT"
Cluster_response_cat_filter <- "Low_Responder"

Condition = "B_cell_depletion_total"

#With Cluster Classification
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  # dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder)) %>%
  left_join(df_Cluster_response_cat, by = "Sample_code") %>%
  dplyr::filter(Cluster_response_cat %in% Cluster_response_cat_filter)

# Construction du modèle
model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)

emm <- emmeans(model, ~ Treatment_reorder)

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

# # Afficher le tableau formaté
# print(contrast_table)

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
  mutate(Treatment = factor(Treatment, levels = unique(Treatment)))

contrast_df <- contrast_df %>% right_join(All_treatment,by = "Treatment") %>% 
  mutate(Treatment = factor(Treatment, levels = All_treatment$Treatment))


# Forest plot

forest_plot_low <- ggplot(contrast_df, aes(x = Treatment, 
                                       y = estimate, 
                                       ymin = lower, 
                                       ymax = upper, 
                                       color = signif)) +
  geom_pointrange(size = 0.9, fatten = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")+
  scale_color_manual(values = c("p < 0.05" = "#D73027", "ns" = "gray50")) +
  labs(
    x = "Treatment",
    color = "Significance",
    y = "Low Responders"
  ) +
  theme_custom()+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_blank(),
    plot.margin = margin(t = 0.5,  # Top margin
                         r = 2,  # Right margin
                         b = 1,  # Bottom margin
                         l = 1))+
  ylim(-30,125)
  # Left margin 

}

#####Medium Responder Forest Plot####
{Control <- "UT"
Cluster_response_cat_filter <- "Medium_Responder"

Condition = "B_cell_depletion_total"

#With Cluster Classification
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  # dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder)) %>%
  left_join(df_Cluster_response_cat, by = "Sample_code") %>%
  dplyr::filter(Cluster_response_cat %in% Cluster_response_cat_filter)

# Construction du modèle
model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)

emm <- emmeans(model, ~ Treatment_reorder)

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

# # Afficher le tableau formaté
# print(contrast_table)

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
  mutate(Treatment = factor(Treatment, levels = unique(Treatment)))

contrast_df <- contrast_df %>% right_join(All_treatment) %>% 
  mutate(Treatment = factor(Treatment, levels = All_treatment$Treatment))



# Forest plot

forest_plot_medium <- ggplot(contrast_df, aes(x = Treatment, 
                                       y = estimate, 
                                       ymin = lower, 
                                       ymax = upper, 
                                       color = signif)) +
  geom_pointrange(size = 0.9, fatten = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")+
  scale_color_manual(values = c("p < 0.05" = "#D73027", "ns" = "gray50")) +
  labs(
    x = "Treatment",
    color = "Significance",
    y="Medium Responders"
  ) +
  theme_custom()+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 0.5,  # Top margin
                         r = 2,  # Right margin
                         b = 1,  # Bottom margin
                         l = 1))+ # Left margin 
  ylim(-30,125)

}


#####High Responder Forest Plot####
{Control <- "UT"
Cluster_response_cat_filter<- "High_Responder"

Condition = "B_cell_depletion_total"

#With Cluster Classification

data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  dplyr::filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>%
  # dplyr::filter(grepl("TCB 0,01 nM",Treatment_reorder)) %>%
  left_join(df_Cluster_response_cat, by = "Sample_code") %>%
  dplyr::filter(Cluster_response_cat %in% Cluster_response_cat_filter)

# Construction du modèle
model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)

emm <- emmeans(model, ~ Treatment_reorder)

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

# # Afficher le tableau formaté
# print(contrast_table)

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
  mutate(Treatment = factor(Treatment, levels = unique(Treatment)))

# Forest plot
contrast_df <- contrast_df %>% right_join(All_treatment) %>% 
  mutate(Treatment = factor(Treatment, levels = All_treatment$Treatment))


forest_plot_high <- ggplot(contrast_df, aes(x = Treatment, 
                                       y = estimate, 
                                       ymin = lower, 
                                       ymax = upper, 
                                       color = signif)) +
  geom_pointrange(size = 0.9, fatten = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")+
  scale_color_manual(values = c("p < 0.05" = "#D73027", "ns" = "gray50")) +
  labs(
    x = "",
    y="High Responders",
    color = "Significance"
  ) +
  theme_custom()+
  theme(
    axis.text.x = element_text(angle=90, hjust=1, size = 7,vjust=0.4),
    legend.position = "none",
    title = element_blank(),
  plot.margin = margin(t = 0.5,  # Top margin
                       r = 2,  # Right margin
                       b = 20,  # Bottom margin
                       l = 1))+ # Left margin 
  ylim(-30,125)

}

###Montage ####
{combined_plots <- plot_grid(
  plot_grid(forest_plot_all, ncol = 1, labels = c("A"),  # Deuxième ligne
  plot_grid(forest_plot_low, ncol = 1, labels = c("B")), # 3-4 Lignes
  plot_grid(forest_plot_medium, ncol = 1, labels = c("C")),
  plot_grid(forest_plot_high, ncol = 1, labels = c("D")),
  nrow=4,
  rel_heights = c(0.335,0.335,0.335,0.6)))

Figure_7_total <- ggdraw() +
  draw_label("Estimated Effect ± 95% CI", 
             x = 0.02, y = 0.6, angle = 90, 
             fontface = "bold", size = 12) +
  draw_plot(combined_plots, x = 0.05, y = 0, width = 0.95, height = 1)

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure7/Figure_7_total.png",
  plot = Figure_7_total,
  device = "png",
  width = 31,        # largeur A4 en cm
  height = 40 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
}

