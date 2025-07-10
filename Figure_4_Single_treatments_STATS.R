{# Charger les packages nécessaires
library(lme4)      # Pour les modèles mixtes
library(lmerTest)  # Pour obtenir les p-values
library(emmeans)   # Pour les comparaisons post-hoc

  Control <- "UT"
  # Cluster_response_cat_filter <- "ADC_TCB_Medium_responders"
  
Condition = "B_cell_depletion_total"
data_treatment <- data_talyies_full %>%
  select(Treatment, Disease, Sample_code, Day,B_cell_depletion_total,Condition) %>%
  mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
  filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL") & B_cell_depletion_total > -40) %>%
  right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
  # filter(Treatment_type == "Single" & !Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
  filter(!Sample_code %in% c("DLBCL1_LN1","FL10_PB1"))
  # filter(grepl("TCB 0,01 nM",Treatment_reorder))

length(setdiff(data_talyies_full$Sample_code,data_treatment$Sample_code))


##With Cluster Classification
# data_treatment <- data_talyies_full %>%
#   select(Treatment, Disease, Sample_code, Day, !!sym(Condition),B_cell_depletion_total) %>%
#   mutate(!!sym(Condition) := replace(!!sym(Condition), !!sym(Condition) == "", NA)) %>%
#   filter(Day == "D6", Disease %in% c("FL", "DLBCL", "tFL")) %>%
#   right_join(all_combinations, by = c("Sample_code", "Treatment")) %>% 
#   filter(grepl("TCB 0,01 nM",Treatment_reorder)) %>% 
#   left_join(df_Cluster_response_cat, by = "Sample_code") %>% 
#   filter(Cluster_response_cat %in% Cluster_response_cat_filter)
#   


# valid_samples <- data_treatment %>%
#   group_by(Sample_code) %>%
#   filter(!is.na(!!sym(Condition))) %>%
#   summarize(treatments_count = n_distinct(Treatment_reorder)) %>%
#   filter(treatments_count == 2) %>%
#   pull(Sample_code)

data_treatment_plot <- data_treatment 
# %>% 
#   filter(Sample_code %in% valid_samples)

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


# Construction du modèle
model <- lmer(B_cell_depletion_total ~ Treatment_reorder + (1|Sample_code), 
              data = data_treatment)

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
}
library(ggplot2)
library(dplyr)

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

final_plot <- plot_grid(barplot, forest_plot,labels = c(""),
                        ncol =2,
                        rel_widths = c(1, 0.6))

print(final_plot)
VarCorr(model)
print(paste("ICC =", round(icc, 3)))

print(contrast_table)
