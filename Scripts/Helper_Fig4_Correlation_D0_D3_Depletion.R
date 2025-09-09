##### Correlation Condition WITHOUT normalization and B Cell Depletion ######
# Préparer les données
{Current_day = "D0"
  Condition = "Per_CD8"
  # Current_treatment = "αTIGIT 0.1 µg/mL"
  
  Current_treatment = "αCD20-TCB 0.1 nM"
  
markers_d0 <- data_talyies_full %>%
  filter(Day %in% Current_day & Pop == "CD3_CD8") %>%
  # select(Sample_code, Per_CD3, Per_CD4, Per_CD8, Per_NK,Per_TGD,Per_CD11b) %>%
  # select(Sample_code, Per_CD22_total, Per_CD20, Per_CD19, Per_CD22_CD10_plus) %>%
  # select(Sample_code, Area_mean, Roundness_mean, Viability, Cell_count_PDLS) %>%
  # select(Sample_code, Tfr, Tfh, Treg, CD3_Naive,CD3_Central_Memory,CD3_Effector_Memory,CD3_TEMRA) %>% 
  select(Sample_code,!!sym(Condition)) 

# Préparer les données de déplétion à D6
data_depletion_d6 <- data_correlation %>%
  select(Sample_code,Treatment,B_cell_depletion_total,Treatment_reorder) %>% 
  filter(Treatment_reorder %in% Current_treatment)

# Fusionner les données des marqueurs à D0 avec les données de déplétion à D6
data_combined <- data_depletion_d6 %>%
  left_join(markers_d0, by = "Sample_code") %>% 
  filter(Sample_code %in% Sample_list)

data_unnormalized <- data_combined %>%
  filter(Treatment_reorder == Current_treatment) %>% 
  mutate(value = (!!sym(Condition))) %>%
  filter(B_cell_depletion_total > -50) 

 
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
  labs(title = paste0("Correlation between",Condition," \n and B Cell Depletion at ",Current_day),
       x = paste0("B Cell Depletion (%) \n for ",Current_treatment),
       y = paste0("",Condition," at ",Current_day) )+
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
