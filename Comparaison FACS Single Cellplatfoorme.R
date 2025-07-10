# Tcells_subset <- subset(Tcells, Cell_type == "CD4" & Sample_code %in% Sample_list)
Sample_metadata_subset <- Sample_metadata

colors <- c("lightskyblue2", "lightseagreen")  # Couleurs pour chaque variable

data_sample_plot_population <- data_talyies_flex %>% select(Sample_code_paper,Sample_code,Day,Per_CD22_total,Per_CD11b,Per_CD4,Per_CD8,Per_NK,Per_CD3,Per_CD20,Treatment) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD11b","CD3","CD4", "CD8", "NK", "TGD","NKT"))
population_CD3_CD20 <- data_sample_plot_population %>%
  filter(Day == "D0", variable %in% c("CD3", "CD20")) %>%
  rename(Proportion = value)

###Correlation Pop FACS vs SingleCell
{Population <- "CD8"

###Pop Facs
data_sample_plot_population_filter <-  data_sample_plot_population %>% filter(Day == "D0" & !variable %in% c("CD3","CD20")) %>% 
  mutate(variable = case_when(
    variable == "CD22" ~ "Bcells",
    variable == "CD11b" ~ "Myeloid",
    TRUE ~ variable # Conserve les autres valeurs inchangées
  )) %>% 
  group_by(Sample_code) %>% 
  mutate(Proportion = value / sum(value, na.rm= TRUE) * 100) %>% 
  filter(variable == Population) %>% 
  select(-variable,-Treatment,-Day) %>% 
  rename(Percentages_Facs = Proportion)

#####Repartition Pop SC ######
data_population<- read_csv("Repartition_cells.csv") %>% 
  select(-Response_type_CD10,-B_cell_depletion_CD22_CD10_plus_glofi,Sample_code) 

data_plot_population <- data_population %>%
  mutate(Total = rowSums(across(c(CD8, CD4, Other, NK, Bcells, Myeloid)), na.rm = TRUE)) %>%
  mutate(across(c(CD8, CD4, Other, NK, Bcells, Myeloid)) / Total *100) %>%
  pivot_longer(cols = c(CD8, CD4, Other, NK, Bcells, Myeloid),
               values_to = "Percentage",
               names_to = "Cell_type") %>% 
  group_by(Cell_type,Sample_code,Total) %>% 
  summarize(Percentage_Single_cells = mean(Percentage, na.rm = TRUE),, .groups = "drop")

data_plot_population$Cell_type <- factor(data_plot_population$Cell_type, levels = c("Bcells", "Myeloid","CD4", "CD8", "NK", "TGD","Other"))


grouped_data_cluster <- subset (data_plot_population, Cell_type == Population) %>% 
  left_join(data_sample_plot_population_filter)

# Calculer le modèle linéaire pour obtenir R^2
fit_cluster <- lm(Percentages_Facs ~ Percentage_Single_cells, data = grouped_data_cluster)
summary_fit_cluster <- summary(fit_cluster)
r_squared_cluster <- summary_fit_cluster$r.squared

ggplot(data = grouped_data_cluster, aes(x = Percentages_Facs, y = Percentage_Single_cells, label = Sample_code)) +
  geom_point(aes(color = Sample_code), size = 3, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE,shape = 21, stroke = 0.5)+
  geom_text_repel(aes(label = Sample_code))+
  scale_color_manual(values = rev(colors))+
  geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
  # geom_text_repel(size = 4, box.padding = 0.3, point.padding = 0.2, max.overlaps =50) +  # Annotations ajustées
  annotate("text", x = 25, y = 0, label = sprintf("R² = %.4f", r_squared_cluster), 
           size = 5, color = "black") + 
  labs(
    x = paste0(Population, "% FACS"),
    y = paste0(Population, " % Single Cell"),
    color = "Sample Code"  # Légende des couleurs
  ) +
  theme_classic() +  # Thème épuré
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Titre centré et en gras
    axis.title = element_text(size = 14),  # Titres des axes en gras
    axis.text = element_text(size = 12, color = "black"),  # Textes des axes en gras et lisibles
    legend = NA  # Position de la légende
  )
}
###Correlation Treg TFH FACS vs ScRNA
#### TREG TFh And TEMRA ######
{cluster <- paste0("CD4_Tfh_newcluster")
cluster_facs <- "Tfh"
  Sample_metadata_subset <- Sample_metadata %>% 
  filter(Cell_type == "CD4" & Sample_code %in% Sample_list) %>% 
  mutate(Cluster_name = if_else(True_Tfh == TRUE, "CD4_Tfh_newcluster", Cluster_name))

data_sample_plot_Tfonction <- data_talyies %>% filter(Day %in% c("D0","D3") & Pop == "CD22") %>% 
  select(Sample_code,Condition,Day,Tfh,Tfr,Treg,CD3_Naive,CD3_Central_Memory,CD3_Effector_Memory,CD3_TEMRA,Treatment,B_cell_depletion_total_glofi) %>% 
  filter(!!sym(Condition) %in% List_condition,Treatment == "UT") %>% 
  melt(id.vars = c(Condition,"Day","Sample_code","Treatment","B_cell_depletion_total_glofi")) %>% 
  mutate(variable = gsub("_Per", "", variable))

data_sample_plot_Tfonction_pop <- data_sample_plot_Tfonction %>% filter(variable == cluster_facs & Day =="D0") %>% 
  select(-variable,-Treatment,-Day) %>% 
  rename(Percentages_Facs = value)


grouped_data_cluster <- Sample_metadata_subset %>%
  group_by(Sample_code) %>%
  mutate(counttotal = n()) %>%
  ungroup() %>%
  group_by(B_cell_depletion_total_glofi, Sample_code, Cluster_name) %>%
  mutate(count = n(),
         percentage = count / counttotal *100) %>% 
  select(Sample_code,B_cell_depletion_total_glofi,count,counttotal,percentage,Cluster_name,Response_type_total) %>% 
  distinct()

grouped_data_cluster <- subset (grouped_data_cluster, Cluster_name == cluster) %>%
  left_join(data_sample_plot_Tfonction_pop)

# Calculer le modèle linéaire pour obtenir R^2
fit_cluster <- lm(Percentages_Facs ~ percentage, data = grouped_data_cluster)
summary_fit_cluster <- summary(fit_cluster)
r_squared_cluster <- summary_fit_cluster$r.squared

Plot_correlation_Tfh <- ggplot(data = grouped_data_cluster, aes(x = Percentages_Facs, y = percentage, label = Sample_code)) +
  geom_point(aes(color = Response_type_total), size = 3, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE,shape = 21, stroke = 0.5)+
  geom_text_repel(aes(label = Sample_code))+
  scale_color_manual(values = rev(colors))+
  geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
  # geom_text_repel(size = 4, box.padding = 0.3, point.padding = 0.2, max.overlaps =50) +  # Annotations ajustées
  annotate("text", x = 10, y = 0, label = sprintf("R² = %.4f", r_squared_cluster), 
           size = 5, color = "black") + 
  labs(
    x = "% FACS",
    y = paste0(cluster, " % Single Cell"),
    color = "Sample Code"  # Légende des couleurs
  ) +
  theme_classic() +  # Thème épuré
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Titre centré et en gras
    axis.title = element_text(size = 14),  # Titres des axes en gras
    axis.text = element_text(size = 12, color = "black"),  # Textes des axes en gras et lisibles
    legend = NA  # Position de la légende
  )

print(Plot_correlation_Tfh)
}
##########Correlation Depletion Pop
{Population <- "NK"
  
  ###Pop Facs
  data_sample_plot_population_filter <-  data_sample_plot_population %>% filter(Day == "D0" & !variable %in% c("CD3","CD20") & Response_type_total %in% c("High","Low")) %>% 
    # filter(!Sample_code %in% c("FL10_LN2","tFL5_PB1","FL3_PB1")) %>% 
    mutate(variable = case_when(
      variable == "CD22" ~ "Bcells",
      variable == "CD11b" ~ "Myeloid",
      TRUE ~ variable # Conserve les autres valeurs inchangées
    )) %>% 
    group_by(Sample_code) %>% 
    mutate(Proportion = value / sum(value, na.rm= TRUE) * 100) %>% 
    filter(variable == Population) %>% 
    select(-variable,-Treatment,-Day) 

  grouped_data_cluster <-  data_sample_plot_population_filter %>% left_join(data_response_type_CD10) 

  # Calculer le modèle linéaire pour obtenir R^2
  fit_cluster <- lm(B_cell_depletion_total_glofi ~ Proportion, data = grouped_data_cluster)
  summary_fit_cluster <- summary(fit_cluster)
  r_squared_cluster <- summary_fit_cluster$r.squared
  
  ggplot(data = grouped_data_cluster, aes(x = B_cell_depletion_total_glofi, y = Proportion, label = Sample_code)) +
    geom_point(aes(color = Sample_code), size = 3, show.legend = FALSE) +
    geom_point(size = 3, show.legend = FALSE,shape = 21, stroke = 0.5)+
    geom_text_repel(aes(label = Sample_code))+
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
    # geom_text_repel(size = 4, box.padding = 0.3, point.padding = 0.2, max.overlaps =50) +  # Annotations ajustées
    annotate("text", x = 25, y = -1, label = sprintf("R² = %.4f", r_squared_cluster), 
             size = 5, color = "red2") + 
    labs(title = paste0("Correlation between ", Population, " % by FACS and Glofitamab depletion"),
      x = paste0("B Cell Depletion (%)"),
      y = paste0("% of ", Population, " by Facs" ),
              color = "Sample Code"  # Légende des couleurs
    ) +
    theme_classic() +  # Thème épuré
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Titre centré et en gras
      axis.title = element_text(size = 14),  # Titres des axes en gras
      axis.text = element_text(size = 12, color = "black"),  # Textes des axes en gras et lisibles
      legend = NA  # Position de la légende
    )
  }

##########Correlation plot ICP CD4+
{ICP <- "TIM3"
  data_sample_plot_icp_CD4
  ###Pop Facs
  data_sample_plot_ICP_filter <-  data_sample_plot_icp_CD4 %>% filter(Day == "D0") %>% 
    # filter(!Sample_code %in% c("FL10_LN2","tFL5_PB1","FL3_PB1")) %>% 
    group_by(Sample_code) %>% 
    filter(variable == ICP) %>% 
    select(-variable,-Day) 
  
  grouped_data_cluster <-  data_sample_plot_ICP_filter %>% left_join(data_response_type_CD10) 
  
  # Calculer le modèle linéaire pour obtenir R^2
  fit_cluster <- lm(B_cell_depletion_total_glofi ~ value, data = grouped_data_cluster)
  summary_fit_cluster <- summary(fit_cluster)
  r_squared_cluster <- summary_fit_cluster$r.squared
  
  ggplot(data = grouped_data_cluster, aes(x = B_cell_depletion_total_glofi, y = value, label = Sample_code)) +
    geom_point(aes(color = Sample_code), size = 3, show.legend = FALSE) +
    geom_point(size = 3, show.legend = FALSE,shape = 21, stroke = 0.5)+
    geom_text_repel(aes(label = Sample_code))+
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
    # geom_text_repel(size = 4, box.padding = 0.3, point.padding = 0.2, max.overlaps =50) +  # Annotations ajustées
    annotate("text", x = 25, y = -10, label = sprintf("R² = %.4f", r_squared_cluster), 
             size = 5, color = "red3") + 
    labs(
      x = paste0("B Cell Depletion (%)"),
      y = paste0("% of CD4+ ", ICP, "+ by Facs" ),
      title = paste0("Correlation between % of CD4+ ",ICP,"+ cells and Glofitamab depletion"),
      color = "Sample Code"  # Légende des couleurs
    ) +
    theme_classic() +  # Thème épuré
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Titre centré et en gras
      axis.title = element_text(size = 14),  # Titres des axes en gras
      axis.text = element_text(size = 12, color = "black"),  # Textes des axes en gras et lisibles
      legend = NA  # Position de la légende
    )
}

##########Correlation plot ICP CD8+
{ICP <- "TIM3"
  data_sample_plot_icp_CD8
  ###Pop Facs
  data_sample_plot_ICP_filter <-  data_sample_plot_icp_CD8 %>% filter(Day == "D0") %>% 
    # filter(!Sample_code %in% c("FL10_LN2","tFL5_PB1","FL3_PB1")) %>% 
    group_by(Sample_code) %>% 
    filter(variable == ICP) %>% 
    select(-variable,-Day) 
  
  grouped_data_cluster <-  data_sample_plot_ICP_filter %>% left_join(data_response_type_CD10) 
  
  # Calculer le modèle linéaire pour obtenir R^2
  fit_cluster <- lm(B_cell_depletion_total_glofi ~ value, data = grouped_data_cluster)
  summary_fit_cluster <- summary(fit_cluster)
  r_squared_cluster <- summary_fit_cluster$r.squared
  
  ggplot(data = grouped_data_cluster, aes(x = B_cell_depletion_total_glofi, y = value, label = Sample_code)) +
    geom_point(aes(color = Sample_code), size = 3, show.legend = FALSE) +
    geom_point(size = 3, show.legend = FALSE,shape = 21, stroke = 0.5)+
    geom_text_repel(aes(label = Sample_code))+
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth= 0.8) +  # Ligne de tendance noire avec ombrage gris clair
    # geom_text_repel(size = 4, box.padding = 0.3, point.padding = 0.2, max.overlaps =50) +  # Annotations ajustées
    annotate("text", x = 25, y = -10, label = sprintf("R² = %.4f", r_squared_cluster), 
             size = 5, color = "red3") + 
    labs(
      x = paste0("B Cell Depletion (%)"),
      y = paste0("% of CD8+ ", ICP, "+ by Facs" ),
      title = paste0("Correlation between % of CD8+ ",ICP,"+ cells and Glofitamab depletion"),
      color = "Sample Code"  # Légende des couleurs
    ) +
    theme_classic() +  # Thème épuré
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Titre centré et en gras
      axis.title = element_text(size = 14),  # Titres des axes en gras
      axis.text = element_text(size = 12, color = "black"),  # Textes des axes en gras et lisibles
      legend = NA  # Position de la légende
    )
}

