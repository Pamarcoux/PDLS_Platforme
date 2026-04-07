##### TLS vs BCN ####
source("../GlofiResistance/00_Metadata_Talyies.R")

# Depletion  --------------------------------------------------------------

data_filtered_depletion_BCN_vs_TLS  <- data_talyies_LN_PBMC%>%
  select(Origin, Disease, Treatment, Sample_code, Location,B_cell_depletion_total, Day, Area_mean,Review) %>%
  filter(Day %in% c("D6") & Treatment == "αCD20-TCB 0,1 nM" & Disease == "FL")

Plot_depletion_BCN_vs_TLS <- ggplot(data_filtered_depletion_BCN_vs_TLS, aes(x =Location , y = B_cell_depletion_total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6, fill = "lightgrey") +  # Ajouter des boxplots
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all)+
  stat_compare_means(method = "wilcox", paired = FALSE, aes(group = Location), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  labs(
    title = "BCN vs TLS",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    # y = "B Cell Depletion",  # Légende de l'axe Y
    y = "B Cell Depletion (%)",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  guides(fill = guide_legend(ncol = 6))  # nombre de colonnes souhaité
# theme(legend.position = "none")
  
  
# Population --------------------------------------------------------------
data_filtered_pop_BCN_vs_TLS  <- data_talyies_LN_PBMC%>%
  select(Origin, Disease, Treatment, Sample_code, Location,B_cell_depletion_total, Day,Review,Per_CD4,Per_CD8,Per_CD22_total,Pop) %>%
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD4")


data_filtered_pop_BCN_vs_TLS_long <- data_filtered_pop_BCN_vs_TLS %>%
  pivot_longer(cols = c(Per_CD4, Per_CD8, Per_CD22_total), 
               names_to = "Population", 
               values_to = "Percentage") |> 
  mutate(Population = recode(Population, Per_CD4 = "CD4+ T cells", Per_CD8 = "CD8+ T cells", Per_CD22_total = "B cells"))

Plot_population_BCN_vs_TLS <- ggplot(data_filtered_pop_BCN_vs_TLS_long, aes(x =Location , y = Percentage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6, fill = "lightgrey") +  # Ajouter des boxplots
  facet_grid(Day ~ Population, switch = "both" )+
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 

  scale_fill_manual(values = sample_colors_all)+
  stat_compare_means(method = "wilcox", paired = FALSE, aes(group = Location), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  labs(
    title = "Percentage of population in PDLS at D0 and D6 between BCN vs TLS centers",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    # y = "B Cell Depletion",  # Légende de l'axe Y
    y = "% Population",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")


# ICP CD4 ---------------------------------------------------------------------
data_filtered_icp_CD4_BCN_vs_TLS  <- data_talyies_LN_PBMC %>%
  select(Disease,Location,Sample_code,Day,Pop,'41BB_Per',
         PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',
         CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')  |> 
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD4") |> 
  melt(id.vars = c("Disease","Day","Pop","Sample_code","Location")) %>% 
  mutate(variable = gsub("_Per", "", variable),
         variable =recode(variable,
                          PD1 = 'PD-1',
                          LAG3 = 'LAG-3',
                          TIM3 = 'TIM-3',
                          '41BB' = '4-1BB'))


graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("4-1BB","CD28")
  } else {
    variable_i <- c("LAG-3", "PD-1", "TIGIT","TIM-3")
  }
  
  data_sample_plot_population_i <- data_filtered_icp_CD4_BCN_vs_TLS %>% 
    filter(variable %in% variable_i) 
  
  plot_icp_cd4 <- ggplot(data_sample_plot_population_i, aes(x = Location, y = value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6, fill = "lightgrey") +  # Ajouter des boxplots
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    scale_fill_manual(values = sample_colors_all)+
    facet_grid(Day ~ variable, switch = "both") +
    stat_compare_means(method = "wilcox.test", aes(group = Location), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "% of CD4+ T-cells expressing ICP ",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "% of CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
}

# ICP CD8 ---------------------------------------------------------------------
data_filtered_icp_CD8_BCN_vs_TLS  <- data_talyies_LN_PBMC %>%
  select(Disease,Location,Sample_code,Day,Pop,'41BB_Per',
         PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',
         CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')  |> 
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD8") |> 
  melt(id.vars = c("Disease","Day","Pop","Sample_code","Location")) %>% 
  mutate(variable = gsub("_Per", "", variable),
         variable =recode(variable,
                          PD1 = 'PD-1',
                          LAG3 = 'LAG-3',
                          TIM3 = 'TIM-3',
                          '41BB' = '4-1BB'))


graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("4-1BB","CD28")
  } else {
    variable_i <- c("LAG-3", "PD-1", "TIGIT","TIM-3")
  }
  
  data_sample_plot_population_i <- data_filtered_icp_CD8_BCN_vs_TLS %>% 
    filter(variable %in% variable_i) 
  
  plot_icp_cd8 <- ggplot(data_sample_plot_population_i, aes(x = Location, y = value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6, fill = "lightgrey") +  # Ajouter des boxplots
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    scale_fill_manual(values = sample_colors_all)+
    facet_grid(Day ~ variable, switch = "both") +
    stat_compare_means(method = "wilcox.test", aes(group = Location), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "% of CD8+ T-cells expressing ICP ",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "% of CD8+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
}


# Montage -----------------------------------------------------------------

Figure_BCN_vs_TLS<- plot_grid(
  Plot_depletion_BCN_vs_TLS, labels = c("A"),
  plot_grid(Plot_population_BCN_vs_TLS, ncol = 1, labels = c("B")),  # Deuxième ligne
  plot_grid(plot_icp_cd4, ncol = 1, labels = c("C")),
  plot_grid(plot_icp_cd8, ncol = 1, labels = c("D")),
  nrow=4 )

ggsave(
  filename = "~/Postdoc/R/PDLS_Platforme/Outputs/Figures/Figure_BCN_vs_TLS.png",
  plot = Figure_BCN_vs_TLS,
  device = "png",
  width = 32,        # largeur A4 en cm
  height = 40 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
