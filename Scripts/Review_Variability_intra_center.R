##### Rep intra Lab ####
source("../GlofiResistance/00_Metadata_Talyies.R")

# Depletion  --------------------------------------------------------------

List_sample_rep <- c("FL1_LN1_1","FL1_LN1_2","FL5_PB1_1","FL5_PB1_2","FL5_PB1_3",
                     "FL6_PB1_1","FL6_PB1_2","FL16_PB1_1","FL16_PB1_2")

data_filtered_intra_Lab <- data_talyies_LN_PBMC %>% 
  filter(Sample_code %in% List_sample_rep) %>%
  mutate(Sample_origin = case_when(
    str_detect(Sample_code, "FL1_LN1")   ~ "FL1_LN1",
    str_detect(Sample_code, "FL5_PB1")   ~ "FL5_PB1",
    str_detect(Sample_code, "FL6_PB1")   ~ "FL6_PB1",
    str_detect(Sample_code, "FL16_PB1")  ~ "FL16_PB1",
    TRUE ~ NA_character_),
  Rep_number = str_extract(Sample_code, "\\d+$")  # récupère le chiffre final
  ) |> 
  relocate(Sample_code,Sample_origin, Rep_number, .after = Patient_code) |> 
  mutate(Sample_origin = factor(Sample_origin, levels = c("FL1_LN1", "FL5_PB1", "FL6_PB1", "FL16_PB1")))


data_filtered_intra_Lab_depletion <- data_filtered_intra_Lab |> 
  filter(Day %in% c("D6") & Treatment == "αCD20-TCB 0,1 nM")

plot_depletion_intra_center <- ggplot(data_filtered_intra_Lab_depletion, aes(x =Rep_number , y = B_cell_depletion_total)) +
  geom_line(aes(group = Sample_origin),
            color = "black",
            size = 0.5, alpha = 0.7) +
  geom_point(aes(fill = Sample_origin, shape = Location),
             size = 2.5, stroke = 0.5, alpha = 0.9) +
  scale_shape_manual(values = c(21, 24)) +  # rond (21) et triangle (24)
  
  # stat_compare_means(method = "wilcox", paired = FALSE, aes(group = Location), hide.ns = TRUE,
  #                    label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') +
  labs(
    title = "Depletion repetition intra center",  # Titre du graphique
    x = "Rep number",  # Légende de l'axe X
    # y = "B Cell Depletion",  # Légende de l'axe Y
    y = "B Cell Depletion (%)",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey70"))
  )

# Pop ---------------------------------------------------------------------
data_filtered_intra_Lab_pop <- data_filtered_intra_Lab |> 
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD4")


data_filtered_intra_Lab_long <- data_filtered_intra_Lab_pop %>%
  pivot_longer(cols = c(Per_CD4, Per_CD8, Per_CD22_total), 
               names_to = "Population", 
               values_to = "Percentage") |> 
  mutate(Population = recode(Population, Per_CD4 = "CD4+ T cells", Per_CD8 = "CD8+ T cells", Per_CD22_total = "B cells"))

Plot_population_intra_center <- ggplot(data_filtered_intra_Lab_long, aes(x = Rep_number, y = Percentage)) +
  facet_grid(Day ~ Population, switch = "both") +
  geom_line(aes(group = Sample_origin),
            color = "black",
            size = 0.5, alpha = 0.7) +
  geom_point(aes(fill = Sample_origin, shape = Location),
             size = 2.5, stroke = 0.5, alpha = 0.9) +
  scale_shape_manual(values = c(21, 24)) +  # rond (21) et triangle (24)

  labs(
    title = "Percentage of population in PDLS at D0 and D3 variability intra center",
    x = "Rep number",
    y = "% Population",
    fill = "Sample Name",
    shape = "Center location"
  ) +
  
  theme_custom() +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey70"))
  )
# ICP CD4 ---------------------------------------------------------------------

data_filtered_intra_Lab_icp_cd4  <- data_filtered_intra_Lab %>%
  select(Disease,Location,Sample_code,Day,Pop,'41BB_Per',
         PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',
         CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per',"Sample_origin","Rep_number")  |> 
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD4") |> 
  melt(id.vars = c("Disease","Day","Pop","Sample_code","Location","Sample_origin","Rep_number")) %>% 
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
  
  data_sample_plot_population_i <- data_filtered_intra_Lab_icp_cd4 %>% 
    filter(variable %in% variable_i) 
  
ggplot(data_sample_plot_population_i, aes(x = Rep_number, y = value)) +
  geom_line(aes(group = Sample_origin),  # Relier les points du même échantillon
            color = "black",  # Couleur des lignes
            size = 0.5, alpha = 0.7) +
    geom_point(aes(fill = Sample_origin),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    facet_grid(Day ~ variable, switch = "both") +
    labs(
      title = "% of CD4+ T-cells expressing ICP variability intra center",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "% of CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
}

# ICP CD8 ---------------------------------------------------------------------
data_filtered_intra_Lab_icp_cd8  <- data_filtered_intra_Lab %>%
  select(Disease,Location,Sample_code,Day,Pop,'41BB_Per',
         PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',
         CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per',"Sample_origin","Rep_number")  |> 
  filter(Day %in% c("D3","D0")& Disease == "FL" & Pop == "CD3_CD8") |> 
  melt(id.vars = c("Disease","Day","Pop","Sample_code","Location","Sample_origin","Rep_number")) %>% 
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
  
  data_sample_plot_population_i <- data_filtered_intra_Lab_icp_cd8 %>% 
    filter(variable %in% variable_i) 
  
  ggplot(data_sample_plot_population_i, aes(x = Rep_number, y = value)) +
    geom_line(aes(group = Sample_origin),  # Relier les points du même échantillon
              color = "black",  # Couleur des lignes
              size = 0.5, alpha = 0.7) +
    geom_point(aes(fill = Sample_origin),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    facet_grid(Day ~ variable, switch = "both") +
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

empty_plot <- ggplot() + theme_void()

plot_depletion_intra_center
Plot_population_intra_center

plot_grid(
  plot_grid(plot_depletion_intra_center, empty_plot, rel_widths = c(0.67,1)),
  Plot_population_intra_center, labels = c("A", "B"),
  nrow = 2, align = "v"
)
