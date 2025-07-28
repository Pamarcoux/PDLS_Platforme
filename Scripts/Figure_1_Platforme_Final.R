library(writexl)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(grid)
library(png)
library(Seurat)
library(gridExtra)
library(cowplot)
library(paletteer)
library(magick)
library(RColorBrewer)
library(scales)

source("~/Postdoc/R/GlofiResistance/00_Metadata_Talyies.R")
####FL####
#### Info plot ###

data_sample_plot_population <- data_talyies_FL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD19","CD22_CD10_plus"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22_CD10_plus = "CD10",
)

##### Plot FL Area ####
data_filtered <- data_talyies_FL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Area_mean) %>%
  filter(Day %in% c("D3", "D6") & Treatment == "UT")

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Area_mean) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_plot_area <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  select(Sample_code, B_cell_depletion_total, Day, Area_mean, Disease, Origin) %>%
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin")) %>%
  mutate(variable = gsub("_Per", "", variable)) %>%
  distinct()

Plot_PDLS_area_FL <- ggplot(data_sample_plot_area, aes(x = Day, y = value/1e6)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_FL)+
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = FALSE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  labs(
    title = "FL",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    # y = "PDLS Area (µm²)",  # Légende de l'axe Y
    y = "PDLS Area (mm²)",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")

##### Plot FL Via #####
data_filtered <- data_talyies_FL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Viability) %>%
  filter(Day %in% c("D0","D3","D6") & Treatment == "UT" & !is.na(Viability))
  
# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D0_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Viability) == FALSE) %>%
  filter(n_distinct(Day) == 3) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_viability<- data_filtered %>% select(Origin,Disease,Treatment,Sample_code,B_cell_depletion_total,Day,Viability) %>% 
  filter(Sample_code %in% samples_with_D0_D3_D6) %>%
  select(Sample_code,B_cell_depletion_total,Day,Viability,Disease,Origin) %>% 
  melt(id.vars = c("B_cell_depletion_total","Day","Sample_code","Disease","Origin")) %>% 
  mutate(variable = gsub("_Per", "", variable)) %>% 
  distinct() 

comparisons <- list(c("D0", "D3"), c("D0", "D6"),c("D3", "D6"))
Plot_PDLS_viability_FL <- ggplot(data_sample_plot_viability, aes(x = Day, y = value)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_y_continuous(limits = c(0, 110),breaks=c(0,20,40,40,60,80,100)) +
  scale_fill_manual(values = sample_colors_FL)+
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = comparisons,
                     hide.ns = FALSE, label = "p.signif", 
                     size = 3.5, label.y.npc = "top", vjust = 0.5)+
  labs(
    title = "FL",  # Titre du graphique
    x = "",  # Légende de l'axe X
    # y = "PDLS Area (µm²)",  # Légende de l'axe Y
    y = "Viability (%)",  # Légende de l'axe Y
    
    fill = "Glofitamab Responders"  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")

##### Plot FL PDLS #####
# Filtrer les données pour D3 et D6
data_filtered <- data_talyies_FL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Cell_count_PDLS) %>%
  filter(Day %in% c("D3", "D6") & Treatment == "UT")

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Cell_count_PDLS) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_PDLS_number <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  select(Sample_code, B_cell_depletion_total, Day, Cell_count_PDLS, Disease, Origin) %>%
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin")) %>%
  mutate(variable = gsub("_Per", "", variable)) %>%
  distinct()

Plot_PDLS_number_FL<- ggplot(data_sample_PDLS_number, aes(x = Day, y = value)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_FL)+
  scale_y_continuous(labels = scales::comma,limits = c(0, 275000), n.breaks =4) +
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = FALSE, 
                     label = "p.signif", size = 3.5, vjust = 0,label.x.npc = 'middle') + 
  labs(
    title = "FL",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    y = "Number of cells / PDLS",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")

##### Plot FL Pop 1 (Bcells, TCd4,TCD8,TGD)######
data_filtered <- data_talyies_FL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD4,Per_CD8,Per_NK,Per_TGD,Treatment,Disease) %>% 
  filter(Treatment == "UT") 
  
# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D0_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Per_CD22_total) == FALSE) %>%
  filter(n_distinct(Day) == 3) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_plot_population <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD11b","CD3","CD4", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22 = "B-Cells",
  CD11b = "Monocytes",
  CD4 = "CD4+ \nT-Cells",
  CD8 = "CD8+ \nT-Cells"
)

data_sample_plot_population1_D0_D3_D6 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3","D6")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

comparisons <- list(c("D0", "D3"),c("D3", "D6"),c("D0", "D6"))

plot_sample_plot_population1_D0_D3_D6_FL <- ggplot(data_sample_plot_population1_D0_D3_D6, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_FL)+
  scale_y_continuous(n.breaks =4) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom", scales = "free") + 
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = comparisons,
                     hide.ns = FALSE, label = "p.signif", 
                     size = 3.5, label.y.npc = "top", vjust = 0,step.increase = 0.1) +
  labs(
    title = "FL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")

##### Plot Correlation D0/D3 FL Pop1 ####
data_sample_plot_population1_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion)
  
# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population1_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population1_D0_D3$variable)[i]
  
  population_name <- paste0(Population, "+ cells")
  data_correlation_population <- data_sample_plot_population1_D0_D3 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- cor_test$estimate  # Coefficient de corrélation
  p_value <- cor_test$p.value   # Valeur p
  
  # Créer le graphique
  p <- ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_FL) +
    labs(title = "",
         x = "% Cells D0 (Tumor)",
         y="% Cells D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none") # Ajuster la taille du texte de l'axe y
  
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_pop1_FL <- plot_grid(plotlist = graph_list, ncol = 5, align = "v", axis = "tblr")



##### Plot FL Pop 2 (Monocytes,Tfh,Tfr,Treg) #####
data_sample_plot_population <- data_talyies_FL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,Tfr,Treg,Tfh,Per_CD11b,Treatment,Disease) %>% 
  rename(Per_Tfh = Tfh, Per_Tfr = Tfr, Per_Treg = Treg) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD11b","CD3", "Tfh","Tfr","Treg", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD11b = "Monocytes",
)

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- "Monocytes"
  } else {
    variable_i <- c("Tfh", "Tfr", "Treg")
  }

data_sample_plot_population_i <- data_sample_plot_population %>% 
  filter(variable %in% variable_i) 

samples_with_D0_D3 <- data_sample_plot_population_i %>% 
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(if (i == 1) n_distinct(value) == 2 else n_distinct(value) == 6) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_population2_D0_D3 <-  data_sample_plot_population_i %>% filter(Day %in% c("D0","D3")) %>% 
  filter(Sample_code %in% samples_with_D0_D3) %>%
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

 graph <- ggplot(data_sample_plot_population2_D0_D3, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_FL)+
  # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  # Formes plus variées pour différencier les groupes
  labs(
    title = "FL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")
 # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
 if (i %% 5 != 1) {
   graph <- graph + theme(axis.title.y = element_blank())
 }
 graph_list[[i]] <- graph
 }

plot_sample_plot_population2_D0_D3_FL <- plot_grid(
  plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.25, 0.75))


##### Plot Correlation D0/D3 FL Pop2 ####
data_sample_plot_population2_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion)

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population2_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population2_D0_D3$variable)[i]
  
  population_name <- paste0(Population, "+ cells")
  data_correlation_population <- data_sample_plot_population2_D0_D3 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- cor_test$estimate  # Coefficient de corrélation
  p_value <- cor_test$p.value   # Valeur p
  
  # Créer le graphique
  p <- ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_FL) +
    labs(title = "",
         x = "% Cells D0 (Tumor)",
         y="% Cells D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none")
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i == 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph 
  }
  
plot_correlation_D0_D3_pop2_FL <- plot_grid(plotlist = graph_list,ncol = 5, align = "v", axis = "tblr")


##### Plot FL B cells #####
data_sample_plot_population <- data_talyies_FL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD19","CD22_CD10_plus"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22_CD10_plus = "CD10",
)
variable_list <- c("CD20", "CD19","CD22","CD10")
graph_list <- list()
for (i in (1:4)) {
  variable_i <- variable_list[i]
  
  data_sample_plot_population_i <- data_sample_plot_population %>% 
    filter(variable %in% variable_i) 
  
  samples_with_D0_D3 <- data_sample_plot_population_i %>% 
    group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
    filter(if (i %in% c(1,2)) n_distinct(value) == 2 else n_distinct(value) == 3) %>%
    ungroup() %>% distinct(Sample_code) %>% 
    pull(Sample_code)
  
  if (i %in% c(1,2)){
    data_sample_plot_b_cells<-  data_sample_plot_population_i %>%
      filter(Day %in% c("D0","D3")) %>%
      filter(Sample_code %in% samples_with_D0_D3) %>%
      group_by(Sample_code,Day) %>% 
      rename(Proportion = value)
    comparisons <- list(c("D0", "D3"))
  }
  else {
    data_sample_plot_b_cells<-  data_sample_plot_population_i %>%
      filter(Day %in% c("D0","D3","D6")) %>%
      filter(Sample_code %in% samples_with_D0_D3) %>%
      group_by(Sample_code,Day) %>% 
      rename(Proportion = value) 
    comparisons <- list(c("D0", "D3"), c("D0", "D6"))
  }
  
  graph <- ggplot(data_sample_plot_b_cells, aes(x = Day, y = Proportion)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
    facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
    stat_compare_means(method = "t.test", paired = TRUE, 
                       comparisons = comparisons,
                       hide.ns = FALSE, label = "p.signif", 
                       size = 3, label.y.npc = "top", vjust = 0) +
    # Formes plus variées pour différencier les groupes
    labs(
      title = "",
      x = "",
      y = "Populations (%)",
      fill = "Patients",
      shape = "Origin",
    ) +
    theme_custom() +
    theme(legend.position = "none")
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  graph_list[[i]] <- graph
}

plot_sample_plot_population_D0_D3_D6_Bcells_FL <- plot_grid(
  plotlist = graph_list, ncol = 4, align = "v", axis = "tblr")+
  plot_annotation(title = "FL", theme = theme_custom())

#####Plot CD3 Pop #####
data_sample_plot_population <- data_talyies_FL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,CD3_Naive,CD3_Central_Memory,CD3_TEMRA,CD3_Effector_Memory,Treatment,Disease) %>% 
  rename(Per_CD3_Naive = CD3_Naive, Per_CD3_Central_Memory = CD3_Central_Memory, Per_CD3_Effector_Memory = CD3_Effector_Memory, Per_CD3_TEMRA = CD3_TEMRA) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable)) %>% mutate(variable = recode(variable,
                         "CD3"                  = "CD3",
                         "CD3_Naive"            = "CD3\nNaive",
                         "CD3_Central_Memory"   = "CD3\nCentral\nMemory",
                         "CD3_Effector_Memory"  = "CD3\nEffector\nMemory",
                         "CD3_TEMRA"            = "CD3\nTEMRA"
))
# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population <- data_sample_plot_population %>%
  mutate(variable = factor(variable, levels = c(
    "CD3", "CD3\nNaive", "CD3\nCentral\nMemory", "CD3\nEffector\nMemory", "CD3\nTEMRA"
  )))

data_sample_plot_populationCD3_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  filter(!Sample_code %in% c("FL3_PB1","FL13_PB1"))

plot_sample_plot_populationCD3_D0_D3_FL <- ggplot(data_sample_plot_populationCD3_D0_D3, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_FL)+
  # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  # Formes plus variées pour différencier les groupes
  labs(
    title = "FL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")


####DLBCL_tFL####
##### Plot DLBCL Area ####
data_filtered <- data_talyies_tFL_DLBCL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Area_mean) %>%
  filter(Day %in% c("D3", "D6") & Treatment == "UT")

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Area_mean) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_plot_area <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  select(Sample_code, B_cell_depletion_total, Day, Area_mean, Disease, Origin) %>%
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin")) %>%
  mutate(variable = gsub("_Per", "", variable)) %>%
  distinct()

Plot_PDLS_area_DLBCL <- ggplot(data_sample_plot_area, aes(x = Day, y = value/1e6)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                     label = "p.signif", size = 6, vjust = 0.5,label.x.npc = 'middle') + 
  labs(
    title = "tFL/DLBCL",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    # y = "PDLS Area (µm²)",  # Légende de l'axe Y
    y = "PDLS Area (mm²)",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")


##### Plot DLBCL Via #####
data_filtered <- data_talyies_tFL_DLBCL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Viability) %>%
  filter(Day %in% c("D0","D3","D6") & Treatment == "UT" & !is.na(Viability))

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D0_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Viability) == FALSE) %>%
  filter(n_distinct(Day) == 3) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_viability<- data_filtered %>% select(Origin,Disease,Treatment,Sample_code,B_cell_depletion_total,Day,Viability) %>% 
  filter(Sample_code %in% samples_with_D0_D3_D6) %>%
  select(Sample_code,B_cell_depletion_total,Day,Viability,Disease,Origin) %>% 
  melt(id.vars = c("B_cell_depletion_total","Day","Sample_code","Disease","Origin")) %>% 
  mutate(variable = gsub("_Per", "", variable)) %>% 
  distinct() 

comparisons <- list(c("D0", "D3"), c("D0", "D6"),c("D3", "D6"))
Plot_PDLS_viability_DLBCL <- ggplot(data_sample_plot_viability, aes(x = Day, y = value)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  scale_y_continuous(limits = c(0, 110),breaks=c(0,20,40,40,60,80,100)) +
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = comparisons,
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, label.y.npc = "top", vjust = 0)+
  labs(
    title = "tFL/DLBCL",  # Titre du graphique
    x = "",  # Légende de l'axe X
    # y = "PDLS Area (µm²)",  # Légende de l'axe Y
    y = "Viability (%)",  # Légende de l'axe Y
    
    fill = "Glofitamab Responders"  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")
##### Plot DLBCL PDLS #####
# Filtrer les données pour D3 et D6
data_filtered <- data_talyies_tFL_DLBCL %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Cell_count_PDLS) %>%
  filter(Day %in% c("D3", "D6") & Treatment == "UT")

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Cell_count_PDLS) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_PDLS_number <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  select(Sample_code, B_cell_depletion_total, Day, Cell_count_PDLS, Disease, Origin) %>%
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin")) %>%
  mutate(variable = gsub("_Per", "", variable)) %>%
  distinct()

Plot_PDLS_number_DLBCL <- ggplot(data_sample_PDLS_number, aes(x = Day, y = value)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  scale_y_continuous(labels = scales::comma,limits = c(0, 275000), n.breaks =4) +
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = FALSE, 
                     label = "p.signif", size = 3.5, vjust = 0 ,label.x.npc = 'middle') + 
  labs(
    title = "tFL/DLBCL",  # Titre du graphique
    x = "  ",  # Légende de l'axe X
    y = "Number of cells / PDLS",  # Légende de l'axe Y
    
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none")

##### Plot DLBCL Pop 1 (Bcells, TCd4,TCD8,TGD)######
data_filtered <- data_talyies_tFL_DLBCL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD4,Per_CD8,Per_NK,Per_TGD,Treatment,Disease) %>% 
  filter(Treatment == "UT") 

# Identifier les Sample_code qui apparaissent à la fois pour D3 et D6
samples_with_D0_D3_D6 <- data_filtered %>%
  group_by(Sample_code) %>% filter(is.na(Per_CD22_total) == FALSE) %>%
  filter(n_distinct(Day) == 3) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

# Filtrer les données originales pour inclure uniquement ces Sample_code
data_sample_plot_population <- data_filtered %>%
  filter(Sample_code %in% samples_with_D3_D6) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD11b","CD3","CD4", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22 = "B-Cells",
  CD11b = "Monocytes",
  CD4 = "CD4+ \nT-Cells",
  CD8 = "CD8+ \nT-Cells"
)

data_sample_plot_population1_D0_D3_D6 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3","D6")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

comparisons <- list(c("D0", "D3"),c("D3", "D6"),c("D0", "D6"))

plot_sample_plot_population1_D0_D3_D6_DLBCL <- ggplot(data_sample_plot_population1_D0_D3_D6, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  scale_y_continuous(n.breaks =4) +
  # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom",scales = "free") + 
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = comparisons,
                     hide.ns = FALSE, label = "p.signif", 
                     size = 3.5, label.y.npc = "top", vjust = 0) +
  # Formes plus variées pour différencier les groupes
  labs(
    title = "tFL/DLBCL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")

##### Plot Correlation D0/D3 DLBCL Pop1 ####
data_sample_plot_population1_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion) %>% 
  filter(!is.na(D0) & !is.na(D3)) # Supprimer les lignes avec NA

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population1_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population1_D0_D3$variable)[i]
  
  population_name <- paste0(Population, "+ cells")
  data_correlation_population <- data_sample_plot_population1_D0_D3 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- cor_test$estimate  # Coefficient de corrélation
  p_value <- cor_test$p.value   # Valeur p
  
  # Créer le graphique
  p <- ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_tFL_DLBCL) +
    labs(title = "",
         x = "% Cells D0 (Tumor)",
         y="% Cells D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none")  # Ajuster la taille du texte de l'axe y
  
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_pop1_DLBCL <- plot_grid(plotlist = graph_list, ncol = 5, align = "v", axis = "tblr")




##### Plot DLBCL Pop 2 (Monocytes,Tfh,Tfr,Treg) #####
data_sample_plot_population <- data_talyies_tFL_DLBCL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,Tfr,Treg,Tfh,Per_CD11b,Treatment,Disease) %>% 
  rename(Per_Tfh = Tfh, Per_Tfr = Tfr, Per_Treg = Treg) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD11b","CD3", "Tfh","Tfr","Treg", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD11b = "Monocytes",
)

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- "Monocytes"
  } else {
    variable_i <- c("Tfh", "Tfr", "Treg")
  }
  
  data_sample_plot_population_i <- data_sample_plot_population %>% 
    filter(variable %in% variable_i) 
  
  samples_with_D0_D3 <- data_sample_plot_population_i %>% 
    group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
    filter(if (i == 1) n_distinct(value) == 2 else n_distinct(value) == 6) %>%
    ungroup() %>% distinct(Sample_code) %>% 
    pull(Sample_code)
  
  data_sample_plot_population2_D0_D3 <-  data_sample_plot_population_i %>% filter(Day %in% c("D0","D3")) %>% 
    filter(Sample_code %in% samples_with_D0_D3) %>%
    group_by(Sample_code,Day) %>% 
    rename(Proportion = value) 
  
graph <- ggplot(data_sample_plot_population2_D0_D3, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") + 
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  # Formes plus variées pour différencier les groupes
  labs(
    title = "tFL/DLBCL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")

# Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
if (i == 1) {
  graph <- graph + theme(axis.title.y = element_blank())
}
graph_list[[i]] <- graph
}

plot_sample_plot_population2_D0_D3_DLBCL <- plot_grid(
  plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.25, 0.75))


##### Plot Correlation D0/D3 tFL/DLBCL Pop2 ####
data_sample_plot_population2_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion) %>% 
  filter(!is.na(D0) & !is.na(D3)) # Supprimer les lignes avec NA

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population2_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population2_D0_D3$variable)[i]
  
  population_name <- paste0(Population, "+ cells")
  data_correlation_population <- data_sample_plot_population2_D0_D3 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- cor_test$estimate  # Coefficient de corrélation
  p_value <- cor_test$p.value   # Valeur p
  
  # Créer le graphique
  p <- ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_tFL_DLBCL) +
    labs(title = "",
         x = "% Cells D0 (Tumor)",
         y="% Cells D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none")  # Ajuster la taille du texte de l'axe y
  
  # Ajouter l'annotation
  graph <- p + annotation_custom(
    grob = textGrob(
      label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
      x = unit(0.05, "npc"), y = unit(0.95, "npc"),
      hjust = 0, vjust = 1,
      gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_pop2_DLBCL <- plot_grid(plotlist = graph_list, ncol = 5, align = "v", axis = "tblr")


#### PLot DLBCL B cells ###
data_sample_plot_population <- data_talyies_tFL_DLBCL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD19","CD22_CD10_plus"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22_CD10_plus = "CD10",
)

##### Plot DLBCL B cells#####
data_sample_plot_population <- data_talyies_tFL_DLBCL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease) %>% 
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD22","CD20", "CD19","CD22_CD10_plus"))
data_sample_plot_population$variable <- dplyr::recode(
  data_sample_plot_population$variable,
  CD22_CD10_plus = "CD10",
)
variable_list <- c("CD20", "CD19","CD22","CD10")
graph_list <- list()
for (i in (1:4)) {
  variable_i <- variable_list[i]
  
  data_sample_plot_population_i <- data_sample_plot_population %>% 
    filter(variable %in% variable_i) 
  
  samples_with_D0_D3 <- data_sample_plot_population_i %>% 
    group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
    filter(if (i %in% c(1,2)) n_distinct(value) == 2 else n_distinct(value) == 3) %>%
    ungroup() %>% distinct(Sample_code) %>% 
    pull(Sample_code)
  
  if (i %in% c(1,2)){
    data_sample_plot_b_cells<-  data_sample_plot_population_i %>%
      filter(Day %in% c("D0","D3")) %>%
      filter(Sample_code %in% samples_with_D0_D3) %>%
      group_by(Sample_code,Day) %>% 
      rename(Proportion = value)
    comparisons <- list(c("D0", "D3"))
  }
  else {
    data_sample_plot_b_cells<-  data_sample_plot_population_i %>%
      filter(Day %in% c("D0","D3","D6")) %>%
      filter(Sample_code %in% samples_with_D0_D3) %>%
      group_by(Sample_code,Day) %>% 
      rename(Proportion = value) 
    comparisons <- list(c("D0", "D3"), c("D0", "D6"))
  }
  
  graph <- ggplot(data_sample_plot_b_cells, aes(x = Day, y = Proportion)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
    stat_compare_means(method = "t.test", paired = TRUE, 
                       comparisons = comparisons,
                       hide.ns = FALSE, label = "p.signif", 
                       size = 3, label.y.npc = "top", vjust = 0) +
    labs(
      title = "",
      x = "",
      y = "Populations (%)",
      fill = "Patients",
      shape = "Origin",
    ) +
    theme_custom() +
    theme(legend.position = "none")
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  graph_list[[i]] <- graph
}

plot_sample_plot_population_D0_D3_D6_Bcells_DLBCL <- plot_grid(
  plotlist = graph_list, ncol = 4, align = "v", axis = "tblr")+
  plot_annotation(title = "tFL/DLBCL", theme = theme_custom())
#####Plot CD3 Pop #####
data_sample_plot_population <- data_talyies_tFL_DLBCL %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,CD3_Naive,CD3_Central_Memory,CD3_TEMRA,CD3_Effector_Memory,Treatment,Disease) %>% 
  rename(Per_CD3_Naive = CD3_Naive, Per_CD3_Central_Memory = CD3_Central_Memory, Per_CD3_Effector_Memory = CD3_Effector_Memory, Per_CD3_TEMRA = CD3_TEMRA) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable)) %>% mutate(variable = recode(variable,
                                                                               "CD3"                  = "CD3",
                                                                               "CD3_Naive"            = "CD3\nNaive",
                                                                               "CD3_Central_Memory"   = "CD3\nCentral\nMemory",
                                                                               "CD3_Effector_Memory"  = "CD3\nEffector\nMemory",
                                                                               "CD3_TEMRA"            = "CD3\nTEMRA"
         ))
# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population <- data_sample_plot_population %>%
  mutate(variable = factor(variable, levels = c(
    "CD3", "CD3\nNaive", "CD3\nCentral\nMemory", "CD3\nEffector\nMemory", "CD3\nTEMRA"
  )))



data_sample_plot_populationCD3_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  filter(!Sample_code %in% c("tFL3_LN1","FL13_PB1"))

plot_sample_plot_populationCD3_D0_D3_DLBCL <- ggplot(data_sample_plot_populationCD3_D0_D3, aes(x = Day, y = Proportion)) +
  geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
             shape = 21,  # Utiliser un shape qui supporte le remplissage
             color = "black",  # Utiliser 'color' pour la bordure
             position = position_dodge(0),
             size = 2.5, stroke = 0.5, alpha = 0.9) + 
  geom_line(aes(group = Sample_code), 
            position = position_dodge(0), 
            color = "grey60", 
            size = 0.2) +  
  scale_fill_manual(values = sample_colors_tFL_DLBCL)+
  # scale_shape_manual(values = c(21, 24))# Ligne plus fine et élégante
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") + 
  stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                     label = "p.signif", size = 3.5, vjust = 0.5,label.x.npc = 'middle') + 
  # Formes plus variées pour différencier les groupes
  labs(
    title = "tFL/DLBCL",
    x = "",
    y = "Populations (%)",
    fill = "Patients",
    shape = "Origin",
  ) +
  theme_custom() +
  theme(legend.position = "none")


####Correlation D0/D3 All Samples ####
data_sample_plot_population_all <- data_talyies_LN_PBMC %>% 
  filter(Treatment == "UT" & Day %in% c("D0","D3"))  %>%
  select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD4,Per_CD8,Per_NK,Per_TGD,Treatment,Disease) %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population_all$variable <- factor(data_sample_plot_population_all$variable, levels = c("CD22","CD20", "CD11b","CD3","CD4", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population_all$variable <- dplyr::recode(
  data_sample_plot_population_all$variable,
  CD22 = "B-Cells",
  CD11b = "Monocytes",
  CD4 = "CD4+ T-Cells",
  CD8 = "CD8+ T-Cells"
)
####
data_sample_plot_population1_D0_D3 <-  data_sample_plot_population_all %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion) %>% 
  mutate(D0 = ifelse((variable == "CD8+ T-Cells" & Sample_code == "tFL3_LN1"),NA,D0))

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population1_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population1_D0_D3$variable)[i]
  
  
  population_name <- paste0(Population)
  data_correlation_population <- data_sample_plot_population1_D0_D3 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- round(cor_test$estimate,2)  # Coefficient de corrélation
  p_value <- formatC(cor_test$p.value, format = "e", digits = 2)
  
  # Créer le graphique
  p <- ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_all) +
    labs(title = paste0(population_name),
         subtitle = paste0("R = ",r_value,", 𝑝 = ",p_value),
         x = "% Cells D0 (Tumor)",
         y="% Cells D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0)  # hjust = 0 → gauche ; 0.5 → centre ; 1 → droite
    )
  
  # Ajouter l'annotation
  graph <- p 
  # + annotation_custom(
  #   grob = textGrob(
  #     label = sprintf("%s\nR = %.4f\np = %.2e", population_name, r_value, p_value),
  #     x = unit(0.05, "npc"), y = unit(0.95, "npc"),
  #     hjust = 0, vjust = 1,
  #     gp = gpar(fontsize = 8, col = "black")))
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_pop1_all <- plot_grid(plotlist = graph_list, ncol = 5, align = "v", axis = "tblr")
###
data_sample_plot_population_all <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,Tfr,Treg,Tfh,Per_CD11b,Treatment,Disease) %>% 
  rename(Per_Tfh = Tfh, Per_Tfr = Tfr, Per_Treg = Treg) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population_all$variable <- factor(data_sample_plot_population_all$variable, levels = c("CD11b","CD3", "Tfh","Tfr","Treg", "CD8", "NK", "TGD","NKT"))
data_sample_plot_population_all$variable <- dplyr::recode(
  data_sample_plot_population_all$variable,
  CD11b = "Monocytes",
)

data_sample_plot_population2_D0_D3 <-  data_sample_plot_population_all %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion)

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_sample_plot_population2_D0_D3$variable))) {
  Population <- unique(data_sample_plot_population2_D0_D3$variable)[i]
  
  population_name <- paste0(Population, "")
  data_correlation_population <- data_sample_plot_population2_D0_D3 %>%
    filter(variable %in% Population) %>% 
    filter(!is.na(D0 & D3))
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_population$D0, data_correlation_population$D3)
  r_value <- round(cor_test$estimate,2)  # Coefficient de corrélation
  p_value <- formatC(cor_test$p.value, format = "e", digits = 2)
  
  # Créer le graphique
graph <-ggplot(data = data_correlation_population, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.5, alpha = 0.9) +
  scale_fill_manual(values = sample_colors_all) +
  labs(title = paste0(population_name),
       subtitle = paste0("R = ",r_value,", 𝑝 = ",p_value),
       x = "% Cells D0 (Tumor)",
       y="% Cells D3 (PDLS)")+
  theme_custom() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)  # hjust = 0 → gauche ; 0.5 → centre ; 1 → droite
  )
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (Population != c("CD3")) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph 
}

plot_correlation_D0_D3_pop2_all <- plot_grid(plotlist = graph_list,ncol = 5, align = "v", axis = "tblr")


#### Montage ####
B <- plot_grid((Plot_PDLS_area_FL ),(Plot_PDLS_area_DLBCL +
                 labs(title = "tFL/DLBCL",
                      y ="")), ncol = 2, labels = c(""))

C <- plot_grid((Plot_PDLS_viability_FL ),(Plot_PDLS_viability_DLBCL +
                                                           labs(title = "tFL/DLBCL",
                                                                y="")), ncol = 2, labels = c(""))


D <- plot_grid((Plot_PDLS_number_FL ),(Plot_PDLS_number_DLBCL +
                                                                labs(title = "tFL/DLBCL",
                                                                     y="")), ncol = 2, labels = c(""))

E <- plot_grid((plot_sample_plot_population1_D0_D3_D6_FL ),(plot_sample_plot_population1_D0_D3_D6_DLBCL +
                                                                                   labs(title = "tFL/DLBCL")), ncol = 2, labels = c(""))
F <- plot_grid((plot_sample_plot_population2_D0_D3_FL ),(plot_sample_plot_population2_D0_D3_DLBCL +
                                                             labs(title = "tFL/DLBCL")), ncol = 2, labels = c(""))
G <- plot_grid((plot_sample_plot_populationCD3_D0_D3_FL ),(plot_sample_plot_populationCD3_D0_D3_DLBCL +
                                                                labs(title = "tFL/DLBCL")), ncol = 2, labels = c(""))


Figure_1_total <- plot_grid(
  plot_grid(B,C, D, ncol = 3, labels = c("B","C", "D")),  # Deuxième ligne
  plot_grid(E, ncol = 1, labels = c("E")), # 3-4 Lignes
  plot_grid(F, ncol = 1, labels = c("F")),
  plot_grid(G, ncol = 1, labels = c("G")),
  nrow=4 )

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure1/Figure_1_total.png",
  plot = Figure_1_total,
  device = "png",
  width = 32,        # largeur A4 en cm
  height = 40 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)


# ggsave(
#   filename = "Paper_platforme/Figure/Figure_1_total.svg",
#   plot = Figure_1_total,
#   device = "svg",
#   width = 29.7,        # largeur A4 en cm
#   height = 21 ,     # hauteur A4 en cm
#   units = "cm",
#   dpi = 300
# )
####FIG 2####

Figure_2_total <- plot_grid(plot_correlation_D0_D3_pop1_all,plot_correlation_D0_D3_pop2_all, 
                            labels = c("A","", ""),
                            nrow=2 )

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure1/Figure_2_total.png",
  plot = Figure_2_total,
  device = "png",
  width = 32,        # largeur A4 en cm
  height = 18 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)



######Suppl Fig ####
ASupl <- plot_grid((plot_sample_plot_population_D0_D3_D6_Bcells_FL ),(plot_sample_plot_population_D0_D3_D6_Bcells_DLBCL +
                                                               labs(title = "tFL/DLBCL",
                                                                    y ="")), ncol = 2, labels = c(""))
  
BSupl1 <- plot_correlation_D0_D3_pop1_DLBCL 
BSupl2 <- plot_correlation_D0_D3_pop2_DLBCL                                                                         
                                                                          
Figure_1_Suplr <- plot_grid(
  plot_grid(ASupl, ncol = 1, labels = c("A")),  
  plot_grid(BSupl1, ncol = 1, labels = c("B")), 
  plot_grid(BSupl2, ncol = 1, labels = c("")),
  nrow =3,
  rel_heights = c(1, 0.8,0.8))

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure1/Figure_1_Suplr.png",
  plot = Figure_1_Suplr,
  device = "png",
  width = 29.7,        # largeur A4 en cm
  height = 25 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
