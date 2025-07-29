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


###FL####
##### Plot ICP FL CD4 ####
data_sample_plot_icp <- data_talyies_FL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per,'41BB_Per') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD4 %>% filter(variable == "PD1") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_CD4 <- data_sample_plot_icp_CD4 %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD4 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "FL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+ CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  #Separer les graphs icp et activateur
  if (i == 1) {
    plot_act_CD4_FL <- graph
  } else
  {
    plot_icp_CD4_FL <- graph
  }

}

##### Plot ICP FL CD8 #####
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per,'41BB_Per') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD8 %>% filter(variable == "PD1") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_CD8 <- data_sample_plot_icp_CD8 %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD8 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "FL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+ CD8+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i == 1) {
    plot_act_CD8_FL <- graph
  } else
  {
    plot_icp_CD8_FL <- graph
  }
}

# plot_icp_CD8_FL <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))
##### Plot ICP RFI FL CD4 ####
data_sample_plot_icp <- data_talyies_FL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_RFI',PDL1_RFI,PDL2_RFI,CD79b_RFI,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,'41BB_RFI',CD28_RFI)
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,CD28_RFI,'41BB_RFI') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_RFI", "", variable))

# Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD4 %>% filter(variable == "PD1") %>% 
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_icp_CD4 <- data_sample_plot_icp_CD4 %>% 
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD4 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "FL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "RFI ICP+ CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  if (i == 1) {
    plot_act_CD4_FL_RFI <- graph
  } else
  {
    plot_icp_CD4_FL_RFI <- graph
  }
}

# plot_icp_CD4_FL_RFI <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))

##### Plot ICP RFI CD8 FL #####
data_sample_plot_icp <- data_talyies_FL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_RFI',PDL1_RFI,PDL2_RFI,CD79b_RFI,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,'41BB_RFI',CD28_RFI)
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,CD28_RFI,'41BB_RFI') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_RFI", "", variable))

# Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD8 %>% filter(variable == "PD1") %>% 
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_icp_CD8 <- data_sample_plot_icp_CD8 %>% 
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD8 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "FL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "RFI ICP+ CD8+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  if (i == 1) {
    plot_act_CD8_FL_RFI <- graph
  } else
  {
    plot_icp_CD8_FL_RFI <- graph
  }
}

# plot_icp_CD8_FL_RFI <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))

##### Plot ICP FL B-Cells ####
data_sample_plot_icp <- data_talyies_FL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_Bcells<- data_sample_plot_icp %>% 
  filter(Pop == "CD22", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,CD79b_Per,PDL1_Per,PDL2_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_Bcells %>% filter(variable == "CD79b") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_Bcells <- data_sample_plot_icp_Bcells %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

variable <- unique(data_sample_plot_icp_Bcells$variable)
graph_list <- list()
for (i in (1:length(variable))) {
 variable_i <- variable[i]
  
  data_sample_plot_population_i <- data_sample_plot_icp_Bcells %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_FL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "FL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+\n CD22+ B-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  graph_list[[i]] <- graph
}

plot_icp_B_FL <- plot_grid(
  plotlist = graph_list, ncol = 3, align = "v", axis = "tblr")

###DLBCL_tFL####
##### Plot ICP DLBCL CD4####
data_sample_plot_icp <- data_talyies_tFL_DLBCL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per,'41BB_Per') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD4 %>% filter(variable == "PD1") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_CD4 <- data_sample_plot_icp_CD4 %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD4 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "tFL/DLBCL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+ CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  #Separer les graphs icp et activateur
  if (i == 1) {
    plot_act_CD4_DLBCL <- graph + 
      theme(axis.title.y = element_blank()) 
  } else
  {
    plot_icp_CD4_DLBCL <- graph
  }
}

# plot_icp_CD4_DLBCL <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))

##### Plot ICP DLBCL CD8 #####
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per,'41BB_Per') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD8 %>% filter(variable == "PD1") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_CD8 <- data_sample_plot_icp_CD8 %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD8 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "tFL/DLBCL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+ CD8+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  #Separer les graphs icp et activateur
  if (i == 1) {
    plot_act_CD8_DLBCL <- graph + 
      theme(axis.title.y = element_blank()) 
  } else
  {
    plot_icp_CD8_DLBCL <- graph
  }
}

# plot_icp_CD8_DLBCL <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))
##### Plot ICP RFI FL CD4 ####
data_sample_plot_icp <- data_talyies_tFL_DLBCL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_RFI',PDL1_RFI,PDL2_RFI,CD79b_RFI,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,'41BB_RFI',CD28_RFI)
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,CD28_RFI,'41BB_RFI') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_RFI", "", variable))

# Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD4 %>% filter(variable == "PD1") %>% 
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_icp_CD4 <- data_sample_plot_icp_CD4 %>% 
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD4 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "tFL/DLBCL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "RFI ICP+ CD4+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  #Separer les graphs icp et activateur
  if (i == 1) {
    plot_act_CD4_DLBCL_RFI <- graph +
      theme(axis.title.y = element_blank()) 
  } else
  {
    plot_icp_CD4_DLBCL_RFI <- graph
  }
}


# plot_icp_CD4_DLBCL_RFI <- plot_grid(
#   plotlist = graph_list, ncol = 2, align = "v", axis = "tblr", rel_widths = c(0.33, 0.66))

##### Plot ICP RFI CD8 FL #####
data_sample_plot_icp <- data_talyies_tFL_DLBCL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_RFI',PDL1_RFI,PDL2_RFI,CD79b_RFI,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,'41BB_RFI',CD28_RFI)
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_RFI,TIM3_RFI,PD1_RFI,TIGIT_RFI,CD28_RFI,'41BB_RFI') %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_RFI", "", variable))

# Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_CD8 %>% filter(variable == "PD1") %>% 
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>% 
  pull(Sample_code)

data_sample_plot_icp_CD8 <- data_sample_plot_icp_CD8 %>% 
  filter(Sample_code %in% samples_with_D0_D3_icp )

graph_list <- list()
for (i in (1:2)) {
  if (i == 1) {
    variable_i <- c("41BB","CD28")
  } else {
    variable_i <- c("LAG3", "PD1", "TIGIT","TIM3")
  }
  
  data_sample_plot_population_i <- data_sample_plot_icp_CD8 %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "tFL/DLBCL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "RFI ICP+ CD8+ T-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  
  #Separer les graphs icp et activateur
  if (i == 1) {
    plot_act_CD8_DLBCL_RFI <- graph + 
      theme(axis.title.y = element_blank()) 
  } else
  {
    plot_icp_CD8_DLBCL_RFI <- graph
  }
}

##### Plot ICP FL B-Cells ####
data_sample_plot_icp <- data_talyies_tFL_DLBCL %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_Bcells<- data_sample_plot_icp %>% 
  filter(Pop == "CD22", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,CD79b_Per,PDL1_Per,PDL2_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

# # Identifier les Sample_code qui apparaissent à la fois pour D0 et D3
samples_with_D0_D3_icp <- data_sample_plot_icp_Bcells %>% filter(variable == "CD79b") %>%
  group_by(Sample_code) %>% filter(is.na(value) == FALSE) %>%
  filter(n_distinct(Day) == 2) %>%
  ungroup() %>% distinct(Sample_code) %>%
  pull(Sample_code)

data_sample_plot_icp_Bcells <- data_sample_plot_icp_Bcells %>%
  filter(Sample_code %in% samples_with_D0_D3_icp )

variable <- unique(data_sample_plot_icp_Bcells$variable)
graph_list <- list()
for (i in (1:length(variable))) {
  variable_i <- variable[i]
  
  data_sample_plot_population_i <- data_sample_plot_icp_Bcells %>% 
    filter(variable %in% variable_i) 
  
  graph <- ggplot(data_sample_plot_population_i, aes(x = Day, y = value)) +
    geom_point(aes(fill = Sample_code),  # Utiliser 'fill' pour la couleur de remplissage
               shape = 21,  # Utiliser un shape qui supporte le remplissage
               color = "black",  # Utiliser 'color' pour la bordure
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) + 
    geom_line(aes(group = Sample_code), 
              position = position_dodge(0), 
              color = "grey60", 
              size = 0.2) +  
    scale_fill_manual(values = sample_colors_tFL_DLBCL)+
    facet_wrap(vars(variable), ncol = 6, strip.position = "bottom") + 
    stat_compare_means(method = "t.test", paired = TRUE, aes(group = Day), hide.ns = TRUE, 
                       label = "p.signif", size = 3, vjust = 1,label.x.npc = 'middle') + 
    labs(
      title = "tFL/DLBCL",  # Titre du graphique
      x = "",  # Légende de l'axe X
      y = "%ICP+\n CD22+ B-cells",  # Légende de l'axe Y
      fill = ""  # Titre de la légende
    ) +
    theme_custom()+
    theme(legend.position = "none")
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i %% 5 != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  graph_list[[i]] <- graph
}

plot_icp_B_DLBCL <- plot_grid(
  plotlist = graph_list, ncol = 3, align = "v", axis = "tblr")

###Correlation ICP CD4 FL + DLBCL ####
sample_colors_all <- c(sample_colors_FL,sample_colors_tFL_DLBCL)

data_sample_plot_icp <- data_talyies_LN_PBMC %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',CD28_Per,LAG3_Per,PD1_Per,TIGIT_Per,TIM3_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

data_correlation_icp_CD4 <- data_sample_plot_icp_CD4 %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion)

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_correlation_icp_CD4$variable))) {
  Population <- unique(data_correlation_icp_CD4$variable)[i]
  
  population_name <- paste0(Population, "")
  data_correlation_graph <- data_correlation_icp_CD4 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_graph$D0, data_correlation_graph$D3)
  r_value <- round(cor_test$estimate,2)  # Coefficient de corrélation
  p_value <- formatC(cor_test$p.value, format = "e", digits = 2)
  
  # Créer le graphique
  graph <- ggplot(data = data_correlation_graph, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_all) +
    labs(title = paste0(population_name),
         subtitle = paste0("R = ",r_value,"\n𝑝 = ",p_value),
         x = "%ICP+ \n CD4+ T-cells  \n D0 (Tumor)",
         y="%ICP+ CD4+ T-cells \n D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0)) # Ajuster la taille du texte de l'axe y
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_ICP_CD4 <- plot_grid(plotlist = graph_list, ncol = 6, align = "v", axis = "tblr")+ labs(title = "FL")

###Correlation ICP CD8 FL + DLBCL ####
sample_colors_all <- c(sample_colors_FL,sample_colors_tFL_DLBCL)

data_sample_plot_icp <- data_talyies_LN_PBMC %>% select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BB_Per',CD28_Per,LAG3_Per,PD1_Per,TIGIT_Per,TIM3_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code")) %>% 
  mutate(variable = gsub("_Per", "", variable))

data_correlation_icp_CD8 <- data_sample_plot_icp_CD8 %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) %>% 
  pivot_wider(names_from = Day, values_from = Proportion)

# Initialiser une liste pour stocker les graphiques
graph_list <- list()

# Boucle pour générer les graphiques
for (i in seq_along(unique(data_correlation_icp_CD8$variable))) {
  Population <- unique(data_correlation_icp_CD8$variable)[i]
  
  population_name <- paste0(Population, "")
  data_correlation_graph <- data_correlation_icp_CD8 %>%
    filter(variable %in% Population)
  
  # Calculer le modèle linéaire pour obtenir R et p
  cor_test <- cor.test(data_correlation_graph$D0, data_correlation_graph$D3)
  r_value <- round(cor_test$estimate,2)  # Coefficient de corrélation
  p_value <- formatC(cor_test$p.value, format = "e", digits = 2)
  
  # Créer le graphique
  graph <- ggplot(data = data_correlation_graph, aes(x = D0, y = D3, label = Sample_code)) +
    geom_smooth(method = "lm", color = "black", fill = "gray80", linewidth = 0.8) +
    geom_point(aes(fill = Sample_code),
               shape = 21,
               color = "black",
               position = position_dodge(0),
               size = 2.5, stroke = 0.8, alpha = 0.9) +
    scale_fill_manual(values = sample_colors_all) +
    labs(title = paste0(population_name),
         subtitle = paste0("R = ",r_value,"\n𝑝 = ",p_value),
         x = "%ICP+ \n CD8+ T-cells  \n D0 (Tumor)",
         y="%ICP+ CD8+ T-cells \n D3 (PDLS)")+
    theme_custom() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0)) 
  
  # Supprimer l'axe y pour les graphiques qui ne sont pas les premiers de la ligne
  if (i != 1) {
    graph <- graph + theme(axis.title.y = element_blank())
  }
  
  # Ajouter le graphique à la liste
  graph_list[[i]] <- graph
}

plot_correlation_D0_D3_ICP_CD8 <- plot_grid(plotlist = graph_list, ncol = 6, align = "v", axis = "tblr")+ labs(title = "FL")
###Montage 2 with Suppl ####
AB <- plot_grid((plot_act_CD4_FL),plot_act_CD4_DLBCL,plot_act_CD8_FL,plot_act_CD8_DLBCL, ncol = 4, labels = c("A","","B",""))
AB_suppl <- plot_grid((plot_act_CD4_FL_RFI),plot_act_CD4_DLBCL_RFI,plot_act_CD8_FL_RFI,plot_act_CD8_DLBCL_RFI, ncol = 4, labels = c("A","","B",""))

C <-  plot_grid(plot_icp_CD4_FL,plot_icp_CD4_DLBCL, ncol = 2, labels = c("C",""))
C_suppl <-  plot_grid(plot_icp_CD4_FL_RFI,plot_icp_CD4_DLBCL_RFI, ncol = 2, labels = c("C",""))

D <-  plot_grid(plot_icp_CD8_FL,plot_icp_CD8_DLBCL, ncol = 2, labels = c("D",""))
D_suppl <-  plot_grid(plot_icp_CD8_FL_RFI,plot_icp_CD8_DLBCL_RFI, ncol = 2, labels = c("D",""))

E <- plot_correlation_D0_D3_ICP_CD4
F <- plot_correlation_D0_D3_ICP_CD8

Figure_3_total <- plot_grid(
  plot_grid(AB, ncol = 1, labels = c("")),  # Deuxième ligne
  plot_grid(C, ncol = 1, labels = c("")),
  plot_grid(D, ncol = 1, labels = c("")),
  plot_grid(E, ncol = 1, labels = c("E")),
  plot_grid(F, ncol = 1, labels = c("F")),

  nrow=5 )

print(Figure_3_total)

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure3/Figure_3_total.png",
  plot = Figure_3_total,
  device = "png",
  width = 32,        # largeur A4 en cm
  height = 50 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)

#### Montage Old####
# A <- plot_grid((plot_act_CD4_FL),plot_act_CD4_DLBCL,plot_act_CD4_FL_RFI,plot_act_CD4_DLBCL_RFI, ncol = 4, labels = c("A","B","",""))
# 
# B <- plot_grid((plot_icp_CD8_FL),(plot_icp_CD8_DLBCL), ncol = 2, labels = c("C","D"))
# 
# C <- plot_grid((plot_icp_CD4_FL_RFI),(plot_icp_CD4_DLBCL_RFI +labs(y ="")), ncol = 2, labels = c("E","F"))
# 
# D <- plot_grid((plot_icp_CD8_FL_RFI),(plot_icp_CD8_DLBCL_RFI +labs(y="")), ncol = 2, labels = c("G","H"))
# 
# E <- plot_correlation_D0_D3_ICP_CD4
# F <- plot_correlation_D0_D3_ICP_CD8
# 
# Figure_3_total <- plot_grid(
#   plot_grid(A, ncol = 1, labels = c("")),  # Deuxième ligne
#   plot_grid(B, ncol = 1, labels = c("")),
#   plot_grid(C, ncol = 1, labels = c("")),
#   plot_grid(D, ncol = 1, labels = c("")),
#   plot_grid(E, ncol = 1, labels = c("I")),
#   plot_grid(F, ncol = 1, labels = c("J")),
#   
#   nrow=6 )
# 
# print(Figure_3_total)

###Suppl Figure LN vs PBMC #####
data_talyies_LN_PBMC <- data_talyies_full %>% 
  filter(Disease %in% c("FL","DLBCL","tFL") & Origin %in% c("LN","PBMC")) %>% 
  mutate(Origin = recode(Origin,
                         "LN" = "LN",
                         "PBMC" = "PB"))

# Associer chaque échantillon à une couleur
sample_colors_LN_PBMC <- setNames(c("#a00000","#1a80bb"), c("LN","PB"))

##### Plot ICP CD4 LNvsPBMC ####
data_sample_plot_icp <- data_talyies_LN_PBMC %>% select(Origin,Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD4<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD4", Day %in% c("D0","D3")) %>% 
  select(Origin,Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code","Origin")) %>% 
  mutate(variable = gsub("_Per", "", variable))

plot_icp_CD4_PBvsLN <- ggplot(data_sample_plot_icp_CD4, aes(x = Day, y = value,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "%ICP+ CD4+ T-cells",  # Légende de l'axe Y
    fill = "Origin"  # Titre de la légende
  ) +
  theme_custom()

##### Plot ICP CD8 PBvsLN #####
data_sample_plot_icp <- data_talyies_LN_PBMC %>% select(Origin,Disease,Sample_code,B_cell_depletion_total,Day,Pop,'41BBL_Per',PDL1_Per,PDL2_Per,CD79b_Per,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,'41BB_Per',CD28_Per,'PD1/LAG3_Per','PD1/TIM3_Per')
data_sample_plot_icp_CD8<- data_sample_plot_icp %>% 
  filter(Pop == "CD3_CD8", Day %in% c("D0","D3")) %>% 
  select(Origin,Disease,Sample_code,B_cell_depletion_total,Day,Pop,LAG3_Per,TIM3_Per,PD1_Per,TIGIT_Per,CD28_Per) %>% 
  melt(id.vars = c("Disease","B_cell_depletion_total","Day","Pop","Sample_code","Origin")) %>% 
  mutate(variable = gsub("_Per", "", variable))

plot_icp_CD8_PBvsLN <- ggplot(data_sample_plot_icp_CD8, aes(x = Day, y = value,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "%ICP+ CD8+ T-cells",  # Légende de l'axe Y
    fill = "Origin"  # Titre de la légende
  ) +
  theme_custom()

#####Montage #####


E_suppl <- plot_grid((plot_icp_CD4_PBvsLN), ncol = 1, labels = c(""))
F_suppl <- plot_grid((plot_icp_CD8_PBvsLN), ncol = 1, labels = c(""))

# E <- plot_grid((plot_icp_B_FL),plot_icp_B_DLBCL, ncol = 2, labels = c(""))
Figure_3_Suppl_RFI<- plot_grid(
  plot_grid(AB_suppl, ncol = 1, labels = c("")),  # Deuxième ligne
  plot_grid(C_suppl, ncol = 1, labels = c("")),
  plot_grid(D_suppl, ncol = 1, labels = c("")),
  nrow=3 )



Figure_3_Suppl<- plot_grid(
  plot_grid(E_suppl, ncol = 1, labels = c("A")),
  plot_grid(F_suppl, ncol = 1, labels = c("B")),
  
  nrow=2 )

print(Figure_3_Suppl)

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure3/Figure_3_Suppl.png",
  plot = Figure_3_Suppl,
  device = "png",
  width = 29,        # largeur A4 en cm
  height = 20 ,     # hauteur A4 en cm
  units = "cm",
  dpi = 300
)
