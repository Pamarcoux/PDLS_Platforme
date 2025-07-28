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


####LNvPBMC####
# setwd("//crct-share.inserm.lan/CRCT09/Team 3D/TALYIES/Data/")
data_talyies_LN_PBMC <- data_talyies_full %>% 
  filter(Disease %in% c("FL","DLBCL","tFL") & Origin %in% c("LN","PBMC")) %>% 
  mutate(Origin = recode(Origin,
                         "LN" = "LN",
                         "PBMC" = "PB"))

# Associer chaque échantillon à une couleur
sample_colors_LN_PBMC <- setNames(c("#a00000","#1a80bb"), c("LN","PB"))



#### Info plot ###
data_sample_plot_population <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease,Area_mean) %>% 
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

#####Plot Area #### 
data_sample_plot_area <- data_sample_plot_population %>%
  filter(Day %in% c("D3", "D6")) %>%
  select(Sample_code, B_cell_depletion_total, Day, Area_mean, Disease, Origin) %>%
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin")) %>%
  mutate(variable = gsub("_Per", "", variable)) %>%
  distinct()

Plot_PDLS_area_LNvPBMC <- ggplot(data_sample_plot_area, aes(x = Day, y = value/1000,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "PDLS Area (mm²)",  # Légende de l'axe Y
    fill = ""  # Titre de la légende
  ) +
  theme_custom()+
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    strip.text.x = element_blank())

##### Plot LNvPBMC Via #####
data_sample_plot_viability<- data_talyies_LN_PBMC %>% select(Origin,Disease,Treatment,Sample_code,B_cell_depletion_total,Day,Viability) %>% 
  filter(Day %in% c("D0","D3","D6") & Treatment == "UT") %>% 
  select(Sample_code,B_cell_depletion_total,Day,Viability,Disease,Origin) %>% 
  melt(id.vars = c("B_cell_depletion_total","Day","Sample_code","Disease","Origin")) %>% 
  mutate(variable = gsub("_Per", "", variable)) %>% 
  distinct() 

comparisons <- list(c("D0", "D3"), c("D0", "D6"))
Plot_PDLS_viability_LNvPBMC <- ggplot(data_sample_plot_viability, aes(x = Day, y = value, fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "Viability (%)",  # Légende de l'axe Y
  ) +
  theme_custom()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.x = element_blank())

##### Plot LNvPBMC PDLS #####
# Filtrer les données pour D3 et D6
data_sample_PDLS_number <- data_talyies_LN_PBMC %>%
  select(Origin, Disease, Treatment, Sample_code, B_cell_depletion_total, Day, Cell_count_PDLS) %>%
  filter(Day %in% c("D3", "D6") & Treatment == "UT") %>% 
  melt(id.vars = c("B_cell_depletion_total", "Day", "Sample_code", "Disease", "Origin","Treatment")) %>%
  distinct()

Plot_PDLS_number_LNvPBMC<- ggplot(data_sample_PDLS_number, aes(x = Day, y = value,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",  # Titre du graphique
    x = "",  # Légende de l'axe X
    y = "Number of cells / PDLS",  # Légende de l'axe Y
    fill = "Origin"  # Titre de la légende
  ) +
  theme_custom()+
  theme( axis.title.x = element_blank(),
        strip.text.x = element_blank())


##### Plot LNvPBMC Pop 1 (Bcells, TCd4,TCD8,TGD)######
data_sample_plot_population <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD4,Per_CD8,Per_NK,Per_TGD,Treatment,Disease) %>% 
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
  CD4 = "CD4+ T-cells",
  CD8 = "CD8+ T-cells"
)

data_sample_plot_population1_D0_D3_D6 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3","D6")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

comparisons <- list(c("D0", "D3"),c("D3", "D6"),c("D0", "D6"))

plot_sample_plot_population1_D0_D3_D6_LNvPBMC <- ggplot(data_sample_plot_population1_D0_D3_D6, aes(x = Day, y = Proportion,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom",scales = "free") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",
    x = "",
    y = "Populations (%)",
    fill = "Origin",
    shape = "Origin",
  ) +
  theme_custom()+
  theme(legend.position = "right",
        axis.title.x = element_blank())


##### Plot LNvPBMC Pop 2 (Monocytes,Tfh,Tfr,Treg) #####
data_sample_plot_population <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,Tfr,Treg,Tfh,Per_CD11b,Treatment,Disease) %>% 
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

data_sample_plot_population2_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

plot_sample_plot_population2_D0_D3_LNvPBMC <- ggplot(data_sample_plot_population2_D0_D3, aes(x = Day, y = Proportion,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",
    x = "",
    y = "Populations (%)",
    fill = "Origin",
    shape = "Origin",
  ) +
  theme_custom()+
  theme(legend.position = "right",
        axis.title.x = element_blank())


##### Plot LNvPBMC B cells #####
data_sample_plot_population <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD22_total,Per_CD19,Per_CD20,Per_CD22_CD10_plus,Treatment,Disease) %>% 
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

data_sample_plot_population_D0_D3_D6_Bcells <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3","D6")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

plot_sample_plot_population_D0_D3_D6_Bcells_LNvPBMC <- ggplot(data_sample_plot_population_D0_D3_D6_Bcells, aes(x = Day, y = Proportion,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",
    x = "",
    y = "Percentage (%)",
    fill = "Origin",
    shape = "Origin",
  ) +
  theme_custom()+
  theme(legend.position = "right",
        axis.title.x = element_blank())

#### PLot LNvPBMC CD3 #####
			

data_sample_plot_population <- data_talyies_LN_PBMC %>% select(Origin,Sample_code,B_cell_depletion_total,Day,Per_CD3,CD3_Naive,CD3_Central_Memory,CD3_TEMRA,CD3_Effector_Memory,Treatment,Disease) %>% 
  rename(Per_CD3_Naive = CD3_Naive, Per_CD3_Central_Memory = CD3_Central_Memory, Per_CD3_Effector_Memory = CD3_Effector_Memory, Per_CD3_TEMRA = CD3_TEMRA) %>%
  filter(Treatment == "UT") %>% 
  distinct() %>%  
  pivot_longer(cols = starts_with("Per_"), 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(variable = gsub("^Per_", "", variable),
         variable = gsub("_total", "", variable))

# Définir l'ordre des variables et leurs couleurs
data_sample_plot_population$variable <- factor(data_sample_plot_population$variable, levels = c("CD3","CD3_Naive", "CD3_Central_Memory","CD3_Effector_Memory","CD3_TEMRA"))

data_sample_plot_populationCD3_D0_D3 <-  data_sample_plot_population %>% filter(Day %in% c("D0","D3")) %>% 
  group_by(Sample_code,Day) %>% 
  rename(Proportion = value) 

plot_sample_plot_populationCD3_D0_D3_LNvPBMC <- ggplot(data_sample_plot_populationCD3_D0_D3, aes(x = Day, y = Proportion,fill = Origin)) +
  geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5) +
  facet_wrap(vars(variable), ncol = 5, strip.position = "bottom") +
  stat_compare_means(method = "t.test", paired = FALSE, aes(group = Origin), 
                     hide.ns = TRUE, label = "p.signif", 
                     size = 3.5, vjust = 0.1, label.x.npc = 'middle') +  
  # scale_fill_manual(values = sample_colors_LN_PBMC)+
  labs(
    title = "",
    x = "",
    y = "Populations (%)",
    fill = "Origin",
    shape = "Origin",
  ) +
  theme_custom()+
  theme(legend.position = "right",
        axis.title.x = element_blank())

#### Montage ####
ABC <- plot_grid(Plot_PDLS_area_LNvPBMC,Plot_PDLS_viability_LNvPBMC,Plot_PDLS_number_LNvPBMC,ncol = 3,labels = c("A","B", "C"))

D <- plot_sample_plot_population1_D0_D3_D6_LNvPBMC


E <- plot_sample_plot_population2_D0_D3_LNvPBMC

F <- plot_sample_plot_population_D0_D3_D6_Bcells_LNvPBMC

G <- plot_sample_plot_populationCD3_D0_D3_LNvPBMC

Figure_2_Suppl <- plot_grid(
  plot_grid(ABC, ncol = 1, labels = c("")),  # Deuxième ligne
  plot_grid(D, ncol = 1, labels = c("D")), # 3-4 Lignes
  plot_grid(E, ncol = 1, labels = c("E")),
  plot_grid(F, ncol = 1, labels = c("F")),
  plot_grid(G, ncol = 1, labels = c("G")),
  nrow=5 )

ggsave(
  filename = "/run/user/309223354/gvfs/smb-share:server=crct-share.inserm.lan,share=crct09/Paul/Platforme_paper/Figure/Figure1/Figure_2_Suplr.png",
  plot = Figure_2_Suppl,
  device = "png",
  width = 29.7,        # largeur A4 en cm
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

