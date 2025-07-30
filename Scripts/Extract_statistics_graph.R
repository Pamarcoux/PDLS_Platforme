extract_statistics_day_columns <- function(data, x_var = "Origin", y_var = "value", scale_factor = 1) {
  
  require(dplyr)
  require(ggpubr)
  
  # 1. Calculer les moyennes et écarts-types
  stats_summary <- data %>%
    group_by(!!sym(x_var)) %>%
    summarise(
      mean = mean(!!sym(y_var)/scale_factor, na.rm = TRUE),
      sd = sd(!!sym(y_var)/scale_factor, na.rm = TRUE),
      n = n(),
      se = sd/sqrt(n),
      .groups = "drop"
    )
  
  # 2. Calculer les p-values entre jours
  p_values <- compare_means(
    formula = as.formula(paste(y_var, "~", x_var)),
    data = data,
    method = "t.test",
    paired = TRUE
  )
  
  # 3. Créer une table qui contient les jours côte à côte et la p-value
  # Obtenir les jours uniques dans l'ordre
  days <- unique(data[[x_var]])
  # Trier les jours si possible (en les convertissant d'abord en caractères)
  days <- sort(as.character(days))
  
  # Initialiser un dataframe de résultats vide
  result_table <- data.frame(matrix(ncol = length(days) + 2, nrow = 1))
  
  # Nommer les colonnes
  colnames(result_table) <- c(days, "p-value","n")
  
  # Pour chaque jour, ajouter la statistique "mean ± sd"
  for (day in days) {
    day_stats <- stats_summary %>% filter(as.character(!!sym(x_var)) == day)
    result_table[1, day] <- paste0(
      round(day_stats$mean, 3), " ± ", round(day_stats$sd, 3)
    )
  }
  
  # Si nous avons exactement 2 jours, ajouter la p-value directement
  if (length(days) == 2) {
    result_table[1, "n"] <- unique(stats_summary$n)
    # Trouver la ligne de p-value pour ces deux jours
    p_val_row <- p_values %>% 
      filter(group1 == days[1] & group2 == days[2]) %>%
      head(1)
    
    if (nrow(p_val_row) > 0) {
      result_table[1, "p-value"] <- paste0(
        round(p_val_row$p, 8), " (", p_val_row$p.signif, ")"
      )
    } else {
      # Essayer l'ordre inverse
      p_val_row <- p_values %>% 
        filter(group1 == days[2] & group2 == days[1]) %>%
        head(1)
      
      if (nrow(p_val_row) > 0) {
        result_table[1, "p-value"] <- paste0(
          round(p_val_row$p, 8), " (", p_val_row$p.signif, ")"
        )
      } else {
        result_table[1, "p-value"] <- "NA"
      }
    }
  } else {
    # Si nous avons plus de 2 jours, créer une ligne pour chaque comparaison de jours
    result_table <- result_table[0,]  # Vider la table
    
    for (i in 1:nrow(p_values)) {
      day1 <- p_values$group1[i]
      day2 <- p_values$group2[i]
      p_val <- p_values$p[i]
      signif <- p_values$p.signif[i]
      
      new_row <- data.frame(matrix(NA, ncol = length(days) + 1, nrow = 1))
      colnames(new_row) <- c(days, "p-value")
      
      # Obtenir les stats pour day1 et day2
      day1_stats <- stats_summary %>% filter(as.character(!!sym(x_var)) == day1)
      day2_stats <- stats_summary %>% filter(as.character(!!sym(x_var)) == day2)
      
      # Remplir les valeurs
      new_row[1, as.character(day1)] <- paste0(
        round(day1_stats$mean, 3), " ± ", round(day1_stats$sd, 3)
      )
      new_row[1, as.character(day2)] <- paste0(
        round(day2_stats$mean, 3), " ± ", round(day2_stats$sd, 3)
      )
      new_row[1, "p-value"] <- paste0(round(p_val, 8), " (", signif, ")")
      new_row[1, "n"] <- unique(stats_summary$n)
      
      # Ajouter à la table principale
      result_table <- rbind(result_table, new_row)
      
    }
  }
  
  return(result_table)
}
data_sample_plot_population1_D0_D3_D6_stat <- data_sample_plot_population2_D0_D3 %>% 
  filter(variable == "Monocytes")
results_table <- extract_statistics_day_columns(data_sample_plot_population1_D0_D3_D6_stat,y_var = "Proportion")

# Afficher les résultats
print(results_table)
data <- data_sample_plot_population1_D0_D3_D6_stat
y_var = "value"
