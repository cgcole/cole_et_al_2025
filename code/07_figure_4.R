
pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(concentration_mM = if_else(concentration_mM <= 0, NA, concentration_mM)) %>% 
  group_by(experiment, compound) %>% 
  mutate(min_concentration_adjustment = min(concentration_mM, na.rm = TRUE)/10) %>% 
  ungroup() %>% 
  mutate(concentration_mM = if_else(is.na(concentration_mM), min_concentration_adjustment, concentration_mM)) 


pfbbr_stats_results <- data.frame()
metabolites <- unique(pfbbr_quant_subset$compound)
days <- unique(pfbbr_quant_subset$day)

for (met in metabolites) {
  for (d in days) {
    sub_data <- pfbbr_quant_subset %>%
      filter((compound == met & day == d) | (group == "SPF" & compound == met))
    
    
    # Dunn's test with control comparisons only
    dunn_result <- dunn_test(sub_data, concentration_mM ~ group, p.adjust.method = "none")
    
    # Filter for comparisons to C
    tidy_dunn <- data.frame(
      compound = met,
      day = d,
      contrast = paste0(dunn_result$group1, " + ", dunn_result$group2),
      p_value = dunn_result$p
    )
    
    pfbbr_stats_results <- bind_rows(pfbbr_stats_results, tidy_dunn)
    
  }
}

# Adjust p-values across all tests (global FDR)
pfbbr_stats_results <- pfbbr_stats_results %>%
  filter(!(contrast %in% c("BpKH6 + BpSCSK", "PBS + BpSCSK", "PBS + BpKH6"))) %>% 
  replace_na(list(p_value = 1)) %>% 
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>% 
  mutate(group = str_replace(contrast, "SPF \\+ ", "")) %>% 
  mutate(significance = case_when(
    adj_p_value < .0001 ~ "****",
    adj_p_value < .001 ~ "***",
    adj_p_value < .01 ~ "**",
    adj_p_value < .05 ~ "*",
    adj_p_value >= .05 ~ NA,
  ))








bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(concentration_mM = if_else(concentration_mM <= 0, NA, concentration_mM)) %>% 
  group_by(experiment, compound) %>% 
  mutate(min_concentration_adjustment = min(concentration_mM, na.rm = TRUE)/10) %>% 
  ungroup() %>% 
  mutate(concentration_mM = if_else(is.na(concentration_mM), min_concentration_adjustment, concentration_mM)) 


bile_acid_stats_results <- data.frame()
metabolites <- unique(bile_acid_quant_subset$compound)
days <- unique(bile_acid_quant_subset$day)

for (met in metabolites) {
  for (d in days) {
    sub_data <- bile_acid_quant_subset %>%
      filter((compound == met & day == d) | (group == "SPF" & compound == met))
    
    
    # Dunn's test with control comparisons only
    dunn_result <- dunn_test(sub_data, concentration_mM ~ group, p.adjust.method = "none")
    
    # Filter for comparisons to C
    tidy_dunn <- data.frame(
      compound = met,
      day = d,
      contrast = paste0(dunn_result$group1, " + ", dunn_result$group2),
      p_value = dunn_result$p
    )
    
    bile_acid_stats_results <- bind_rows(bile_acid_stats_results, tidy_dunn)
    
  }
}

# Adjust p-values across all tests (global FDR)
bile_acid_stats_results <- bile_acid_stats_results %>%
  filter(!(contrast %in% c("BpKH6 + BpSCSK", "PBS + BpSCSK", "PBS + BpKH6"))) %>% 
  replace_na(list(p_value = 1)) %>% 
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>% 
  mutate(group = str_replace(contrast, "SPF \\+ ", "")) %>% 
  mutate(significance = case_when(
    adj_p_value < .0001 ~ "****",
    adj_p_value < .001 ~ "***",
    adj_p_value < .01 ~ "**",
    adj_p_value < .05 ~ "*",
    adj_p_value >= .05 ~ NA,
  ))









##Subset and read in data
pfbbr_spf_quant_subset <- pfbbr_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF")) %>% 
  group_by(group, compound, compound_class) %>% 
  summarise(mean_spf_concentration = mean(concentration_mM),
            n = n()) %>% 
  ungroup() %>% 
  select(compound, compound_class, mean_spf_concentration) 


pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(concentration_mM = if_else(concentration_mM <= 0, NA, concentration_mM)) %>% 
  group_by(experiment, compound) %>% 
  mutate(min_concentration_adjustment = min(concentration_mM, na.rm = TRUE)/10) %>% 
  ungroup() %>% 
  mutate(concentration_mM = if_else(is.na(concentration_mM), min_concentration_adjustment, concentration_mM)) %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  group_by(day, group, compound, compound_class) %>% 
  summarise(mean_concentration = mean(concentration_mM),
            n = n()) %>% 
  ungroup() %>% 
  left_join(pfbbr_spf_quant_subset) %>% 
  filter(!is.na(mean_spf_concentration)) %>% 
  mutate(fold_change = log2(mean_concentration / mean_spf_concentration)) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day)))  %>% 
  arrange(desc(compound)) %>% 
  mutate(compound = factor(compound, levels =  unique(compound)))



bile_acid_spf_quant_subset <- bile_acid_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(concentration_mM = if_else(concentration_mM <= 0, NA, concentration_mM)) %>% 
  group_by(experiment, compound) %>% 
  mutate(min_concentration_adjustment = min(concentration_mM, na.rm = TRUE)/10) %>% 
  ungroup() %>% 
  mutate(concentration_mM = if_else(is.na(concentration_mM), min_concentration_adjustment, concentration_mM)) %>% 
  filter(group %in% c("SPF")) %>% 
  group_by(group, compound, compound_class) %>% 
  summarise(mean_spf_concentration = mean(concentration_mM),
            n = n()) %>% 
  ungroup() %>% 
  select(compound, compound_class, mean_spf_concentration)  


bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(concentration_mM = if_else(concentration_mM <= 0, NA, concentration_mM)) %>% 
  group_by(experiment, compound) %>% 
  mutate(min_concentration_adjustment = min(concentration_mM, na.rm = TRUE)/10) %>% 
  ungroup() %>% 
  mutate(concentration_mM = if_else(is.na(concentration_mM), min_concentration_adjustment, concentration_mM)) %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  group_by(day, group, compound, compound_class) %>% 
  summarise(mean_concentration = mean(concentration_mM),
            n = n()) %>% 
  ungroup() %>% 
  left_join(bile_acid_spf_quant_subset) %>% 
  filter(!is.na(mean_spf_concentration)) %>% 
  mutate(fold_change = log2(mean_concentration / mean_spf_concentration)) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  arrange(desc(compound)) %>% 
  mutate(compound = factor(compound, levels =  unique(compound))) 



#Figure 4A
legend_limit <- ceiling(max(abs(min(pfbbr_quant_subset$fold_change)), abs(max(pfbbr_quant_subset$fold_change))))

f4a <- pfbbr_quant_subset %>% 
  left_join(pfbbr_stats_results %>% 
              mutate(day = factor(day, levels = unique(day)))) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
  ggplot(aes(x = day, y = compound, fill = fold_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance), size = 1, color = "black", na.rm = TRUE) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    name = expression(bold(Log[2] ~ "(Fold Change)")),
    limit = c(-10, 10),
    guide = guide_colorbar(ticks.colour = "black", 
                           frame.colour = "black",
                           frame.linewidth = .25,
                           ticks.outside = TRUE),
  ) +
  labs(x = "Day") +
  facet_grid(compound_class ~ group, scales = "free_y", space = "free") + 
  theme(legend.key.height = unit(1.5, "mm"), 
        legend.key.width = unit(1.5, "mm"),
        axis.title.y = element_blank(),
        legend.ticks.length = unit(-.65, "mm"),
        legend.spacing = unit(0, "pt"),
        legend.justification = "top")   
ggsave("./plots/figure4a.pdf", 
       height= 3.8, 
       width = 10.3,
       units = "cm")



#Figure 4B
legend_limit <- ceiling(max(abs(min(bile_acid_quant_subset$fold_change)), abs(max(bile_acid_quant_subset$fold_change))))

f4b <- bile_acid_quant_subset %>% 
  left_join(bile_acid_stats_results %>% 
              mutate(day = factor(day, levels = unique(day)))) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
  ggplot(aes(x = day, y = compound, fill = fold_change)) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance), size = 1, color = "black", na.rm = TRUE) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, 
    name = expression(bold(Log[2] ~ "(Fold Change)")),
    limit = c(-legend_limit, legend_limit),
    guide = guide_colorbar(ticks.colour = "black", 
                           frame.colour = "black",
                           frame.linewidth = .25,
                           ticks.outside = TRUE),
  ) +
  labs(x = "Day") +
  facet_grid(compound_class ~ group, scales = "free_y", space = "free") + 
  theme(legend.key.height = unit(1.95, "mm"),  
        legend.key.width = unit(1.5, "mm"),
        axis.title.y = element_blank(),
        legend.ticks.length = unit(-.65, "mm"),
        legend.spacing = unit(0, "pt"),
        legend.justification = "top") 
ggsave("./plots/figure4b.pdf", 
       height= 5.7, 
       width = 12,
       units = "cm")


rm(f4a, f4b, pfbbr_spf_quant_subset, pfbbr_quant_subset,
   bile_acid_spf_quant_subset, bile_acid_quant_subset, legend_limit)



