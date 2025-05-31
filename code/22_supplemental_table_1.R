

#Average with standard deviation

pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day == 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  group_by(group, compound, compound_class) %>% 
  summarise(mean = format(mean(concentration_mM), digits = 3),
            SD = format(sd(concentration_mM), digits = 3), 
            N = n(), 
            .groups = "keep") %>% 
  ungroup() %>% 
  mutate(mean_sd = paste0(mean, " ± ", SD, " mM")) %>% 
  select(group, compound, compound_class, mean_sd) %>% 
  pivot_wider(names_from = group, values_from = mean_sd) %>% 
  arrange(compound_class) %>% 
  rename(compound = "Compound",
  compound_class = "Compound Class")


bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day == 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  group_by(group, compound, compound_class) %>% 
  summarise(mean = format(mean(concentration_µM), digits = 3),
            SD = format(sd(concentration_µM), digits = 3), 
            N = n(), 
            .groups = "keep") %>% 
  ungroup() %>% 
  mutate(mean_sd = paste0(mean, " ± ", SD, " µM")) %>% 
  select(group, compound, compound_class, mean_sd) %>% 
  pivot_wider(names_from = group, values_from = mean_sd) %>% 
  arrange(compound_class)  %>% 
  rename(compound = "Compound",
         compound_class = "Compound Class")


st1 <- pfbbr_quant_subset %>% 
  full_join(bile_acid_quant_subset)

write_xlsx(st1, "./tables/supplementary_table_1.xlsx")

rm(st1, pfbbr_quant_subset, bile_acid_quant_subset)

