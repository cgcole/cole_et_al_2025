
#Subsetting the data

pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment == "CCME15") %>% 
  filter(day >= 28) %>% 
  filter(day <= 35) %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK"))) %>% 
  mutate(day = day - 28) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup() %>% 
  arrange(compound) %>% 
  mutate(compound = factor(compound, levels =  unique(compound)))


bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment == "CCME15") %>% 
  filter(day >= 28) %>% 
  filter(day <= 35) %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK"))) %>% 
  mutate(day = day - 28) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup() %>% 
  arrange(desc(compound)) %>% 
  mutate(compound = factor(compound, levels =  unique(compound))) 


# Figure S7A
pfbbr_quant_subset %>%
  ggplot(aes(x=day, y=concentration_mM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = .3)+
  facet_grid(compound_class + compound ~ group, scales = "free", space = "free_x",
             labeller = labeller(group = as_labeller(c(
               "SPF + NVMC" = "SPF +\nNVMC",
               "PBS" = "PBS",
               "BpSCSK" = "BpSCSK")))) +
  # facet_grid(compound_class ~ group, scales = "free_y", space = "free") + 
  ylab("Concentration (mM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure7a.pdf", 
       height=2.75, 
       width=2.75,
       units = "in")


# Figure S7B
bile_acid_quant_subset %>%
  ggplot(aes(x=day, y=concentration_µM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = .3)+
  facet_grid(compound_class + compound ~ group, scales = "free", space = "free_x",
             labeller = labeller(group = as_labeller(c(
               "SPF + NVMC" = "SPF +\nNVMC",
               "PBS" = "PBS",
               "BpSCSK" = "BpSCSK")))) +
  ylab("Concentration (µM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure7b.pdf",
       height=4.85, 
       width=3.75, 
       units = "in")


rm(pfbbr_quant_subset, bile_acid_quant_subset)
