
#Subsetting the data

pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment %in% c("CCME11")) %>% 
  filter(day >= 28) %>%
  filter(day <= 44) %>% 
  filter(!(group == "SPF" & day == 44)) %>% 
  filter(group %in% c("SPF",  "FMT Control", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "FMT Control", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup() %>% 
  arrange(compound) %>% 
  mutate(compound = factor(compound, levels =  unique(compound)))


bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment %in% c("CCME11")) %>% 
  filter(day >= 28) %>%
  filter(day <= 44) %>% 
  filter(!(group == "SPF" & day == 44)) %>% 
  filter(group %in% c("SPF",  "FMT Control", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "FMT Control", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup() %>% 
  arrange(desc(compound)) %>% 
  mutate(compound = factor(compound, levels =  unique(compound))) 


# Figure S5A
pfbbr_quant_subset %>%
  ggplot(aes(x=day, y=concentration_mM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = point_size)+
  facet_grid(compound_class+compound ~ group, scales = "free", space = "free_x") +
  # facet_grid(compound_class ~ group, scales = "free_y", space = "free") + 
  ylab("Concentration (mM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure5a.pdf", 
       height=2.75, 
       width=4.75,
       units = "in")


#Supplemental Figure S5B
bile_acid_quant_subset %>%
  ggplot(aes(x=day, y=concentration_µM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = point_size)+
  facet_grid(compound_class + compound ~ group, scales = "free", space = "free_x") +
  ylab("Concentration (µM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure5b.pdf", height=4.75, width=6, units = "in")


rm(pfbbr_quant_subset, bile_acid_quant_subset)

