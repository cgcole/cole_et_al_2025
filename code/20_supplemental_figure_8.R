

#Subsetting the data

pfbbr_quant_subset <- pfbbr_quant %>% 
  filter(experiment == "CCME6") %>% 
  mutate(group = factor(group, levels = c("ARM", "CBBP"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  arrange(compound) %>% 
  mutate(compound = factor(compound, levels =  unique(compound)))


bile_acid_quant_subset <- bile_acid_quant %>% 
  filter(experiment == "CCME6") %>% 
  mutate(group = factor(group, levels = c("ARM", "CBBP"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day)))%>% 
  arrange(desc(compound)) %>% 
  mutate(compound = factor(compound, levels =  unique(compound))) 


# Figure S8A
pfbbr_quant_subset %>%
  ggplot(aes(x=day, y=concentration_mM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = .3)+
  facet_grid(compound_class+compound ~ group, scales = "free", space = "free_x") +
 # facet_grid(compound_class ~ group, scales = "free_y", space = "free") + 
  ylab("Concentration (mM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure8a.pdf", 
       height=2.75, 
       width=3.5,
       units = "in")


# Figure S8B
bile_acid_quant_subset %>%
  ggplot(aes(x=day, y=concentration_µM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = .3)+
  facet_grid(compound_class + compound ~ group, scales = "free", space = "free_x") +
  ylab("Concentration (µM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure8b.pdf", 
       height=4.75, 
       width=4.25, 
       units = "in")

rm(pfbbr_quant_subset, bile_acid_quant_subset)