
#Subsetting the data
t_subset <- seq_table_filtered %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day <= 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()



#Figure S3A
sf3a <- t_subset %>% 
  select(experiment, day, group, treatment, sampleid, samplename, copy_per_g, copy_per_mg, total_16s_copies, total_16s_copy_concentration) %>% 
  distinct() %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  group_by(day, group) %>% 
  mutate(mean_copy_per_g = mean(copy_per_g)) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  ggplot(aes(x=day,y=copy_per_g, color = group)) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = .35), size = point_size) +
  geom_point(aes(y=mean_copy_per_g, x=day), shape="\U005F", size= 3, colour="black", show.legend = FALSE) +
  ylab("16S rRNA Copy Number/g of Feces")+
  xlab("Day") +
  guides(color = "none")  +
  theme(panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(margin = margin(r = 3)),
        panel.border = element_rect(color = "black", linewidth = .15),
        axis.ticks = element_line(linewidth = .15)) +
  facet_grid( . ~ group, scales = "free", space = "free") +
  scale_color_manual(values = figure_colors) +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = trans_breaks("log10", function(x) 10^x, n = 5), 
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits = c(100000000,200000000000)) +
  annotation_logticks(sides = "l", outside = TRUE, linewidth = .15, 
                      long = unit(4, "pt"),
                      mid = unit(2.5, "pt"),
                      short = unit(1.75, "pt")) +
  coord_cartesian(clip = "off")
ggsave("./plots/supplmentary_figure3a.pdf", 
       device = "pdf",
       width = 10.5,
       height = 3.5,
       unit = "cm")



rm(sf3a, t_subset)


