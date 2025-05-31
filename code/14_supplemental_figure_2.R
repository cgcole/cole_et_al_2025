
#Read in timeline
timeline <- read_csv("./mouse_timelines/supplemental_figure_2_timeline.csv") %>% 
  mutate(size = if_else(is.na(important), small_circle_size,
                        if_else(important == "yes", large_circle_size, NA))) %>% 
  group_by(group) %>% 
  mutate(radius = size/2,
         distance = radius + lag(radius, default = 0),
         point_location = cumsum(distance) - radius[1])  %>% 
  ungroup() %>% 
  mutate(size = case_when(
    sample_collection == "yes" ~ size - (point_size * 2),
    TRUE ~ size))


#Subsetting the data
t_subset <- seq_table_filtered %>% 
  filter(experiment == "CCME1") %>% 
  filter(sample.type == "Fecal") %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()



#Figure S2A 

#Group 1 PBS, BpKH6, and BpSCSK
group_1_timeline <- timeline %>% 
  filter(group == "PBS,BpKH6,BpSCSK")

group_1_range <- max(group_1_timeline$day)-min(group_1_timeline$day)

group_1_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_1_timeline, sample_collection == "yes"), 
             aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_1_timeline, important == "yes"), 
            aes(label = day), size = 5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.3))
ggsave("./plots/supplemental_figure2a.pdf", 
       device = "pdf",
       width = group_1_range * .215,
       height = 2.5,
       units = "cm")






##Figure S2B
sf2b_1 <- t_subset %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, group,
           day) %>%
  summarise(rel_abundance=sum(rel_abundance/mouse.total)) %>%
  ungroup() %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>% 
  group_by(day, group) %>% 
  arrange(Genus) %>% 
  mutate(cum.pct = cumsum(rel_abundance), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label= if_else(grepl("unclassified",Genus), 
                            "unclassified",
                            as.character(Genus))) %>%
  mutate(tax.label = if_else(rel_abundance >= .1, tax.label, "")) %>%
  ggplot(aes(x=day,y=rel_abundance)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  #  geom_text(aes(y=1-y.text,x=day,label=tax.label), angle=90,
  #            lineheight=0.6, size.unit = "pt", size = 3) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ group, scales = "free",space = "free") +
  ylab("16S Relative Abundance") +
  xlab("") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0.0015,0.0015))



sf2b_2 <- t_subset %>%
  mutate(Strain = ifelse(seq %in% BpSCSK, "BpSCSK",
                         ifelse(seq %in% BpKH6, "BpKH6", NA))) %>% 
  filter(!is.na(Strain)) %>% 
  filter(Strain %in% c("BpSCSK", "BpKH6")) %>% 
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, group,
           day, mouse.number, Strain) %>%
  summarise(rel_abundance=sum(rel_abundance)) %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, group,
           day, Strain) %>%
  summarise(mean_rel_abundance = mean(rel_abundance),
            SD = sd(rel_abundance), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>%
  mutate_at(vars(SD, CI), ~ ifelse(. == 0, NA, .)) %>% 
  mutate(SD_ymin = mean_rel_abundance - SD,
         CI_ymin = mean_rel_abundance - CI) %>% 
  mutate(SD_ymin = ifelse(SD_ymin < 0, 0, SD_ymin)) %>% 
  mutate(CI_ymin = ifelse(CI_ymin < 0, 0, CI_ymin)) %>% 
  ggplot(aes(x=day,y=mean_rel_abundance)) +
  geom_bar(aes(fill=Strain), stat="identity", na.rm = TRUE) +
 # geom_errorbar(aes(ymin = SD_ymin, 
 #                   ymax = mean_rel_abundance + SD),
 #               width = 0.4, linewidth = .5) +
  scale_fill_manual(values=figure_colors) +
  facet_grid(. ~ group, scales = "free",space = "free", labeller = labeller(.multi_line = FALSE), drop = TRUE) +
  ylab("16S Relative Abundance") +
  xlab("Day") +
  guides(fill=guide_legend(title="Strain:")) +
  theme(strip.text.x = element_blank(),
        legend.spacing = unit(0, "pt")) +
  scale_y_continuous(expand = c(0.0015,0.0015), limits = c(0, 1), breaks = c(0.00, 0.25,0.50, 0.75, 1.00))


pdf(file = "./plots/supplemental_figure2b.pdf", height = 2.1, width = 2.5)
ggstack <- gg.stack(sf2b_1,sf2b_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()

rm(t_subset, sf2b_1, sf2b_2, timeline, group_1_timeline, group_1_range)
