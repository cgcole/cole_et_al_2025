
#Read in timeline
timeline <- read_csv("./mouse_timelines/supplemental_figure_7_timeline.csv") %>% 
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
  filter(experiment == "CCME6") %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day >= 0) %>% 
  mutate(group = factor(group, levels = c("ARM", "CBBP"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()


#Figure S7A 
#Group CBBP
group_1_timeline <- timeline %>% 
  filter(group == "CBBP")

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
  scale_x_continuous(expand = expansion(mult = 0.2))
ggsave("./plots/supplemental_figure7a.pdf", 
       device = "pdf",
       width = group_1_range * .2,
       height = 2.5,
       units = "cm")




##Figure S7B
sf7b_1 <- t_subset %>%
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
#            lineheight=0.6,size=3, size.unit = "pt") +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ group, scales = "free",space = "free") +
  ylab("16S Relative Abundance") +
  xlab("") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0.001,0.001))



sf7b_2 <- t_subset %>%
  mutate(Strain = ifelse(seq %in% BpSCSK, "B. pseudococcoides SCSK", 
                         ifelse(seq %in% PdCBBP, "P. distasonis CBBP", 
                                ifelse(seq %in% CbCBBP, "C. bolteae CBBP",
                                       ifelse(seq %in% BsCBBP, "B. sartorii CBBP", NA))))) %>% 
  mutate(Strain = factor(Strain, levels = c("C. bolteae CBBP", "B. pseudococcoides SCSK", 
                                            "B. sartorii CBBP", "P. distasonis CBBP"))) %>% 
  filter(!is.na(Strain)) %>% 
  group_by(group, day, mouse.number, Strain) %>%
  summarise(rel_abundance=sum(rel_abundance)) %>%
  group_by(group, day, Strain) %>%
  summarise(mean_rel_abundance = mean(rel_abundance),
            SD = sd(rel_abundance), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>%
  mutate_at(vars(mean_rel_abundance, SD, CI), ~ ifelse(. == 0, NA, .)) %>% 
  mutate(SD_ymin = mean_rel_abundance - SD,
         CI_ymin = mean_rel_abundance - CI) %>% 
  mutate(SD_ymin = ifelse(SD_ymin < 0, 0, SD_ymin)) %>% 
  mutate(CI_ymin = ifelse(CI_ymin < 0, 0, CI_ymin)) %>% 
  ggplot(aes(x=day,y=mean_rel_abundance)) +
  geom_bar(aes(fill=Strain), stat="identity", na.rm = TRUE) +
  # geom_errorbar(aes(ymin = SD_ymin, 
  #                 ymax = mean_rel_abundance + SD,),
  #               width = .4, linewidth = .5) +
  scale_fill_manual(values= figure_colors,
                    labels = c(expression(italic("E. bolteae") ~ " CBBP"), 
                               expression(italic("B. pseudococcoides") ~ " SCSK"), 
                               expression(italic("P. sartorii") ~ " CBBP"),
                               expression(italic("P. distasonis") ~ " CBBP"))) +
  facet_grid(. ~ group, scales = "free",space = "free", labeller = labeller(.multi_line = FALSE), drop = TRUE) +
  ylab("16S Relative Abundance") +
  xlab("Day") +
  theme(strip.text.x = element_blank()) +
  guides(fill=guide_legend(title="Strain:")) +
  scale_y_continuous(expand = c(0.001,0.001), limits = c(0, 1))


pdf(file = "./plots/supplemental_figure7b.pdf", height = 2.1, width = 2.7559)
ggstack <- gg.stack(sf7b_1,sf7b_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()




#Figure S7C
sf7c <- alpha_diversity %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  mutate(group = factor(group, levels = c("ARM", "CBBP"))) %>% 
  group_by(day, group) %>% 
  summarise(mean_observed = mean(Observed),
            SD = sd(Observed), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>% 
  mutate(SD_ymin = mean_observed - SD,
         CI_ymin = mean_observed - CI) %>% 
  mutate(SD_ymin = ifelse(SD_ymin < 0, 0, SD_ymin)) %>% 
  mutate(CI_ymin = ifelse(CI_ymin < 0, 0, CI_ymin)) %>% 
  ggplot(aes(y = mean_observed, x = day, color = group)) +
  geom_line(stat= "identity", linewidth = line_linewidth, alpha = .75) +
  geom_point(stat= "identity", size = line_point_size, shape = 15, alpha = .75) +
  geom_errorbar(aes(ymin = CI_ymin, 
                    ymax = mean_observed + CI),
                width = 1, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Unique ASVs") +
  scale_x_continuous(limits = c(-1, 15), breaks = c(0,6,10,14)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 350, 50)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/supplemental_figure7c.pdf", 
       device = "pdf",
       width = 4.8,
       height = 3.1,
       units = "cm")


#Figure S7D
sf7d <- alpha_diversity %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  mutate(group = factor(group, levels = c("ARM", "CBBP"))) %>% 
  group_by(day, group) %>% 
  summarise(mean_shannon = mean(Shannon),
            SD = sd(Shannon), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>% 
  mutate(SD_ymin = mean_shannon - SD,
         CI_ymin = mean_shannon - CI) %>% 
  mutate(SD_ymin = ifelse(SD_ymin < 0, 0, SD_ymin)) %>% 
  mutate(CI_ymin = ifelse(CI_ymin < 0, 0, CI_ymin)) %>% 
  ggplot(aes(y = mean_shannon, x = day, color = group)) +
  geom_line(stat= "identity", linewidth = line_linewidth, alpha = .75) +
  geom_point(stat= "identity", size = line_point_size, shape = 15, alpha = .75) +
  geom_errorbar(aes(ymin = CI_ymin, 
                    ymax = mean_shannon + CI),
                width = 1, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Shannon Diversity") +
  scale_x_continuous(limits = c(-1, 15), breaks = c(0,6,10,14)) +
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4, 1)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/supplemental_figure7d.pdf", 
       device = "pdf",
       width = 4.7,
       height = 3.1,
       units = "cm")



rm(sf7b_1, sf7b_2, sf7c, sf7d, t_subset, ggstack,
   timeline, group_1_timeline, group_1_range)
