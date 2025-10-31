

#Read in timeline
timeline <- read_csv("./mouse_timelines/supplemental_figure_4_timeline.csv") %>% 
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
  filter(experiment %in% c("CCME11")) %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day >= 28) %>%
  filter(day <= 44) %>% 
  filter(!(group == "SPF" & day == 44)) %>% 
  mutate(day = day - 30) %>%  
  filter(group %in% c("SPF", "FMT Control")) %>% 
  mutate(group = factor(group, levels = c("SPF", "FMT Control"))) %>%
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup() 



## Figure S4A

#Group 1 PBS, BpKH6, and BpSCSK
group_1_timeline <- timeline %>% 
  filter(group == "fmt_control")

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
ggsave("./plots/supplemental_figure4a.pdf", 
       device = "pdf",
       width = group_1_range * .21,
       height = 2.5,
       units = "cm")



## Figure S4B
sf4b <- t_subset %>%
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
  xlab("Day") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0.0015,0.0015))
ggsave("./plots/supplemental_figure4b.pdf", 
       device = "pdf",
       width = 3.5,
       height = 4,
       units = "cm")





## Figure S4C
sf4c <- alpha_diversity %>% 
  mutate(day = day - 30) %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  filter(group %in% c("FMT Control")) %>% 
  group_by(day, group) %>% 
  summarise(mean_observed = mean(Observed),
            SD = sd(Observed), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>% 
  ggplot(aes(y = mean_observed, x = day, color = group)) +
  geom_line(stat= "identity", linewidth = line_linewidth, alpha = .75) +
  geom_point(stat= "identity", size = line_point_size, shape = 15, alpha = .75) +
  geom_errorbar(aes(ymin = mean_observed - CI, 
                    ymax = mean_observed + CI),
                width = 1.2, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Unique ASVs") +
  scale_x_continuous(limits = c(-3, 15), breaks = c(-2,1,3,7,14)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/supplemental_figure4c.pdf", 
       device = "pdf",
       width = 5.1,
       height = 3.1,
       units = "cm")



#Figure S4D
sf4d <- alpha_diversity %>% 
  mutate(day = day - 30) %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  filter(group %in% c("FMT Control")) %>% 
  group_by(day, group) %>% 
  summarise(mean_shannon = mean(Shannon),
            SD = sd(Shannon), 
            N = n(), 
            SEM = SD/sqrt(N),
            CI = SEM * qt(0.975, N-1),
            .groups = "keep") %>% 
  ungroup() %>% 
  ggplot(aes(y = mean_shannon, x = day, color = group)) +
  geom_line(stat= "identity", linewidth = line_linewidth, alpha = .75) +
  geom_point(stat= "identity", size = line_point_size, shape = 15, alpha = .75) +
  geom_errorbar(aes(ymin = mean_shannon - CI, 
                    ymax = mean_shannon + CI),
                width = 1.2, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Shannon Diversity") +
  scale_x_continuous(limits = c(-3, 15), breaks = c(-2,1,3,7,14)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 4.5, 1)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/supplemental_figure4d.pdf", 
       device = "pdf",
       width = 5,
       height = 3.1,
       units = "cm")



##Supplemental Figure S4E
phy_subset <- prune_samples(phy@sam_data$samplename %in% t_subset$samplename, 
                            phy)
phy_subset <- prune_samples(phy_subset@sam_data$group %in% c("SPF", "FMT Control"), 
                            phy_subset)
phy_subset <- prune_samples(phy_subset@sam_data$day != 28 | phy_subset@sam_data$group != "FMT Control" , 
                            phy_subset)
phy_subset <- transform_sample_counts(phy_subset, function(x) 1E6 * x/sum(x))


ord <- ordinate(phy_subset, "PCoA", "bray")
ord_data <- plot_ordination(phy_subset, 
                            ordinate(phy_subset, "PCoA", "bray"), 
                            color="group")$data
eigvec <- ord$values$Eigenvalues
fracvar <- eigvec[1:2] / sum(eigvec)  # Get first two axes' variance
percvar <- round(100 * fracvar, 1)  # Convert to percentage


sf4e <- ord_data %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color=group)) +  
  geom_point(size=point_size) +
  scale_color_manual(values = figure_colors) +
  stat_ellipse(geom = "polygon", alpha=0.3, aes(fill=group), linewidth = .25) +
  scale_fill_manual(values = figure_colors) +
  guides(fill=guide_legend(title="Group"), color = guide_legend(title="Group"))  +
  labs(x = paste0("Axis.1 [", percvar[1], "%]"),
       y = paste0("Axis.2 [", percvar[2], "%]")) +
  ggtitle("PCoA: Bray-Curtis Dissimilarity")
ggsave("./plots/supplemental_figure4e.pdf", 
       device = "pdf",
       width = 5.2,
       height = 3.5,
       units = "cm")


rm(sf4b, sf4c, sf4d, sf4e, t_subset, phy_subset,
   ord, ord_data, eigvec, fracvar, percvar, timeline, group_1_range, group_1_timeline)





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


# Figure S4F
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
ggsave("./plots/supplemental_figure4f.pdf", 
       height=2.75, 
       width=4.75,
       units = "in")


#Supplemental Figure S4G
bile_acid_quant_subset %>%
  ggplot(aes(x=day, y=concentration_µM)) +
  geom_bar(aes(fill = compound),stat="summary", alpha = .9, fun = "mean") +
  geom_point(size = point_size)+
  facet_grid(compound_class + compound ~ group, scales = "free", space = "free_x") +
  ylab("Concentration (µM)") +
  xlab("Day") +
  theme(legend.position = "none") +
  scale_fill_manual(values = metabolite_figure_colors)
ggsave("./plots/supplemental_figure4g.pdf", height=4.75, width=6, units = "in")


rm(pfbbr_quant_subset, bile_acid_quant_subset)


