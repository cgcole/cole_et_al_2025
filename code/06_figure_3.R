

#Read in timeline
timeline <- read_csv("./mouse_timelines/figure_3_timeline.csv") %>% 
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




#Figure 3A 

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
  scale_size_identity()
ggsave("./plots/figure3a_group1_timeline.pdf", 
       device = "pdf",
       width = group_1_range * point_spacing,
       height = 2.5,
       units = "cm")


#Group 2 SPF
group_2_timeline <- timeline %>% 
  filter(group == "SPF")

group_2_range <- max(group_2_timeline$day)-min(group_2_timeline$day)

group_2_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_2_timeline, sample_collection == "yes"), aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_2_timeline, important == "yes"), aes(label = day), size=5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity()
ggsave("./plots/figure3a_group2_timeline.pdf", 
       device = "pdf",
       width = group_2_range * point_spacing,
       height = 2.5,
       units = "cm")






##Figure 3B
f3b_1 <- t_subset %>%
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
  scale_y_continuous(expand = c(0.003,0.003))



f3b_2 <- t_subset %>%
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
  mutate_at(vars(mean_rel_abundance, SD, CI), ~ ifelse(. == 0, NA, .)) %>% 
  mutate(SD_ymin = mean_rel_abundance - SD,
         CI_ymin = mean_rel_abundance - CI) %>% 
  mutate(SD_ymin = ifelse(SD_ymin < 0, 0, SD_ymin)) %>% 
  mutate(CI_ymin = ifelse(CI_ymin < 0, 0, CI_ymin)) %>% 
  ggplot(aes(x=day,y=mean_rel_abundance)) +
  geom_bar(aes(fill=Strain), stat="identity", na.rm = TRUE) +
  geom_errorbar(aes(ymin = SD_ymin, 
                    ymax = mean_rel_abundance + SD),
                width = .5, linewidth = error_bar_linewidth) +
  scale_fill_manual(values=figure_colors) +
  facet_grid(. ~ group, scales = "free",space = "free", labeller = labeller(.multi_line = FALSE), drop = TRUE) +
  ylab("16S Relative Abundance") +
  xlab("Day") +
  guides(fill=guide_legend(title="Strain:")) +
  theme(strip.text.x = element_blank()) +
  scale_y_continuous(expand = c(0.003,0.003), limits = c(0, 1))


pdf(file = "./plots/figure3b.pdf", height = 2.1, width = 3.5)
ggstack <- gg.stack(f3b_1,f3b_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()



##Figure 3C
f3c <- alpha_diversity %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  filter(samplename %in% t_subset$samplename) %>% 
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
                width = 2, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Unique ASVs") +
  scale_x_continuous(limits = c(-3, 40), breaks = c(-2, 1, 3, 7, 14, 28)) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/figure3c.pdf", 
       device = "pdf",
       width = 5.2,
       height = 3.1,
       units = "cm")


#Figure 3D
f3d <- alpha_diversity %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  filter(samplename %in% t_subset$samplename) %>% 
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
                width = 2, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Shannon Diversity") +
  scale_x_continuous(limits = c(-3, 40), breaks = c(-2, 1, 3, 7, 14, 28)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 4.5, 1)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/figure3d.pdf", 
       device = "pdf",
       width = 5.2,
       height = 3.1,
       units = "cm")

##Figure 3E
phy_subset <- prune_samples(phy@sam_data$samplename %in% t_subset$samplename, 
                            phy)
phy_subset <- prune_samples(phy_subset@sam_data$day %in% c("1", "3", "7", "14", "28"), phy_subset)
phy_subset <- transform_sample_counts(phy_subset, function(x) x/sum(x))


ord <- ordinate(phy_subset, "PCoA", "bray")
ord_data <- plot_ordination(phy_subset, 
                            ordinate(phy_subset, "PCoA", "bray"), 
                            color="group")$data
eigvec <- ord$values$Eigenvalues
fracvar <- eigvec[1:2] / sum(eigvec)  # Get first two axes' variance
percvar <- round(100 * fracvar, 1)  # Convert to percentage


f3e <- ord_data %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color=group)) +  
  geom_point(size=point_size) +
  scale_color_manual(values = figure_colors) +
  stat_ellipse(geom = "polygon", alpha=0.3, aes(fill=group), linewidth = .25) +
  scale_fill_manual(values = figure_colors) +
  guides(fill=guide_legend(title="Group"), color = guide_legend(title="Group"))  +
  labs(x = paste0("Axis.1 [", percvar[1], "%]"),
       y = paste0("Axis.2 [", percvar[2], "%]")) +
  ggtitle("PCoA: Bray-Curtis Dissimilarity")
ggsave("./plots/figure3e.pdf", 
       device = "pdf",
       width = 5.4,
       height = 4,
       units = "cm")



rm(f3b_1, f3b_2, f3c, f3d, f3e, t_subset, phy_subset, ggstack,
   ord, ord_data, eigvec, fracvar, percvar, timeline, group_1_timeline,
   group_2_timeline, group_1_range, group_2_range)




#Subsetting the data
t_subset <- seq_table_filtered %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day == 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()


figure3_alpha_diversity <- alpha_diversity %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  filter(samplename %in% t_subset$samplename) 


#Figure 3C stats

lm_model <- lm(Observed ~ group, figure3_alpha_diversity)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Was not normal

levene_test(figure3_alpha_diversity, Observed ~ group, center = median)

#not normal and unequal variance
kruskal_test(figure3_alpha_diversity, Observed ~ group)
figure3_observed_dunn <- dunn_test(figure3_alpha_diversity, Observed ~ group, p.adjust.method = "none")

figure3_observed_dunn$p.adj <- p.adjust(figure3_observed_dunn$p, method = "BH")


#Figure 3D stats

lm_model <- lm(Shannon ~ group, figure3_alpha_diversity)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Is normal

bartlett.test(Shannon ~ group, figure3_alpha_diversity)

#Normal and Equal variance
aov_model <- aov(figure3_alpha_diversity$Shannon ~ figure3_alpha_diversity$group)
summary(aov_model)
TukeyHSD(aov_model)








