
#Read in timeline
timeline <- read_csv("./mouse_timelines/figure_5_timeline.csv") %>% 
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
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()

#Figure 5A 
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
  scale_x_continuous(expand = expansion(mult = 0.2))
ggsave("./plots/figure5a_timeline.pdf", 
       device = "pdf",
       width = group_1_range * .19,
       height = 2.5,
       units = "cm")






##Figure 5B
f5b_1 <- t_subset %>%
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



f5b_2 <- t_subset %>%
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


pdf(file = "./plots/figure5b.pdf", height = 2.1, width = 3)
ggstack <- gg.stack(f5b_1,f5b_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()



##Figure 5C
f5c <- alpha_diversity %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
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
                width = 1.25, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Unique ASVs") +
  scale_x_continuous(limits = c(27, 48), breaks = c(28,31,33,37,44)) +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 125, 50)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/figure5c.pdf", 
       device = "pdf",
       width = 5.2,
       height = 3.1,
       units = "cm")



#Figure 5D
f5d <- alpha_diversity %>% 
  filter(samplename %in% t_subset$samplename) %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
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
                width = 1.25, linewidth = error_bar_linewidth, alpha = .75) +
  guides(color = guide_legend(title="Group"))  +
  labs(x = "Day", y = "Shannon Diversity") +
  scale_x_continuous(limits = c(27, 48), breaks = c(28,31,33,37,44)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 4.5, 1)) +
  scale_color_manual(values=figure_colors) 
ggsave("./plots/figure5d.pdf", 
       device = "pdf",
       width = 5.2,
       height = 3.1,
       units = "cm")



##Figure 5E
phy_subset <- prune_samples(phy@sam_data$samplename %in% t_subset$samplename, 
                            phy)
phy_subset <- prune_samples(phy_subset@sam_data$group %in% c("SPF", "PBS", "BpKH6", "BpSCSK"), 
                            phy_subset)
phy_subset <- transform_sample_counts(phy_subset, function(x) 1E6 * x/sum(x))


ord <- ordinate(phy_subset, "PCoA", "bray")
ord_data <- plot_ordination(phy_subset, 
                            ordinate(phy_subset, "PCoA", "bray"), 
                            color="group")$data
eigvec <- ord$values$Eigenvalues
fracvar <- eigvec[1:2] / sum(eigvec)  # Get first two axes' variance
percvar <- round(100 * fracvar, 1)  # Convert to percentage

f5e <- ord_data %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color=group)) +  
  geom_point(size=point_size) +
  scale_color_manual(values = figure_colors) +
  stat_ellipse(geom = "polygon", alpha=0.3, aes(fill=group), linewidth = .25) +
  scale_fill_manual(values = figure_colors) +
  guides(fill=guide_legend(title="Group"), color = guide_legend(title="Group"))  +
  labs(x = paste0("Axis.1 [", percvar[1], "%]"),
       y = paste0("Axis.2 [", percvar[2], "%]")) +
  ggtitle("PCoA: Bray-Curtis Dissimilarity")
ggsave("./plots/figure5e.pdf", 
       device = "pdf",
       width = 5.4,
       height = 4,
       units = "cm")

rm(f5b_1, f5b_2, f5c, f5d, f5e, t_subset, phy_subset, ggstack, 
   ord, ord_data, eigvec, fracvar, percvar, group_1_timeline, group_1_range, timeline)




#Subsetting the data
t_subset <- seq_table_filtered %>% 
  filter(experiment %in% c("CCME11")) %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day == 44) %>% 
  filter(!(group == "SPF" & day == 44)) %>% 
  filter(group %in% c("PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()


figure5_alpha_diversity <- alpha_diversity %>% 
  mutate(group = factor(group, levels = c("PBS", "BpKH6", "BpSCSK"))) %>% 
  filter(samplename %in% t_subset$samplename) 


#Figure 5C stats

lm_model <- lm(Observed ~ group, figure5_alpha_diversity)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Is normal

bartlett.test(Observed ~ group, figure5_alpha_diversity)

#Normal but unequal variance
oneway.test(Observed ~ group, data = figure5_alpha_diversity, var.equal = FALSE)
figure_5c_welch <- pairwise.t.test(figure5_alpha_diversity$Observed, figure5_alpha_diversity$group, p.adjust.method = "none", pool.sd = FALSE)


p.adjust(c(0.39520717, 0.03632503, 0.01372334), method = "BH")


#Figure 5D stats

lm_model <- lm(Shannon ~ group, figure5_alpha_diversity)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Is normal

bartlett.test(Shannon ~ group, figure5_alpha_diversity)

#Normal and equal variance
aov_model <- aov(figure5_alpha_diversity$Shannon ~ figure5_alpha_diversity$group)
summary(aov_model)
TukeyHSD(aov_model)


oneway.test(Shannon ~ group, data = figure5_alpha_diversity, var.equal = FALSE)
figure_5d_welch <- pairwise.t.test(figure5_alpha_diversity$Shannon, figure5_alpha_diversity$group, p.adjust.method = "none", pool.sd = FALSE)


p.adjust(c(0.198013789, 0.002374019, 1.599e-05), method = "BH")