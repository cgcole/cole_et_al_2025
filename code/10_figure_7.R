
#Read in timeline
timeline <- read_csv("./mouse_timelines/figure_7_timeline.csv") %>% 
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



kleb_t_subset <- seq_table_filtered %>% 
  filter(experiment == "CCME16") %>% 
  filter(sample.type == "Fecal") %>% 
  filter(day > 28) %>% 
  mutate(group = str_replace(group, "SPF", "SPF + Amp")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + Amp", "PBS", "BpSCSK"))) %>% 
  mutate(day = day - 34) %>% 
  arrange(day) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()



kleb_cfu_t_subset <- kleb_cfu_t %>% 
  mutate(group = str_replace(group, "SPF", "SPF + Amp")) %>% 
  group_by(experiment, mouse.number, day, days.post.infection,
           sampleid, group, treatment, challenge) %>% 
  summarise(CFU.g = mean(CFU.g),
            n = n()) %>% 
  ungroup() %>% 
  mutate(CFU.g = ifelse(CFU.g <= 1000, 1000, as.numeric(CFU.g))) %>% 
  group_by(day, group) %>% 
  mutate(mean_CFU.g = mean(CFU.g)) %>% 
  ungroup() %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + Amp", "PBS", "BpSCSK"))) 



c_diff_t_subset <- seq_table_filtered %>% 
  filter(experiment == "CCME15") %>% 
  filter(sample.type == "Fecal") %>% 
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
  ungroup()

c_diff_symptom_t_subset <- c_diff_symptom_t %>% 
  filter(mouse.number != 175) %>% 
  group_by(mouse.number, group) %>% 
  summarise(weight_loss_percentage = 100 - min(weight_percentage)) %>% 
  ungroup() %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK")))


c_diff_cfu_t_subset <- c_diff_cfu_t %>% 
  filter(mouse.number != 175) %>% 
  group_by(mouse.number, experiment, day, days.post.infection,
           sampleid, group, treatment, challenge) %>% 
  summarise(CFU.g = mean(CFU.g),
            n = n()) %>% 
  ungroup() %>% 
  mutate(CFU.g = ifelse(CFU.g <= 1000, 1000, as.numeric(CFU.g))) %>% 
  group_by(day, group) %>% 
  mutate(mean_CFU.g = mean(CFU.g)) %>% 
  ungroup() %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK")))



#Figure 7A 
#Group Before Kleb
group_1_timeline <- timeline %>% 
  filter(group == "kleb_before")

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
ggsave("./plots/figure7a_1.pdf", 
       device = "pdf",
       width = group_1_range * .23,
       height = 2.5,
       units = "cm")


#Group After Kleb
group_2_timeline <- timeline %>% 
  filter(group == "kleb_after")

group_2_range <- max(group_2_timeline$day)-min(group_2_timeline$day)

group_2_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_2_timeline, sample_collection == "yes"), 
             aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_2_timeline, important == "yes"), 
            aes(label = day), size = 5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.4))
ggsave("./plots/figure7a_2.pdf", 
       device = "pdf",
       width = group_2_range * .2,
       height = 2.5,
       units = "cm")


#Group Kleb Control
group_3_timeline <- timeline %>% 
  filter(group == "kleb_control")

group_3_range <- max(group_3_timeline$day)-min(group_3_timeline$day)

group_3_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_3_timeline, sample_collection == "yes"), 
             aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_3_timeline, important == "yes"), 
            aes(label = day), size = 5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.2))
ggsave("./plots/figure7a_3.pdf", 
       device = "pdf",
       width = group_3_range * .16,
       height = 2.5,
       units = "cm")


##Figure 7B
f7b_1 <- kleb_t_subset %>%
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
  #            lineheight=0.6,size=3, size.unit = 'pt') +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ group, scales = "free",space = "free") +
  ylab("16S Relative Abundance") +
  xlab("Days P.I.") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0.003,0.003))



f7b_2 <- kleb_t_subset %>%
  mutate(Strain = ifelse(seq %in% KpMH258, "K. pneumoniae MH258", NA)) %>% 
  filter(Strain %in% c("K. pneumoniae MH258")) %>% 
  group_by(group,
           day, mouse.number, Strain) %>%
  summarise(rel_abundance=sum(rel_abundance)) %>%
  group_by( group,
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
                width = .4, linewidth = .5) +
  scale_fill_manual(values= c("#bb4d4c"),
                    labels = c(expression(italic("K. pneumoniae") ~ " MH258"))) +
  facet_grid(. ~ group, scales = "free",space = "free", labeller = labeller(.multi_line = FALSE), drop = TRUE) +
  ylab("16S Relative Abundance") +
  xlab("Day") +
  theme(strip.text.x = element_blank()) +
  guides(fill=guide_legend(title="Strain")) +
  scale_y_continuous(expand = c(0.002,0.002), limits = c(0, 1))

pdf(file = "./plots/figure7b.pdf", height = 3, width = 3.5)
ggstack <- gg.stack(f7b_1,f7b_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()



##Figure 7C
f7c <- kleb_cfu_t_subset %>% 
  ggplot(aes(x = days.post.infection, y = CFU.g, color = group)) +
  geom_point(position = position_jitterdodge(dodge.width = .7, jitter.width = .3), size = point_size) +
  geom_point(aes(x = days.post.infection, y = mean_CFU.g, color = group, fill = group),
             shape="\U005F", size=1.75, colour="black", 
             position = position_dodge(width = .7), show.legend = FALSE) +
  geom_hline(yintercept=1000,linetype="dashed", color = "black", alpha = .5, linewidth = .1) +
  guides(color = guide_legend(title="Group"))  +
  theme(legend.position = "right",
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(margin = margin(r = 3)),
        panel.border = element_rect(color = "black", linewidth = .15),
        axis.ticks = element_line(linewidth = .15),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(linewidth = .05)) +
  scale_color_manual(values=figure_colors) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x, n = 10),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(100,100000000000)) +
  annotation_logticks(sides = "l", outside = TRUE, linewidth = .15, 
                      long = unit(4, "pt"),
                      mid = unit(2.5, "pt"),
                      short = unit(1.75, "pt")) +
  coord_cartesian(clip = "off") +
  labs(y = "CFU/g feces", x = "Days P.I.") 
ggsave("./plots/figure7c.pdf", 
       device = "pdf",
       width = 5.6,
       height = 3.1,
       units = "cm")





#Figure 7D
#Group Before C. diff
group_1_timeline <- timeline %>% 
  filter(group == "cdiff_before")

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
ggsave("./plots/figure7d_1.pdf", 
       device = "pdf",
       width = group_1_range * .23,
       height = 2.5,
       units = "cm")


#Group After C. diff
group_2_timeline <- timeline %>% 
  filter(group == "cdiff_after")

group_2_range <- max(group_2_timeline$day)-min(group_2_timeline$day)

group_2_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_2_timeline, sample_collection == "yes"), 
             aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_2_timeline, important == "yes"), 
            aes(label = day), size = 5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.4))
ggsave("./plots/figure7d_2.pdf", 
       device = "pdf",
       width = group_2_range * .2,
       height = 2.5,
       units = "cm")


#Group C diff control
group_3_timeline <- timeline %>% 
  filter(group == "cdiff_control")

group_3_range <- max(group_3_timeline$day)-min(group_3_timeline$day)

group_3_timeline %>% 
  ggplot(aes(x = point_location, y = group, fill = treatment, color = treatment)) +
  geom_point(aes(size = size), shape = 21, stroke = 0) +
  geom_point(data = subset(group_3_timeline, sample_collection == "yes"), 
             aes(size = size), shape = 21, color = "black", stroke = sample_stroke) +
  geom_text(data = subset(group_3_timeline, important == "yes"), 
            aes(label = day), size = 5, size.unit = "pt", color = "black", fontface = "bold")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = timeline_colors) +
  scale_color_manual(values = timeline_colors) +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.2))
ggsave("./plots/figure7d_3.pdf", 
       device = "pdf",
       width = group_3_range * .16,
       height = 2.5,
       units = "cm")




##Figure 7E
f7e_1 <- c_diff_t_subset %>%
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
  facet_grid(. ~ group, scales = "free",space = "free",
             labeller = as_labeller(c(
               "SPF + NVMC" = "SPF +\nNVMC",
               "PBS" = "PBS",
               "BpSCSK" = "BpSCSK"))) +
  ylab("16S Relative Abundance") +
  xlab("Days P.I.") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0.003,0.003))


f7e_2 <- c_diff_t_subset %>%
  mutate(Strain = ifelse(seq %in% CdR20291, "C. difficile R20291", NA)) %>% 
  filter(Strain %in% c("C. difficile R20291")) %>% 
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
                width = 0.4, linewidth = .5) +
  scale_fill_manual(values= c("#1dc5aa"),
                    labels = c(expression(italic("C. difficile") ~ " R20291"))) +
  facet_grid(. ~ group, scales = "free",space = "free", labeller = labeller(.multi_line = FALSE), drop = TRUE) +
  ylab("16S Relative Abundance") +
  xlab("Day") +
  theme(strip.text.x = element_blank()) +
  guides(fill=guide_legend(title="Strain:")) +
  scale_y_continuous(expand = c(0.001,0.001), limits = c(0, 1))


pdf(file = "./plots/figure7e.pdf", height = 3, width = 3)
ggstack <- gg.stack(f7e_1,f7e_2,heights=c(4,2.5), newpage = F, gap = 2)
dev.off()




##Figure 7F
f7f <- c_diff_cfu_t_subset %>% 
  ggplot(aes(x = days.post.infection, y = CFU.g, color = group)) +
  geom_point(position = position_jitterdodge(dodge.width = .7, jitter.width = .3), size = point_size) +
  geom_point(aes(x = days.post.infection, y = mean_CFU.g, color = group, fill = group),
             shape="\U005F", size=1.75, colour="black", 
             position = position_dodge(width = .7), show.legend = FALSE) +
  geom_hline(yintercept=1000,linetype="dashed", color = "black", alpha = .5, linewidth = .1) +
  guides(color = guide_legend(title="Group"))  +
  theme(legend.position = "right",
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(margin = margin(r = 3)),
        panel.border = element_rect(color = "black", linewidth = .15),
        axis.ticks = element_line(linewidth = .15),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(linewidth = .05)) +
  scale_color_manual(values=figure_colors) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x, n = 10),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(100,1000000000)) +
  annotation_logticks(sides = "l", outside = TRUE, linewidth = .15, 
                      long = unit(4, "pt"),
                      mid = unit(2.5, "pt"),
                      short = unit(1.75, "pt")) +
  coord_cartesian(clip = "off") +
  labs(y = "CFU/g feces", x = "Days P.I.") 
ggsave("./plots/figure7f.pdf", 
       device = "pdf",
       width = 5.1,
       height = 3.1,
       units = "cm")



#Figure 7G
f7g <- c_diff_symptom_t_subset %>% 
  ggplot(aes(y = weight_loss_percentage, x = group, fill = group)) +
  geom_bar(stat= "summary", fun = "mean", alpha = .9) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = 1.75), size = point_size) +
  guides(fill = guide_legend(title="Group"))  +
  labs(x = "", y = "Percent Weightloss (8 days)") +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(linewidth = .05)) +
  scale_fill_manual(values=figure_colors) +
  scale_y_continuous(expand = c(0.005,0.005), limits = c(0,25)) +
scale_x_discrete(labels = c(
  "SPF + NVMC" = "SPF +\nNVMC",
  "PBS" = "PBS",
  "BpSCSK" = "BpSCSK")) 
ggsave("./plots/figure7g.pdf", 
       device = "pdf",
       width = 3.7,
       height = 3.7,
       units = "cm")


rm(f7b_1, f7b_2, f7c, f7e_1, f7e_2, f7f, f7g, ggstack,
   c_diff_cfu_t_subset, c_diff_symptom_t_subset, kleb_t_subset, c_diff_t_subset, kleb_cfu_t_subset,
   timeline, group_1_timeline,group_2_timeline, group_3_timeline, group_1_range, group_2_range, group_3_range)
