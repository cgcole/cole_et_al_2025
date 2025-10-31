

#Tax Hits for the different lan genes
gene_calls %>% 
  filter(include == "yes") %>% 
  group_by(highest_tax_hit) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  arrange(desc(total)) %>% 
  slice(1:100) %>% 
  mutate(highest_tax_hit = factor(highest_tax_hit, levels = highest_tax_hit)) %>% 
  ggplot(aes(x = highest_tax_hit, y = total)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_bw() +
  ylab("# of lanA Genes Detected") +
  xlab("Taxonomy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"), face = "bold", size = 15),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"), face = "bold", size = 15),
        axis.text = element_text(size = 10),
        plot.margin = margin(.5,.5,.5,3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
ggsave("../plots/analysis/tax_hits.pdf", 
       device = "pdf",
       width = 60,
       height = 20,
       units = "cm")


#Tax Hits for the different lan genes top 20
gene_calls %>% 
  filter(include == "yes") %>% 
  group_by(highest_tax_hit) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  arrange(desc(total)) %>% 
  slice(1:20) %>% 
  mutate(highest_tax_hit = factor(highest_tax_hit, levels = highest_tax_hit)) %>% 
  ggplot(aes(x = highest_tax_hit, y = total)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_bw() +
  ylab("# of lanA Genes Detected") +
  xlab("Taxonomy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"), face = "bold", size = 15),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"), face = "bold", size = 15),
        axis.text = element_text(size = 10),
        plot.margin = margin(.5,.5,.5,3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
ggsave("../plots/analysis/tax_hits_top_20.pdf", 
       device = "pdf",
       width = 25,
       height = 20,
       units = "cm")

  
#Number of times that a lan is detected top 100
gene_calls %>% 
  filter(include == "yes") %>% 
  group_by(lan_gene) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  arrange(desc(total)) %>% 
  slice(1:100) %>% 
  mutate(lan_gene = factor(lan_gene, levels = lan_gene)) %>% 
  ggplot(aes(x = lan_gene, y = total)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_bw() +
  ylab("# of Samples") +
  xlab("Lan Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"), face = "bold", size = 15),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"), face = "bold", size = 15),
        axis.text = element_text(size = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
ggsave("../plots/analysis/lan_hits.pdf", 
       device = "pdf",
       width = 60,
       height = 20,
       units = "cm")


# Number of times that a lan is detected top 50
gene_calls %>% 
  filter(include == "yes") %>% 
  filter(lan_gene != "lanA") %>% 
  group_by(lan_gene) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  arrange(desc(total)) %>% 
  slice(1:50) %>% 
  mutate(lan_gene = factor(lan_gene, levels = lan_gene)) %>% 
  ggplot(aes(x = lan_gene, y = total)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_bw() +
  ylab("# of Samples") +
  xlab("Lan Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"), face = "bold", size = 15),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"), face = "bold", size = 15),
        axis.text = element_text(size = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
ggsave("../plots/analysis/lan_hits_top50.pdf", 
       device = "pdf",
       width = 40,
       height = 20,
       units = "cm")




# Number of species that the lan is found in
gene_calls %>% 
  filter(str_detect(lan_gene, "lanA")) %>%
  filter(lan_gene != "lanA") %>% 
  group_by(lan_gene) %>% 
  summarise(total = length(unique(highest_tax_hit))) %>% 
  ungroup() %>% 
  arrange(desc(total)) %>% 
  slice(1:100) %>% 
  mutate(lan_gene = factor(lan_gene, levels = lan_gene)) %>% 
  ggplot(aes(x = lan_gene, y = total)) +
  geom_bar(stat = "identity", fill = "black") +
  theme_bw() +
  ylab("# of Species") +
  xlab("Lan Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"), face = "bold", size = 15),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"), face = "bold", size = 15),
        axis.text = element_text(size = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
ggsave("../plots/analysis/lan_species_hits.pdf", 
       device = "pdf",
       width = 60,
       height = 20,
       units = "cm")



# Boxplot of lan lengths
gene_calls %>% 
  filter(str_detect(annotated_lan_genes, "lan")) %>% 
  filter(lan_gene != "lanA") %>% 
  mutate(annotated_lan_genes = factor(annotated_lan_genes, levels = c("lanA", "lanB", "lanC", "lanM", "lanKC",
                                                                      "lanF", "lanE", "lanG", "lanI",
                                                                      "lanT", "lanR", "lanK", "lanP"))) %>% 
  ggplot(aes(y = aa_length, x = annotated_lan_genes, fill = annotated_lan_genes)) +
  scale_fill_manual(values = gene_colors) +
  geom_boxplot() +
  labs( y = "AA Length") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        axis.title.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1100, 100))
ggsave("../plots/analysis/lan_gene_lengths.pdf", 
       device = "pdf",
       width = 20,
       height = 10,
       units = "cm")



