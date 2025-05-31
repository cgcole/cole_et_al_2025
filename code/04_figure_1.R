
#Table of lanA data
lan_data <- clinical_donor_t %>% 
  select(patient_ID, lan_gene, predicted_class, blast_genus_species, blast_taxid, blast_genus, blast_species, lan_patient_count,
         lanthipeptide, label) %>% 
  distinct() %>% 
  filter(!is.na(lan_gene)) %>% 
  mutate(predicted_class = factor(predicted_class, levels = c("I", "II", "III", "Unknown"))) %>% 
  arrange(desc(lan_patient_count)) %>% 
  mutate(rank = dense_rank(desc(lan_patient_count))) 




#Figure 1A Patients that have a lanA genes
f1a <- clinical_donor_t %>% 
  group_by(patient_ID) %>%
  distinct(lan_gene) %>%  
  summarize(lanA_gene_count = sum(!is.na(lan_gene))) %>% 
  ungroup() %>% 
  mutate(lan_presence = if_else(lanA_gene_count >= 1, "Yes", "No")) %>% 
  group_by(lan_presence) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(aes(y = lan_presence, x = n)) +
  geom_bar(stat = "identity", fill = "black") +
  xlab("# of Patients") +
  ylab("lanA Genes Detected") +
  scale_x_continuous(expand = expansion(mult = c(0,.1))) 
ggsave("./plots/figure1a.pdf", 
       device = "pdf",
       width = 4,
       height = 2.75,
       units = "cm")


#Figure 1B
f1b <- clinical_donor_t %>% 
  select(patient_ID, lan_gene, predicted_class) %>% 
  distinct() %>% 
  filter(!is.na(lan_gene)) %>% 
  group_by(predicted_class) %>% 
  summarise(n = length(unique(lan_gene))) %>% 
  ungroup() %>% 
  mutate(predicted_class = factor(predicted_class, levels = c("Unknown", "III", "II", "I"))) %>% 
  ggplot(aes(y = predicted_class, x = n, fill = predicted_class)) +
  geom_bar(stat = "identity") +
  xlab("# of Unique lanA Genes") +
  ylab("Lan Class") +
  theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = c(0,.1))) +
  scale_fill_manual(values = lantibiotic_class_colors) 
ggsave("./plots/figure1b.pdf", 
       device = "pdf",
       width = 4,
       height = 2.75,
       units = "cm")


#Figure 1C
f1c <- clinical_donor_t %>% 
  filter(!is.na(lan_gene)) %>% 
  select(patient_ID, lan_gene, blast_genus_species, blast_genus, blast_species, predicted_class) %>% 
  distinct() %>% 
  group_by(blast_genus, predicted_class) %>% 
  summarise(total = n()) %>% 
  group_by(blast_genus) %>% 
  mutate(total_sum = sum(total)) %>% 
  ungroup() %>% 
  arrange(desc(total_sum)) %>% 
  mutate(rank = dense_rank(desc(total_sum))) %>% 
  filter(rank <= 21) %>% 
  filter(blast_genus != "Unclassified", blast_genus != "Eisenbergiella") %>% #This is because there was a tie! 
  mutate(blast_genus = factor(blast_genus, levels = unique(blast_genus))) %>%
  ggplot(aes(x = blast_genus, y = total, fill = predicted_class)) +
  geom_bar(stat = "identity") +
  ylab("# of lanA Genes Detected") +
  guides(fill = guide_legend(title="Lan Class"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.003,.1)))　+
  scale_fill_manual(values = lantibiotic_class_colors)
ggsave("./plots/figure1c.pdf", 
       device = "pdf",
       width = 7.5,
       height = 4.5,
       units = "cm")



#Figure 1D species for top hits
f1d_1 <- lan_data %>% 
  select(lan_gene, rank, predicted_class, lan_patient_count) %>% 
  distinct() %>% 
  filter(rank <= 16) %>% 
  arrange(desc(lan_patient_count)) %>% 
  mutate(lan_gene = factor(lan_gene, levels = unique(.$lan_gene))) %>% 
  ggplot(aes(x = lan_gene, y = lan_patient_count, fill = predicted_class)) +
  geom_bar(stat = "identity") +
  ylab("# of Patients") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(.005,.1))) +
  scale_fill_manual(values = lantibiotic_class_colors) +
  labs(fill='Lan Class')



x_axis_order <- lan_data %>% 
  select(lan_gene, rank, predicted_class, lan_patient_count) %>% 
  distinct() %>% 
  filter(rank <= 16) %>% 
  arrange(desc(lan_patient_count)) %>% 
  pull(lan_gene)

y_axis_order <- lan_data %>% 
  select(lan_gene, rank, predicted_class, blast_genus_species) %>% 
  filter(rank <= 16) %>% 
  group_by(blast_genus_species) %>% 
  summarise(n= n()) %>% 
  arrange(n) %>% 
  pull(blast_genus_species)


f1d_2 <- lan_data %>% 
  select(lan_gene, rank, predicted_class, blast_genus_species, lan_patient_count) %>% 
  filter(rank <= 16) %>% 
  group_by(lan_gene, blast_genus_species) %>% 
  summarise(n = n()) %>% 
  mutate(lan_gene = factor(lan_gene, levels = x_axis_order)) %>% 
  mutate(blast_genus_species = factor(blast_genus_species, levels = y_axis_order)) %>% 
  ggplot(aes(x = lan_gene, y = blast_genus_species)) +
  geom_point(aes(size = n)) +  # Dots scaled by counts
  scale_size_continuous(range = c(.05, 1.5)) +  # Adjust dot sizes
  ylab("Taxonomic Classification") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(face = "italic")) +
  scale_y_discrete(expand = c(.03,.03))



pdf(file.path("./plots/figure1d_1.pdf"), height = 3.75, width = 4)
gg.stack(f1d_1, f1d_2, 
         heights = c(2, 4),
         newpage = F)
dev.off()



f1d_3 <- lan_data %>% 
  select(lan_gene, rank, predicted_class, blast_genus_species, blast_taxid) %>% 
  filter(rank <= 16) %>% 
  group_by(blast_genus_species) %>% 
  summarise(n= n()) %>% 
  ungroup() %>% 
  arrange(n) %>% 
  mutate(blast_genus_species = factor(blast_genus_species, levels = unique(.$blast_genus_species))) %>% 
  ggplot(aes(x = n, y = blast_genus_species)) +
  geom_bar(stat = "identity", fill = "black") +
  xlab("# of Genes") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_discrete(expand = c(.03,.03))
ggsave("./plots/figure1d_2.pdf", 
       device = "pdf",
       width = 1,
       height = 2.58,
       units = "in")



rm(lan_data, f1a, f1b, f1c, f1d_1, f1d_2, f1d_3, y_axis_order, x_axis_order)





