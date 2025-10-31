
#Table of lanA data
lan_data <- clinical_donor_t %>% 
  select(patient_ID, lan_patient_count, lan_gene, predicted_class, blast_genus_species, blast_taxid, blast_genus, blast_species, 
         lanthipeptide, label) %>% 
  distinct() %>% 
  filter(!is.na(lanthipeptide)) %>% 
  mutate(predicted_class = factor(predicted_class, levels = c("I", "II", "III", "Unknown"))) %>% 
  group_by(lanthipeptide) %>% 
  mutate(lan_patient_count = length(unique(patient_ID))) %>% 
  ungroup() %>% 
  arrange(desc(lan_patient_count)) %>% 
  mutate(rank = dense_rank(desc(lan_patient_count))) 



# Figure S1A
sf1a <- clinical_donor_t %>% 
  filter(!is.na(lan_gene)) %>% 
  select(patient_ID, lan_gene, blast_genus_species, blast_genus, blast_species, predicted_class, blast_taxid) %>% 
  distinct() %>% 
  mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                         blast_taxid), blast_genus_species)) %>% 
  group_by(blast_genus_species, predicted_class) %>% 
  summarise(total = n()) %>% 
  group_by(blast_genus_species) %>% 
  mutate(total_sum = sum(total)) %>% 
  ungroup() %>% 
  arrange(desc(total_sum)) %>% 
  mutate(rank = dense_rank(desc(total_sum))) %>% 
  filter(rank <= 21) %>% 
  filter(blast_genus_species != "Unclassified") %>% 
  mutate(blast_genus_species = factor(blast_genus_species, levels = unique(blast_genus_species))) %>%
  ggplot(aes(x = blast_genus_species, y = total, fill = predicted_class)) +
  geom_bar(stat = "identity") +
  ylab("# of lanA Genes Detected") +
  guides(fill = guide_legend(title="Lan Class"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.003,.1)))　+
  scale_fill_manual(values = lantibiotic_class_colors)
ggsave("./plots/supplemental_figure1a.pdf", 
       device = "pdf",
       width = 9.1,
       height = 4.75,
       units = "cm")






# Figure S1B species for top hits
sf1b_1 <- lan_data %>% 
  select(rank, predicted_class, lan_patient_count, lanthipeptide) %>% 
  distinct() %>% 
  arrange(desc(lan_patient_count)) %>% 
  mutate(lanthipeptide = factor(lanthipeptide, levels = unique(.$lanthipeptide))) %>% 
  ggplot(aes(x = lanthipeptide, y = lan_patient_count, fill = predicted_class)) +
  geom_bar(stat = "identity") +
  ylab("# of Patients") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(.005,.1))) +
  scale_fill_manual(values = lantibiotic_class_colors) +
  labs(fill='Lan Class')



x_axis_order <- lan_data %>% 
  select(lanthipeptide, rank, predicted_class, lan_patient_count) %>% 
  distinct() %>% 
  arrange(desc(lan_patient_count)) %>% 
  pull(lanthipeptide)

y_axis_order <- lan_data %>% 
  select(lanthipeptide, rank, predicted_class, blast_genus_species) %>% 
  group_by(blast_genus_species) %>% 
  summarise(n= n()) %>% 
  arrange(n) %>% 
  pull(blast_genus_species)


sf1b_2 <- lan_data %>% 
  select(lanthipeptide, rank, predicted_class, blast_genus_species, lan_patient_count) %>% 
  group_by(lanthipeptide, blast_genus_species) %>% 
  summarise(n = n()) %>% 
  mutate(lanthipeptide = factor(lanthipeptide, levels = x_axis_order)) %>% 
  mutate(blast_genus_species = factor(blast_genus_species, levels = y_axis_order)) %>% 
  ggplot(aes(x = lanthipeptide, y = blast_genus_species)) +
  geom_point(aes(size = n)) +  # Dots scaled by counts
  scale_size_continuous(range = c(.05, 1.5)) +  # Adjust dot sizes
  ylab("Taxonomic Classification") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic")) +
  scale_y_discrete(expand = c(.03,.03))



pdf(file.path("./plots/supplemental_figure1b_1.pdf"), height = 3.75, width = 5)
gg.stack(sf1b_1, sf1b_2, 
         heights = c(2, 4),
         newpage = F)
dev.off()



sf1b_3 <- lan_data %>% 
  select(lanthipeptide, rank, predicted_class, blast_genus_species, blast_taxid) %>% 
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
ggsave("./plots/supplemental_figure1b_2.pdf", 
       device = "pdf",
       width = 1,
       height = 2.41,
       units = "in")


rm(sf1a, sf1b_1, sf1b_2, sf1b_3, x_axis_order, y_axis_order, lan_data)

