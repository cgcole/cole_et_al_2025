
# Determine fitness ratio
producing_strains <- clinical_donor_t %>% 
  mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                        blast_taxid), blast_genus_species)) %>% 
  mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                             blast_taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
  filter(blast_genus_species == genus_species) %>% 
  group_by(lan_gene, blast_genus_species, genus_species) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5)


nonproducing_strains <- clinical_donor_t %>% 
  mutate(blast_genus_species = if_else(blast_species == " sp\\.", paste0(blast_genus_species, " ", 
                                                                         blast_taxid), blast_genus_species)) %>% 
  mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                             taxid), genus_species)) %>%
  filter(blast_genus_species != genus_species | is.na(blast_genus_species) | is.na(genus_species)) %>%
  select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs) %>% 
  distinct() %>% 
  group_by(genus_species) %>% 
  summarise(mean_nonproducer_pctseqs = mean(pctseqs),
            n_nonproducer = n()) %>% 
  ungroup() %>% 
  filter(n_nonproducer >= 5)


fitness_data <- producing_strains %>% 
  left_join(nonproducing_strains) %>% 
  filter(!is.na(mean_nonproducer_pctseqs)) %>% 
  mutate(fold_change = pctseqs/mean_nonproducer_pctseqs) %>% 
  group_by(lan_gene, predicted_class,label, lanthipeptide) %>% 
  summarise(mean_fold_change = mean(fold_change),
            n_lan_gene = n()) %>% 
  ungroup() %>% 
  mutate(lan_gene = if_else(!is.na(label), label, lan_gene)) %>% 
  mutate(log2_mean_fold_change = log2(mean_fold_change)) %>% 
  arrange(desc(log2_mean_fold_change)) %>% 
  mutate(lan_gene = factor(lan_gene, levels = unique(lan_gene)))  %>% 
  mutate(predicted_class = if_else(is.na(predicted_class), "Unknown", as.character(predicted_class))) %>% 
  mutate(predicted_class = factor(predicted_class, levels = c("I", "II", "III", "Unknown"))) %>% 
  mutate(advantage = if_else(log2_mean_fold_change > 0, "yes", "no")) 


#Figure 2A
f2a <- fitness_data %>% 
  arrange(log2_mean_fold_change) %>% 
  mutate(lan_gene = factor(lan_gene, levels = unique(lan_gene))) %>% 
  ggplot(aes(x = log2_mean_fold_change, y = lan_gene, fill = predicted_class)) +
  geom_bar(stat= "identity", alpha = .85) +
  guides(fill = guide_legend(title="Lan Class"))  +
  labs(x  = expression(bold(Log[2] ~ "(Fold Change)")),) +
  theme(legend.position = "right",
        axis.title.y = element_blank()) +
  scale_fill_manual(values = lantibiotic_class_colors) 
ggsave("./plots/figure2a.pdf", 
       device = "pdf",
       width = 7,
       height = 18,
       units = "cm")



#Figure 2B
lanA_sequences <- clinical_donor_t %>% 
  select(lan_gene, aa_sequence) %>%
  filter(!is.na(lan_gene)) %>% 
  distinct() %>% 
  mutate(order = str_replace(lan_gene, "lanA_", "")) %>% 
  type_convert() %>% 
  arrange(order) %>% 
  select(-order)

write.fasta(as.list(lanA_sequences$aa_sequence), lanA_sequences$lan_gene, "./fasta_files/clinical_lanA_aa_sequences.fa")





sequences <- readAAStringSet("./fasta_files/clinical_lanA_aa_sequences.fa")
alignment <- msa(sequences, method = "Muscle")
alignment_converted <- msaConvert(alignment, type = "seqinr::alignment")
dist_matrix <- dist.alignment(alignment_converted, "identity")
umap_result <- umap(as.matrix(dist_matrix))


umap_data <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2])

umap_data <- umap_data %>% 
  rownames_to_column(var = "lan_gene") %>% 
  left_join(clinical_donor_t %>% 
              select(lan_gene, predicted_class, aa_sequence) %>% 
              distinct() %>% 
              filter(!is.na(lan_gene))) %>% 
  left_join(clinical_donor_t %>% 
              select(lan_gene, lanthipeptide, label, lan_patient_count) %>% 
              distinct() %>% 
              filter(!is.na(lan_gene)) %>% 
              arrange(lan_patient_count) %>% 
              distinct(lanthipeptide, .keep_all = TRUE)) %>% 
  mutate(predicted_class = if_else(is.na(predicted_class), "Unknown", as.character(predicted_class))) %>% 
  mutate(predicted_class = factor(predicted_class, levels = c("I", "II", "III", "Unknown"))) %>% 
  left_join(fitness_data %>% 
              select(lan_gene, advantage, log2_mean_fold_change)) %>% 
  mutate(advantage = if_else(is.na(advantage), "no", advantage))



f2b <- umap_data %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = predicted_class, fill = predicted_class,
             shape = advantage, size = advantage, label = label)) +
  geom_point(alpha = .75) +
  labs(x = "UMAP Axis 1", y = "UMAP Axis 2") +
  scale_color_manual(values = lantibiotic_class_colors) +
  scale_fill_manual(values = lantibiotic_class_colors) +
  guides(color = guide_legend(title="Lan Class"), 
         shape = "none", 
         size = "none",
         fill = "none")  +
  scale_size_manual(values = c("yes" = .75, "no" = .25)) +
  geom_text_repel(size = 1.5, na.rm = TRUE, show.legend = FALSE, 
                  max.overlaps = 100, max.time = 60, segment.size = .15,
                  force = 5) +
  coord_cartesian(clip = "off", ylim = c(min(umap_data$UMAP2), max(umap_data$UMAP2) * 1.2)) 

ggsave("./plots/figure2b.pdf", 
       device = "pdf",
       width = 9.5,
       height = 7,
       units = "cm")


rm(f2a, f2b,
   producing_strains, nonproducing_strains, fitness_data, 
   lanA_sequences, sequences, alignment, alignment_converted,
   dist_matrix, umap_result, umap_data)