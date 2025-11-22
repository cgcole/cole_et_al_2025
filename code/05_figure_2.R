
fitness_data <- tibble("patient_ID" = NA, "ID"= NA, "shotgunSeq_id"= NA, "genus_species"= NA, 
                       "pctseqs"= NA, "lan_gene"= NA, "n"= NA, "mean_nonproducer_pctseqs"= NA, 
                       "n_nonproducers"= NA, "predicted_class" = NA,"label" = NA, "lanthipeptide" = NA)


figure2A_stats_results <- data.frame()


producing_strains <- clinical_donor_t %>% 
  mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                        blast_taxid), blast_genus_species)) %>% 
  mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                             taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
  mutate(blast_genus_species = if_else(str_detect(blast_genus_species," bacterium"), paste0(blast_genus_species, " ", 
                                                                                            blast_taxid), blast_genus_species)) %>% 
  mutate(genus_species = if_else(str_detect(genus_species, " bacterium"), paste0(genus_species, " ", 
                                                                                 taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
  filter(blast_genus_species == genus_species) %>% 
  select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs, lan_gene) %>% 
  distinct() %>% 
  group_by(lan_gene, genus_species) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  mutate(species_lan = paste0(genus_species, " ", lan_gene))


for (lan in unique(producing_strains$species_lan)) {
  lanthipeptide_producers <- clinical_donor_t %>% 
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                          blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                               taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," bacterium"), paste0(blast_genus_species, " ", 
                                                                                              blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " bacterium"), paste0(genus_species, " ", 
                                                                                   taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    filter(blast_genus_species == genus_species) %>% 
    select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs, lan_gene, predicted_class,label, lanthipeptide) %>% 
    distinct() %>% 
    group_by(lan_gene, genus_species) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    filter(n >= 5) %>% 
    mutate(species_lan = paste0(genus_species, " ", lan_gene)) %>% 
    filter(species_lan == lan) %>% 
    mutate(category = "producer")
  
  nonproducing_strains <- clinical_donor_t %>% 
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                          blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                               taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," bacterium"), paste0(blast_genus_species, " ", 
                                                                                              blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " bacterium"), paste0(genus_species, " ", 
                                                                                   taxid), genus_species)) %>% 
    mutate(lan_gene = if_else(is.na(lan_gene), "none", lan_gene)) %>% 
    mutate(species_lan = paste0(genus_species, " ", lan_gene, ",")) %>% 
    filter(genus_species %in% lanthipeptide_producers$genus_species) %>% 
    select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs, lan_gene, species_lan) %>% 
    group_by(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs) %>% 
    summarise(all_species_lan = paste0(species_lan, collapse = "")) %>% 
    ungroup() %>% 
    filter(!str_detect(all_species_lan, paste0(lan, ","))) %>% 
    mutate(category = "nonproducer")
  
  producer_nonproducer_table <- lanthipeptide_producers %>% 
    full_join(nonproducing_strains) 
  
  # Dunn's test with control comparisons only
  dunn_result <- dunn_test(producer_nonproducer_table, pctseqs ~ category, p.adjust.method = "none") %>% 
    mutate(species_lan = lan)
  
  figure2A_stats_results <- bind_rows(figure2A_stats_results, dunn_result)
}


figure2A_stats_results <- figure2A_stats_results %>%
  mutate(adj_p_value = p.adjust(p, method = "BH")) %>% 
  mutate(significance = case_when(
    adj_p_value < .0001 ~ "****",
    adj_p_value < .001 ~ "***",
    adj_p_value < .01 ~ "**",
    adj_p_value < .05 ~ "*",
    adj_p_value >= .05 ~ NA,
  ))





fitness_data <- data.frame()


for (lan in unique(producing_strains$species_lan)) {
  lanthipeptide_producers <- clinical_donor_t %>% 
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                          blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                               taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," bacterium"), paste0(blast_genus_species, " ", 
                                                                                              blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " bacterium"), paste0(genus_species, " ", 
                                                                                   taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    filter(blast_genus_species == genus_species) %>% 
    select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs, lan_gene, predicted_class,label, lanthipeptide) %>% 
    distinct() %>% 
    group_by(lan_gene, genus_species) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    filter(n >= 5) %>% 
    mutate(species_lan = paste0(genus_species, " ", lan_gene)) %>% 
    filter(species_lan == lan) %>% 
    mutate(category = "producer")
  
  nonproducing_strains <- clinical_donor_t %>% 
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," sp\\."), paste0(blast_genus_species, " ", 
                                                                                          blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " sp\\."), paste0(genus_species, " ", 
                                                                               taxid), genus_species)) %>% #add taxids to samples with no specified species for matching
    mutate(blast_genus_species = if_else(str_detect(blast_genus_species," bacterium"), paste0(blast_genus_species, " ", 
                                                                                              blast_taxid), blast_genus_species)) %>% 
    mutate(genus_species = if_else(str_detect(genus_species, " bacterium"), paste0(genus_species, " ", 
                                                                                   taxid), genus_species)) %>% 
    mutate(lan_gene = if_else(is.na(lan_gene), "none", lan_gene)) %>% 
    mutate(species_lan = paste0(genus_species, " ", lan_gene, ",")) %>% 
    filter(genus_species %in% lanthipeptide_producers$genus_species) %>% 
    select(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs, lan_gene, species_lan) %>% 
    group_by(patient_ID, ID, shotgunSeq_id, genus_species, pctseqs) %>% 
    summarise(all_species_lan = paste0(species_lan, collapse = "")) %>% 
    ungroup() %>% 
    filter(!str_detect(all_species_lan, paste0(lan, ","))) %>% 
    mutate(category = "nonproducer") %>% 
    group_by(genus_species) %>% 
    summarise(mean_nonproducer_pctseqs = mean(pctseqs),
              n_nonproducers = n()) %>% 
    ungroup()
  
  fc_table <- lanthipeptide_producers %>% 
    left_join(nonproducing_strains) %>% 
    mutate(fold_change = pctseqs/mean_nonproducer_pctseqs) 
  
  fitness_data <- bind_rows(fitness_data, fc_table)
}


fitness_results <- fitness_data %>% 
  group_by(genus_species, lan_gene, predicted_class, label, lanthipeptide, n, n_nonproducers, mean_nonproducer_pctseqs, species_lan) %>% 
  summarise(mean_fold_change = mean(fold_change)) %>% 
  ungroup() %>% 
  mutate(log2_mean_fold_change = log2(mean_fold_change)) %>% 
  left_join(figure2A_stats_results) %>% 
  mutate(species_lan_label = paste0("<i>", genus_species, "</i> ", lan_gene)) %>% 
  mutate(advantage = if_else(log2_mean_fold_change > 0, "yes", "no")) 






#Figure 2A
f2a <- fitness_results %>% 
  arrange(log2_mean_fold_change) %>% 
  mutate(species_lan_label = factor(species_lan_label, levels = unique(species_lan_label))) %>% 
  ggplot(aes(x = log2_mean_fold_change, y = species_lan_label, fill = predicted_class)) +
  geom_bar(stat= "identity", alpha = .85) +
  geom_text(aes(label = significance,
                hjust = ifelse(log2_mean_fold_change > 0, -0.1, 1.1),
                vjust = .75), 
            size = 3, na.rm = TRUE) +
  guides(fill = guide_legend(title="Lan Class"))  +
  labs(x  = expression(bold(Log[2] ~ "(Fold Change)")),) +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.text.y = element_markdown() ) +
  scale_fill_manual(values = lantibiotic_class_colors)  +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1, 2, 3, 4),
    limits = c(-2.5, 4),
    expand = expansion(mult = c(0.05, .1))  # add 5% space on left, 20% on right
  )
ggsave("./plots/figure2a.pdf", 
       device = "pdf",
       width = 10,
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


set.seed(834344)  

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
  left_join(fitness_results %>% 
              select(lan_gene, log2_mean_fold_change, advantage)) %>% 
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
                  max.overlaps = Inf, max.time = 200, segment.size = .15,
                  force = 100, box.padding = 0.5) +
  coord_cartesian(clip = "off", ylim = c(min(umap_data$UMAP2), max(umap_data$UMAP2) * 1.2)) 

ggsave("./plots/figure2b.pdf", 
       device = "pdf",
       width = 11.5,
       height = 9,
       units = "cm")


rm(f2a, f2b,
   producing_strains, nonproducing_strains, fitness_data, 
   lanA_sequences, sequences, alignment, alignment_converted,
   dist_matrix, umap_result, umap_data)


