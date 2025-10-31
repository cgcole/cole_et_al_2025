


lantibiotic_analysis <- gene_calls %>% 
  group_by(gene_callers_id) %>% 
  filter(any(str_detect(aa_sequence, lanA_genes$aa_sequence))) %>% 
  ungroup()
  
lantibiotic_analysis <- gene_calls %>%
  filter(include == "yes") 


#Get the aa_sequence for partial matches to a rodeo_aa_sequence
rodeo_gene_matches <- rodeo_table %>% 
  left_join(gene_calls %>% 
              select(gene_callers_id, locus_id, contig, start, stop, direction, aa_sequence, aa_length, 
                     COG20_FUNCTION_accession, COG20_FUNCTION_function, KOfam_accession, KOfam_function, Pfam_accession, Pfam_function, 
                     bakta_accession, bakta_function, ghostkoala_accession, ghostkoala_function), 
            relationship = "many-to-many") %>% 
  mutate(rodeo_match = if_else(aa_sequence != rodeo_aa_sequence & str_detect(aa_sequence, rodeo_aa_sequence) & rodeo_start >= start & rodeo_stop <= stop, 
                               "Partial Match", 
                               if_else(aa_sequence == rodeo_aa_sequence & rodeo_start == start & rodeo_stop == stop, "Full Match", NA))) %>% 
  filter(!is.na(rodeo_match)) %>% 
  full_join(rodeo_table %>% 
              select(ID, locus_id, Leader, Core, Start, End, LEN, `Class Call`, rodeo_aa_sequence, rodeo_start, rodeo_stop, rodeo_direction, fasta_file)) %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_match), "No Match", rodeo_match)) %>% 
  mutate(lanA_source = "RODEO") %>% 
  mutate(COG20_FUNCTION_accession = if_else(rodeo_match == "No Match", "N/A", COG20_FUNCTION_accession)) %>% 
  mutate(COG20_FUNCTION_function = if_else(rodeo_match == "No Match", "N/A", COG20_FUNCTION_function)) %>% 
  mutate(KOfam_accession = if_else(rodeo_match == "No Match", "N/A", KOfam_accession)) %>% 
  mutate(KOfam_function = if_else(rodeo_match == "No Match", "N/A", KOfam_function)) %>% 
  mutate(Pfam_accession = if_else(rodeo_match == "No Match", "N/A", Pfam_accession)) %>% 
  mutate(Pfam_function = if_else(rodeo_match == "No Match", "N/A", Pfam_function)) %>% 
  mutate(bakta_accession = if_else(rodeo_match == "No Match", "N/A", bakta_accession)) %>% 
  mutate(bakta_function = if_else(rodeo_match == "No Match", "N/A", bakta_function)) %>% 
  mutate(ghostkoala_accession = if_else(rodeo_match == "No Match", "N/A", ghostkoala_accession)) %>% 
  mutate(ghostkoala_function = if_else(rodeo_match == "No Match", "N/A", ghostkoala_function)) #These are to fix none annotated genes
  
  
table(rodeo_gene_matches$rodeo_match)






lanA_gene_calls <- gene_calls %>% 
  filter(str_detect(lan_gene, "lanA_"))


write_csv(lanA_gene_calls, "~/Desktop/lanA_gene_calls.csv")












lanA_gene_calls <- gene_calls %>% 
  filter(annotated_lan_genes == "lanA") %>% 
  mutate(orphan_lanA = if_else(str_detect(all_lan_genes, "lanA") & str_detect(all_lan_genes, "lanC|lanM|lanKC"), "No",
                              if_else(str_detect(all_lan_genes, "lanA") & !str_detect(all_lan_genes, "lanC|lanM|lanKC"), "Yes", NA))) %>% 
  full_join(rodeo_table %>% 
               select(locus_id, rodeo_start, rodeo_stop, rodeo_direction, rodeo_aa_sequence)) %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_match), "Not Detected", rodeo_match)) %>% 
  mutate(lan_gene = if_else(is.na(lan_gene), "" , lan_gene)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes), " " , annotated_lan_genes))






lanA_gene_calls <- gene_calls %>% 
  group_by(gene_callers_id) %>% 
  filter(any(str_detect(aa_sequence, lanA_genes$aa_sequence))) %>% 
  left_join()
  
  
  mutate(lanA = if_else(any(str_detect(aa_sequence, lanA_genes$aa_sequence), "Yes", "No")))
  filter(any(str_detect(aa_sequence, lanA_genes$aa_sequence))) %>% 
  mutate(orphan_lanA = if_else(str_detect(all_lan_genes, "lanA") & str_detect(all_lan_genes, "lanC|lanM|lanKC"), "No",
                               if_else(str_detect(all_lan_genes, "lanA") & !str_detect(all_lan_genes, "lanC|lanM|lanKC"), "Yes", NA))) %>% 
  full_join(rodeo_gene_matches) %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_match), "Not Detected", rodeo_match)) %>% 
  mutate(lan_gene = if_else(is.na(lan_gene), "" , lan_gene)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes), " " , annotated_lan_genes))
  

table(lanA_gene_calls$orphan_lanA)  
table(lanA_gene_calls$lanA_source) 
table(lanA_gene_calls$rodeo_match)

write_csv(lanA_gene_calls, "~/Desktop/lanA_gene_calls.csv")

 



#Plot contigs for all full rodeo matches
lanA_full_matches <- lanA_gene_calls %>% 
  filter(rodeo_match == "Full Match") %>% 
  select(contig) %>% 
  distinct(contig) %>% 
  pull(contig)



for (i in lanA_full_matches) {
  
  plot_height <- gene_calls %>% 
    filter(contig == i) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  gene_calls %>% 
    filter(contig == i) %>%  
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    scale_fill_manual(values = gene_colors) +
    facet_grid(highest_tax_hit~., scales = "free", space = "free_y") +
    theme_bw() +
    labs(y = "Contig", x = "") +
    theme(legend.position = "none",
          axis.ticks.x = element_line(linewidth = .75, color = "black"),
          axis.line.x = element_line(color = "black", linewidth = .75),
          axis.text = element_text(size = 9, color = "black"),
          title = element_text(face = "bold", size = 12),
          axis.ticks.length = unit(.17, "cm"),
          panel.grid.major.y = element_line(size = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 

  
  ggsave(paste0("../plots/rodeo_matches/full_matches/", i, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}




#Plot contigs for all partial rodeo matches
lanA_partial_matches <- lanA_gene_calls %>% 
  filter(rodeo_match == "Partial Match") %>% 
  select(contig) %>% 
  distinct(contig) %>% 
  pull(contig)



for (i in lanA_partial_matches) {
  
  plot_height <- gene_calls %>% 
    filter(contig == i) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  partial_lan_gene_calls <- lanA_gene_calls %>% 
    filter(contig == i, rodeo_match == "Partial Match") %>% 
    select(contig, start, stop, direction, annotated_lan_genes, rodeo_start, rodeo_stop, rodeo_direction, rodeo_aa_sequence) %>% 
    left_join(lanA_genes %>% 
                dplyr::rename(rodeo_aa_sequence = aa_sequence) %>% 
                filter(.$rodeo_aa_sequence %in% rodeo_lanA_genes$rodeo_aa_sequence))
  
  gene_calls %>% 
    filter(contig == i) %>%  
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = partial_lan_gene_calls, aes(xsubmin = rodeo_start, xsubmax = rodeo_stop), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = partial_lan_gene_calls, aes(x = (rodeo_start + rodeo_stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = -0.3) + 
    scale_fill_manual(values = gene_colors) +
    facet_grid(highest_tax_hit~., scales = "free", space = "free_y") +
    theme_bw() +
    labs(y = "Contig", x = "") +
    theme(legend.position = "none",
          axis.ticks.x = element_line(linewidth = .75, color = "black"),
          axis.line.x = element_line(color = "black", linewidth = .75),
          axis.text = element_text(size = 9, color = "black"),
          title = element_text(face = "bold", size = 12),
          axis.ticks.length = unit(.17, "cm"),
          panel.grid.major.y = element_line(size = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  
  ggsave(paste0("../plots/rodeo_matches/partial_matches/", i, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}









lanA_partial_matches <- lanA_gene_calls %>% 
  filter(rodeo_match == "Partial Match") %>% 
  select(contig) %>% 
  distinct(contig) %>% 
  pull(contig)

lanA_no_matches_matches <- lanA_gene_calls %>% 
  filter(rodeo_match == "Partial Match") %>% 
  select(contig) %>% 
  distinct(contig) %>% 
  pull(contig)






















#rodeo analysis
rodeo_analysis <- read_tsv("../metadata/lan_gene_calls.txt") %>%
  mutate(locus_id = str_extract(contig, ".*(?=_)")) %>% 
  filter(locus_id %in% modified_rodeo_table$locus_id) %>% 
  filter(aa_sequence %in% modified_rodeo_table$rodeo_aa_sequence) %>% 
  group_by(locus_id, aa_sequence) %>% 
  mutate(lan_gene_number = seq_along(aa_sequence)) %>% 
  ungroup() %>% 
  select(locus_id, contig, aa_sequence, lan_gene_number) %>% 
  full_join(modified_rodeo_table) %>% 
  mutate(rodeo_match = if_else(rodeo_aa_sequence == aa_sequence, "Match", rodeo_match)) %>% 
  mutate(rodeo_match = if_else(is.na(contig), "RODEO Only", rodeo_match)) %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_aa_sequence), "Prodigal Only", rodeo_match)) %>% 
  left_join(combined_gene_annotations %>%
              mutate(locus_id = str_extract(contig, ".*(?=_)")) %>%
              select(locus_id, aa_sequence, KOfam_function, Pfam_function, ghostkoala_function, bakta_function) %>% 
              distinct())


rodeo_analysis_t1 <- rodeo_analysis %>% 
  select(rodeo_match) %>% 
  group_by(rodeo_match) %>% 
  summarise(n = n())


write_csv(rodeo_analysis_t1, "../tables/rodeo_analysis_table_1.csv")
rm(rodeo_analysis_t1)





rodeo_analysis_t2 <- rodeo_analysis %>% 
  distinct(contig) %>% 
  filter(!is.na(contig)) %>% 
  summarise(n = n())


write_csv(rodeo_analysis_t2, "../tables/rodeo_analysis_table_2.csv")
rm(rodeo_analysis_t2)


rodeo_analysis_t3 <- rodeo_analysis %>% 
  left_join(lan_genes) %>% 
  filter(!is.na(lan_gene)) %>% 
  group_by(lan_gene) %>% 
  summarise(n = n())



rodeo_table_annotations <- modified_rodeo_table %>% 
  left_join(bakta_annotations %>% 
              select(aa_sequence, bakta_function) %>% 
              filter(!is.na(bakta_function)) %>% 
              distinct()) %>% 
  left_join(combined_gene_annotations %>% 
              select(aa_sequence, ghostkoala_function) %>% 
              filter(!is.na(ghostkoala_function)) %>% 
              distinct())


