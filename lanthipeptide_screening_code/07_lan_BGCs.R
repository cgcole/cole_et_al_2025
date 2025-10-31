library(gggenes)
library(ggrepel)


lanA_class <- read_csv("../metadata/lanA_class.csv")

# Make a table of all lanA genes with names and class
lanA_genes <- rodeo_lanA_genes %>% 
  dplyr::rename(aa_sequence = rodeo_aa_sequence) %>% 
  full_join(read_tsv("../metadata/alternate_lanA_genes.txt", col_names = "aa_sequence")) %>% 
  full_join(lanA_class) %>% 
  mutate(order = as.numeric(str_replace(lan_gene, "lanA_", ""))) %>% 
  arrange(order)%>% 
  select(-order) %>% 
  rownames_to_column(var = "row") %>% 
  mutate(lan_gene = if_else(is.na(lan_gene), paste0("lanA_",row), lan_gene)) %>% 
  select(-row) %>% 
  mutate(detected = if_else(aa_sequence %in% combined_gene_annotations$aa_sequence, "yes", "no"))

new_lanA_genes <- lanA_genes %>% 
  filter(!(lan_gene %in% lanA_class$lan_gene)) %>% 
  filter(!(aa_sequence %in% blast_lanA_genes$aa_sequence))

#Use this to remake the lanA_class file if lanA_gene numbers change
#write_csv(lanA_genes, "../metadata/lanA_class.csv")

# Generate fasta file of all lanA contigs to run blastp on
#write.fasta(as.list(lanA_genes$aa_sequence), lanA_genes$lan_gene, "../fasta_files/lanA_genes.fa")
#write.fasta(as.list(new_lanA_genes$aa_sequence), new_lanA_genes$lan_gene, "../fasta_files/lanA_genes_02.fa")
#write.fasta(as.list(new_lanA_genes$aa_sequence), new_lanA_genes$lan_gene, "../fasta_files/lanA_genes_03.fa")
#write.fasta(as.list(new_lanA_genes$aa_sequence), new_lanA_genes$lan_gene, "../fasta_files/lanA_genes_04.fa")
#write.fasta(as.list(new_lanA_genes$aa_sequence), new_lanA_genes$lan_gene, "../fasta_files/lanA_genes_05.fa")
#write.fasta(as.list(new_lanA_genes$aa_sequence), new_lanA_genes$lan_gene, "../fasta_files/lanA_genes_06.fa")

# Identify the aa_sequence for partial and full matches to a rodeo_aa_sequence
rodeo_gene_calls <- combined_gene_annotations %>% 
  select(gene_callers_id, contig, start, stop, direction, aa_sequence, aa_length) %>% 
  mutate(locus_id = str_extract(contig, ".*(?=_)")) %>%
  right_join(rodeo_table %>% 
               select(locus_id, rodeo_start, rodeo_stop, 
                      rodeo_direction, rodeo_aa_sequence), 
             relationship = "many-to-many") %>% 
  mutate(rodeo_match = if_else(aa_sequence != rodeo_aa_sequence & str_detect(aa_sequence, 
                                                                             rodeo_aa_sequence) & rodeo_start >= start & rodeo_stop <= stop, 
                               "Partial Match", 
                               if_else(aa_sequence == rodeo_aa_sequence & rodeo_start == start & rodeo_stop == stop,
                                       "Full Match", NA))) %>% 
  filter(!is.na(rodeo_match)) %>% 
  full_join(rodeo_table %>% 
              select(locus_id, rodeo_start, rodeo_stop, rodeo_direction, rodeo_aa_sequence)) %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_match), "No Match", rodeo_match)) %>% 
  mutate(aa_sequence = if_else(is.na(aa_sequence), rodeo_aa_sequence, aa_sequence)) %>% 
  left_join(lanA_genes %>% 
              select(-class) %>% 
              dplyr::rename(rodeo_lanA_gene = lan_gene,
                            rodeo_aa_sequence = aa_sequence)) %>% 
  left_join(lanA_genes)
  

table(rodeo_gene_calls$rodeo_match)


# Generate the full gene call table with all information and identify lan genes
gene_calls <- combined_gene_annotations %>% 
  arrange(gene_callers_id) %>% 
  mutate(locus_id = str_extract(contig, ".*(?=_)")) %>% 
  left_join(lanA_genes) %>% 
  left_join(rodeo_gene_calls %>% 
              filter(rodeo_match != "No Match")) %>% 
  mutate(lan_gene = if_else(is.na(lan_gene) & !is.na(rodeo_lanA_gene), "lanA", lan_gene)) %>% 
  mutate(annotated_lan_genes = case_when(
    ghostkoala_accession == "K20484" ~ "lanC",
    ghostkoala_accession == "K27862" ~ "lanM",
    ghostkoala_accession == "K27863" ~ "lanL",
    ghostkoala_accession == "K24914" ~ "lanKC",
    ghostkoala_accession == "K20483" ~ "lanB",
    ghostkoala_accession == "K20485" ~ "lanT",
    ghostkoala_accession %in% c("K20487", "K14988") ~ "lanK",
    ghostkoala_accession %in% c("K20488", "K14989") ~ "lanR",
    ghostkoala_accession %in% c("K20490", "K20459") ~ "lanF",
    ghostkoala_accession %in% c("K20491", "K20460") ~ "lanE",
    ghostkoala_accession %in% c("K20492", "K20461") ~ "lanG",
    ghostkoala_accession == "K20486" ~ "lanP")) %>% 
  left_join(read_csv("../metadata/bakta_lan_annotations.csv")) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes) & 
                                         Pfam_accession %in% c("PF05147,PF13575", "PF00330,PF05147",
                                                               "PF01532,PF05147,PF13575", "PF05147,PF07944,PF13575"),
                                       "lanM", annotated_lan_genes)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes) & 
                                         Pfam_accession == "PF18218", "lanI", annotated_lan_genes)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes) & 
                                         Pfam_accession == "PF05147", "lanC", annotated_lan_genes)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes) & !is.na(bakta_lan_gene), 
                                       bakta_lan_gene, annotated_lan_genes)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes) & str_detect(lan_gene, "lanA"), 
                                       "lanA", annotated_lan_genes)) %>% 
  mutate(lan_gene = if_else(is.na(lan_gene), annotated_lan_genes, lan_gene)) %>% 
  mutate(lan_gene = if_else(is.na(lan_gene), "" , lan_gene)) %>% 
  mutate(rodeo_lanA_gene = if_else(is.na(rodeo_lanA_gene), "" , rodeo_lanA_gene)) %>% 
  mutate(annotated_lan_genes = if_else(is.na(annotated_lan_genes), " " , annotated_lan_genes)) %>% 
  group_by(contig) %>% 
  mutate(lanA_count = sum(annotated_lan_genes == "lanA")) %>%
  mutate(all_lan_genes = paste0(lan_gene, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(all_lan_genes = paste0(all_lan_genes, ",")) %>% #This will be used to search for the lan otherwise things like lanA_1 would also match lanA_11
  # group_by(gene_callers_id, aa_sequence) %>% 
  # mutate(lanA_source = case_when(any(str_detect(aa_sequence, rodeo_lanA_genes$rodeo_aa_sequence)) ~ "RODEO", #Will label partial matches as RODEO DOESN'T work I need a new method
  #                                aa_sequence %in% interpro_lanA_genes$aa_sequence ~ "Interpro",
  #                                aa_sequence %in% biobank_lanA_genes$aa_sequence ~ "Biobank",
  #                                aa_sequence %in% annotated_lanA_genes$aa_sequence ~ "Annotations",
  #                                aa_sequence %in% manual_called_lanA_genes$aa_sequence ~ "Manual")) %>% 
  # ungroup() %>% 
  mutate(rodeo_match = if_else(is.na(rodeo_match) & annotated_lan_genes == "lanA", "Not Detected", rodeo_match)) %>% 
  select(gene_callers_id, locus_id, contig, start, stop, direction, aa_sequence, aa_length, 
         COG20_FUNCTION_accession, COG20_FUNCTION_function, #KOfam_accession, KOfam_function, 
         Pfam_accession, Pfam_function, 
         bakta_accession, bakta_function, ghostkoala_accession, ghostkoala_function, 
         annotated_lan_genes, lan_gene, class, predicted_class, rodeo_lanA_gene, all_lan_genes, lanA_count, highest_tax_hit,
         rodeo_start, rodeo_stop, rodeo_direction, rodeo_aa_sequence, rodeo_match, include, predicted_lanA) %>% 
  group_by(gene_callers_id, contig, start, stop) %>% 
  mutate(bakta_function = paste(bakta_function, collapse = ",")) %>%
  ungroup() %>% 
  distinct()


lan_gene_ranges <- gene_calls %>% 
  group_by(contig) %>% 
  mutate(gene_number = row_number()) %>% 
  ungroup() %>% 
  filter(str_detect(annotated_lan_genes, "lan")) %>% 
  group_by(contig) %>% 
  mutate(max_gene_call_id = max(gene_number) + 7,
         min_gene_call_id =min(gene_number) - 7) %>% 
  ungroup() %>% 
  select(contig, max_gene_call_id, min_gene_call_id) %>% 
  distinct()

relevant_genes <- gene_calls %>% 
  group_by(contig) %>% 
  mutate(gene_number = row_number()) %>% 
  ungroup() %>% 
  left_join(lan_gene_ranges) %>% 
  filter(min_gene_call_id <= gene_number & gene_number <= max_gene_call_id) %>% 
  group_by(contig) %>% 
  mutate(min_start = min(start) - 1) %>% 
  ungroup() %>% 
  mutate(start = start - min_start,
         stop = stop - min_start)


write_csv(relevant_genes, "../data/gene_call_table.csv")


# Samples and lanA genes used in the plot for loop
sample_ids <- relevant_genes %>% 
  select(locus_id) %>% 
  distinct() %>% 
  pull()


lan_ids <- lanA_genes %>% 
  mutate(lan_gene = paste0(lan_gene, ",")) %>% 
  pull(lan_gene)


confirmed_lan_ids <- lanA_genes %>% 
  filter(include == "yes") %>% 
  mutate(lan_gene = paste0(lan_gene, ",")) %>% 
  pull(lan_gene)

unconfirmed_lan_ids <- lanA_genes %>% 
  filter(is.na(include)) %>% 
  mutate(lan_gene = paste0(lan_gene, ",")) %>% 
  pull(lan_gene)


# Color scheme for lan genes
gene_colors <- c("lanA" = "#FB8072",
                 "lanF" = "#8DD3C7",
                 "lanE" = "#A6CEE3",
                 "lanG" = "#BEBADA",
                 "lanI" = "#80B1D3",
                 "lanNSR" = "#80B1D3",
                 "lanB" = "#FDB462",
                 "lanT" = "#B3DE69",
                 "lanC" = "#FCCDE5",
                 "lanM" = "#33A02C",
                 "lanR" = "#BC80BD",
                 "lanK" = "#CCEBC5",
                 "lanKC" = "#FFED6F",
                 "lanP" = "#1F78B4",
                 " " = "darkgrey")






# Plot contigs for all the confirmed lantibiotics
for (i in confirmed_lan_ids) {
  
  plot_height <- relevant_genes %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  rodeo_lan_annotation <- relevant_genes %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    mutate(rodeo_start = start - min_start,
           rodeo_stop = stop - min_start)
  
  relevant_genes %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                        xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  name <- str_replace(i, ",", "")
  
  ggsave(paste0("../plots/confirmed_lantibiotics/", name, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}




# Plot contigs for all the unconfirmed lantibiotics
for (i in unconfirmed_lan_ids) {
  
  plot_height <- relevant_genes %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  rodeo_lan_annotation <- relevant_genes %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    mutate(rodeo_start = start - min_start,
           rodeo_stop = stop - min_start)
  
  relevant_genes %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                        xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  name <- str_replace(i, ",", "")
  
  ggsave(paste0("../plots/unconfirmed_lantibiotics/", name, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}





# Plot contigs for each lantibiotic
for (i in lan_ids) {
  
  plot_height <- gene_calls %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  rodeo_lan_annotation <- gene_calls %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(str_detect(all_lan_genes, i))
  
  gene_calls %>% 
    filter(str_detect(all_lan_genes, i)) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                          xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  name <- str_replace(i, ",", "")
  
  ggsave(paste0("../plots/lantibiotics/", name, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}


# Plot each contig
for (i in unique(gene_calls$contig)) {
  
  plot_height <- gene_calls %>% 
    filter(contig == i) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  rodeo_lan_annotation <- gene_calls %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(str_detect(all_lan_genes, i))
  
  gene_calls %>% 
    filter(contig == i) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                          xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  
  ggsave(paste0("../plots/contigs/", i, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}





partial_rodeo_lan_ids <- gene_calls %>% 
  select(rodeo_lanA_gene, rodeo_match) %>% 
  filter(str_detect(rodeo_lanA_gene, "lanA"), rodeo_match == "Partial Match") %>% 
  distinct() %>% 
  arrange(rodeo_lanA_gene) %>% 
  pull(rodeo_lanA_gene)



# Plot contigs for each lantibiotic
for (i in partial_rodeo_lan_ids) {
  
  plot_height <- gene_calls %>% 
    filter(rodeo_lanA_gene == i) %>% 
    summarise(n = length(unique(contig))) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  
  contigs <- gene_calls %>% 
    filter(rodeo_lanA_gene == i) %>% 
    select(contig)
  
  rodeo_lan_annotation <- gene_calls %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(contig %in% contigs$contig)
  
  gene_calls %>% 
    filter(contig %in% contigs$contig) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                        xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0,8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 

  
  ggsave(paste0("../plots/partial_rodeo_lantibiotics/", i, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
}


rm(plot_height, rodeo_lan_annotation, contigs, i, name)




# Plot contigs for each sample
for (i in sample_ids) {
  
  plot_height <- gene_calls %>% 
    filter(locus_id == i) %>% 
    distinct(contig) %>% 
    summarise(n = n()) %>% 
    pull(n)
  
  plot_height <- plot_height * 2
  
  if (plot_height == 0) {
    next
  }
  
  rodeo_lan_annotation <- gene_calls %>% 
    filter(str_detect(rodeo_lanA_gene, "lanA")) %>% 
    filter(locus_id == i)
  
  gene_calls %>% 
    filter(locus_id == i) %>% 
    ggplot(aes(xmin = start, xmax = stop, y = contig, fill = annotated_lan_genes, forward = direction, label = lan_gene)) +
    geom_gene_arrow(arrowhead_height = unit(6, "mm"),
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(6, "mm")) +
    geom_text_repel(aes(x = (start + stop) / 2, label = lan_gene),
                    size = 5,
                    nudge_y = 0.3) + 
    geom_subgene_arrow(data = rodeo_lan_annotation, aes(xmin = start, xmax = stop, y = contig, xsubmin = rodeo_start, 
                                                        xsubmax = rodeo_stop, forward = direction), 
                       fill= "purple", alpha=.7,
                       arrowhead_height = unit(6, "mm"),
                       arrowhead_width = unit(1, "mm"),
                       arrow_body_height = unit(6, "mm")) +
    geom_text_repel(data = rodeo_lan_annotation, aes(x = (rodeo_start + rodeo_stop) / 2, label = rodeo_lanA_gene),
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
          panel.grid.major.y = element_line(linewidth = .75, color = "black"),
          axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"), face = "bold"),
          strip.background = element_blank(),
          strip.text.y=element_text(size=15, angle = 0, face = "bold")) 
  
  
  ggsave(paste0("../plots/samples/", i, ".pdf"), limitsize = FALSE, width = 25, height = plot_height)
  
}

rm(plot_height, partial_lan_annotation, i)





 