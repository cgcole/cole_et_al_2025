

# Read in all gene calls and rename the gene_callers_id to unique names for anvi-gene-contig-database
# This are additional contigs to add to previously annotated contigs
lan_contig_gene_calls_5 <- read_tsv("../metadata/lan_contig_gene_calls.txt") %>% 
  dplyr::rename(original_gene_callers_id = gene_callers_id) %>% 
  full_join(read_tsv("../metadata/gene_callers_id.txt") %>% 
              filter(subset %in% c(1,2,3,4))) %>% 
  arrange(gene_callers_id) %>% 
  rowid_to_column(var = "new_gene_callers_id") %>% 
  filter(is.na(gene_callers_id)) %>%
  mutate(gene_callers_id = new_gene_callers_id) %>% 
  select(-new_gene_callers_id) %>% 
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id)) %>% 
  mutate(subset = 5)


# Add new names to gene_callers_id list
gene_callers_id <- read_tsv("../metadata/gene_callers_id.txt") %>% 
  filter(subset != 5) %>% 
  full_join(lan_contig_gene_calls_5 %>% 
              select(original_gene_callers_id, gene_callers_id, contig, aa_sequence, subset))

#write_tsv(gene_callers_id, "../metadata/gene_callers_id.txt")

rm(gene_callers_id)

# Read in all gene calls and rename the gene_callers_id to unique names for anvi-gene-contig-database
# This are additional contigs to add to previously annotated contigs
# Run this after the previous code has already been run
lan_contig_gene_calls_5 <- read_tsv("../metadata/lan_contig_gene_calls.txt", num_threads = 5) %>%
  dplyr::rename(original_gene_callers_id = gene_callers_id) %>%
  right_join(read_tsv("../metadata/gene_callers_id.txt", num_threads = 5)) %>%
  arrange(gene_callers_id) %>%
  filter(subset == 5) %>%
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id))


# Create list of contigs in this subset to specify which ones to run through bakta
lan_contig_id_subset_05 <- lan_contig_gene_calls_5 %>% 
  distinct(contig)

write_tsv(lan_contig_id_subset_05, "../metadata/lan_contig_id_subset_05.txt", col_names = FALSE)

rm(lan_contig_id_subset_05)

# Save fasta file of all lan contigs genes for GhostKOALA annotation
write.fasta(as.list(lan_contig_gene_calls_5$aa_sequence), lan_contig_gene_calls_5$contig_gene_id, "../fasta_files/lan_contig_genes_5.fa")
write.fasta(as.list(lan_contig_gene_calls_5$aa_sequence[1:500000]), lan_contig_gene_calls_5$contig_gene_id[1:500000], "../fasta_files/lan_contig_genes_subset_14.fa")
write.fasta(as.list(lan_contig_gene_calls_5$aa_sequence[500001:518475]), lan_contig_gene_calls_5$contig_gene_id[500001:518475], "../fasta_files/lan_contig_genes_subset_15.fa")

# Save gene calls to input into anvi-gen-contig-database to make contig database for annotation in Anvio
lan_contig_annotation_gene_calls_5 <- lan_contig_gene_calls_5 %>% 
  select(-contig_gene_id, -original_gene_callers_id, -subset)%>% 
  relocate(gene_callers_id)

write_tsv(lan_contig_annotation_gene_calls_5, "../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_5.txt")

rm(lan_contig_annotation_gene_calls_5)
