

# Read in all gene calls and rename the gene_callers_id to unique names for anvi-gene-contig-database
# This are additional contigs to add to previously annotated contigs
lan_contig_gene_calls_2 <- read_tsv("../metadata/lan_contig_gene_calls.txt") %>% 
  dplyr::rename(original_gene_callers_id = gene_callers_id) %>% 
  left_join(read_tsv("../metadata/gene_callers_id.txt")) %>% 
  arrange(gene_callers_id) %>% 
  rowid_to_column(var = "new_gene_callers_id") %>% 
  filter(is.na(gene_callers_id)) %>% 
  mutate(gene_callers_id = new_gene_callers_id) %>% 
  select(-new_gene_callers_id) %>% 
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id)) 


# New code now that this code had been run previously
lan_contig_gene_calls_2 <- read_tsv("../metadata/lan_contig_gene_calls.txt", num_threads = 5) %>% 
  dplyr::rename(original_gene_callers_id = gene_callers_id) %>% 
  right_join(read_tsv("../metadata/gene_callers_id.txt", num_threads = 5)) %>% 
  filter(subset == 2) %>% 
  arrange(gene_callers_id) %>% 
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id)) 

# Create list of contigs in this subset to specify which ones to run through bakta
lan_contig_id_subset_02 <- lan_contig_gene_calls_2 %>% 
  distinct(contig)

write_tsv(lan_contig_id_subset_02, "../metadata/lan_contig_id_subset_02.txt", col_names = FALSE)

rm(lan_contig_id_subset_02)

# Save fasta file of all lan contigs genes for GhostKOALA annotation
write.fasta(as.list(lan_contig_gene_calls_2$aa_sequence), lan_contig_gene_calls_2$contig_gene_id, "../fasta_files/lan_contig_genes_2.fa")
write.fasta(as.list(lan_contig_gene_calls_2$aa_sequence[1:301328]), lan_contig_gene_calls_2$contig_gene_id[1:301328], "../fasta_files/lan_contig_genes_subset_11.fa")

# Save gene calls to input into anvi-gen-contig-database to make contig database for annotation in Anvio
lan_contig_annotation_gene_calls_2 <- lan_contig_gene_calls_2 %>% 
  select(-contig_gene_id, -original_gene_callers_id, -subset)%>% 
  relocate(gene_callers_id)

write_tsv(lan_contig_annotation_gene_calls_2, "../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_2.txt")

rm(lan_contig_annotation_gene_calls_2)


