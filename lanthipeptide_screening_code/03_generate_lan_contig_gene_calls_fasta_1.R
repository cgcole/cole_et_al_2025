
# Read in all gene calls and rename the gene_callers_id to unique names for anvi-gene-contig-database
# This should be run if this is the first time running. Otherwise, use the later code.
lan_contig_gene_calls <- read_tsv("../metadata/lan_contig_gene_calls.txt", num_threads = 5) %>% 
  mutate(original_gene_callers_id = gene_callers_id) %>% 
  select(-gene_callers_id) %>% 
  rowid_to_column(var = "gene_callers_id") %>% 
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id)) 

gene_callers_id <- lan_contig_gene_calls %>% 
  select(original_gene_callers_id, gene_callers_id, contig, aa_sequence) %>% 
  mutate(subset = 1)

#write_tsv(gene_callers_id, "../metadata/gene_callers_id.txt")

rm(gene_callers_id)

# Read in previously changed gene_callers_id info if available. Ignore this if running for the first time.
lan_contig_gene_calls <- read_tsv("../metadata/lan_contig_gene_calls.txt", num_threads = 5) %>% 
  dplyr::rename(original_gene_callers_id = gene_callers_id) %>% 
  right_join(read_tsv("../metadata/gene_callers_id.txt", num_threads = 5)) %>% 
  filter(subset == 1) %>% 
  arrange(gene_callers_id) %>% 
  mutate(contig_gene_id = paste0(contig, "_gene_", gene_callers_id)) 


# Create list of contigs in this subset to specify which ones to run through bakta
lan_contig_id_subset_01 <- lan_contig_gene_calls %>%
  distinct(contig)

write_tsv(lan_contig_id_subset_01, "../metadata/lan_contig_id_subset_01.txt", col_names = FALSE)

rm(lan_contig_id_subset_01)


# Save fasta file of all lan contigs genes for GhostKOALA annotation
write.fasta(as.list(lan_contig_gene_calls$aa_sequence), lan_contig_gene_calls$contig_gene_id, "../fasta_files/lan_contig_genes.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[1:500000]), lan_contig_gene_calls$contig_gene_id[1:500000], "../fasta_files/lan_contig_genes_subset_01.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[500001:1000000]), lan_contig_gene_calls$contig_gene_id[500001:1000000], "../fasta_files/lan_contig_genes_subset_02.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[1000001:1500000]), lan_contig_gene_calls$contig_gene_id[1000001:1500000], "../fasta_files/lan_contig_genes_subset_03.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[1500001:2000000]), lan_contig_gene_calls$contig_gene_id[1500001:2000000], "../fasta_files/lan_contig_genes_subset_04.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[2000001:2500000]), lan_contig_gene_calls$contig_gene_id[2000001:2500000], "../fasta_files/lan_contig_genes_subset_05.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[2500001:3000000]), lan_contig_gene_calls$contig_gene_id[2500001:3000000], "../fasta_files/lan_contig_genes_subset_06.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[3000001:3500000]), lan_contig_gene_calls$contig_gene_id[3000001:3500000], "../fasta_files/lan_contig_genes_subset_07.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[3500001:4000000]), lan_contig_gene_calls$contig_gene_id[3500001:4000000], "../fasta_files/lan_contig_genes_subset_08.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[4000001:4500000]), lan_contig_gene_calls$contig_gene_id[4000001:4500000], "../fasta_files/lan_contig_genes_subset_09.fa")
write.fasta(as.list(lan_contig_gene_calls$aa_sequence[4500001:4733739]), lan_contig_gene_calls$contig_gene_id[4500001:4733739], "../fasta_files/lan_contig_genes_subset_10.fa")

# Save gene calls to input into anvi-gen-contig-database to make contig database for annotation in Anvio
lan_contig_annotation_gene_calls <- lan_contig_gene_calls %>% 
  select(-contig_gene_id, -original_gene_callers_id, -subset) %>% 
  relocate(gene_callers_id)

write_tsv(lan_contig_annotation_gene_calls, "../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls.txt")

rm(lan_contig_annotation_gene_calls)


