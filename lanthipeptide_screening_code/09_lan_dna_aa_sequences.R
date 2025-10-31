

library(Biostrings)
library(stringi)
library(seqinr)


confirmed_lanA_genes <- gene_calls %>%
  filter(include == "yes") 

# This is to make a fasta file of contigs with a confirmed lanA gene to run on blast and get species level information
lan_contig_fasta_file <- readDNAStringSet("../fasta_files/lan_contigs.fa")

lanA_contig_dna_sequences <- tibble(
  contig = names(lan_contig_fasta_file),
  dna_sequence = as.character(lan_contig_fasta_file)) %>% 
  right_join(confirmed_lanA_genes %>% 
              select(gene_callers_id, contig, start, stop, direction, aa_sequence, aa_length, lan_gene)) %>% 
  distinct(contig, .keep_all = TRUE)

write.fasta(as.list(lanA_contig_dna_sequences$dna_sequence), lanA_contig_dna_sequences$contig, "../fasta_files/lanA_full_contigs.fa")




lanA_blastn_results_01 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_01.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) 

lanA_blastn_results_02 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_02.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) 

lanA_blastn_results_03 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_03.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))


lanA_blastn_results_04 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_04.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))

lanA_blastn_results_05 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_05.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))

lanA_blastn_results_06 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_06.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))

lanA_blastn_results_07 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_07.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))


lanA_blastn_results_08 <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_08.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp"))



lanA_blastn_new_contigs <- lanA_contig_dna_sequences %>% 
  filter(!(contig %in% lanA_blastn_results_01$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_02$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_03$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_04$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_05$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_06$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_07$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_08$qseqid)) %>% 
  filter(!(contig %in% lanA_blastn_results_09$qseqid))

write.fasta(as.list(lanA_blastn_new_contigs$dna_sequence), lanA_blastn_new_contigs$contig, "../fasta_files/lanA_full_contigs_10.fa")




# This is to get the dna sequence for the lanA genes specifically
reverse_lan_contig_fasta_file <- complement(readDNAStringSet("../fasta_files/lan_contigs.fa"))

lanA_dna_sequences <- tibble(
  contig = names(lan_contig_fasta_file),
  dna_sequence = as.character(lan_contig_fasta_file)) %>% 
  full_join(tibble(
    contig = names(reverse_lan_contig_fasta_file),
    reverse_dna_sequence = as.character(reverse_lan_contig_fasta_file))) %>% 
  right_join(confirmed_lanA_genes %>% 
               select(gene_callers_id, contig, start, stop, direction, aa_sequence, aa_length, lan_gene)) %>% 
  mutate(lanA_dna_sequence = str_sub(dna_sequence, start = start, end = stop)) %>% 
  mutate(lanA_reverse_dna_sequence = str_sub(reverse_dna_sequence, start = -stop, end = -start)) %>% 
  mutate(lanA_dna_sequence = if_else(direction == 0, lanA_reverse_dna_sequence, lanA_dna_sequence)) %>% 
  select(lan_gene, lanA_dna_sequence, aa_sequence) %>% 
  distinct() %>% 
  mutate(lan_number = as.numeric(str_replace(lan_gene, "lanA_", ""))) %>% 
  arrange(lan_number) %>% 
  select(-lan_number) %>% 
  group_by(lan_gene) %>% 
  mutate(lan_gene = paste0(lan_gene, "_seq", row_number())) %>% 
  ungroup()

write.fasta(as.list(lanA_dna_sequences$lanA_dna_sequence), lanA_dna_sequences$lan_gene, "../fasta_files/lanA_dna_sequences.fa")



# Generate a fasta file of lanA AA sequences including ones not detected for comparison
lanA_aa_sequences <- lanA_genes %>% 
  filter(predicted_lanA == "yes") %>%
  select(aa_sequence) %>% 
  full_join(published_lanA_genes %>% 
              select(aa_sequence)) %>% 
  full_join(interpro_lanA_genes) %>% 
  full_join(refseq_lanA_genes) %>% 
  full_join(walker_lanA_genes, relationship = "many-to-many") %>% 
  distinct() %>% 
  left_join(lanA_class %>% 
              select(lan_gene, aa_sequence)) %>% 
  mutate(core = str_sub(aa_sequence, -18),
         core_label = paste0(lan_gene, "_core")) 

write.fasta(as.list(lanA_aa_sequences$aa_sequence), lanA_aa_sequences$lan_gene, "../fasta_files/lanA_aa_sequences.fa")
write.fasta(as.list(lanA_aa_sequences$core), lanA_aa_sequences$core_label, "../fasta_files/lanA_core_aa_sequences.fa")
























#Using blastn on non lanA genes
confirmed_lan_genes <- gene_calls %>%
  filter(str_detect(annotated_lan_genes, "lan")) %>% 
  filter(!str_detect(all_lan_genes, "lanA"))

# This is to make a fasta file of contigs with a confirmed lanA gene to run on blast and get species level information
lan_contig_fasta_file <- readDNAStringSet("../fasta_files/lan_contigs.fa")

lan_contig_dna_sequences <- tibble(
  contig = names(lan_contig_fasta_file),
  dna_sequence = as.character(lan_contig_fasta_file)) %>% 
  right_join(confirmed_lan_genes %>% 
               select(gene_callers_id, contig, start, stop, direction, aa_sequence, aa_length, lan_gene)) %>% 
  distinct(contig, .keep_all = TRUE)

write.fasta(as.list(lan_contig_dna_sequences$dna_sequence), lan_contig_dna_sequences$contig, "../fasta_files/lan_gene_full_contigs.fa")


lan_blastn_results_01 <- read_tsv("../lan_contig_annotations/lan_contig_blast_taxonomy_annotations/lan_contig_blast_results_01.txt", col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                                                                                                                                      "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                                                                                                                                      "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                                                                                                                                      "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                                                                                                                                      "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) 

                                                                                                                                                   
                                                                                                                                                   
lan_blastn_new_contigs <- lan_contig_dna_sequences %>% 
  filter(!(contig %in% lanA_blastn_results_01$qseqid))

write.fasta(as.list(lan_blastn_new_contigs$dna_sequence), lan_blastn_new_contigs$contig, "../fasta_files/lan_gene_full_contigs_02.fa")












