install.packages("tidyverse")
library(tidyverse)

lanA_contig_annotations <- read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_01.txt", 
        col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
           "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
           "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
           "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
           "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>% 
  mutate(staxids = as.character(staxids)) %>% 
  mutate(sallgi = as.character(sallgi)) %>% 
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_02.txt", 
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send", 
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", 
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids", 
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", 
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_03.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_04.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_05.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_06.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_07.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_08.txt",
                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
              mutate(staxids = as.character(staxids)) %>%
              mutate(sallgi = as.character(sallgi))) %>%
#  full_join(read_tsv("../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_09.txt",
#                     col_names = c("qseqid", "sseqid", "qacc", "sacc", "qlen", "slen", "qstart", "qend", "sstart", "send",
#                                   "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
#                                   "gapopen", "gaps", "positive", "ppos", "frames", "qframe", "sframe", "btop", "staxids",
#                                   "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles",
#                                   "sallseqid", "sallgi", "sallacc", "qcovs", "qcovhsp")) %>%
#              mutate(staxids = as.character(staxids)) %>%
#              mutate(sallgi = as.character(sallgi))) %>%
  distinct() %>% 
  mutate(meets_threshold = if_else(evalue < 1e-3 & pident > 75 & qcovs > 60, "yes", "no")) %>% 
  filter(meets_threshold == "yes") %>% 
  group_by(qseqid) %>% 
  slice_max(bitscore, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(clean_title = str_replace(stitle, "MAG: ", "")) %>% 
  mutate(clean_title = str_replace(clean_title, "MAG TPA_asm: ", "")) %>% 
  mutate(clean_title = str_replace(clean_title, "uncultured ", "")) %>% 
  mutate(clean_title = str_replace(clean_title, "Uncultured human fecal virus clone ", "")) %>%
  mutate(clean_title = if_else(str_detect(clean_title, "R.gnavus"), "Mediterraneibacter gnavus", clean_title)) %>% 
  mutate(clean_title = if_else(str_detect(clean_title, "BlautiaA"), "Blautia sp.", clean_title)) %>% 
  mutate(clean_title = if_else(str_detect(clean_title, "B.wexlerae"), "Blautia wexlerae", clean_title)) %>% 
  mutate(clean_title = if_else(str_detect(clean_title, "D.formicigenerans"), "Dorea formicigenerans", clean_title)) %>% 
  mutate(clean_title = if_else(str_detect(clean_title, "A.hadrus"), "Anaerostipes hadrus", clean_title)) %>% 
  mutate(blast_genus_species = word(clean_title, 1, 2)) %>% 
  mutate(blast_genus_species = str_replace(blast_genus_species, "\\[Ruminococcus\\] torques", "Mediterraneibacter torques"),
         blast_genus_species = str_replace(blast_genus_species, "Ruminococcus torques", "Mediterraneibacter torques"),
         blast_genus_species = str_replace(blast_genus_species, "Ruminococcus gnavus", "Mediterraneibacter gnavus"),
         blast_genus_species = str_replace(blast_genus_species, "Ruminococcus obeum", "Blautia obeum"),
         blast_genus_species = str_replace(blast_genus_species, "\\[Clostridium\\] hylemonae", "Lachnoclostridium  hylemonae"),
         blast_genus_species = str_replace(blast_genus_species, "\\[Clostridium\\] innocuum", "Erysipelatoclostridium innocuum"),
         blast_genus_species = str_replace(blast_genus_species, "\\[Clostridium\\] scindens", "Lachnoclostridium scindens"),
         blast_genus_species = str_replace(blast_genus_species, "\\[Eubacterium\\] hallii", "Anaerobutyricum hallii")) %>% 
  mutate(blast_genus = word(blast_genus_species, 1, 1)) %>% 
  mutate(blast_species = word(blast_genus_species, 2, 2)) %>% 
  select(qseqid, staxids, stitle, blast_genus_species, blast_genus, blast_species, clean_title) %>% 
  dplyr::rename(contig = qseqid) %>% 
  right_join(read_csv("../data/gene_call_table.csv") %>% 
               select(locus_id, contig, lan_gene, predicted_class, all_lan_genes, highest_tax_hit, include, aa_sequence) %>% 
               filter(include == "yes")) %>% 
  mutate(blast_genus_species = if_else(is.na(blast_genus_species), "Unclassified", blast_genus_species)) %>% 
  mutate(blast_genus = if_else(is.na(blast_genus), "Unclassified", blast_genus)) %>% 
  mutate(blast_species = if_else(is.na(blast_species), "Unclassified", blast_species))
  
write_csv(lanA_contig_annotations, "../data/lanA_gene_calls.csv")

