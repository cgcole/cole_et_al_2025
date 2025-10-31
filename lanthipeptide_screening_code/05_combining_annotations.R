
# Anvio gene annotations
anvio_gene_annotations_01 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))

anvio_gene_annotations_02 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_2.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_2) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))

anvio_gene_annotations_03 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_3.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_3) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))


anvio_gene_annotations_04 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_4.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_4) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))


anvio_gene_annotations_05 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_5.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_5) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))


anvio_gene_annotations_06 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_6.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_6) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))

anvio_gene_annotations_07 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_7.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_7) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))


anvio_gene_annotations_08 <- read_tsv("../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_8.txt", num_threads = 5) %>% 
  distinct() %>% 
  mutate(accession = if_else(source == "Pfam", str_extract(accession, "PF.*(?=\\.)"), accession)) %>% 
  arrange(accession) %>% 
  group_by(gene_callers_id, source) %>% 
  mutate(accession = paste(accession, collapse = ","),
         `function` = paste(`function`, collapse = ","),
         e_value = paste(e_value, collapse = ",")) %>%
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(id_cols = gene_callers_id, 
              names_from = source, 
              values_from = c("accession", "function", "e_value"),
              names_sep = "_",
              values_fill = list(values = NA),
              names_glue = "{source}_{.value}") %>% 
  right_join(lan_contig_gene_calls_8) %>%
  mutate(start = start + 1) %>% #anvio starts at 1 lower but ends correctly for whatever reason
  mutate(direction = if_else(direction == "f", 1,
                             if_else(direction == "r", 0, NA)))



anvio_gene_annotations <- anvio_gene_annotations_01 %>% 
  full_join(anvio_gene_annotations_02) %>% 
  full_join(anvio_gene_annotations_03) %>% 
  full_join(anvio_gene_annotations_04) %>% 
  full_join(anvio_gene_annotations_05) %>% 
  full_join(anvio_gene_annotations_06) %>% 
  full_join(anvio_gene_annotations_07) %>% 
  full_join(anvio_gene_annotations_08)


# Read in GhostKOALA Annotations
ghostkoala_annotations <- read_tsv("../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_annotations.txt", num_threads = 5) %>%
  fill(contig_gene_id) %>% 
  filter(!is.na(ghostkoala_accession)) %>% 
  select(-second_best, -second_best_score) %>% 
  arrange(ghostkoala_accession) %>% 
  group_by(contig_gene_id) %>% 
  mutate(ghostkoala_accession = paste(ghostkoala_accession, collapse = ", "),
         ghostkoala_function = paste(ghostkoala_function, collapse = ", "),
         ghostkoala_score = paste(ghostkoala_score, collapse = ", ")) %>% 
  ungroup() %>% 
  distinct()

#Read in GhostKOALA gene taxonomy
ghostkoala_taxonomy <- read_tsv("../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_taxonomy.txt", num_threads = 5) %>% 
  mutate(contig_gene_id = str_replace(contig_gene_id, "^user:", "")) 

contig_taxonomy <- ghostkoala_taxonomy %>% 
  mutate(contig = str_replace(contig_gene_id, "_gene_.*", "")) %>% 
  filter(tax_score >= 100) %>% 
  group_by(contig, kingdom, phylum, genus) %>% 
  summarise(n = n()) %>% 
  filter(n >= 2) %>% 
  group_by(contig) %>% 
  arrange(desc(n)) %>% 
  summarise(tax_hits = paste0("k_", kingdom, "_p_", phylum, "_g_", genus, "_", n, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(highest_tax_hit = word(tax_hits, 1, sep = ",")) %>% 
  mutate(highest_tax_hit = str_replace(highest_tax_hit, "_[:digit:]+$", ""))


# Read in Bakta Annotations
bakta_fasta_file_01 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_01/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_01 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_01/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_01) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_01/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_01) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()

bakta_fasta_file_02 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_02/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_02 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_02/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_02) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_02/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_02) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_fasta_file_03 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_03/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_03 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_03/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_03) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_03/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_03) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_fasta_file_04 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_04/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_04 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_04/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_04) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_04/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_04) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()



bakta_fasta_file_05 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_05/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_05 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_05/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_05) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_05/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_05) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_fasta_file_06 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_06/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_06 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_06/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_06) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_06/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_06) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()



bakta_fasta_file_07 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_07/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_07 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_07/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_07) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_07/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_07) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_fasta_file_08 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_08/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_08 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_08/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_08) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_08/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_08) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_fasta_file_09 <- read.fasta("../lan_contig_annotations/bakta_annotations/bakta_annotations_09/lan_contigs-bakta.faa", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1)

bakta_annotations_09 <- read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_09/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% 
  dplyr::rename(bakta_accession = Gene,
                bakta_function = Product,
                contig = `#Sequence Id`,
                bakta_start = Start,
                bakta_stop = Stop,
                direction = Strand) %>% 
  right_join(bakta_fasta_file_09) %>% 
  dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
  mutate(direction = if_else(direction == "+", 1,
                             if_else(direction == "-", 0, NA))) %>% 
  select(-Type, -DbXrefs, -`Locus Tag`) %>% 
  left_join(anvio_gene_annotations %>% 
              select(contig, start, stop, direction, aa_sequence), relationship = "many-to-many") %>% 
  mutate(match = if_else((start >= bakta_start & stop <= bakta_stop) | start == bakta_start | stop == bakta_stop, "Yes", "No")) %>% 
  filter(match == "Yes") %>% 
  full_join(read_tsv("../lan_contig_annotations/bakta_annotations/bakta_annotations_09/lan_contigs-bakta.tsv", skip = 5, num_threads = 5) %>% #This brings back bakta gene calls not found by Anvi'o
              dplyr::rename(bakta_accession = Gene,
                            bakta_function = Product,
                            contig = `#Sequence Id`,
                            bakta_start = Start,
                            bakta_stop = Stop,
                            direction = Strand) %>% 
              right_join(bakta_fasta_file_09) %>% 
              dplyr::rename(bakta_aa_sequence = aa_sequence) %>% 
              mutate(direction = if_else(direction == "+", 1,
                                         if_else(direction == "-", 0, NA))) %>% 
              select(-Type, -DbXrefs, -`Locus Tag`)) %>% 
  select(-match) %>% 
  distinct()


bakta_annotations <- bakta_annotations_01 %>% 
  full_join(bakta_annotations_02) %>% 
  full_join(bakta_annotations_03) %>% 
  full_join(bakta_annotations_04) %>% 
  full_join(bakta_annotations_05) %>% 
  full_join(bakta_annotations_06) %>% 
  full_join(bakta_annotations_07) %>% 
  full_join(bakta_annotations_08) %>% 
  full_join(bakta_annotations_09) %>% 
  mutate(start = if_else(is.na(start), bakta_start, start)) %>% 
  mutate(stop = if_else(is.na(stop), bakta_stop, stop)) %>%
  mutate(aa_sequence = if_else(is.na(aa_sequence), bakta_aa_sequence, aa_sequence))


# Combine all gene calls
combined_gene_annotations <- anvio_gene_annotations %>% 
  left_join(ghostkoala_annotations) %>% 
  full_join(bakta_annotations) %>% 
  left_join(contig_taxonomy) %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  mutate(highest_tax_hit = if_else(is.na(highest_tax_hit), "Unknown", highest_tax_hit))

write_tsv(combined_gene_annotations, "../data/combined_gene_annotations.txt")

rm(bakta_fasta_file_01, bakta_fasta_file_02, bakta_annotations_01, bakta_annotations_02, 
   anvio_gene_annotations_01, anvio_gene_annotations_02, anvio_gene_annotations_03, anvio_gene_annotations_04, anvio_gene_annotations_05,
   bakta_fasta_file_03, bakta_annotations_03, bakta_fasta_file_04, bakta_annotations_04, bakta_annotations_05, bakta_annotations_06,
   bakta_fasta_file_06, bakta_fasta_file_05, bakta_annotations_07, bakta_fasta_file_07, anvio_gene_annotations_06, anvio_gene_annotations_07,
   anvio_gene_annotations, bakta_annotations_08, bakta_fasta_file_08, bakta_annotations)











