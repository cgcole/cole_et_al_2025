
# Lantibiotic genes pulled from interpro database
interpro_lanA_genes <- read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR006079.fasta", seqtype = "AA", as.string = TRUE) %>% 
  tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
  select(-1) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR012519.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR027632.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR031031.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR046016.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR048275.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR007682.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  full_join(read.fasta("../fasta_files/interpro_lanA_genes/protein-matching-IPR029243.fasta", seqtype = "AA", as.string = TRUE) %>% 
              tibble(`Locus Tag` = names(.), aa_sequence = unlist(.)) %>% 
              select(-1)) %>% 
  select(aa_sequence) %>% 
  distinct() 


# KEGG lanA genes for K20482
kegg_lanA_genes <- read_tsv("../metadata/kegg_lanA_genes.txt", col_names = "aa_sequence") %>% 
  distinct()


# These are lanA genes pulled from refseq that were used in Jerry's paper
refseq_lanA_genes_file <- read.fasta("../fasta_files/ED.File.1.lanA_in_refseq.fasta", as.string = TRUE)

refseq_lanA_genes <- as.data.frame(getSequence(refseq_lanA_genes_file, as.string = TRUE)) %>% 
  pivot_longer(1:128, names_to = "name", values_to = "aa_sequence") %>% 
  select(-name) %>% 
  mutate(aa_sequence = str_to_upper(aa_sequence))

rm(refseq_lanA_genes_file)


# LanA genes that were identified using something like RODEO on our biobank samples (Eric L. and Jerry)
biobank_lanA_genes <- read_tsv("../data/all_hits_per_family_per_nisin_blauticin.txt") %>% 
  dplyr::rename(aa_sequence = peptide) %>% 
  full_join(read_tsv("../metadata/biobank_lanA_genes.txt", col_names = "aa_sequence")) %>% 
  distinct(aa_sequence) 



# LanA genes identified through annotations 
# Bakta annotations are separated due to the different gene calling it contains which produces some that don't line up with Anvi'o gene calls
bakta_lanA_functions <- read_tsv("../metadata/bakta_lanA_functions.txt", col_names = c("bakta_function"))


bakta_annotated_lanA_genes <- combined_gene_annotations %>% 
  filter(bakta_function %in% bakta_lanA_functions$bakta_function) %>% 
  select(bakta_aa_sequence) %>%
  dplyr::rename(aa_sequence = bakta_aa_sequence) %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  filter(aa_length < 150) %>% 
  distinct(aa_sequence)



annotated_lanA_genes <- combined_gene_annotations %>% 
  filter(bakta_function %in% bakta_lanA_functions$bakta_function |
          # KOfam_accession == "K20482" | KOfam_accession == "K20383" | KOfam_accession == "K20384" | KOfam_accession == "K20383,K20384" |
           ghostkoala_accession == "K20482" | ghostkoala_accession == "K20383" | ghostkoala_accession == "K20384" | ghostkoala_accession == "K20383,K20384" |
           str_detect(Pfam_accession, "PF14867") | str_detect(Pfam_accession, "PF02052") | str_detect(Pfam_accession, "PF16934") | 
           str_detect(Pfam_accession, "PF04604") | str_detect(Pfam_accession, "PF08130") | str_detect(Pfam_accession, "PF19402")) %>%
  full_join(bakta_annotated_lanA_genes) %>% 
  select(aa_sequence) %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  filter(aa_length < 150) %>% 
  distinct(aa_sequence)



# Manually called lanA genes
manual_called_lanA_genes <- read_tsv("../metadata/manual_called_lanA_genes.txt", col_names = "aa_sequence")


# lanA genes identifed in walker et al from the van der Donk lab
walker_lanA_genes <- read_tsv("../metadata/walker_et_al_lanA_genes.txt") %>% 
  select(aa_sequence)


#These are lanA genes that I have curated myself from previously published papers
published_lanA_genes <- read_tsv("../metadata/known_lan_genes2.txt") %>%
  mutate(leader = replace_na(leader, "")) %>% 
  mutate(aa_sequence = paste0(leader, core))


#lanA genes from NCBI
blast_lanA_genes <- readAAStringSet("../fasta_files/ncbi_lanA_genes/ncbi_lanA_aa_genes.fa") %>% 
  {tibble(stitle = names(.), aa_sequence = as.character(.)) %>% 
      mutate(qseqid = word(stitle, start = 1, end = 1, sep = fixed(" "))) %>% 
      mutate(clean_title = str_replace(stitle, " \\[.*\\]", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "MULTISPECIES: ", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "LOW QUALITY PROTEIN: ", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "RecName: Full=", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "MAG: ", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "UNVERIFIED_ORG: ", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "MAG TPA: ", "")) %>% 
      mutate(clean_title = str_replace(clean_title, "; Flags: Precursor", "")) %>% 
      mutate(clean_title = str_extract(clean_title, "(?<=\\s).*")) %>% 
      filter(!str_detect(clean_title, "hypothetical protein")) %>% 
      mutate(potential_lanA = str_detect(clean_title, 
                                         "lantipeptide|lanthipeptide|lantibiotic|columbicin|mutA'|nisA|nisin|pldA|salivaricin|sboA|gallidermin|labyrinthopeptin|lactocin|lanthonine-containing|lichenicidin|nsuA|pneA")) %>% 
      filter(potential_lanA == TRUE) %>% 
      filter(!str_detect(clean_title, "transport|immunity|resistance|response|biosynthesis|protease|dehydratase|kinase|modifying|synthetase|efflux|modification|production|exporter|dehydrogenase|permease|cyclase|membrane|protection|regulator|processing|decarboxylase")) %>% 
      select(aa_sequence) %>% 
      distinct()}



# Addition of adjusted lanA genes that contain methionine as the start amino acid by both replaceing
# and by deleting amino acids up to the first methionine
adjusted_lanA_genes <- interpro_lanA_genes %>% 
  full_join(kegg_lanA_genes) %>% 
  full_join(refseq_lanA_genes) %>% 
  full_join(biobank_lanA_genes) %>% 
  full_join(annotated_lanA_genes) %>%
  full_join(manual_called_lanA_genes) %>% 
  full_join(rodeo_lanA_genes %>% 
              dplyr::rename(aa_sequence = rodeo_aa_sequence)) %>% 
  full_join(walker_lanA_genes) %>% 
  full_join(published_lanA_genes %>% 
              select(aa_sequence), relationship = "many-to-many") %>% 
  full_join(blast_lanA_genes) %>% 
  distinct() %>% 
  mutate(start_aa = str_sub(aa_sequence, 1, 1)) %>% 
  filter(start_aa != "M") %>% 
  mutate(aa_sequence_1 = str_replace(aa_sequence, "^.{1}","M")) %>% 
  mutate(aa_sequence_2 = str_extract(aa_sequence, "M.*")) %>% 
  select(-start_aa) %>% 
  pivot_longer(everything(), names_to = "id", values_to = "aa_sequence") %>% 
  select(aa_sequence) %>% 
  filter(!is.na(aa_sequence))


# Combine all lanA genes and remove genes already found by RODEO
alternate_lanA_genes <- interpro_lanA_genes %>% 
  full_join(kegg_lanA_genes) %>% 
  full_join(refseq_lanA_genes) %>% 
  full_join(biobank_lanA_genes) %>% 
  full_join(annotated_lanA_genes) %>%
  full_join(manual_called_lanA_genes) %>% 
  full_join(walker_lanA_genes) %>% 
  full_join(published_lanA_genes %>% 
              select(aa_sequence), relationship = "many-to-many") %>% 
  full_join(blast_lanA_genes) %>% 
  full_join(adjusted_lanA_genes, relationship = "many-to-many") %>% 
  distinct() %>% 
  filter(!(aa_sequence %in% rodeo_lanA_genes$rodeo_aa_sequence)) %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  filter(aa_length >= 17) %>% 
  select(-aa_length)


write_tsv(alternate_lanA_genes, "../metadata/alternate_lanA_genes.txt", col_names = FALSE)


# Previously published lanA genes without the leader sequence
known_lanA_genes <- read_tsv("../metadata/known_lan_genes.txt", col_names = "aa_sequence") %>% 
  distinct()

#Blast matches to peptides given to Che to run on the clinical samples. These may be lanA genes. 
lanA_blast_search <- read_tsv("../metadata/lanA_blast_sequences.txt", col_names = "aa_sequence") %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  filter(aa_length >= 17) 

write_tsv(lanA_blast_search, "../metadata/lanA_blast_sequences.txt")


core_lanA_sequences <- alternate_lanA_genes %>% 
  full_join(known_lanA_genes) %>% 
  full_join(rodeo_lanA_genes %>% 
              rename(aa_sequence = rodeo_aa_sequence)) %>% 
  full_join(read_tsv("../metadata/blastp_lanA_genes.txt", col_names = "aa_sequence")) %>% 
  mutate(aa_sequence = str_sub(aa_sequence, -10)) %>% 
  distinct()

write_tsv(core_lanA_sequences, "../metadata/core_lanA_genes.txt", col_names = FALSE)


leader_lanA_sequences <- alternate_lanA_genes %>% 
  full_join(known_lanA_genes) %>% 
  full_join(rodeo_lanA_genes %>% 
              rename(aa_sequence = rodeo_aa_sequence)) %>% 
  full_join(read_tsv("../metadata/blastp_lanA_genes.txt", col_names = "aa_sequence")) %>% 
  mutate(aa_sequence = str_sub(aa_sequence, 1, 10)) %>% 
  distinct()

write_tsv(leader_lanA_sequences, "../metadata/leader_lanA_genes.txt", col_names = FALSE)
  
