
blast_database_lanA_annotations <- readAAStringSet("../fasta_files/ncbi_lanA_genes/ncbi_lanA_aa_genes.fa") %>% 
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
      filter(!str_detect(clean_title, "transport|immunity|resistance|response|biosynthesis|protease|dehydratase|kinase|modifying|synthetase|efflux|modification|production|exporter|dehydrogenase|permease|cyclase|membrane|protection|regulator|processing|decarboxylase|export|peptidase|peptidase")) %>% 
      select(aa_sequence, clean_title) %>% 
      distinct()}

# Get annotations and check
blast_lanA_annotations <- read_tsv("../lan_contig_annotations/blast_lanA_annotations/combined_lanA_gene_blastp_results.txt", num_threads = 5) %>% 
  select(qseqid, stitle) %>% 
  rename(lan_gene = "qseqid",
         blast_function = "stitle") %>% 
  distinct() 



lanA_annotations <- lanA_genes %>% 
  left_join(combined_gene_annotations %>% 
              select(aa_sequence, bakta_accession, bakta_function, Pfam_accession, Pfam_function, ghostkoala_accession, ghostkoala_function,
                     COG20_FUNCTION_accession, COG20_FUNCTION_function) %>% 
              distinct()) %>% 
  mutate(bakta_function = replace_na(bakta_function, ""),
         bakta_accession = replace_na(bakta_accession, ""),
         Pfam_function = replace_na(Pfam_function, ""),
         Pfam_accession = replace_na(Pfam_accession, ""),
         ghostkoala_function = replace_na(ghostkoala_function, ""),
         ghostkoala_accession = replace_na(ghostkoala_accession, ""),
         COG20_FUNCTION_function = replace_na(COG20_FUNCTION_function, ""),
         COG20_FUNCTION_accession = replace_na(COG20_FUNCTION_accession, "")) %>% 
  group_by(lan_gene) %>% 
  mutate(bakta_function = paste(bakta_function, collapse = ","),
         bakta_accession = paste(bakta_accession, collapse = ","),
         Pfam_function = paste(Pfam_function, collapse = ","),
         Pfam_accession = paste(Pfam_accession, collapse = ","),
         ghostkoala_function = paste(ghostkoala_function, collapse = ","),
         ghostkoala_accession = paste(ghostkoala_accession, collapse = ","),
         COG20_FUNCTION_function = paste(COG20_FUNCTION_function, collapse = ","),
         COG20_FUNCTION_accession = paste(COG20_FUNCTION_accession, collapse = ",")) %>% 
  ungroup() %>% 
  distinct() %>% 
  left_join(blast_lanA_annotations, relationship = "many-to-many")   %>% 
  group_by(lan_gene) %>% 
  mutate(blast_function = paste(blast_function, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(check = str_detect(blast_function, "lantipeptide|lanthipeptide|lantibiotic|columbicin|mutA'|nisA|nisin|pldA|salivaricin|sboA|gallidermin|labyrinthopeptin|lactocin|lanthonine-containing|lichenicidin|nsuA|pneA")) %>% 
  left_join(blast_database_lanA_annotations, relationship = "many-to-many") %>% 
  mutate(check = if_else(!is.na(clean_title), TRUE, check)) %>% 
  arrange(desc(check)) %>% 
  distinct() 


write_csv(lanA_annotations, "../metadata/lanA_annotation_check.csv")


#Update the lanA_class file
processed_lanA_annotations <- read_csv("../metadata/lanA_annotation_check.csv") %>% 
  select(lan_gene, aa_sequence, class, predicted_class, include, predicted_lanA, detected) %>% 
  distinct() %>% 
  mutate(order = as.numeric(str_replace(lan_gene, "lanA_", ""))) %>% 
  arrange(order)%>% 
  select(-order)

write_csv(processed_lanA_annotations, "../metadata/lanA_class.csv")




blastp_lanA_annotations <- lanA_genes %>% 
  left_join(blast_lanA_annotations, relationship = "many-to-many")  %>% 
  filter(predicted_lanA == "yes") %>% 
  filter(!is.na(blast_function), !str_detect(blast_function, "hypothetical protein")) %>% 
  distinct(blast_function) %>% 
  mutate(blast_function = str_replace(blast_function, " \\[.*\\]", "")) %>% 
  mutate(blast_function = str_replace(blast_function, "MULTISPECIES: ", "")) %>% 
  mutate(blast_function = str_replace(blast_function, "LOW QUALITY PROTEIN: ", "")) %>% 
  mutate(blast_function = str_replace(blast_function, "RecName: Full=", "")) %>% 
  mutate(blast_function = str_replace(blast_function, "MAG: ", "")) %>% 
  mutate(blast_function = str_replace(blast_function, "; Flags: Precursor", "")) %>% 
  distinct() %>% 
  arrange(blast_function)

#write_tsv(blastp_lanA_annotations, "../metadata/blastp_lanA_annotations.txt", col_names = FALSE)

#blastp_lanA_annotations <- read_tsv("../metadata/blastp_lanA_annotations.txt", col_names = "blast_function") %>% 
#  mutate

lantibiotic, plantaricin C family
















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
      filter(!str_detect(clean_title, "transport|immunity|resistance|response|biosynthesis|protease|dehydratase|kinase|modifying|synthetase|efflux|modification|production|exporter|dehydrogenase|permease|cyclase|membrane|protection|regulator|processing|decarboxylase"))}

unique_blast_annotations <- blast_lanA_genes %>% 
  distinct(clean_title)

write_tsv(unique_blast_annotations, "../metadata/blast_lanA_annotations.txt", col_names = FALSE)


x <- gene_calls %>% 
  filter(aa_sequence %in% blast_lanA_genes$aa_sequence)





blast_lan_genes <- readAAStringSet("../fasta_files/ncbi_lanA_genes/ncbi_lanA_aa_genes.fa") %>% 
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
      filter(potential_lanA == TRUE)}




blast_identified_lanA_genes <- gene_calls %>% 
  left_join(blast_lan_genes %>% 
              select(aa_sequence, clean_title) %>% 
              rename(blast_annotation = "clean_title") %>% 
              distinct() %>% 
              group_by(aa_sequence) %>% 
              mutate(blast_annotation = paste0(blast_annotation, collapse = ",")) %>% 
              ungroup() %>% 
              distinct()) %>% 
  filter(!is.na(blast_annotation)) %>% 
  filter(!str_detect(annotated_lan_genes, "lan")) %>% 
  filter(blast_annotation %in% c("class I lanthipeptide", "class III lanthipeptide", "SapB/AmfS family lantipeptide")) %>% 
  distinct(aa_sequence)














blast_lan_genes <- readAAStringSet("../fasta_files/ncbi_lanA_genes/ncbi_lanA_aa_genes.fa") %>% 
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
      distinct(aa_sequence)}





#This code pertains to the blast searches that were run by Che on the lanA genes I gave him
biobank_lanA_blast_search <- read_tsv("../metadata/lanA_blast_sequences.txt", col_names = "aa_sequence") %>% 
  mutate(aa_length = str_length(aa_sequence)) %>% 
  filter(aa_length >= 17) 


