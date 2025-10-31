
#Generate list of lan associated pfam domains
lan_gene_pfams <- gene_calls %>% 
  filter(str_detect(lan_gene, "lan")) %>%
  filter(!is.na(Pfam_accession)) %>%
  distinct(Pfam_accession) %>% 
  separate_rows(Pfam_accession, sep = ",") %>% 
  distinct(Pfam_accession) %>% 
  arrange(Pfam_accession) %>% 
  mutate(Pfam_accession = str_replace(Pfam_accession, "\\..*", ""))
  
write_tsv(lan_gene_pfams, "../metadata/lan_gene_pfam_accessions.txt", 
          col_names = FALSE)

