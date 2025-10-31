library(readxl)
library(seqinr)
library(tidyverse)


# Read in all RODEO screening results
rodeo_table <- read_xlsx("../data/240201_metaRODEO_lanthis_Pamer.xlsx") %>% 
  type_convert() %>% 
  mutate(ID = str_replace(Accession_id, ".*DFIclinical_contigs/", "")) %>% 
  mutate(ID = str_replace(ID, "\\.megahit.contigs.fasta", "")) %>% 
  mutate(fasta_file = str_extract(Accession_id, "(?<=/).*")) %>% 
  mutate(rodeo_aa_sequence = str_extract(FASTA, "(?<=\\n).*")) %>% 
  mutate(locus_id = str_replace_all(ID, "-", "_")) %>% 
  mutate(locus_id = str_replace_all(locus_id, "\\.", "_")) %>% 
  filter(`SRD calls` == 1) %>%
  group_by(FASTA) %>% 
  mutate(rodeo_start = min(Start, End)) %>% 
  mutate(rodeo_stop = max(Start, End)) %>% 
  mutate(rodeo_direction = if_else(Start < End, 1, 
                             if_else(Start > End, 0, NA))) %>% 
  ungroup() %>% 
  mutate(rodeo_start = rodeo_start + 1) %>% #starts at 0 and you can't have a 0th nucleotide
  select(ID, locus_id, Leader, Core, Start, End, LEN, `Class Call`, 
         rodeo_aa_sequence, fasta_file, rodeo_start, rodeo_stop, rodeo_direction) 


# Identify unique lanA genes from the RODEO results
rodeo_lanA_genes <- rodeo_table %>% 
  distinct(rodeo_aa_sequence) 

write_tsv(rodeo_lanA_genes, "../metadata/rodeo_lanA_genes.txt", col_names = FALSE)

