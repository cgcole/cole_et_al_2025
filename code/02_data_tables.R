
## read sequencing data in
seq_table <- read_csv("./data/mouse_16S_sequencing.csv") 

  
## filter for bacteria and sequences with above .0001 relative abundance
seq_table_filtered <- seq_table %>% 
  select(samplename, sampleid, experiment, group, treatment, day, mouse.number, sample.type, sample.weight, reads.in, 
         asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, seq, numseqs, 
         total_16s_copies, total_16s_copy_concentration, copy_per_g, copy_per_mg) %>% 
  filter(Kingdom == "Bacteria") %>% #Remove any Archaea and Eukaryote sequences
  filter(Class != "Chloroplast") %>% #Remove Chloroplast sequences
  group_by(samplename) %>% 
  mutate(total_seqs = sum(numseqs)) %>% 
  ungroup() %>% 
  mutate(rel_abundance = numseqs / total_seqs) %>% 
  filter(rel_abundance > .0001) %>% #Filter for ASVs greater than .0001 abundance
  group_by(samplename) %>% 
  mutate(total_seqs = sum(numseqs)) %>% 
  ungroup() %>% 
  mutate(rel_abundance = numseqs / total_seqs) %>% 
  complete(nesting(samplename, sampleid, experiment, group, treatment, day, mouse.number, sample.type, sample.weight, reads.in, total_seqs, 
                   total_16s_copies, total_16s_copy_concentration, copy_per_g, copy_per_mg), 
           nesting(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, seq), fill = list(rel_abundance = 0, numseqs = 0))



## create otu table
otu_df <- seq_table_filtered %>% 
  filter(numseqs > 0) %>% 
  distinct(samplename, asv, numseqs) %>% 
  pivot_wider(names_from = asv, values_from = numseqs, values_fill = list(numseqs = 0)) %>% 
  column_to_rownames(var = "samplename")

## create taxonomy table
tax_df <- seq_table_filtered %>%  
  distinct(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%  
  column_to_rownames(var = "asv") %>%  
  as.matrix()

## create meta data
meta_df <- seq_table_filtered %>%  
  distinct(samplename, sampleid, group, treatment, experiment, day,
           sample.weight, mouse.number) %>%  
  mutate(samplename_row_id = samplename) %>% 
  column_to_rownames(var = "samplename_row_id")

## read in ASV sequences
asv_df <- seq_table_filtered %>%  
  distinct(asv, seq) %>% 
  group_by(seq) %>% 
  mutate(n = n())
asv_fa <- DNAStringSet(asv_df$seq)
names(asv_fa) <- asv_df$asv

## optional: perform alignment of ASVs so you can calculate UniFrac distance
#aln_16s_muscle <- msaMuscle(asv_fa)
#aln_16s <- msaConvert(aln_16s_muscle, type = "seqinr::alignment")
#dist <- seqinr::dist.alignment(aln_16s, "identity")
#tree <- ape::nj(dist)

# assemble the phyloseq object
phy <- phyloseq(otu_table(otu_df, taxa_are_rows = F), 
                asv_fa, 
                tax_table(tax_df),
                # phy_tree(tree), # optional: feel free to comment it out
                sample_data(meta_df))

rm(otu_df, tax_df, meta_df, asv_df, asv_fa)


## Table of all samples
sample_table <- seq_table %>% 
  distinct(samplename, experiment, sampleid, group, treatment, day, mouse.number)


## Calculate alpha diversity metrics
alpha_diversity <- estimate_richness(phy, measures = c("Observed", "Chao1", "ACE", 
                                                       "Shannon", "Simpson", "InvSimpson")) %>% 
  mutate(samplename = rownames(.)) %>% 
  as_tibble() %>% 
  left_join(sample_table)


## PFBBr quant table
pfbbr_quant <- read_csv("./data/mouse_pfbbr_quant_metabolomics.csv") 


## Bile acid quant table
bile_acid_quant <- read_csv("./data/mouse_bile_acid_quant_metabolomics.csv") 


## K. pneumoniae CFU data
kleb_cfu_t <- read_csv("./data/mouse_kleb_cfu.csv") 


## C. difficile CFU data
c_diff_cfu_t <- read_csv("./data/mouse_c_diff_cfu.csv") 


## C. difficile Symptom data
c_diff_symptom_t <- read_csv("./data/mouse_c_diff_weights.csv") %>% 
  group_by(mouse.number, group, treatment) %>% 
  mutate(weight_percentage = (mouse_weight_g/mouse_weight_g[days.post.infection == 0])*100) 


## Table of all the previously published lanthipeptides that matched to the identified lanthipeptides
published_lan_sequences <- read_csv("./data/unique_identified_lanthipeptides.csv") %>% 
  cross_join(read_csv("./data/previously_published_lanthipeptides.csv") %>% 
               select(lanthipeptide, core)) %>% 
  filter(str_detect(aa_sequence, core)) %>% 
  mutate(label = paste0(lan_gene, " (",lanthipeptide,")"))


## Clinical donor data
clinical_donor_t <- read_csv("./data/human_clinical_donor_table.csv") %>% 
  group_by(patient_ID) %>% 
  mutate(lan_within_patient_count = n_distinct(lan_gene, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(lan_gene) %>% 
  mutate(lan_patient_count = n_distinct(patient_ID, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(lan_patient_count = if_else(is.na(lan_gene), NA, lan_patient_count)) %>% 
  left_join(published_lan_sequences %>% 
              select(-core))


clinical_lanA_sequences <- clinical_donor_t %>% 
  select(lan_gene, aa_sequence) %>%
  filter(!is.na(lan_gene)) %>% 
  distinct() %>% 
  mutate(order = str_replace(lan_gene, "lanA_", "")) %>% 
  type_convert() %>% 
  arrange(order) %>% 
  select(-order)

write.fasta(as.list(clinical_lanA_sequences$aa_sequence), clinical_lanA_sequences$lan_gene, "./fasta_files/clinical_lanA_aa_sequences.fa")



## Healthy donor data
healthy_donor_t <- read_csv("./data/human_healthy_donor_table.csv") %>% 
  group_by(patient_ID) %>% 
  mutate(lan_within_patient_count = n_distinct(lan_gene, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(lan_gene) %>% 
  mutate(lan_patient_count = n_distinct(patient_ID, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(lan_patient_count = if_else(is.na(lan_gene), NA, lan_patient_count)) %>% 
  left_join(published_lan_sequences %>% 
              select(-core))


healthy_lanA_sequences <- healthy_donor_t %>% 
  select(lan_gene, aa_sequence) %>%
  filter(!is.na(lan_gene)) %>% 
  distinct() %>% 
  mutate(order = str_replace(lan_gene, "lanA_", "")) %>% 
  type_convert() %>% 
  arrange(order) %>% 
  select(-order)

write.fasta(as.list(healthy_lanA_sequences$aa_sequence), healthy_lanA_sequences$lan_gene, "./fasta_files/healthy_lanA_aa_sequences.fa")







#
n_distinct(clinical_donor_t$patient_ID)
n_distinct(clinical_donor_t$lan_gene)
n_distinct(clinical_donor_t$shotgunSeq_id) + n_distinct(healthy_donor_t$shotgunSeq_id)
n_distinct(healthy_donor_t$lan_gene)
n_distinct(clinical_donor_t$contig)
n_distinct(clinical_donor_t$lanthipeptide, na.rm = TRUE)


clinical_donor_t %>% 
  distinct(lan_gene, predicted_class) %>% 
  filter(!is.na(lan_gene)) %>% 
  group_by(predicted_class) %>% 
  summarize(n = n())


clinical_donor_t %>% 
  select(patient_ID, lan_gene) %>% 
  filter(!is.na(lan_gene)) %>% 
  distinct(patient_ID) %>% 
  summarize(n = n())


healthy_donor_t %>% 
  select(patient_ID, lan_gene) %>% 
  filter(!is.na(lan_gene)) %>% 
  distinct(patient_ID) %>% 
  summarize(n = n())



clinical_donor_t %>% 
  distinct(contig, blast_genus_species) %>% 
  filter(!is.na(contig)) %>% 
  group_by(blast_genus_species) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n))


clinical_donor_t %>% 
  distinct(label) %>% 
  filter(!is.na(label)) %>% 
  summarize(n = n()) 

