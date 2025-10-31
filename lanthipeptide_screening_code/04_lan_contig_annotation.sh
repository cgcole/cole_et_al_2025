
# Split the contigs into smaller fasta files for bakta to process
seqkit grep -f ../metadata/lan_contig_id_subset_01.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_01.fa
seqkit grep -f ../metadata/lan_contig_id_subset_02.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_02.fa
seqkit grep -f ../metadata/lan_contig_id_subset_03.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_03.fa
seqkit grep -f ../metadata/lan_contig_id_subset_04.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_04.fa
seqkit grep -f ../metadata/lan_contig_id_subset_05.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_05.fa
seqkit grep -f ../metadata/lan_contig_id_subset_06.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_06.fa
seqkit grep -f ../metadata/lan_contig_id_subset_07.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_07.fa
seqkit grep -f ../metadata/lan_contig_id_subset_08.txt ../fasta_files/lan_contigs.fa > ../fasta_files/lan_contig_subset_08.fa

# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_01/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_01 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contigs_01.fa
fi

# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_02/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_02 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contigs_02.fa
fi


# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_01.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
            -T 31 
            
# Kegg Annoations
anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions.txt
            
            
            
            
## Additional contigs to annotate for subset 2   
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_03/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_03 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_02.fa
fi

# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_02.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_2.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_2-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_2.txt
            
            
            
## Additional contigs to annotate for subset 3 
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_04/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_04 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_03.fa
fi


# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_03.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_3.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_3.txt
            

      
            
## Additional contigs to annotate for subset 4
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_05/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_05 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_04.fa
fi


# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_04.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_4-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_4.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_4-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_4-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_4-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_4-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_4.txt
            


## Additional contigs to annotate for subset 5
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_06/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_06 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_05.fa
fi




# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_05.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_5-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_5.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_5-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_5-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_5-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_5-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_5.txt
            




## Additional contigs to annotate for subset 6
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_07/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_07 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_06.fa
fi




# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_06.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_6-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_6.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_6-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_6-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_6-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_6-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_6.txt
            
                            



## Additional contigs to annotate for subset 7
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_08/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_08 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_07.fa
fi




# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_07.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_7-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_7.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_7-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_7-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_7-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_7-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_7.txt
            
            
            
## Additional contigs to annotate for subset 8
# Run bakta on lan contig fasta file

if [ -f ../lan_contig_annotations/bakta_annotations/bakta_annotations_09/lan_contigs-bakta.gbk ]
then
continue
else
  bakta --db ../../../bakta_db \
      --prefix lan_contigs-bakta \
      --output ../lan_contig_annotations/bakta_annotations/bakta_annotations_09 \
      --meta \
      --threads 15 \
      --keep-contig-headers \
      --skip-trna \
      --skip-crispr \
      --skip-rrna \
      --skip-tmrna \
      --skip-ncrna \
      --skip-ncrna-region \
      ../fasta_files/lan_contig_subset_08.fa
fi




# Generate a contigs database for Anvio annotation
anvi-gen-contigs-database -f ../fasta_files/lan_contig_subset_08.fa \
                              --project-name lan_contig_annotation \
                              -o ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_8-contigs.db \
                              --external-gene-calls ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_gene_calls_8.txt \
                              -T 31


# NCBI Cog annotations
anvi-run-ncbi-cogs -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_8-contigs.db \
            -T 31 
            
# Kegg Annoations
#anvi-run-kegg-kofams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_3-contigs.db \
#            -T 31 
            
# Pfam Annoations
anvi-run-pfams -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_8-contigs.db \
            -T 31 
            
# Cazyme Annoations
anvi-run-cazymes -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_8-contigs.db \
            -T 31 

# Export gene annotations from contigs database
anvi-export-functions -c ../lan_contig_annotations/anvio_annotations/lan_contig_annotation_8-contigs.db \
            -o ../lan_contig_annotations/anvio_annotations/lan_contig_gene_functions_8.txt
            
                                           
                      
    
#Combine all annotation files into one tsv file
printf "contig_gene_id\tghostkoala_accession\tghostkoala_function\tghostkoala_score\tsecond_best\tsecond_best_score\n" > ../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_annotations.txt

for annotation_file in $(find ../lan_contig_annotations/ghostkoala_annotations/* -iname "user_ko_definition.txt")
do
  cat $annotation_file >> ../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_annotations.txt
done


#Combine all taxonomy files into one tsv file
printf "contig_gene_id\tghostkoala_accession\tkingdom\tphylum\tgenus\ttax_entry\ttax_score\n" > ../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_taxonomy.txt

for taxonomy_file in $(find ../lan_contig_annotations/ghostkoala_annotations/* -iname "user.out.top.gz")
do
  if [ ! -f ${taxonomy_file%.gz} ]  
  then
    gzip -dk $taxonomy_file
  fi
  cat ${taxonomy_file%.gz} >> ../lan_contig_annotations/ghostkoala_annotations/combined_ghostkoala_taxonomy.txt
done


            
            
            
