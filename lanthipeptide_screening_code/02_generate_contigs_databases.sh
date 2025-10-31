
# Make necessary directories
mkdir -p ../01_FASTA
mkdir -p ../02_CONTIGS
mkdir -p ../lan_contig_annotations


# Create contig databases and extract gene calls
for fasta_file in $(find ../DFIclinical_contigs/* -iname "*.fasta.gz")
do
    
    name=$(basename -- "${fasta_file%.megahit.contigs.fasta.gz}") 
    fixedname=$(echo ${name} | sed -e 's/\./_/g' -e 's/-/_/g' -e 's/ /_/g' -e 's/\//_/g')

    if [ -f ../02_CONTIGS/$fixedname-contigs.db ]
        then
            continue
        else
            mkdir -p ../01_FASTA/$fixedname 

            gunzip -c ${fasta_file} > /tmp/$fixedname.fasta &

            anvi-script-reformat-fasta /tmp/$fixedname.fasta \
                               -o ../01_FASTA/$fixedname/$fixedname-contigs.fa \
                               --simplify-names \
                               --prefix $fixedname \
                               --seq-type NT

            rm /tmp/$fixedname.fasta

            anvi-gen-contigs-database -f ../01_FASTA/$fixedname/$fixedname-contigs.fa \
                              --project-name $fixedname \
                              -o ../02_CONTIGS/$fixedname-contigs.db \
                              -T 31

            anvi-export-gene-calls -c ../02_CONTIGS/$fixedname-contigs.db \
                           --gene-caller prodigal \
                           -o ../01_FASTA/$fixedname/${fixedname}_gene_calls.txt
    fi
done


# Generate the directory containing only lan gene pfams
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-file ../metadata/lan_gene_pfam_accessions.txt \
                                              --output-directory ../lan_gene_pfams \
                                              --force-overwrite



# Run hmm search on contig databases to identify pfam domains associated with lan genes
for contig in $(find ../02_CONTIGS/* -iname "*-contigs.db")
do 
    
    anvi-delete-functions -c $contig \
                          --annotation-source lan_gene_pfams
    
    anvi-run-hmms -c $contig \
                  --hmm-profile-dir ../lan_gene_pfams \
                  --add-to-functions-table \
                  --hmmer-program hmmsearch \
                  -T 8
    
    echo "Finished processing $contig"
done



# Extract all the lan gene calls that pfam identified and move them to a single table
printf "gene_callers_id\tsource\taccession\tfunction\te_value\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_gene_pfam_gene_calls.txt    

for contig in $(find ../02_CONTIGS/* -iname "*-contigs.db")
do 
    
    samplename=$(basename "$contig" -contigs.db)
    
    anvi-export-functions -c $contig \
            -o ../metadata/${samplename}_functions.tmp \
            --annotation-sources lan_gene_pfams
    
    tail -n +2 ../metadata/${samplename}_functions.tmp | sort -t $'\t' -k1,1 > ../metadata/${samplename}_sorted_functions.tmp
    tail -n +2 ../01_FASTA/${samplename}/${samplename}_gene_calls.txt | sort -t $'\t' -k1,1 > ../metadata/${samplename}_modified_gene_calls.tmp

    join -a 1 -1 1 -2 1 -t $'\t' ../metadata/${samplename}_sorted_functions.tmp ../metadata/${samplename}_modified_gene_calls.tmp > ../metadata/joined_table.tmp
 
    cat ../metadata/joined_table.tmp >> ../metadata/lan_gene_pfam_gene_calls.txt
 
    rm ../metadata/*.tmp
    
    echo "Finished extracting lan_gene_pfams from $contig"
                  
done


# Remove unnessary pfam domains to reduce number of contigs to investigate
printf "gene_callers_id\tsource\taccession\tfunction\te_value\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_gene_pfam_gene_calls_subset.txt    

awk -F'\t' '$3 == "PF04738.18" || $3 == "PF02052.20" || $3 == "PF04604.19" || $3 == "PF16934.10" || $3 == "PF13817.11" || $3 == "PF18923.6" || $3 == "PF03572.23" || $3 == "PF01964.24" || $3 == "PF14867.11" || $3 == "PF20693.3" || $3 == "PF14028.11" || $3 == "PF05147.19" || $3 == "PF13575.12" || $3 == "PF00082.27" || $3 == "PF12730.13" || $3 == "PF18218.6"' ../metadata/lan_gene_pfam_gene_calls.txt >> ../metadata/lan_gene_pfam_gene_calls_subset.txt

printf "gene_callers_id\tsource\taccession\tfunction\te_value\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_gene_pfam_gene_calls_subset_2.txt    

awk -F'\t' '$3 == "PF04738.18" || $3 == "PF02052.20" || $3 == "PF04604.19" || $3 == "PF16934.10" || $3 == "PF18923.6" || $3 == "PF14867.11" || $3 == "PF14028.11" || $3 == "PF05147.1" || $3 == "PF13575.11" || $3 == "PF18218.6" || $3 == "PF01532.25" || $3 == "PF08130.16"' ../metadata/lan_gene_pfam_gene_calls.txt >> ../metadata/lan_gene_pfam_gene_calls_subset_2.txt

printf "gene_callers_id\tsource\taccession\tfunction\te_value\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_gene_pfam_gene_calls_subset_3.txt    

awk -F'\t' '$3 == "PF10439.14" || $3 == "PF03902.18"' ../metadata/lan_gene_pfam_gene_calls.txt >> ../metadata/lan_gene_pfam_gene_calls_subset_3.txt

printf "gene_callers_id\tsource\taccession\tfunction\te_value\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_gene_pfam_gene_calls_subset_4.txt    

awk -F'\t' '$3 == "PF16934.11" || $3 == "PF17914.6" || $3 == "PF04369.19" || $3 == "PF19402.4"' ../metadata/lan_gene_pfam_gene_calls.txt >> ../metadata/lan_gene_pfam_gene_calls_subset_4.txt




# Extract all of the gene calls that contain a lanA aa sequence identified by RODEO and other sources
printf "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lanA_gene_calls.txt

cat ../metadata/rodeo_lanA_genes.txt ../metadata/alternate_lanA_genes.txt ../metadata/core_lanA_genes.txt ../metadata/leader_lanA_genes.txt ../metadata/core_lanA_genes.txt ../metadata/lanA_blast_sequences.txt > ../metadata/all_lanA_genes.txt

for file in $(find ../01_FASTA/*/* -iname "*gene_calls.txt")
do 
    echo "Processing $file"
    
    grep -Ff ../metadata/all_lanA_genes.txt ${file} >> ../metadata/lanA_gene_calls.txt
    
done



# Make a list of all contig names that contain an lan aa sequences
awk 'NR > 1 {print $2}' ../metadata/lanA_gene_calls.txt | sort -u > ../metadata/lanA_contig_ids.txt
awk 'NR > 1 {print $6}' ../metadata/lan_gene_pfam_gene_calls_subset.txt | sort -u > ../metadata/lan_pfam_contig_ids.txt
awk 'NR > 1 {print $6}' ../metadata/lan_gene_pfam_gene_calls_subset_2.txt | sort -u > ../metadata/lan_pfam_contig_ids_2.txt
awk 'NR > 1 {print $6}' ../metadata/lan_gene_pfam_gene_calls_subset_3.txt | sort -u > ../metadata/lan_pfam_contig_ids_3.txt
awk 'NR > 1 {print $6}' ../metadata/lan_gene_pfam_gene_calls_subset_4.txt | sort -u > ../metadata/lan_pfam_contig_ids_4.txt

cat ../metadata/manual_contig_ids.txt ../metadata/lanA_contig_ids.txt ../metadata/lan_pfam_contig_ids.txt ../metadata/lan_pfam_contig_ids_2.txt ../metadata/lan_pfam_contig_ids_3.txt ../metadata/lan_pfam_contig_ids_4.txt | sort -u > ../metadata/lan_contig_ids.txt

printf "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n" > ../metadata/lan_contig_gene_calls.txt

for file in $(find ../01_FASTA/*/* -iname "*gene_calls.txt")
do 
    echo "Processing $file"
    grep -Ff ../metadata/lan_contig_ids.txt ${file} >> ../metadata/lan_contig_gene_calls.txt
done



# Extract the nucleotide sequences for lan contigs and compile into a fasta file for bakta annotation
rm ../fasta_files/lan_contigs.fa #In case you need to rerun

for file in $(find ../01_FASTA/*/* -iname "*-contigs.fa")
do 
    echo "Processing $file"
    grep -A 1 -f ../metadata/lan_contig_ids.txt ${file} | grep -v "^--$" >> ../fasta_files/lan_contigs.fa
done






