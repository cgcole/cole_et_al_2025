#!/bin/bash

# Generate the directory containing only lan gene pfams
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-file ../metadata/lan_gene_pfam_accessions.txt \
                                              --output-directory ../lan_gene_pfams


mkdir -p ../00_LOGS
touch ../00_LOGS/lan_gene_pfam_search_log.txt

# Run hmm on search 
for contig in $(find ../02_CONTIGS/* -iname "*-contigs.db")
do 
    
    anvi-run-hmms -c $contig \
                  --hmm-profile-dir ../lan_gene_pfams \
                  --add-to-functions-table \
                  --hmmer-program hmmsearch \
                  -T 31 
    
    echo "Finished processing $contig"
    
    echo "$contig\n" >> ../00_LOGS/lan_gene_pfam_search_log.txt
    
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

