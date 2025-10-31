
mkdir -p ../lan_contig_annotations/blast_lanA_annotations
mkdir -p ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations


#Run blastp on all lanA_genes
blastp -query ../fasta_files/lanA_genes.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 16
      
       
       
#Run blastp on all lanA_genes
blastp -query ../fasta_files/lanA_genes_02.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results_02.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25       
       
       
#Run blastp on all lanA_genes
blastp -query ../fasta_files/lanA_genes_03.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results_03.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25       
       

#Run blastp on all lanA_genes
blastp -query ../fasta_files/lanA_genes_04.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results_04.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25    
       
       
blastp -query ../fasta_files/lanA_genes_05.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results_05.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25       
       
       
       
blastp -query ../fasta_files/lanA_genes_06.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/lanA_gene_blastp_results_06.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25       
       
 
 #Combine all annotation files into one tsv file
printf "qseqid\tsseqid\tqacc\tsacc\tqlen\tslen\tqstart\tqend\tsstart\tsend\tqseq\tsseq\tevalue\tbitscore\tscore\tlength\tpident\tnident\tmismatch\tgapopen\tgaps\tpositive\tppos\tframes\tqframe\tsframe\tbtop\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\tsallseqid\tsallgi\tsallacc\tqcovs\tqcovhsp\n" > ../lan_contig_annotations/blast_lanA_annotations/combined_lanA_gene_blastp_results.txt

for annotation_file in $(find ../lan_contig_annotations/blast_lanA_annotations/* -iname "lanA_gene_blastp_results_*")
do
  cat $annotation_file >> ../lan_contig_annotations/blast_lanA_annotations/combined_lanA_gene_blastp_results.txt
done

 
 
#Run blastp on all the potential lanA_genes
blastp -query ../fasta_files/potential_lanA_genes.fa \
       -db ../../../blast_nr_db/nr \
       -out ../lan_contig_annotations/blast_lanA_annotations/potential_lanA_gene_blastp_results.txt \
       -outfmt '6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp' \
       -max_target_seqs 5 \
       -num_threads 25






#Run blastn to get taxonomy for contigs
blastn -query ../fasta_files/lanA_full_contigs.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results.txt \
       -num_threads 25 \
       -max_target_seqs 5



blastn -query ../fasta_files/lanA_full_contigs_02.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_02.txt \
       -num_threads 25 \
       -max_target_seqs 5


blastn -query ../fasta_files/lanA_full_contigs_04.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_04.txt \
       -num_threads 25 \
       -max_target_seqs 5
       

blastn -query ../fasta_files/lanA_full_contigs_05.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_05.txt \
       -num_threads 25 \
       -max_target_seqs 5       
       
       
blastn -query ../fasta_files/lanA_full_contigs_06.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_06.txt \
       -num_threads 25 \
       -max_target_seqs 5       
       
       
blastn -query ../fasta_files/lanA_full_contigs_07.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_07.txt \
       -num_threads 25 \
       -max_target_seqs 5     
       
       
blastn -query ../fasta_files/lanA_full_contigs_08.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_08.txt \
       -num_threads 25 \
       -max_target_seqs 5     
       
blastn -query ../fasta_files/lanA_full_contigs_09.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lanA_contig_blast_taxonomy_annotations/lanA_contig_blast_results_09.txt \
       -num_threads 25 \
       -max_target_seqs 5     
       


#blastn for contigs that have other lan genes without a lanA      
blastn -query ../fasta_files/lan_gene_full_contigs.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lan_contig_blast_taxonomy_annotations/lan_contig_blast_results_01.txt \
       -num_threads 25 \
       -max_target_seqs 5       


blastn -query ../fasta_files/lan_gene_full_contigs_02.fa \
       -db ../../../blast_nt_db/nt \
       -outfmt "6 qseqid sseqid qacc sacc qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch gapopen gaps positive ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sallseqid sallgi sallacc qcovs qcovhsp" \
       -out ../lan_contig_annotations/lan_contig_blast_taxonomy_annotations/lan_contig_blast_results_02.txt \
       -num_threads 25 \
       -max_target_seqs 5       







