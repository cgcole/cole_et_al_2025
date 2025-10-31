#!/bin/bash

mkdir -p ../fasta_files/ncbi_lanA_genes

# Iterate over each annotation in the input file
while IFS= read -r annotation; do
  echo "Querying NCBI for: ${annotation}"
  
  annotation_encoded=$(echo "$annotation" | sed 's/ /+/g')
  
  # Use esearch to search for the annotation in the protein database
  esearch -db protein -query "${annotation_encoded}" < /dev/null | \
  efetch -format fasta >> "../fasta_files/ncbi_lanA_genes/ncbi_lanA_aa_genes.fa"

  # Optional: Add a delay to avoid overwhelming NCBI's servers
  sleep 1  # Adjust sleep time as needed
done < ../metadata/blastp_lanA_annotations.txt


# Iterate over each annotation in the input file
while IFS= read -r annotation; do
  echo "Querying NCBI for: ${annotation}"
  
  annotation_encoded=$(echo "$annotation" | sed 's/ /+/g')

  # Use esearch to search for the annotation in the protein database
  esearch -db protein -query "${annotation_encoded}" < /dev/null | \
  elink -target nuccore | \
  efetch -format fasta >> "../fasta_files/ncbi_lanA_genes/ncbi_lanA_dna_genes.fa"

  # Optional: Add a delay to avoid overwhelming NCBI's servers
  sleep 1  # Adjust sleep time as needed
done < ../metadata/blastp_lanA_annotations.txt



ERRORS:
  Lantibiotic duramycin; AltName: Full=Antibiotic PA48009; AltName: Full=Leucopeptin
