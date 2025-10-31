
# Make directory to hold all the fasta files from BV-BRC
mkdir -p ../fasta_files

# Download files from BV-BRC database
for file in $(p3-ls /dfi_uchicago@patricbrc.org/DFIclinical_contigs)
do
    p3-cp ws:/dfi_uchicago@patricbrc.org/DFIclinical_contigs/${file} ../fasta_files
done

# gzip all files (some are already gzipped)
gzip ../fasta_files/*.fasta 


