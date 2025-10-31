# Data Repo for Cole et al., 2025

Title: Lantibiotic-producing bacteria impact microbiome resilience and colonization resistance

For any inquires regarding the contents of this repo, please reach out to Cody Cole (cgcole@uchicago.edu).

## Folder structure

1. **code**: contains all R scripts for generating figures. 
2. **data**: data used to generate figures
- `Blautia_producta_KH6-ghostkoala.txt`: GhostKoala annotations for Blautia producta KH6.
- `Blautia_pseudococcoides_SCSK-ghostkoala.txt`: GhostKoala annotations for Blautia pseudococcoides SCSK.
- `human_clinical_donor_table.csv`: Shotgun sequencing results of clinical donor fecal samples paired with lanthipeptides detected.
- `human_healthy_donor_table.csv`: Shotgun sequencing results of healthy human donor fecal samples paired with lanthipeptides detected.
- `mouse_16S_sequencing.csv`: 16S rRNA gene sequencing data for all mouse experiments.
- `mouse_bile_acid_quant_metabolomics.csv`: Quantitative metabolomics for bile acids in mouse fecal pellets.
- `mouse_c_diff_cfu.csv`: CFU counts for C. difficile in the fecal pellets of challenged mice.
- `mouse_c_diff_weights.csv`: Mouse weights for mice challenged with C. difficile.
- `mouse_kleb_cfu.csv`: CFU counts for K. pneumoniae in the fecal pellets of challenged mice.
- `mouse_pfbbr_quant_metaboloimcs.csv`: Quantitative metabolomics from PFBBr panel in mouse fecal pellets.
- `previously_published_lanthipeptides.csv`: Lanthipeptides that have been characterized and previously published on.
- `unique_identified_lanthipeptides.csv`: All unique lanthipeptides identified in this study.
3. **fasta_files**: This folder contains the assembled genomes of B. producta KH6 and B. pseudococcoides SCSK that were used for generating the gut metabolic modules (GMMs) for Figure S2A.
4. **lanthipeptide_screening_code**: This folder contains the code that was used to screen for lanthipeptides in the human fecal shotgun metagenomic samples. The code is provided as is for transparency. It is not executable in this repository. 
5. **mouse_timelines**: Timelines for mouse experiments.


## How to Run

1. Open `sandbox.Rproj` in RStudio.
2. Run the `code/00_master.R` file, or run the scripts individually.
3. Output will generate the `plots/` and `tables/` folders that contain all of the plots and tables used in the paper.

