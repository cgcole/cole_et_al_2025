## Packages

# First ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages if not already present
bioc_pkgs <- c("phyloseq", "Biostrings", "msa")  # Add any others
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}


if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, 
  readxl, 
  phyloseq,
  Biostrings, 
  ggpubr,
  scales,
  writexl,
  rstatix,
  ggrepel,
  VennDiagram,
  lubridate,
  ComplexUpset,
  seqinr,
  umap,
  msa,
  pheatmap
)
pacman::p_load_gh("ying14/yingtools2")

rm(bioc_pkgs, pkg)

# Visit https://github.com/omixer/omixer-rpmR for instructions to install Omixer


## Generate folders to store generated tables and plots
if (!file.exists("./plots")) {
  dir.create("./plots")
}

if (!file.exists("./tables")) {
  dir.create("./tables")
} 

if (!file.exists("./fasta_files")) {
  dir.create("./fasta_files")
} 




## Functions
getShades <- function(spdf){
  
  hexcol <- unique(spdf$color)
  
  if (nrow(spdf) > 3){
    resdf <- spdf %>%
      mutate(cols = rep(shades(hexcol, variation = 0.25),
                        length.out = nrow(spdf)
      )
      )
  } else {
    resdf <- spdf %>%
      mutate(cols = shades(hexcol, variation = 0.25, ncolor = nrow(spdf)))
  }
  
  return(resdf)
}



getRdpPal <- function(tax) {
  
  require(tidyverse)
  require(yingtools2)
  
  tax <- tax %>%
    ungroup()
  
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax.obj)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family",
             "Genus")
  
  if (!all(ranks %in% names(tax))) {
    stop("Error: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  
  tax.dict <- tax %>%
    dplyr::select(all_of(ranks)) %>%
    distinct()
  
  # set all color to gray as base
  tax.dict <- tax.dict %>%
    mutate(color = rep(shades("gray", variation = 0.25),
                       length.out = nrow(tax.dict)))
  
  # color for each level ------------------------------------------------------------
  phypal <- tibble(Phylum = c("Proteobacteria",
                              "Thermodesulfobacteriota",
                              "Pseudomonadota",
                              "Actinobacteria",
                              "Actinomycetota",
                              "Bacteroidetes",
                              "Bacteroidota"),
                   phycol = c("red","red","red",
                              "#A77097","#A77097","#51AB9B","#51AB9B"))
  ordpal <- tibble(Order = c("Clostridiales", "Eubacteriales"),
                   ordcol = c("#9C854E", "#9C854E"))
  fampal <- tibble(Family = c("Lachnospiraceae","Ruminococcaceae","Oscillospiraceae","Erysipelotrichaceae", "Lactobacillaceae", "Staphylococcaceae"),
                   famcol = c("#EC9B96","#9AAE73","#9AAE73","orange","#3b51a3", "#f1eb25"))
  genpal <- tibble(Genus = c("Enterococcus","Streptococcus", "Turicibacter", "Akkermansia"),
                   gencol = c("#129246","#9FB846", "#bf5700", "#83AAD7"))
  
  tax.split <- tax.dict %>%
    left_join(phypal) %>%
    left_join(ordpal) %>%
    left_join(fampal) %>%
    # ambiguous genus match
    # mutate(gencol = case_when(
    #   grepl("Enterococcus$", Genus) ~ "#129246",
    #   grepl("Streptococcus$", Genus) ~ "#9FB846",
    #   grepl("Staphylococcus$", Genus) ~ "#f1eb25",
    #   TRUE ~ NA_character_
    # )) %>%
    left_join(genpal) %>%
    mutate(color = case_when(
      !is.na(gencol) ~ gencol,
      !is.na(famcol) ~ famcol,
      !is.na(ordcol) ~ ordcol,
      !is.na(phycol) ~ phycol,
      TRUE ~ color)
    ) %>%
    dplyr::select(Kingdom:Genus, color) %>%
    group_split(color)
  
  tax.color <- bind_rows(lapply(tax.split,getShades))
  tax.palette <- structure(tax.color$cols, names = as.character(tax.color$Genus))
  return(tax.palette)
}



