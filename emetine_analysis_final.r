#Emetine Analysis - mostly pulled from 1.2 fold change script
library(dplyr)
library(biomaRt)
library(Biostrings)
library(WebGestaltR)

#Load and split starfish names
fp <- "C:/Users/Brennan/Documents/Kettenbach Lab/Starfish Gene Ontology Project/"
sf_names <- read.delim(paste0(fp, "emetine proteins names.txt"), stringsAsFactors = F, header = F)
comma_splitter <- function(input){
  split_name <- unlist(strsplit(input, "; ")) #Returns character vectors that can have multiple values
  return(split_name)
}
colnames(sf_names) <- "name"
sf_names$split_names <- lapply(sf_names$name, comma_splitter)
all_sf_names <- unique(unlist(sf_names$split_names)) #169 names

#We're just going to do a fresh BLAST on these
sf_proteins <- readAAStringSet(paste0(fp, "MITEchinobaseassembly_Patiriaminiata_Combined.fasta")) #loads 84328 proteins...
sf_proteins_names <- names(sf_proteins) #May be necessary to allow searching of FASTA names

#Sub function that does the work of extracting the starfish sequences
#Now handles sf_name in row 1...
extract_sf_seqs <- function(input){
  sf_protein <- grep(input, sf_proteins_names, fixed = T)
  if(length(sf_protein) > 1){
    warning(sf_protein_name, " Caused more than one protein to be found")
  }
  return(sf_proteins[sf_protein[1]])
}
sf_seqs <- lapply(all_sf_names, extract_sf_seqs)
sf_seqs <- do.call(c, sf_seqs)
writeXStringSet(sf_seqs, paste0(fp, "emetine proteins.fasta"))

#BLASTed these on Discovery
blast_fp <- "C:/Users/Brennan/Documents/Kettenbach Lab/Starfish Gene Ontology Project/BLAST Info/"
blast_results <- read.csv(paste0(blast_fp, "emetine_proteins.csv"), header = F, stringsAsFactors = F)
colnames(blast_results) <- c("sf_prot",	"ncbi_accession",	"pcnt_identity",
                                 "alignment_length",	"mismatches",	"gap_opens",	"sf_start",
                                 "sf_end",	"ncbi_start",	"ncbi_end",	"E_value", "BLAST_score")

filter_prefer <- function(input_df, target_col, preference){
  filter_prefer2 <- function(input, preference){
    return(grepl(preference, input))
  }
  input_df$check <- lapply(input_df[[target_col]], filter_prefer2, preference)
  if(all(input_df$check == F)){
    return(input_df[ ,which(colnames(input_df) != "check")])
  }
  else{
    output_df <- input_df[input_df$check == T, ]
    return(output_df[ ,which(colnames(output_df) != "check")])
  }
}

blast_results <- blast_results %>% dplyr::select(sf_prot, ncbi_accession, E_value, BLAST_score, pcnt_identity) %>%
  filter(E_value < .01) %>% group_by(sf_prot) %>% filter_prefer("ncbi_accession", "NP_") %>% top_n(-1, E_value) %>% 
  slice(1) %>% ungroup()

name_fetcher <- function(input_vector){
  result <- NA
  for(i in c(1:length(input_vector))){
    if(length(grep(input_vector[i], blast_results$sf_prot)) >= 1){
      match_row <- grep(input_vector[i], blast_results$sf_prot)
    }
    else if(length(grep(input_vector[i], blast_results$sf_prot)) == 0){
      result[i] <- NA
      next
    }
    result[i] <- blast_results$ncbi_accession[match_row[1]]
  }
  return(result)
}
sf_names$ncbi_accession <- lapply(sf_names$split_names, name_fetcher)

all_match <- function(input_vector){
  if(length(unlist(input_vector)) == 1){
    if(is.na(input_vector)){
      return(F)
    }
  }
  if(length(-which(is.na(input_vector))) >0)
    input_vector <- input_vector[-which(is.na(input_vector))]
  if(length(unique(input_vector)) == 1)
    return(TRUE)
  else
    return(FALSE)
}
sf_names$all_match <- lapply(sf_names$ncbi_accession, all_match)

name_condenser <- function(input_row, target_col){
  if(input_row$all_match == T){
    result <- unlist(input_row[target_col])
    if(length(-which(is.na(result))) >0)
      result <- target_names[-which(is.na(result))]
    return(result[1])
  }
  else{
    return(NA)
  }
}
sf_names$refseq_accession <- apply(sf_names, 1, name_condenser, "ncbi_accession")
sf_names$refseq_accession <- sub("\\..*", "", sf_names$refseq_accession)

#Fetch UniProt names and numbers

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
biomart_result <- getBM(attributes = c("refseq_peptide", "uniprot_gn_symbol", "uniprot_gn_id"),
                        filters = "refseq_peptide",
                        values = sf_names$refseq_accession,
                        mart = ensembl)
biomart_result <- biomart_result %>% group_by(refseq_peptide) %>% slice(1) %>% ungroup #Only 61 unique entries...

match_name_rows <- function(match_name, search_table, search_row, keep_rows){
  match_row <- grep(match_name, search_table[[search_row]])
  if(length(match_row >= 1)){
    result_row <- search_table[match_row[1], ]
    return(list(result_row[keep_rows]))
  }
  else{
    return(rep(NA, length(keep_rows)))
  }
}
matched_name_rows <- sapply(sf_names$refseq_accession, match_name_rows, biomart_result, "refseq_peptide", c(2:3))
matched_name_rows <- data.frame(matrix(unlist(matched_name_rows), ncol = 2, byrow = T))
colnames(matched_name_rows) <- colnames(biomart_result[ ,c(2:3)])
sf_names <- cbind(sf_names, matched_name_rows)

blast_results$ncbi_accession <- sub("\\..*", "", blast_results$ncbi_accession)
matched_name_rows <- sapply(sf_names$refseq_accession, match_name_rows, blast_results, "ncbi_accession", c(3:5))
matched_name_rows <- data.frame(matrix(unlist(matched_name_rows), ncol = 3, byrow = T))
colnames(matched_name_rows) <- colnames(blast_results[ ,c(3:5)])
sf_names <- cbind(sf_names, matched_name_rows)
#Honestly match_name_rows could probably be wrapped up in another function layer to fully automate up to the cbind stuff

ont_list <- c("Biological_Process", "Cellular_Component", "Molecular_Function")
output_fp <- "C:/Users/Brennan/Documents/Kettenbach Lab/Starfish Gene Ontology Project/GSEA Output"
anal_name <- "emetine_proteins" #Change project name every time.
gene_set <- as.character(unlist(sf_names$refseq_accession)) #Change gene_set depending on what you're doing

for(i in c(1:3)){
  WebGestaltR(enrihMethod = "ORA", organism = "hsapiens",
              enrichDatabase = paste0("geneontology_", ont_list[i], "_noRedundant"),
              interestGene = gene_set,
              interestGeneType = "refseq_peptide", #This may need to change depending on what you want to use
              referenceSet = "genome_protein-coding", #refseq_peptide and uniprotswissprot are the most likely two
              sigMethod = "top",
              outputDirectory = output_fp,
              projectName = paste0(anal_name, ont_list[i])
  )
}

sf_names$all_refseq_accessions <- sapply(sf_names$ncbi_accession, function(input) paste(unlist(input), collapse = "; "))

output_csv <- sf_names %>% dplyr::select(name, all_refseq_accessions, refseq_accession, uniprot_gn_symbol, uniprot_gn_id,
                                  E_value, BLAST_score, pcnt_identity)
write.csv(output_csv, paste(fp, "emetine_proteins_analysis.csv"))
