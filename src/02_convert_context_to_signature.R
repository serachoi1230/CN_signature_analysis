### Convert CN context into mutational signatures 
### Sera Choi, Jan 2024

### Introduction
# Aim is to decompose CN context file into signatures using reference from 
# COSMIC (https://cancer.sanger.ac.uk/signatures/cn/)

### Context


get_signature_names <- function(ref){
  # get CN signature names (e.g. CN1, CN2 ...)
  label <- colnames(ref)[2:ncol(ref)]
  return(label)
}

prepare_reference <- function(ref){
  # create matrix from reference
  ref <- ref[,2:ncol(ref)] %>% as.matrix()
  return(ref)
}

signature_decomposition <- function(ref, context, signature_names){
  require(nnls)
  result <- nnls::nnls(ref, context[,2])$x
  signature_df <- data.frame("signature" = signature_names,
                             "count" = result)
  return(signature_df)
}

calculate_and_add_proportion <- function(df){
  # add proportion 
  df$proportion <- df$count / sum(df$count)
  return(df)
}


convert_context_to_signature <- function(ref_path, context){
  ref <- read.csv(ref_path)
  signature_names <- get_signature_names(ref)
  ref <- prepare_reference(ref)
  signature_df <- signature_decomposition(ref, context, signature_names)
  signature_df <- calculate_and_add_proportion(signature_df)
  return(signature_df)
}


########
source("/Users/schoi/Github_repos/CN_signature_analysis/src/01_convert_txt_to_context.R")
cna_path <- "/Users/schoi/Github_repos/CN_signature_analysis/data/0b7da678-0837-4fc0-8523-62df1a964939/6a556615-d80f-4b48-94af-a6f25089b40c.wgs.ASCAT.copy_number_variation.seg.txt" 
context <- convert_txt_to_context(cna_path)
ref_path <- "/Users/schoi/Github_repos/CN_signature_analysis/data/COSMIC_v3.3_CN_GRCh37.csv"
result <- convert_context_to_signature(ref_path, context)
result





