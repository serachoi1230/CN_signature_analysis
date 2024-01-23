### Convert CN context into mutational signatures 
### Sera Choi, Jan 2024

### Introduction
# Aim is to decompose CN context file into signatures using reference from 
# COSMIC (https://cancer.sanger.ac.uk/signatures/cn/)


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
  # fit signatures using non-negative least squares
  # output data frame with signature contribution count
  require(nnls)
  result <- nnls::nnls(ref, context[,2])$x
  signature_df <- data.frame("signature" = signature_names,
                             "count" = result)
  return(signature_df)
}

calculate_and_add_proportion <- function(df){
  # calculate proportion from count
  df$proportion <- df$count / sum(df$count)
  return(df)
}

### Master function
convert_context_to_signature <- function(ref_path, context){
  ref <- read.csv(ref_path)
  signature_names <- get_signature_names(ref)
  ref <- prepare_reference(ref)
  signature_df <- signature_decomposition(ref, context, signature_names)
  signature_df <- calculate_and_add_proportion(signature_df)
  return(signature_df)
}

