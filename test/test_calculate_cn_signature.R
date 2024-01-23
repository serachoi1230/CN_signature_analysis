### Unit tests for functions to calculate CN signatures 
### Sera Choi, Jan 2024

### setup
library(testthat)
source("/Users/schoi/Github_repos/CN_signature_analysis/src/01_convert_txt_to_context.R")
source("/Users/schoi/Github_repos/CN_signature_analysis/src/02_convert_context_to_signature.R")


### testing read_and_classify_cna()
test_read_and_classify_cna <- function(){
  
  path <- "/Users/schoi/Github_repos/CN_signature_analysis/test/test_copy_number_variation.seg.txt"
  df <- read_and_classify_cna(path)
  
  expect_equal(dim(df), c(10,10))
  expect_equal(df[10,]$cn_class, 1)
  expect_equal(df[10,]$event_class, "LOH")
  
}


### testing make_and_count_catalogue_label
test_make_and_count_catalogue_label <- function(){
  
  cn_class <- c("3-4","1")
  event_class <- c("het","LOH")
  size_class <- c("10Mb-40Mb", "0-100kb")
  
  df <- data.frame(cn_class, event_class, size_class)
  df <- make_and_count_catalogue_label(df)
  
  expect_equal(dim(df), c(48,2))
  expect_equal(sum(df$count), 2)
  
}


### testing convert_context_to_signature
test_convert_context_to_signature <- function(){
  
  ref_path <- "/Users/schoi/Github_repos/CN_signature_analysis/data/COSMIC_v3.3_CN_GRCh37.csv"
  context_path <- "/Users/schoi/Github_repos/CN_signature_analysis/test/test_copy_number_variation.seg.txt"
  
  context <- convert_txt_to_context(context_path)
  df <- convert_context_to_signature(ref_path, context)
  
  expect_equal(dim(df), c(24,3))
  expect_equal(df[1,]$count, 0)
  expect_equal(round(df[15,]$count), 2)

}


### Master test function
test_calculate_cn_signature <- function(){
  
  test_read_and_classify_cna()
  test_make_and_count_catalogue_label()
  test_convert_context_to_signature()
  
}
