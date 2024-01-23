### Convert .txt into copy-number context 
### Sera Choi, Jan 2024

### Introduction
# Aim is to convert .txt file containing CNAs into 48-channel CN classification
# context data frame (https://cancer.sanger.ac.uk/signatures/cn/)

### Content 
# 1. read_and_classify_cna() 
# 2. make_and_count_catalogue_label() 
# 3. convert_txt_to_context()

### setup
library(tidyverse)

#############################

read_cna_txt <- function(path){
  # read CNA.txt file 
  # create data frame 
  df <- read_delim(path)
  return(df)
}

classify_cn <- function(df){
  # add "cn_class" column
  # classify copy numbers into categories 
  df <- mutate(df, cn_class = case_when(
    Copy_Number == 0 ~ "0",
    Copy_Number == 1 ~ "1",
    Copy_Number == 2 ~ "2",
    Copy_Number == 3 | Copy_Number == 4 ~ "3-4",
    Copy_Number >=5 & Copy_Number < 9 ~ "5-8",
    Copy_Number >=9 ~ "9+"
  ))
  return(df)
}

classify_event <- function(df){
  # add "event_class" column
  # classify mutation event into categories
  df <- mutate(df, event_class = case_when(
    Minor_Copy_Number == 0 & Major_Copy_Number == 0 ~ "homdel",
    Minor_Copy_Number == 0 & Major_Copy_Number > 0 ~ "LOH",
    Minor_Copy_Number > 0 & Major_Copy_Number > 0 ~ "het",
  ))
  return(df)
}

classify_size <- function(df){
  # add "size_class" column
  # classify size into categories 
  df <- mutate(df, size_class = case_when(
    Copy_Number == 0 & (End - Start) >= 1000000 ~ ">1Mb",
    (End - Start) < 100000 ~ "0-100kb",
    (End - Start) >= 100000 & (End - Start) < 1000000 ~ "100kb-1Mb",
    Copy_Number != 0 & (End - Start) >= 1000000 & (End - Start) < 10000000 ~ "1Mb-10Mb",
    Copy_Number != 0 & (End - Start) >= 10000000 & (End - Start) < 40000000 ~ "10Mb-40Mb",
    Copy_Number != 0 & (End - Start) >= 40000000 ~ ">40Mb"   
  ))
  return(df)
}

### Sub-function 1
read_and_classify_cna <- function(path){
  cna <- read_cna_txt(path)
  cna <- classify_cn(cna)
  cna <- classify_event(cna)
  cna <- classify_size(cna)
  return(cna)
}

add_catalogue_labels <- function(df){
  # make labels with cn_class, event_class and size_class
  # add "catalogue_label" column
  df$catalogue_label <- paste0(df$cn_class, ":", df$event_class, ":", df$size_class)
  return(df)
}

load_reference_catalogue_label <- function(){
  catalogue_labels <- c("1:LOH:0-100kb","1:LOH:100kb-1Mb","1:LOH:1Mb-10Mb","1:LOH:10Mb-40Mb",'1:LOH:>40Mb',
                        '2:het:0-100kb','2:het:100kb-1Mb','2:het:1Mb-10Mb','2:het:10Mb-40Mb','2:het:>40Mb',
                        '2:LOH:0-100kb','2:LOH:100kb-1Mb','2:LOH:1Mb-10Mb','2:LOH:10Mb-40Mb','2:LOH:>40Mb',
                        '3-4:het:0-100kb','3-4:het:100kb-1Mb','3-4:het:1Mb-10Mb','3-4:het:10Mb-40Mb','3-4:het:>40Mb',
                        '3-4:LOH:0-100kb','3-4:LOH:100kb-1Mb','3-4:LOH:1Mb-10Mb','3-4:LOH:10Mb-40Mb','3-4:LOH:>40Mb',
                        '5-8:het:0-100kb','5-8:het:100kb-1Mb','5-8:het:1Mb-10Mb','5-8:het:10Mb-40Mb','5-8:het:>40Mb',
                        '5-8:LOH:0-100kb','5-8:LOH:100kb-1Mb','5-8:LOH:1Mb-10Mb','5-8:LOH:10Mb-40Mb','5-8:LOH:>40Mb',
                        '9+:het:0-100kb','9+:het:100kb-1Mb','9+:het:1Mb-10Mb','9+:het:10Mb-40Mb','9+:het:>40Mb',
                        '9+:LOH:0-100kb','9+:LOH:100kb-1Mb','9+:LOH:1Mb-10Mb','9+:LOH:10Mb-40Mb','9+:LOH:>40Mb',
                        '0:homdel:0-100kb','0:homdel:100kb-1Mb','0:homdel:>1Mb')
  return(catalogue_labels)
}

count_catalogue_label <- function(df, reference_catalogue_label){
  # count catalogue label
  # output count data frame
  count_df <- factor(df$catalogue_label, 
                     levels = reference_catalogue_label) %>% 
    table() %>% 
    data.frame()
  
  colnames(count_df) <- c("mutation_type", "count")
  return(count_df)
}

### Sub-function 2
make_and_count_catalogue_label <- function(cna){
  cna <- add_catalogue_labels(cna)
  context_count <- count_catalogue_label(cna, load_reference_catalogue_label())
  return(context_count)
}

### Master function
convert_txt_to_context <- function(path){
  cna <- read_and_classify_cna(path)
  context_count <- make_and_count_catalogue_label(cna)
  return(context_count)
}
