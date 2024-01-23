### Master function to calculate CN signatures 
### Sera Choi, Jan 2024

source("/Users/schoi/Github_repos/CN_signature_analysis/src/01_convert_txt_to_context.R")
source("/Users/schoi/Github_repos/CN_signature_analysis/src/02_convert_context_to_signature.R")

### Master function to calculate CN signatures from txt file 
calculate_cn_signature <- function(cna_path, ref_path){
  context <- convert_txt_to_context(cna_path)
  signatures <- convert_context_to_signature(ref_path, context)
  return(signatures)
}
