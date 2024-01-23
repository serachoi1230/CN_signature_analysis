### Plot CN signature decomposition plot 
### Sera Choi, Jan 2024


signature_cols <- function(){
  cols <- c("#CAFF70", "#8FBC8F", "#228B22", "#FFA07A", "#FF6A6A", "#FF3030", "#CD3700", "#8B0000",
            "#BBFFFF", "#7EC0EE", "#3A5FCD", "#36648B", "#EEAEEE", "#9F79EE", "#BA55D3", "#551A8B",
            "#1B9E77", "#FFC125", "#FFDAB9", "#8B0A50", "#D2691E", "#CDC8B1", "#C1CDCD", "#8B7355",
            "#F1B6DA")
  return(cols)
}


plot_signature_decomposition <- function(df){
  
  df$label <- "signatures"
  df$signature <- factor(df$signature, levels = df$signature)
  
  p <- ggplot(df, aes(x=label, y=proportion, fill=signature)) +
    geom_bar(stat = "identity", position = "fill", alpha=0.85) +
    scale_fill_manual(values=signature_cols(), drop=FALSE) +
    theme_bw() +
    xlab("")
  
  return(p)

}


