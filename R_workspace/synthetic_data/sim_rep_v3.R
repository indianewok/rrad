generate_synthetic_sequences <- function(whitelist, n_cells = 10){
  bases <- c("A", "C", "T", "G")
  forw_primer<-"CTACACGACGCTCTTCCGATCT"
  rev_primer<-"ACTCTGCGTTGATACCACTGCTT"
  bcr<-immuneSIM(number_of_seqs = 10000, species = "hs", receptor = "ig", chain = "h")
  bcr_seqs<-bcr$sequence %>% str_to_upper(.)
  random_barcodes <- bits_to_barcodes(sample(whitelist$whitelist_bcs, size = n_cells, replace = FALSE))
  generate_cell_data<-function(barcode){
    n_seqs<-sample(1:5000, 1)
    bcrseqs<-sample(bcr_seqs, size = n_seqs, replace = TRUE)
    generated_umis<-replicate(10, paste(sample(bases, size = 10, replace = TRUE), collapse = ""))
    cell_data<-data.table(
      id = sapply(sample(generated_umis, n_seqs, replace = TRUE), function(umi) paste(barcode, umi, sep = "_")),
      sequence = paste0(
        forw_primer,
        barcode,
        sample(generated_umis, n_seqs, replace = TRUE),
        paste(rep("T", 35), collapse = ""),
        bcrseqs,
        rev_primer
      ),
      stringsAsFactors = FALSE
    )
    # Randomly reverse complement some sequences
    revcomp_indices<-sample(1:n_seqs, size = round(runif(1, 0.1, 0.3) * n_seqs))
    cell_data[revcomp_indices, "sequence"] <- sapply(cell_data[revcomp_indices, "sequence"], function(seq){
      bajranger::revcomp(seq)
    })
    return(cell_data)
  }
  synthetic_sequences<-pblapply(random_barcodes, generate_cell_data)
  return(synthetic_sequences)
}
