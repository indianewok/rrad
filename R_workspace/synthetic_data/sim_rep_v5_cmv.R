generate_variable_sequence <- function(unique_sequences = 10, sequence_length = 16) {
  bases <- c("A", "C", "T", "G")
  replicate(unique_sequences, paste(sample(bases, size = sequence_length, replace = TRUE), collapse = ""))
}
# Function to generate the synthetic sequences
generate_synth_reads<-function(read_layout_csv, n_cells) {
  # Load the read layout CSV
  read_layout<-data.table::fread(file = read_layout_csv, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  synth_layout[.(type = "static"), class := "forw_primer", on = "type", mult = "first"][.(type = "static"), 
    class := "rev_primer", on = "type", mult = "last"]  
  synth_layout$expected_length[which(!is.na(synth_layout$seq))]<-stringr::str_length(synth_layout$seq[which(!is.na(synth_layout$seq))])
  
  synth_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length),expected_length := 12]
  synth_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
  synth_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
  
  if(any(duplicated(synth_layout$id))){
    synth_layout$id[which(duplicated(synth_layout$id)|duplicated(synth_layout$id, 
      fromLast = TRUE))]<-which(duplicated(synth_layout$id)|duplicated(synth_layout$id,
        fromLast = TRUE)) %>%
      synth_layout$id[.] %>% paste0(., "_", seq(length(.)))
  }
  
  if(any(duplicated(synth_layout$class))){
    synth_layout$class[which(duplicated(synth_layout$class)|duplicated(synth_layout$class, fromLast = TRUE))] %>% unique %>%
      sapply(., FUN = function(x){
        out<-which(synth_layout$class == x) %>%
          synth_layout[., class_id := paste0(class, "_", seq(.))]
        return(out)
      }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
    synth_layout$class_id[which(is.na(synth_layout$class_id))]<-synth_layout$class[which(is.na(synth_layout$class_id))]
  } else {
    synth_layout$class_id<-synth_layout$class
  }
  # Define the default read length
  default_read_length <- 250
  # Initialize lists to collect barcodes and UMIs
  all_barcodes<-vector("list", n_cells)
  all_umis<-vector("list", n_cells)
  # Initialize the list of synthetic sequences
  synth_sequences<-vector("list", n_cells)
  for (i in 1:n_cells) {
    synth_sequences[[i]] <- ""
    all_barcodes[[i]] <- vector("list")
    all_umis[[i]]<-vector("list")
  }
  # Process each type of component based on its class
  for (row in 1:nrow(read_layout)) {
    class<-read_layout$class[row]
    type<-read_layout$type[row]
    print(type)
    seq_length<-ifelse(is.na(read_layout$expected_length[row]), default_read_length, read_layout$expected_length[row])
    if(type == "static "){
      static_seq<-read_layout$seq[row]
      print(paste0("static ", class))
      if(class == "poly_a") {
        static_seq<-paste0(rep("A", 15), collapse = "")
      } else {
        if(class == "poly_t"){
          static_seq<-paste0(rep("T", 15), collapse = "")
        }
      }
      for (i in 1:n_cells) {
        synth_sequences[[i]]<-paste0(synth_sequences[[i]], static_seq)
      }
    }
    if (type == "variable") {
      print(paste0("variable ", class))
      variable_seq <- generate_variable_sequence(unique_sequences = n_cells, sequence_length = seq_length)
      if (class == "barcode") {
        # Sample barcodes from the whitelist if available
        if (!is.na(read_layout$source[row]) && file.exists(read_layout$source[row])) {
          whitelist <- fread(read_layout$source[row], header = FALSE)$V1
          barcodes <- sample(whitelist, size = n_cells, replace = TRUE)
        } else {
          barcodes <- variable_seq
        }
        # Append barcodes to each synthetic sequence and collect them
        for (i in 1:n_cells) {
          synth_sequences[[i]] <- paste0(synth_sequences[[i]], barcodes[i])
          all_barcodes[[i]][[class]] <- barcodes[i]
        }
      } else if (class == "umi") {
        # Append UMIs to each synthetic sequence and collect them
        for (i in 1:n_cells) {
          synth_sequences[[i]] <- paste0(synth_sequences[[i]], variable_seq[i])
          all_umis[[i]][[class]] <- variable_seq[i]
        }
      } else {
        # Append other variable sequences to each synthetic sequence
        for (i in 1:n_cells) {
          synth_sequences[[i]] <- paste0(synth_sequences[[i]], variable_seq[i])
        }
      }
    }
  }
  # Return a list containing synthetic reads, all used barcodes, and all used UMIs
  return(list(synthetic_reads = synth_sequences, barcodes = all_barcodes, umis = all_umis))
}