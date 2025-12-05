variable_sequence_generation <- function(unique_sequences = 10, sequence_length = 16) {
  bases <- c("A", "C", "T", "G")
  variable_seq <- replicate(unique_sequences, paste(sample(bases, size = sequence_length, replace = TRUE), collapse = ""))
  while (length(unique(variable_seq)) < unique_sequences) {
    variable_seq <- unique(variable_seq)
    new_barcodes <- variable_sequence_generation(sequence_length = sequence_length, 
      unique_sequences = unique_sequences - length(unique(variable_seq)))
    variable_seq <- c(variable_seq, new_barcodes)
  }
  return(variable_seq)
}

synth_gen <- function(synthetic_read_layout, n_cells, prep = FALSE) {
  print("Importing read layout and figuring out its order!")
  tidytable::setDTthreads(threads = 1)
  
  synth_layout<-data.table::fread(file = synthetic_read_layout, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  
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
  if(prep == TRUE){
  synth_layout<-synth_layout %>%
    .[,direction := "forward"] %>%
    copy(.) %>%
    .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
    .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"),paste0("rc_", id))] %>%
    .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
    .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
    .[, direction := "reverse"] %>%
    .[order(seq(nrow(.)), decreasing = TRUE),] %>%
    #This part does the reversing
    {rbind(synth_layout, .)} %>%
    # This code annotates the reverse variable elements with the "rc" in the type column.
    {.[,class_id := ifelse(test = {.$direction == "reverse" &! (.$class_id == "poly_a" | .$class_id == "poly_t")},
      yes = paste0("rc_", .$class_id),
      no = .$class_id)]} %>%
    {.[,"order" := seq(1:nrow(.))]} %>%
    #{.[,"direction" := NULL]} %>%
    {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
      yes = "poly_tail",
      no = .$class)]}  
  synth_adapters<-synth_layout$seq[!is.na(synth_layout$seq)]
  names(synth_adapters)<-synth_layout$id[!is.na(synth_layout$seq)]
  }
  
  # synth_template<-paste0(ifelse(!is.na(synth_layout$seq[synth_layout$direction == "forward"]), 
  #   ifelse(synth_layout$class != "poly_tail", synth_layout$seq, paste0("{", synth_layout$class_id, "}")),
  #   paste0("{", synth_layout$class_id, "}")), collapse = "")
  
  output_list<-list()
  for (element in synth_layout$class_id) {
    if(synth_layout[class_id == element, type =="static"]){
      if(!is.na(synth_layout[class_id == element, "seq"])){
        output_list[[element]]<-synth_layout[class_id == element, "seq"] %>% as.character
      }
      if(synth_layout[class_id == element, class == "poly_a"]){
          output_list[[element]]<-paste0(rep("A", 18), collapse = "")
      } else 
        {
          if(synth_layout[class_id == element, class == "poly_t"]){
            output_list[[element]]<-paste0(rep("T", 18), collapse = "")
          }
      }
    }
    if(synth_layout[class_id == element, class == "read"]){
      output_list[[element]]<-variable_sequence_generation(unique_sequences = 1, sequence_length = 250)
    }
    if (synth_layout[class_id == element, class == "barcode"]) {
      if (!is.na(synth_layout[synth_layout$class_id == element, source])) {
        whitelist <- data.table::fread(as.character(synth_layout[class_id == element, source]), header = FALSE, 
          col.names = "whitelist_bcs")
        output_list[[element]] <- sample(whitelist$whitelist_bcs, size = 10000, replace = FALSE)
      } else if (!is.na(synth_layout[synth_layout$class_id == element, expected_length])) {
        output_list[[element]] <- variable_sequence_generation(unique_sequences = 10000,
          sequence_length = synth_layout[synth_layout$class_id == element, expected_length])
      }
    } else {
      if(synth_layout[class_id == element, class == "umi"]){
        output_list[[element]]<-variable_sequence_generation(unique_sequences = 100000, 
          sequence_length = synth_layout[class_id == element, expected_length])
      }
    }
  }
  new_output<-ifelse(test = sapply(output_list, function(x){length(x) == 1}), 
    yes = output_list, no = paste0( "!",names(output_list), "!")) %>%
    unlist(., use.names = FALSE) %>% paste0(collapse = "") %>% 
    str_split(., "!") %>% 
    unlist %>% 
    .[which(sapply(., FUN = function(x){str_length(x)}) > 1)]
}



