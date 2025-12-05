prep_seq<-function(read_layout_form, external_path_form, create_output_dir = TRUE,
  data_table_nthreads = 1, ...){
  {
    print("Importing read layout and figuring out its order!")
    tidytable::setDTthreads(threads = 1)
    read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)

    if(any(read_layout$class != "forw_primer") | any(read_layout$class != "rev_primer")){
      which(read_layout$type == "static") %>% {read_layout[min(.),class := "forw_primer"]}
      which(read_layout$type == "static") %>% {read_layout[max(.),class := "rev_primer"]}
    }

    read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])

    read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length),
      expected_length := 12]
    read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
    read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]

    if(any(duplicated(read_layout$id))){
      read_layout$id[which(duplicated(read_layout$id)|duplicated(read_layout$id, fromLast = TRUE))]<-which(duplicated(read_layout$id)|duplicated(read_layout$id,
         fromLast = TRUE)) %>%
        read_layout$id[.] %>% paste0(., "_", seq(length(.)))
    }
    if(any(duplicated(read_layout$class))){
      read_layout$class[which(duplicated(read_layout$class)|duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
        sapply(., FUN = function(x){
          out<-which(read_layout$class == x) %>%
            read_layout[., class_id := paste0(class, "_", seq(.))]
          return(out)
        }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
      read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
    } else {
      read_layout$class_id<-read_layout$class
    }
    read_layout<-read_layout %>%
      .[,direction := "forward"] %>%
       copy(.) %>%
      .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
      .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"),paste0("rc_", id))] %>%
      .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
      .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
      .[, direction := "reverse"] %>%
      .[order(seq(nrow(.)), decreasing = TRUE),] %>%
      #This part does the reversing
      {rbind(read_layout, .)} %>%
      # This code annotates the reverse variable elements with the "rc" in the type column.
      {.[,class_id := ifelse(test = {.$direction == "reverse" &! (.$class_id == "poly_a" | .$class_id == "poly_t")},
       yes = paste0("rc_", .$class_id),
       no = .$class_id)]} %>%
      {.[,"order" := seq(1:nrow(.))]} %>%
      #{.[,"direction" := NULL]} %>%
      {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
        yes = "poly_tail",
        no = .$class)]}
    #adapters<-read_layout$seq[grep(pattern = "adapter", x = read_layout$type)]
    #names(adapters)<-read_layout$id[grep(pattern = "adapter", x = read_layout$type)] #list2env this one

    adapters<-read_layout$seq[-which(is.na(read_layout$seq))]
    names(adapters)<-read_layout$class_id[-which(is.na(read_layout$seq))]
    #super quick fix to make the alignment adapter stuff work bc of hardcoded adapter name assumptions downstream, fix

  }#preparing read layout and whatnot
  {
    print("Checking executable paths!")
    path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
      key = "path_type", data.table = TRUE)
    if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
      dir.create(path = path_layout[path_type == "output_dir", actual_path])
    }
    whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
  }#preparing external paths for input
  return(list2env(
    x = list(read_layout = read_layout,
      path_layout = path_layout,
      adapters = adapters,
      whitelist = whitelist), envir = .GlobalEnv))
}
