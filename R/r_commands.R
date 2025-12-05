.datatable.aware = TRUE

whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  cat("Importing a little bit of your whitelist to figure out its type...\n")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1)
  if(ncol(whitelist) == 1){
    data.table::setnames(whitelist, new = "whitelist_bcs")
  } else {
    cat("More than one column detected! Assuming the first one is the one with the barcodes.\n")
    data.table::setnames(whitelist, old = "V1", new = "whitelist_bcs")
  }
  if(class(whitelist$whitelist_bcs) == "character" || 
      !grepl(pattern = "bitlist", x = whitelist_path, fixed = TRUE)){
    cat("Character type detected! Importing the whitelist.\n")
    if(ncol(whitelist) == 1){
      whitelist<-data.table::fread(
        input = whitelist_path,
        header = FALSE,
        data.table = TRUE,
        nThread = 1,
        col.names = "whitelist_bcs"
        )
    } else {
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1)
      whitelist<-whitelist[,1]
      data.table::setnames(whitelist, "whitelist_bcs")
    }
    if(convert != FALSE){
      cat("Converting the whitelist into a more efficient format!\n")
      whitelist$whitelist_bcs<-sequence_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-data.table::setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",tools::file_path_sans_ext(whitelist_path)
          %>% basename(.),
          "_bitlist.csv.gz")
        cat(paste0("Now saving the integer whitelist to ",new_path,"!\n"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      cat("Integer64-type whitelist detected! Now importing...\n")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
        col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      cat(paste0("Your whitelist is ", nrow(whitelist), " barcodes long!\n"))
      print(head(whitelist))
      return(whitelist)
    }
  }
}

write_fastqas<-function(df, fn, type = "fq", append = FALSE, nthreads = 1){
  if(type == "fq"){
    fq_list<-df %>%
      dplyr::mutate(dummy = "+") %>%
      dplyr::select(id, seq, dummy, qual) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fq_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
  }
  if(type == "fa"){
    if(any(grepl("^>", df$id)) == FALSE){
      df$id<-paste0(">",df$id)
      fa_list<-df %>%
        dplyr::select(id, seq) %>%
        as.matrix() %>%
        t() %>%
        as.character() %>%
        tibble::as_tibble()
      data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    } else {
    fa_list<-df %>%
      dplyr::select(id, seq) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    }
  }
}

read_fastqas<-function(fn, type = "fq", full_id = FALSE, ...){
  rfq_single<-function(fn, type, full_id = FALSE, ...){
    ext<-tools::file_ext(fn)
    sysname<-Sys.info()[["sysname"]]
    if (ext == "gz"){
      if(sysname == "Linux"){
        cat_cmd<-"zcat"
      }else if(sysname =="Darwin"){
        cat_cmd<-"gunzip -c"
      }else{
        stop("Unsupported OS")
      }
    } else {
      cat_cmd<-"cat"
    }
    if (type == "fq"){
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | paste - - - - | cut -f1,2,4"),
        col.names = c("id", "seq", "qual"),
        sep = "\t", header = FALSE, quote = "", data.table = TRUE, ...
      )
      res<-res[, .(id = data.table::tstrsplit(id, " ", fixed = TRUE)[[1]], qual, seq)]
      if (full_id) {
        res[, full_id := id]
      }
      return(res)
    }
    if(type == "fa"){
      awk_cmd<-'awk \'/^>/{if (NR>1) printf("\\n"); printf("%s\\t",$0); next} {printf("%s",$0);} END {printf("\\n");}\''
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | {awk_cmd}"),
        col.names = c("id", "seq"),
        sep = "\t", data.table = TRUE, ...
      )
      return(res)
    }
  }
  if (is.list(fn) || length(fn) > 1) {
    out<-pbapply::pblapply(fn, FUN = function(X) {
      rfq_single(fn = X, type = type, full_id = full_id, ...)
    }, ...)
    dt<-data.table::rbindlist(out)
    return(dt)
  } else {
    return(rfq_single(fn = fn, type = type, full_id = full_id, ...))
  }
}

prep_read_layout<-function(read_layout_form, return_env = .GlobalEnv) {
  cat("Importing read layout and figuring out its order!\n")
  tidytable::setDTthreads(threads = 1)
  read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
  
  read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length), expected_length := 12]
  read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
  read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
  
  new_row_top<-data.table::data.table(id = "seq_start", seq = NA, expected_length = 0,
    type = "static", class = "start", direction = "forward", class_id = "seq_start")
  new_row_bottom<-data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, 
    type = "static", class = "stop", direction = "forward", class_id = "seq_stop")
  
  read_layout<-data.table::rbindlist(list(new_row_top, read_layout, new_row_bottom), fill = TRUE)
  
  if(any(read_layout$class != "forw_primer")){
    start_index<-which(read_layout$class == "start")
    if(read_layout$type[start_index+1] == "static"|is.na(read_layout$class[start_index+1])){
      read_layout[start_index+1, class := "forw_primer"]
    }
  }
  
  if(any(read_layout$class != "rev_primer")){
    stop_index<-which(read_layout$class == "stop")
    if(read_layout$type[stop_index-1] == "static"|is.na(read_layout$class[stop_index-1])){
      read_layout[stop_index-1, class := "rev_primer"]
    }
  }
  
  if (any(duplicated(read_layout$id))) {
    read_layout$id[which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE))]<-
      which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE)) %>%
      read_layout$id[.] %>% paste0(., "_", seq(length(.)))
  }
  if (any(duplicated(read_layout$class))) {
    read_layout$class[which(duplicated(read_layout$class) | duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
      sapply(., FUN = function(x){
        out<-which(read_layout$class == x) %>%
          read_layout[., class_id := paste0(class, "_", seq(.))]
        return(out)
      }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
    read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
  } else {
    read_layout$class_id <- ifelse(test = (read_layout$class == "start" | read_layout$class == "stop"), 
      yes = read_layout$class_id, 
      no = read_layout$class)
  }
  
  read_layout<-read_layout %>%
    .[, direction := "forward"] %>%
    data.table::copy(.) %>%
    .[!class %in% c("start", "stop")] %>%
    .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
    .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"), paste0("rc_", id))] %>%
    .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
    .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
    .[, direction := "reverse"] %>%
    .[order(seq(nrow(.)), decreasing = TRUE), ] %>%
    {data.table::rbindlist(list(data.table::data.table(id = "seq_start", seq = NA, expected_length = 0, type = "static",
      class = "start", direction = "reverse", class_id = "seq_start"),
      ., data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, type = "static",
        class = "stop", direction = "reverse", class_id = "seq_stop")),
      fill = TRUE)} %>%
    {rbind(read_layout, .)} %>%
    {.[, class_id := ifelse(test = {.$direction == "reverse" &!
        (.$class_id == "poly_a" | .$class_id == "poly_t")},# | .$class_id == "seq_start" | .$class_id == "seq_stop")},
      yes = paste0("rc_", .$class_id),
      no = .$class_id)]} %>%
    {.[, "order" := seq(1:nrow(.))]} %>%
    {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
      yes = "poly_tail",
      no = .$class)]}
  
  adapters <- read_layout$seq[-which(is.na(read_layout$seq))]
  names(adapters) <- read_layout$class_id[-which(is.na(read_layout$seq))]
  return(list2env(x = list(read_layout = read_layout, adapters = adapters), envir = return_env))
}

stat_collector<-function(df, read_layout, mode = "stats"){
  forward_adapters<-read_layout[type %in% c("static") & direction == "forward" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  reverse_adapters<-read_layout[type %in% c("static") & direction == "reverse" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  df<-tidytable::as_tidytable(df)
  #changed best_edit_distance from == 0 to <= 1. 
  #not the best fix, as it leaves room for perfect concatenates with no errors or something weird like that
  #but it also adds for some more sequences that might otherwise be discarded.
  forward_zero_ids<-df %>%
    tidytable::filter(query_id %in% forward_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
    tidytable::pull(id)
  reverse_zero_ids<-df %>%
    tidytable::filter(query_id %in% reverse_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
    tidytable::pull(id)
  if(mode == "stats"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    misalignment_thresholds<-rbind(stats_forward, stats_reverse)
    return(misalignment_thresholds)
  }
  if(mode == "graphics"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id)
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id)
    total_stats<-data.table::rbindlist(list(stats_forward, stats_reverse))
    return(total_stats)
  }
}

prep_external_path<-function(external_path_form){
  print("Checking executable paths!")
  path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
    key = "path_type", data.table = TRUE)
  if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
    dir.create(path = path_layout[path_type == "output_dir", actual_path])
  }
  whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
  return(list2env(x = list(path_layout = path_layout, whitelist = whitelist), envir = .GlobalEnv))
}

read_paf<-function(path, ...){
  df<-data.table::fread(input = path, 
    fill = Inf, sep = "\t", strip.white = FALSE, header = FALSE, na.strings = "")
  df[df == '']<-NA
  ids<-c("id","ql","qs","qe","qr","sn","sl","ss","se","mm","mt","mq")
  data.table::setnames(df, c(ids,
    df[,ncol(df):ncol(df),] %>% {
      which(!is.na(.))} %>% 
      df[.,13:ncol(df),] %>% 
      {substr(unlist(.), start = 1, stop = 2)} %>% 
      unique))
  for (col in (length(ids) + 1):ncol(df)) {
    tag_prefixes<-substr(df[[col]], start = 1, stop = 2)
    expected_tag<-colnames(df)[col]
    mismatch_indices<-which(tag_prefixes != expected_tag & !is.na(tag_prefixes))
    if (length(mismatch_indices) > 0) {
      # Shifting columns to the right for mismatched rows
      df[mismatch_indices, (col + 1):ncol(df) := .SD, .SDcols = col:(ncol(df) - 1)]
      df[mismatch_indices, (col) := NA]  # Set current column to NA for mismatched rows
    }
    df[, (expected_tag) := gsub(pattern = paste0(expected_tag,":.:"), replacement = "", x = get(expected_tag))]
  }
  for (col in names(df)) {
    numeric_conversion <- suppressWarnings(as.numeric(df[[col]]))
    if (!any(is.na(numeric_conversion)) || all(is.na(numeric_conversion) == is.na(df[[col]]))) {
      df[[col]] <- numeric_conversion
    }
  }
  return(df)
}

aggc<-function(){
  gc()
  sysname<-Sys.info()[["sysname"]]
  if(sysname == "Linux"){
    mallinfo::malloc.trim()
  }
}

aggc <- function() {
  sysname <- Sys.info()[["sysname"]]
  if (sysname == "Linux") {
    mallinfo::malloc.trim()
  } 
  if (sysname == "Darwin") {  
    # macOS
    system("purge", ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
}


process_barcodes<-function(barcode_path, read_layout) {
  barcodes<-data.table::fread(barcode_path, col.names = c("seq", "count","concatenate_count"))
  data.table::setkey(read_layout, "class_id")
  
  dir_path<-paste0(dirname(barcode_path),"/")
  barcode_id<-basename(tools::file_path_sans_ext(barcode_path, compression = TRUE))
  
  expected_length<-read_layout[barcode_id, expected_length]
  barcodes[, filtered := NA_character_]
  barcodes[(stringr::str_length(seq) != expected_length), filtered := "length_filter"]
  rescue_filter<-ceiling(expected_length * 0.5)
  patterns<-sprintf("%s{%d,}", c("A", "C", "T", "G"), rescue_filter + 1)
  poly_runs<-data.table::rbindlist(lapply(patterns, function(pattern) {
    barcodes[is.na(filtered) & grepl(pattern, seq), .(seq, count, filtered = "poly_run_filtered")]
  }))
  if (nrow(poly_runs) > 0){
    barcodes[poly_runs, on = .(seq), filtered := i.filtered]
  }
  static_seqs<-read_layout[type == "static" & !is.na(expected_length) & expected_length > 0 & class != "poly_tail", seq]
  static_qgrams<-colnames(stringdist::qgrams(static_seqs, q = expected_length, useNames = TRUE))
  static_qgrams_mutations<-static_qgrams %>% 
    sequence_to_bits(.) %>% 
    generate_recursive_quaternary_mutations(., mutation_rounds =  2, sequence_length = 16) %>%
    bits_to_sequence(.)
  static_qgrams<-append(static_qgrams, static_qgrams_mutations)
  
  barcodes[is.na(filtered) & seq %in% static_qgrams, filtered := "adapter_segment_filtered"]
  whitelist_path<-read_layout[barcode_id, whitelist]
  cat(paste0("The whitelist path is ", whitelist_path,"\n"))
  whitelist<-NULL
  whitelist_path<-switch(
    whitelist_path,
    "10x_3v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"), #don't know if true
    "10x_3v2" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
    "10x_3v3" = system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"),
    "10x_3v3.1" = system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"),
    "10x_3HTv3.1" = system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"),
    "10x_3v4" = system.file(package = "rad", "extdata", "3M-3pgex-may-2023.txt_bitlist.csv.gz"),
    "10x_5v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
    "10x_5v2" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
    "10x_5HTv2" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
    "10x_5v3" = system.file(package = "rad", "extdata", "3M-5pgex-jan-2023.txt_bitlist.csv.gz"),
    "10x_Vis_V1" = system.file(package = "rad", "extdata", "visium-v1_v2_bitlist.csv.gz"),
    "10x_Vis_V2" = system.file(package = "rad", "extdata", "visium-v1_v2_bitlist.csv.gz"),
    "10x_Vis_V3" = system.file(package = "rad", "extdata", "visium-v3_v4_bitlist.csv.gz"),
    "10x_Vis_V4" = system.file(package = "rad", "extdata", "visium-v3_v4_bitlist.csv.gz"),
    "10x_Vis_V5" = system.file(package = "rad", "extdata", "visium-v5_bitlist.csv.gz"),
    default = whitelist_path
  )
  cat(paste0("The whitelist path is ", whitelist_path, "\n"))
  
  if(!is.na(whitelist_path) && file.exists(whitelist_path)){
    ext<-tools::file_ext(whitelist_path)
    convert<-ext == "csv" || ext == "txt"
    whitelist<-whitelist_importer(whitelist_path, convert = convert)
  }
  
  barcodes[, pois_dist := stats::ppois(q = count, lambda = mean(count))]
  
  barcodes[is.na(filtered), ':='(
    int64_seq = sequence_to_bits(seq),
    int64_rcseq = sequence_to_bits(revcomp(seq))
  )]
  if (!is.null(whitelist)){
    barcodes[(is.na(filtered) & stringr::str_length(seq) == expected_length), filtered := ifelse(
      int64_seq %in% whitelist$whitelist_bcs | int64_rcseq %in% whitelist$whitelist_bcs,
      "whitelist_barcode", "barcode_to_correct"
    )]
  }
  whitelist_output<-density_estimator_v13(barcode_df = barcodes)
  if(exists("kde_density", where = whitelist_output)){
    kde_density<-ggpubr::ggpar(
      whitelist_output$kde_density,
      title = paste0(basename(dirname(dir_path)), " KDE Density with Threshold",
        subtitle = paste0("\nThreshold: ", 
          whitelist_output$threshold,
          "\nPredicted Count: ",
          length(which(whitelist_output$results$kde_threshold)))),
      legend = "none", font.x = 14, font.y = 14, font.main = 14)
    
    ggplot2::ggsave(
      filename = paste0(dir_path, barcode_id, "_kde_density.pdf"), 
      plot =  kde_density,
      width = 11, height = 8.5, 
      unit = "in"
    )
  }
  data.table::fwrite(x = barcodes, file = barcode_path, col.names = TRUE)
}

correct_barcodes<-function(input_file,
                           verbose = FALSE,
                           nthreads = 1,
                           mutation_rounds = 3,
                           max_shift = 3,
                           chunk_size = 50000,
                           read_layout = read_layout
                           ){
  barcode<-data.table::fread(
    input_file, key = "filtered")
  data.table::setkey(read_layout, "class_id")
  compressed<-grepl(
    pattern = ".gz", 
    x = input_file, fixed = TRUE)
  barcode_id<-basename(tools::file_path_sans_ext(input_file, compression = compressed))
  
  dir_path<-paste0(dirname(input_file),"/")
  input_fastq_path<-paste0(dirname(dir_path),"/demuxed_reads.fastq.gz")
  filtered_fastq_path<-paste0(dirname(dir_path),"/filtered_reads.fastq.gz")
  unfiltered_fastq_path<-paste0(dirname(dir_path),"/unfiltered_reads.fastq.gz")
  barcode_counts_output_path<-paste0(dir_path, barcode_id, "_counts.csv")
  verbose_barcode_output_path<-paste0(dir_path, barcode_id, "_correction_counts.csv")
  barcode_mutation_data_path<-paste0(dir_path, barcode_id, "_mutation_counts.csv.gz")
  expected_length<-read_layout[barcode_id, expected_length]
  
  true_barcodes<- bit64::as.integer64(barcode["pois_validated_barcode"]$int64_seq)
  invalid_barcodes<-bit64::as.integer64(barcode["barcode_to_correct"]$int64_seq)
  true_counts<-as.vector(barcode["pois_validated_barcode"]$count)
  invalid_counts<-as.vector(barcode["barcode_to_correct"]$count)
  
  correction_map<-barcode_correction_v3(
    true_barcodes = true_barcodes,
    invalid_barcodes = invalid_barcodes,
    true_counts = true_counts,
    invalid_counts = invalid_counts,
    mutation_rounds = mutation_rounds,
    sequence_length = expected_length,
    nthread = nthreads,
    max_shift = max_shift, 
    verbose = FALSE)
  data.table::fwrite(correction_map, file = barcode_mutation_data_path)
  correction_map<-data.table::fread(input = barcode_mutation_data_path, na.strings = "")
  
  correction_map<-correction_map[which(correction_map$resolved_status == TRUE),]
  
  #this guy really runs slow, but i think it's chunk size dependent
  fastq_correction(
    input_file = input_fastq_path,
    output_file = filtered_fastq_path,
    incorrect_bcs = correction_map$incorrect_barcodes,
    correct_bcs = correction_map$corrected_barcodes,
    barcode_id = barcode_id,
    nthreads = nthreads,
    chunk_size = chunk_size,
    downsample = -1
  )
}

process_sig<-function(file_path,
                      output_prefix,
                      chunk_size = 50000,
                      max_sequences = -1,
                      compress = TRUE,
                      output_type = "summary",
                      nthreads = 1){
  tabulate_sigs(file_path = file_path, 
                output_prefix = output_prefix, 
                chunk_size = chunk_size,
                max_sequences = max_sequences, 
                compress = compress, 
                nthreads = nthreads,
                output_type = output_type)
  summary_results<-data.table::fread(file = paste0(output_prefix,"/summary_table.csv",  
                                                   ifelse(test = compress, yes = ".gz", no = "")),
                                     header = TRUE)
  valid_reads<-sum(summary_results$unique_id_count[which(summary_results$direction == "F"|summary_results$direction == "R")])
  total_reads<-sum(summary_results$unique_id_count)
  summary_plot<-ggpubr::ggbarplot(data = summary_results, 
                                  x = "direction", 
                                  y = "unique_id_count",
                                  color = "concatenate",
                                  fill = "concatenate", 
                                  title = "SigString Summary", 
                                  legend = "right", 
                                  legend.title = "Concatenate Status", 
                                  xlab = "Read (or Filter) Type",
                                  ylab = "Number of Reads", 
                                  order = c("F", 
                                            "R", 
                                            "missing_read_or_barcode_undecided", 
                                            "invalid_barcode_length_undecided",
                                            "invalid_read_length_undecided"), 
                                  subtitle = paste0(
                                    valid_reads, 
                                    " valid reads out of a total ", 
                                    total_reads, " total reads."), 
                                  x.text.angle = 45, palette = "npj")
  ggplot2::ggsave(filename = paste0(output_prefix,"summary_barplot.pdf"), 
                  plot = summary_plot, 
                  device = "pdf", 
                  width = 8, 
                  height = 11, 
                  units = "in", 
                  dpi = "retina")
  results<-data.table::data.table(unique_id_count = paste0(valid_reads, " valid reads captured from ", total_reads, " total reads."))
  summary_results<-data.table::rbindlist(l = list(summary_results, results), fill = TRUE)
  data.table::fwrite(summary_results, paste0(output_prefix,"/summary_table.csv", ifelse(test = compress, yes = ".gz", no = "")))
}

rad_run<-function(
    fastq_file_or_directory_path, 
    read_layout_path, 
    output_directory_path, 
    compress = TRUE,
    generate_sigstring_diagnostics = FALSE,
    sigstring_diagnostic_verbose = FALSE,
    misalignment_threshold_path = NULL,
    tabulated_sigstring_count = -1,
    chunk_size = 50000,
    mutation_rounds = 3,
    nthreads = 1){
  start_time = Sys.time()
  if(!dir.exists(output_directory_path)){
    dir.create(output_directory_path)
  }
  prep_read_layout(read_layout_form = read_layout_path)
  if(!is.null(misalignment_threshold_path)){
    misalignment_threshold<-data.table::fread(misalignment_threshold_path, header = FALSE)
  } else {
    out<-generate_misalignment_info(read_input_path = fastq_file_or_directory_path, read_layout_path = read_layout_path)
    data.table::fwrite(out$misalignment_threshold, file = paste0(output_directory_path,"/misalignment_threshold.csv"))
    data.table::fwrite(out$masked_misalignment_threshold, file = paste0(output_directory_path,"/masked_misalignment_threshold.csv"))
  }
  
  sigstream(input_path = fastq_file_or_directory_path, 
            nthreads = nthreads, 
            compress = compress, 
            write_fastq = TRUE, 
            min_length = 100,
            chunkSize = chunk_size,
            outputPath = paste0(output_directory_path, "/demuxed_reads.fastq"),
            adapters = out$masked_adapters,
            misalignment_threshold = out$masked_misalignment_threshold,
            read_layout = read_layout,
            sigstringsFilePath = paste0(output_directory_path, "/sigstrings.txt")
  )
  gc()
  if(generate_sigstring_diagnostics){
    if(sigstring_diagnostic_verbose){
      output_type<-"diagnostic"
    } else {
      output_type<-"summary"
    }
    process_sig(
      file_path = paste0(output_directory_path, "/sigstrings.txt", ifelse(test = compress, yes = ".gz", no = "")),
      output_prefix = output_directory_path,
      max_sequences = tabulated_sigstring_count,
      chunk_size = chunk_size,
      compress = compress,
      nthreads = nthreads,
      output_type = output_type
    )
    gc()
  } else {
    cat("Generate sigstring diagnostics set to FALSE--skipping sigstring tabulation.\n")
  }
  dir.create(path = paste0(output_directory_path,"/variable_seqs/"))
  cat("Tabulating barcodes...\n")
  extractable_elements<-read_layout[type == "variable" & class != "read" & direction == "forward", class_id]
  tabulate_variable_sequences(
    input_files = paste0(output_directory_path, "/demuxed_reads.fastq", 
                         ifelse(test = compress, yes = ".gz", no = "")), 
    identity_elements = c(extractable_elements),
    output_prefix = paste0(output_directory_path,"/variable_seqs/"), 
    compress = TRUE)
  barcode_files<-list.files(path = paste0(output_directory_path,"/variable_seqs/"), full.names = TRUE, pattern = "barcode")
  #need to figure out how to loop over multiple barcodes and rewrite the output file
  out<-lapply(X = barcode_files, FUN = function(X){
    process_barcodes(barcode_path = X, read_layout = read_layout)
    #think I need to do it here, because the input file gets chosen internally in correct_barcodes
    correct_barcodes(input_file = X, read_layout = read_layout, nthreads = nthreads, chunk_size = chunk_size)
  })
  end_time = Sys.time()
  cat(paste0("Finished running! Total run time: ", difftime(end_time, start_time, units='mins'), " minutes.\n"))
}

generate_synthetic_reads<-function(
  read_layout_path,
  num_cells,
  read_length,
  num_reads = c(1,50)){
  # Prep layout and set key for class_id
  prep_read_layout(read_layout_path, return_env = parent.frame())
  data.table::setkey(read_layout, "class_id")
  read_layout <- read_layout[direction == "forward" & (expected_length > 0 | class == "read"), ]
  data.table::setkey(read_layout, "class_id")
  
  # Initialize layout as a data.table with num_cells rows
  layout <- data.table(cell_id = as.character(1:num_cells))
  
  # Pre-generate the number of reads per cell
  num_reads_per_cell <- sample(seq(num_reads[1], num_reads[2]), num_cells, replace = TRUE)
  
  # Now we generate the static elements (primers, poly_t, etc.)
  for (unique_class in unique(read_layout$class)) {
    # Handle static elements (e.g., primers, poly_t)
    if (unique_class %in% c("forw_primer", "poly_tail","tso", "rev_primer")) {
      print(paste("Dealing with static element:", unique_class))
      for (static_id in unique(read_layout[class == unique_class, class_id])) {
        static_seq<-read_layout[class_id == static_id, seq]
        if(unique_class == "poly_tail"){
          print(static_seq)
          min_reps<-as.numeric(DescTools::StrExtractBetween(x = static_seq, left = "\\{", right = ","))
          print(min_reps)
          poly_seq<-ifelse(
            grepl("A", static_seq, fixed = TRUE), "A",
            ifelse(grepl("C", static_seq, fixed = TRUE), "C",
              ifelse(grepl("T", static_seq, fixed = TRUE), "T",
                ifelse(grepl("G", static_seq, fixed = TRUE), "G",
                  stop("Invalid base detected")))))
          static_seq<-paste0(rep(x = poly_seq, times = ceiling(min_reps*1.5)), collapse = "")
          print(static_seq)
        }
        layout[, paste0(static_id) := static_seq]
      }
    }
    
    # Handle barcodes (per cell)
    if (unique_class == "barcode") {
      print("Dealing with barcode")
      for (barcode_id in unique(read_layout[class == "barcode", class_id])) {
        layout[, paste0(barcode_id) := {
          whitelist_path <- switch(
            read_layout[barcode_id, whitelist],
            "10x_3v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
            "10x_3v3" = system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"),
            "10x_3v4" = system.file(package = "rad", "extdata", "3M-3pgex-may-2023.txt_bitlist.csv.gz"),
            "10x_5v1" = system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"),
            "10x_5v3" = system.file(package = "rad", "extdata", "3M-5pgex-jan-2023.txt_bitlist.csv.gz"),
            NA_character_
          )
          whitelist <- whitelist_importer(whitelist_path = whitelist_path)
          barcode_sample <- sample(whitelist$whitelist_bcs, size = num_cells, replace = FALSE)
          barcode_length <- read_layout[barcode_id, expected_length]
          bits_to_sequence(barcode_sample, sequence_length = barcode_length)
        }]
      }
      print("done with barcode")
    }
  }
  # Now clone the rows for each barcode according to num_reads_per_cell
  layout<-layout[rep(seq_len(num_cells), num_reads_per_cell)]
  # Now generate UMIs and reads for each read (post-cloning of rows)
  for (unique_class in unique(read_layout$class)) {
    
    # Handle UMIs
    if (unique_class == "umi") {
      print("Dealing with UMI")
      for (umi_id in unique(read_layout[class == "umi", class_id])) {
        umi_length <- as.numeric(read_layout[class_id == umi_id, expected_length])
        if (is.na(umi_length) | umi_length <= 0) stop(paste("Invalid UMI length for", umi_id))
        
        # Generate a unique UMI for each read within each cell
        layout[, paste0(umi_id) := unlist(lapply(1:num_cells, function(cell) {
          replicate(num_reads_per_cell[cell], paste0(sample(c("A", "T", "C", "G"), umi_length, replace = TRUE), collapse = ""))
        }))]
      }
    }
    
    # Handle reads
    if (unique_class == "read") {
      print("Dealing with reads")
      for (read_id in unique(read_layout[class == "read", class_id])) {
        # Generate a unique read for each read within each cell
        layout[, paste0(read_id) := unlist(lapply(1:num_cells, function(cell) {
          replicate(num_reads_per_cell[cell], paste0(sample(c("A", "T", "C", "G"), read_length, replace = TRUE), collapse = ""))
        }))]
      }
    }
  }
  
  #   # Combine multiple barcode and UMI columns if needed
  layout[, combined_barcode := do.call(paste, c(.SD, sep = "_")), .SDcols = grep("barcode", names(layout), value = TRUE)]
  layout[, combined_umi := do.call(paste, c(.SD, sep = "_")), .SDcols = grep("umi", names(layout), value = TRUE)]
  # Create cell_id by combining barcode, UMI, and read number
  layout[, cell_id := paste0(
    combined_barcode, "_",
    combined_umi, "_",
    sequence(.N)  # Reset numbering within each group
  ), by = .(combined_barcode)]  # Group by combined_barcode
  
  # Reshape and return the final read layout
  ordered_class_ids<-read_layout[order(order), class_id]
  reshaped_layout<-tidytable::unite(layout, final_read, all_of(ordered_class_ids), sep = "")
  reshaped_layout <- reshaped_layout[, .(id = cell_id, seq = final_read)]
  return(reshaped_layout)
}

density_estimator_v12<-function(barcode_df){
  data.table::setkey(barcode_df, "filtered")
  # Calculate ncpm and log1p_ncpm
  barcode_df[, ncpm := (count / sum(count)) * 1e6] 
  barcode_df[, log1p_ncpm := log1p(ncpm)]
  barcode_df[, ncpm_pois := stats::ppois(q = ncpm, lambda = mean(ncpm))]
  # Subset the data based on conditions
  data_subset<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, log1p_ncpm]
  # Perform GMM clustering
  gmm_scan<-mclust::densityMclust(data = data_subset, plot = FALSE, modelNames = c("V", "E"))
  density_output<-density(data_subset)
  plot_df<-data.table::data.table(x = density_output$x, y = density_output$y)
  # Find peaks and valleys
  peaks<-pracma::findpeaks(density_output$y)
  valleys<-pracma::findpeaks(-density_output$y)
  # If more than 3 peaks, drop the lowest peak
  if (nrow(peaks) > 3) {
    min_peak_index<-which.min(peaks[, 1])
    peaks<-peaks[-min_peak_index, , drop = FALSE]
  }
  # If still more than 3 peaks, iterate over bandwidths to reduce peaks to 2
  bw<-density_output$bw
  iter<-1
  while(nrow(peaks) > 2 && iter <= 50){
    bw<-bw + 0.05  # Increment bandwidth
    density_output<-density(data_subset, bw = bw)
    peaks<-pracma::findpeaks(density_output$y)
    valleys<-pracma::findpeaks(-density_output$y)
    iter<-iter + 1
  }
  if(nrow(peaks) == 2){
    peak1_position<-peaks[1, 2]
    peak2_position<-peaks[2, 2]
    valley_in_between<-valleys[valleys[, 2] > peak1_position & valleys[, 2] < peak2_position, ]
    if (!is.matrix(valley_in_between)) {
      valley_in_between<-matrix(valley_in_between, nrow = 1)
    }
    if (nrow(valley_in_between) > 0) {
      valid_valley<-valley_in_between[1, 2]
    } else {
      stop("No valley found between the two peaks.")
    }
    # Threshold is the x-position of the valid valley
    threshold<-density_output$x[valid_valley]
    # Calculate the x positions of the peaks
    peak_one<-density_output$x[peaks[1, 2]]
    peak_two<-density_output$x[peaks[2, 2]]
    # Classify points based on threshold
    plot_df[, cluster := ifelse(x >= threshold, "Above Threshold", "Below Threshold")]
    # Extract GMM results
    best_gmm<-gmm_scan
    uncertainty_percentage<-round(best_gmm$uncertainty, digits = 4) * 100
    best_cluster<-best_gmm$classification
    kde_threshold<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95]$log1p_ncpm >= threshold
    cluster_threshold<-uncertainty_percentage <= 5
    results<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, 
      .(seq, count, concatenate_count, log1p_ncpm, 
        uncertainty_percentage, best_cluster,
        kde_threshold, cluster_threshold)]
    kde_density_plot<-ggplot2::ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "blue") +  
      ggplot2::geom_vline(xintercept = peak_one, linetype = "dashed", color = "red") + 
      ggplot2::geom_vline(xintercept = peak_two, linetype = "dashed", color = "red") +  
      labs(title = "KDE Density with Threshold", x = "log1p_ncpm", y = "Density") +
      ggplot2::scale_color_manual(values = c("Below Threshold" = "black", "Above Threshold" = "green")) +
      ggplot2::theme_minimal()
    
    data.table::setkey(barcode_df, "seq")
    data.table::setkey(results, "kde_threshold")
    barcode_df[results[J(TRUE)]$seq, filtered:="pois_validated_barcode"]
    # Return the results, plots, and GMM model
    return(list(
      threshold = threshold,
      density = density_output,
      results = results,
      kde_density = kde_density_plot,
      best_gmm = best_gmm
    ))
  }
  if(nrow(peaks)  < 2 && is.null(valleys)){
    xrange<-extendrange(gmm_scan$data, f = 0.1)
    eval.points<-seq(from = xrange[1], to = xrange[2], length = 1000)
    gmm_density<-predict(gmm_scan, newdata = eval.points, what = "dens")
    ggdensity_plot<-data.table::data.table(x = eval.points, y = gmm_density)
    
    peaks<-pracma::findpeaks(ggdensity_plot$y)
    valleys<-pracma::findpeaks(-ggdensity_plot$y)
    
    if (nrow(peaks) >= 3) {
      min_peak_index<-which.min(peaks[, 1])
      peaks<-peaks[-min_peak_index, , drop = FALSE]
    }
    
    if(nrow(peaks) == 2){
      peak1_position<-peaks[1, 2]
      peak2_position<-peaks[2, 2]
      valley_in_between<-valleys[valleys[, 2] > peak1_position & valleys[, 2] < peak2_position, ]
      if (!is.matrix(valley_in_between)) {
        valley_in_between<-matrix(valley_in_between, nrow = 1)
      }
      if (nrow(valley_in_between) > 0) {
        valid_valley<-valley_in_between[1, 2]
      } else {
        stop("No valley found between the two peaks.")
      }
      # Threshold is the x-position of the valid valley
      threshold<-ggdensity_plot$x[valid_valley]
      # Calculate the x positions of the peaks
      peak_one<-ggdensity_plot$x[peaks[1, 2]]
      peak_two<-ggdensity_plot$x[peaks[2, 2]]
      # Classify points based on threshold
      ggdensity_plot[, cluster := ifelse(x >= threshold, "Above Threshold", "Below Threshold")]
      # Extract GMM results
      best_gmm<-gmm_scan
      uncertainty_percentage<-round(best_gmm$uncertainty, digits = 4) * 100
      best_cluster<-best_gmm$classification
      kde_threshold<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95]$log1p_ncpm >= threshold
      cluster_threshold<-uncertainty_percentage <= 5
      results<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, 
        .(seq, count, concatenate_count, log1p_ncpm, 
          uncertainty_percentage, best_cluster,
          kde_threshold, cluster_threshold)]
      kde_density_plot<-ggplot2::ggplot(ggdensity_plot, aes(x = x, y = y, color = cluster)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "blue") +  
        ggplot2::geom_vline(xintercept = peak_one, linetype = "dashed", color = "red") + 
        ggplot2::geom_vline(xintercept = peak_two, linetype = "dashed", color = "red") +  
        labs(title = "KDE Density with Threshold", x = "log1p_ncpm", y = "Density") +
        ggplot2::scale_color_manual(values = c("Below Threshold" = "black", "Above Threshold" = "green")) +
        ggplot2::theme_minimal()
      
      data.table::setkey(barcode_df, "seq")
      data.table::setkey(results, "kde_threshold")
      barcode_df[results[J(TRUE)]$seq, filtered:="pois_validated_barcode"]
      # Return the results, plots, and GMM model
      return(list(
        threshold = threshold,
        density = density_output,
        results = results,
        kde_density = kde_density_plot,
        best_gmm = best_gmm
      ))
    }
  }
}


density_estimator_v13<-function(barcode_df){
  data.table::setkey(barcode_df, "filtered")
  # Calculate ncpm and log1p_ncpm
  barcode_df[, ncpm := (count / sum(count)) * 1e6] 
  barcode_df[, log1p_ncpm := log1p(ncpm)]
  barcode_df[, ncpm_pois := stats::ppois(q = ncpm, lambda = mean(ncpm))]
  
  # Subset the data based on conditions
  data_subset<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, log1p_ncpm]
  
  # Perform GMM clustering
  gmm_scan<-mclust::densityMclust(data = data_subset, plot = FALSE, modelNames = c("V", "E"))
  density_output<-density(data_subset)
  plot_df<-data.table::data.table(x = density_output$x, y = density_output$y)
  
  # Find peaks and valleys
  peaks<-pracma::findpeaks(density_output$y)
  valleys<-pracma::findpeaks(-density_output$y)
  
  # Function to classify and return results
  classify_and_return<-function(threshold, peaks, density_output, plot_df, barcode_df) {
    peak_one<-density_output$x[peaks[1, 2]]
    peak_two<-density_output$x[peaks[2, 2]]
    
    # Classify points based on threshold
    plot_df[, cluster := ifelse(x >= threshold, "Above Threshold", "Below Threshold")]
    
    # Extract GMM results
    uncertainty_percentage<-round(gmm_scan$uncertainty, digits = 4) * 100
    best_cluster<-gmm_scan$classification
    
    kde_threshold<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95]$log1p_ncpm >= threshold
    cluster_threshold<-uncertainty_percentage <= 5
    
    results<-barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, 
      .(seq, count, concatenate_count, log1p_ncpm, 
        uncertainty_percentage, best_cluster, kde_threshold, cluster_threshold)]
    
    kde_density_plot<-ggplot2::ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "blue") +  
      ggplot2::geom_vline(xintercept = peak_one, linetype = "dashed", color = "red") + 
      ggplot2::geom_vline(xintercept = peak_two, linetype = "dashed", color = "red") +  
      labs(title = "KDE Density with Threshold", x = "log1p_ncpm", y = "Density") +
      ggplot2::scale_color_manual(values = c("Below Threshold" = "black", "Above Threshold" = "green")) +
      ggplot2::theme_minimal()
    
    data.table::setkey(barcode_df, "seq")
    data.table::setkey(results, "kde_threshold")
    barcode_df[results[J(TRUE)]$seq, filtered := "pois_validated_barcode"]
    
    return(list(
      threshold = threshold,
      density = density_output,
      results = results[J(TRUE),],
      kde_density = kde_density_plot,
      best_gmm = gmm_scan
    ))
  }
  
  # If more than 3 peaks, drop the lowest peak
  if (nrow(peaks) > 3) {
    min_peak_index<-which.min(peaks[, 1])
    peaks<-peaks[-min_peak_index, , drop = FALSE]
  }
  
  # Iterate over bandwidths to reduce peaks to 2 or 1
  bw<-density_output$bw
  iter<-1
  while (nrow(peaks) > 2 && iter <= 50) {
    bw<-bw + 0.05  # Increment bandwidth
    density_output<-density(data_subset, bw = bw)
    peaks<-pracma::findpeaks(density_output$y)
    valleys<-pracma::findpeaks(-density_output$y)
    iter<-iter + 1
  }
  
  # If there are 2 peaks, return threshold and plot
  if (nrow(peaks) == 2) {
    peak1_position<-peaks[1, 2]
    peak2_position<-peaks[2, 2]
    valley_in_between<-valleys[valleys[, 2] > peak1_position & valleys[, 2] < peak2_position, ]
    valley_in_between<-as.matrix(t(valley_in_between))
    if (nrow(valley_in_between) > 0) {
      valid_valley<-valley_in_between[1, 2]
      threshold<-density_output$x[valid_valley]
      return(classify_and_return(threshold, peaks, density_output, plot_df, barcode_df))
    } else {
      stop("No valley found between the two peaks.")
    }
  }
  
  # If valleys are null or fewer than 2 peaks, fall back to GMM
  if (nrow(peaks) < 2 || is.null(valleys)) {
    xrange<-extendrange(gmm_scan$data, f = 0.1)
    eval.points<-seq(from = xrange[1], to = xrange[2], length = 1000)
    gmm_density<-predict(gmm_scan, newdata = eval.points, what = "dens")
    
    ggdensity_plot<-data.table::data.table(x = eval.points, y = gmm_density)
    
    peaks<-pracma::findpeaks(ggdensity_plot$y)
    valleys<-pracma::findpeaks(-ggdensity_plot$y)
    
    # Apply the same logic to the GMM density
    if (nrow(peaks) >= 3) {
      min_peak_index<-which.min(peaks[, 1])
      peaks<-peaks[-min_peak_index, , drop = FALSE]
    }
    
    if (nrow(peaks) == 2) {
      peak1_position<-peaks[1, 2]
      peak2_position<-peaks[2, 2]
      valley_in_between<-valleys[valleys[, 2] > peak1_position & valleys[, 2] < peak2_position, ]
      
      if (nrow(valley_in_between) > 0) {
        valid_valley<-valley_in_between[1, 2]
        threshold<-ggdensity_plot$x[valid_valley]
        return(classify_and_return(threshold, peaks, ggdensity_plot, plot_df, barcode_df))
      } else {
        stop("No valley found between the two peaks.")
      }
    }
  }
  # Default to using all barcodes in filtered conditions if GMM and KDE don't work
  barcode_df[filtered == "whitelist_barcode" & ncpm_pois >= 0.95, filtered := "pois_validated_barcode"]
  data.table::setkey(barcode_df, 'filtered')
  return(list(
    threshold = min(barcodes['pois_validated_barcode', log1p_ncpm]),
    density = density_output,
    results = barcode_df[filtered == "pois_validated_barcode",],
    best_gmm = gmm_scan
  ))
}

blaze_data_processor<-function(
  path_to_matched_fastq, 
  path_to_putative_bcs
  ){
  putative_bcs<-data.table::fread(path_to_putative_bcs, na.strings = "", select = 2) #2 is the putative_bc column
  putative_bcs<-putative_bcs[, .N, by = putative_bc]
  data.table::setnames(putative_bcs, "N", "original")
  final_counts<-extract_blaze_barcode(path_to_matched_fastq)
  final_counts<-as.data.table(final_counts)
  merged_barcodes<-final_counts[putative_bcs, on = c("barcode" = "putative_bc"), nomatch = 0L]
  merged_barcodes$corrected<-merged_barcodes$total-merged_barcodes$original
  return(merged_barcodes)
}

generate_misalignment_info<-function(read_input_path, read_layout_path){
  prep_read_layout(read_layout_form = read_layout_path)
  misalignment_data<-n_masked_misalignment_cigar_string(
    input_path = read_input_path, 
    adapters = adapters, 
    nthreads = 8, 
    max_sequences = 500000, 
    chunk_size = 50000)
  data.table::setDT(misalignment_data)
  misalignment_threshold<-stat_collector_new(df = misalignment_data, read_layout = read_layout, mode = 'stats')
  data.table::setDT(misalignment_threshold, key = 'query_id')
  misalignment_threshold[, lower_misal_bounds := round(misal_threshold - misal_sd)]
  out<-sapply(misalignment_threshold[,query_id], FUN = function(adapter){
    seq<-misalignment_data[query_id == adapter & best_edit_distance < 
                             misalignment_threshold[adapter, lower_misal_bounds], aligned_sequence]
    full_consensus<-cmv_consensusMatrix(seq)
    seq<-full_consensus %>% 
      {.[,c(1:stringr::str_length(adapters[adapter]))]} %>%
      {prop.table(., margin = 2) * 100} %>%
      # Apply function to get consensus with 'N' for non-dominant bases
      {apply(., MARGIN = 2, FUN = function(X) {
        if (max(X) > 50) {
          rownames(full_consensus)[which.max(X)]
        } else {
          "N"
        }
      }, simplify = TRUE)} %>%
      paste0(., collapse = '')
    return(list(seq = seq, full_consensus = full_consensus))
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  masked_adapters<-sapply(out, FUN = function(x){x$seq}, USE.NAMES = TRUE, simplify = TRUE)
  masked_adapters<-append(masked_adapters,adapters[-which(names(adapters) %in% names(masked_adapters))])
  
  masked_misalignment_threshold<-n_masked_misalignment_cigar_string(
    input_path = read_input_path, 
    adapters = masked_adapters, 
    nthreads = 8, 
    max_sequences = 500000, 
    chunk_size = 50000) %>% {stat_collector_new(df = ., read_layout = read_layout, mode = 'stats')}
  
  return(
      list(
        misalignment_threshold = misalignment_threshold,
        masked_misalignment_threshold = masked_misalignment_threshold,
        masked_adapters = masked_adapters,
        adapters = adapters
      )
  )
}

stat_collector_new<-function(df, read_layout, mode = "stats") {
  data.table::setDT(df)
  data.table::setDT(read_layout)
  adapters <- read_layout[type == "static" & 
      direction %in% c("forward", "reverse") & 
      !class %in% c("poly_tail", "start", "stop"),
    .(class_id, direction)]
  
  forward_adapters <- adapters[direction == "forward", class_id]
  reverse_adapters <- adapters[direction == "reverse", class_id]
  
  # Find IDs with good alignments (edit distance <= 1) more efficiently
  forward_zero_ids <- df[query_id %in% forward_adapters, 
    if(all(best_edit_distance == 0)) .SD[1], 
    by = id][, id]
  
  reverse_zero_ids <- df[query_id %in% reverse_adapters,
    if(all(best_edit_distance == 0)) .SD[1],
    by = id][, id]
  
  if(mode == "stats") {
    # Calculate misalignment stats (wrong direction)
    stats_reverse <- df[id %in% forward_zero_ids & 
        query_id %in% reverse_adapters,
      .(misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance),
        misal_mean_start = mean(best_start_pos/sequence_length)*100,
        misal_mean_stop = mean(best_stop_pos/sequence_length)*100,
        misal_sd_start = sd(best_start_pos/sequence_length)*100,
        misal_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    stats_forward <- df[id %in% reverse_zero_ids & 
        query_id %in% forward_adapters,
      .(misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance),
        misal_mean_start = mean(best_start_pos/sequence_length)*100,
        misal_mean_stop = mean(best_stop_pos/sequence_length)*100,
        misal_sd_start = sd(best_start_pos/sequence_length)*100,
        misal_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    # Calculate correct alignment stats (right direction)
    correct_stats_forward <- df[id %in% forward_zero_ids & 
        query_id %in% forward_adapters,
      .(align_mean_start = mean(best_start_pos/sequence_length)*100,
        align_mean_stop = mean(best_stop_pos/sequence_length)*100,
        align_sd_start = sd(best_start_pos/sequence_length)*100,
        align_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    correct_stats_reverse <- df[id %in% reverse_zero_ids & 
        query_id %in% reverse_adapters,
      .(align_mean_start = mean(best_start_pos/sequence_length)*100,
        align_mean_stop = mean(best_stop_pos/sequence_length)*100,
        align_sd_start = sd(best_start_pos/sequence_length)*100,
        align_sd_stop = sd(best_stop_pos/sequence_length)*100),
      by = query_id]
    
    # Combine misalignment results
    misalignment_thresholds <- data.table::rbindlist(list(
      stats_forward[, direction := "forward"],
      stats_reverse[, direction := "reverse"]
    ))
    
    # Add correct alignment stats
    correct_alignments <- data.table::rbindlist(list(
      correct_stats_forward[, direction := "forward"],
      correct_stats_reverse[, direction := "reverse"]
    ))
    
    # Merge misalignment and correct alignment stats
    results <- merge(
      misalignment_thresholds,
      correct_alignments,
      by = c("query_id", "direction"),
      all = TRUE
    )
    
    # Add metadata
    results <- merge(
      results,
      read_layout[, .(class_id, class, type)],
      by.x = "query_id",
      by.y = "class_id",
      all.x = TRUE
    )
    
    return(results)
    
  } else if(mode == "graphics") {
    out<-data.table::rbindlist(list(
      # Misaligned sequences
      df[id %in% reverse_zero_ids & 
          query_id %in% forward_adapters, 
        .SD, 
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "forward", alignment_type = "misaligned")],
      
      df[id %in% forward_zero_ids & 
          query_id %in% reverse_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "reverse", alignment_type = "misaligned")],
      
      # Correctly aligned sequences
      df[id %in% forward_zero_ids & 
          query_id %in% forward_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "forward", alignment_type = "correct")],
      
      df[id %in% reverse_zero_ids & 
          query_id %in% reverse_adapters,
        .SD,
        .SDcols = c("id", "query_id", "best_edit_distance", 
          "best_start_pos", "best_stop_pos", 
          "sequence_length",'cigar_matches')
      ][, `:=`(direction = "reverse", alignment_type = "correct")]
    ))
    out[,'best_start_perc' := (best_start_pos/sequence_length)*100]
          }
}

substring_odds<-function(subsequence_length, sequence_length){
  return((mapply(sequence_length, subsequence_length, FUN = 
                   function(seq_len, sub_seq){
                     1-(1-(0.25)^sub_seq)^(seq_len-sub_seq+1)
  }, USE.NAMES = FALSE, SIMPLIFY = TRUE)))
}

Rcpp::cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Parse CIGAR and append triplets (i=row, j=col) for covered ref bases
inline void emit_cigar_triplets(const std::string& cigar,
                                int aln_start, int /*aln_end*/,
                                int region_start, int region_end,
                                int row1,                     // 1-based row index
                                std::vector<int>& I,
                                std::vector<int>& J) {
  long long ref_pos = aln_start;
  const char* p = cigar.c_str();
  const char* const pend = p + cigar.size();

  while (p < pend) {
    // parse length
    long long len = 0;
    if (*p < \'0\' || *p > \'9\') { ++p; continue; }
    while (p < pend && *p >= \'0\' && *p <= \'9\') {
      len = len * 10 + (*p - \'0\');
      ++p;
    }
    if (p >= pend) break;
    char op = *p++;

    switch (op) {
      case \'M\': case \'=\': {
        long long seg_start = ref_pos;
        long long seg_end   = ref_pos + len;
        long long cs = std::max(seg_start, (long long)region_start);
        long long ce = std::min(seg_end,   (long long)region_end);
        if (ce > cs) {
          // add (row1, col1..)
          int base_col0 = (int)(cs - region_start); // 0-based
          int n = (int)(ce - cs);
          I.reserve(I.size() + n);
          J.reserve(J.size() + n);
          for (int k = 0; k < n; ++k) {
            I.push_back(row1);                  // 1-based
            J.push_back(base_col0 + k + 1);     // 1-based
          }
        }
        ref_pos += len;
        break;
      }
      case \'D\': case \'N\': case \'X\':
        ref_pos += len; // consume reference only
        break;
      // I, S, H, P, B: do not advance reference
      default:
        break;
    }
  }
}

// [[Rcpp::export]]
Rcpp::List cigar_onehot_triplets(Rcpp::List cigars,
                                 Rcpp::IntegerVector aln_starts,
                                 Rcpp::IntegerVector aln_ends,
                                 int region_start,  // 0-based inclusive
                                 int region_end     // 0-based exclusive
) {
  int n = cigars.size();
  if (n != aln_starts.size() || n != aln_ends.size())
    Rcpp::stop("Lengths of cigars, aln_starts, aln_ends must match.");
  if (region_end < region_start)
    Rcpp::stop("region_end must be >= region_start.");

  std::vector<int> I; I.reserve(1024);
  std::vector<int> J; J.reserve(1024);

  for (int r = 0; r < n; ++r) {
    std::string cg = Rcpp::as<std::string>(cigars[r]);
    int s = aln_starts[r];
    int e = aln_ends[r]; // not strictly needed, but kept for symmetry
    (void)e;
    // quick skip if no overlap with region (cheap bound test)
    //if (s >= region_end) continue;
    // we can’t know exact M/=X coverage without parsing, so parse anyway
    emit_cigar_triplets(cg, s, e, region_start, region_end, r + 1, I, J);
  }

  return Rcpp::List::create(
    Rcpp::Named("i")    = Rcpp::IntegerVector(I.begin(), I.end()),
    Rcpp::Named("j")    = Rcpp::IntegerVector(J.begin(), J.end()),
    Rcpp::Named("nrow") = n,
    Rcpp::Named("ncol") = region_end - region_start
  );
}
')

# ---- R helper to build a sparse matrix ----
sparse_onehot_matrix <- function(cigars, aln_starts, aln_ends,
                                 region_start, region_end,
                                 rownames = NULL,
                                 logical_matrix = FALSE,
                                 bin = 1L,
                                 agg = c("any", "sum")) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  agg <- match.arg(agg)
  if (!is.integer(bin)) bin <- as.integer(bin)
  
  trip <- cigar_onehot_triplets(cigars, aln_starts, aln_ends, region_start, region_end)
  i  <- trip$i
  j  <- trip$j
  nr <- trip$nrow
  nc <- trip$ncol
  
  if (length(i) == 0L) {
    # Empty matrix fast-path
    nc_eff <- if (bin > 1L) ceiling(nc / bin) else nc
    M <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                              dims = c(nr, nc_eff), repr = "C")
    if (!is.null(rownames)) rownames(M) <- rownames
    if (logical_matrix) M <- (M != 0)
    attr(M, "bin") <- bin
    attr(M, "region_start") <- region_start
    attr(M, "region_end") <- region_end
    attr(M, "bin_edges0") <- seq.int(from = region_start,
                                     to   = region_start + nc_eff*bin,
                                     by   = bin)[1:nc_eff]
    return(M)
  }
  
  if (bin > 1L) {
    # Map original 1-based columns j to 1-based bin columns jb
    jb     <- ((as.integer(j) - 1L) %/% bin) + 1L
    nc_bin <- as.integer((nc + bin - 1L) %/% bin)
    
    if (agg == "any") {
      # Any coverage in a bin -> 1
      x <- rep.int(1L, length(jb))
      M <- Matrix::sparseMatrix(i = i, j = jb, x = x,
                                dims = c(nr, nc_bin), repr = "C")
      if (logical_matrix) {
        M <- (M != 0)  # lgCMatrix
      } else {
        # Binarize in case duplicates summed > 1
        M@x[] <- 1L
      }
    } else { # agg == "sum" = number of covered bases per bin
      # Each triplet represents one base; summing duplicates yields base count per bin.
      x <- rep.int(1L, length(jb))
      M <- Matrix::sparseMatrix(i = i, j = jb, x = x,
                                dims = c(nr, nc_bin), repr = "C")
      if (logical_matrix) {
        # For logical output, turn any nonzero sum to TRUE
        M <- (M != 0)
      }
    }
    
    if (!is.null(rownames)) rownames(M) <- rownames
    attr(M, "bin") <- bin
    attr(M, "region_start") <- region_start
    attr(M, "region_end") <- region_end
    # 0-based half-open starts for each bin on the reference
    attr(M, "bin_edges0") <- seq.int(from = region_start,
                                     to   = region_start + nc_bin*bin,
                                     by   = bin)[1:nc_bin]
    return(M)
  }
  
  # No binning (bin == 1): original behavior
  x <- rep.int(1L, length(i))
  M <- Matrix::sparseMatrix(i = i, j = j, x = x,
                            dims = c(nr, nc), repr = "C")
  if (!is.null(rownames)) rownames(M) <- rownames
  if (logical_matrix) {
    M <- (M != 0)
  } else {
    M@x[] <- 1L
  }
  attr(M, "bin") <- 1L
  attr(M, "region_start") <- region_start
  attr(M, "region_end") <- region_end
  attr(M, "bin_edges0") <- seq.int(from = region_start,
                                   to   = region_start + nc,
                                   by   = 1L)[1:nc]
  M
}

