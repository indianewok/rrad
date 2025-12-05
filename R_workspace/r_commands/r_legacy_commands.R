#' @importFrom data.table fread setkey := .SD copy
#' @importFrom magrittr %>%
#' @importFrom stringr str_length str_split
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom tidytable filter group_by summarise pull as_tidytable
#' @importFrom glue glue
#' @importFrom tools file_ext file_path_sans_ext
#' @export
# stream_fastqas<-function(fn, type, seqkit_path = NULL, start, stop,...){
#   ext<-tools::file_ext(fn)
#   if (ext == "gz"){
#     cat_cmd<-"zcat"
#   }
#   else{
#     cat_cmd<-"cat"
#   }
#   if(type == "fq"){
#     res<-fread(cmd = glue::glue("{seqkit_path} range -r {start}:{stop} {fn} | paste - - - - | cut -f1,2,4"), 
#       col.names = c("id", "fastq_files","qc"), sep = "\t", quote = "", header = FALSE, ...)
#     res<-res[,c(1,3,2)]
#     res$id<-str_split(res$id, pattern = " ") %>% sapply(., function(x){x[1]}, USE.NAMES = FALSE)
#     return(res)
#   }
#   if(type == "fa"){
#     res<-fread(cmd = glue::glue("{seqkit_path} range -r {start}:{stop} | paste - - | cut -f1,2"),
#       quote ="", header = FALSE,
#       col.names = c("id", "seq"), sep = "\t", ...)
#     return(res)
#   }
# }
#' @export
whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  print("Importing a little bit of sample data...")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1,
    col.names = "whitelist_bcs")
  if(class(whitelist$whitelist_bcs) == "character"){
    print("Character type detected! Importing the whitelist.")
    whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1, col.names = "whitelist_bcs")
    if(convert != FALSE){
      print("Converting the whitelist into a more efficient format!")
      whitelist$whitelist_bcs<-barcodes_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",file_path_sans_ext(whitelist_path)
          %>% basename(.),
          "_bitlist.csv.gz")
        print(paste0("Now saving the integer whitelist to ",new_path,"!"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      print("Integer64-type whitelist detected! Now importing...")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
        col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      return(whitelist)
    }
  }
}
#' @export
# prepare_to_anger<-function(read_layout_form, external_path_form, create_output_dir = TRUE,
#    data_table_nthreads = 1, seqkit_stats = TRUE, ...){
#   {
#     print("Importing read layout and figuring out its order!")
#     tidytable::setDTthreads(threads = 1)
#     read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
#     read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
#     if(any(read_layout$type %in% "poly_a" | read_layout$type %in% "poly_t")){
#       if(("poly_a" %in% read_layout$type) == TRUE){
#         read_layout$seq[which(read_layout$type %in% "poly_a")]<-"A{12,}+"
#       }
#       else{
#         read_layout$seq[which(read_layout$type %in% "poly_t")]<-"T{12,}+"
#       }
#     }
#     read_layout$order<-seq(1:nrow(read_layout))
#     read_layout$direction<-"forward"
#     barcode_position<-read_layout[order == min(order[type == "barcode"])]$order
#     adapters_near_barcode<-read_layout[type == "adapter" & (order - barcode_position) <= 1]
#     if(nrow(adapters_near_barcode) > 0){
#       flanking_adapter_id<-adapters_near_barcode[order(-expected_length), head(.SD, 1)]$id
#       read_layout[id == flanking_adapter_id, type := "flanking_adapter"]
#     }
#     
#     reverses_table<-copy(read_layout) %>%
#       .[, seq := ifelse(type == "poly_t", "A{12,}+", ifelse(type == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>% 
#       .[, id := ifelse(type %in% c("poly_a", "poly_t"), ifelse(type == "poly_a", "poly_t", "poly_a"),paste0("rc_", id))] %>% 
#       .[, type := ifelse(type == "poly_a", "poly_t", ifelse(type == "poly_t", "poly_a", type))] %>%
#       .[, direction := "reverse"]
#     
#     read_layout<-rbind(read_layout, reverses_table) #list2env this one
#     read_layout$order<-seq(1:nrow(read_layout))
#     #adapters<-read_layout$seq[grep(pattern = "adapter", x = read_layout$type)]
#     #names(adapters)<-read_layout$id[grep(pattern = "adapter", x = read_layout$type)] #list2env this one
#     adapters<-read_layout$seq[-which(is.na(read_layout$seq))]
#     names(adapters)<-read_layout$id[-which(is.na(read_layout$seq))]
#   }#preparing read layout and whatnot
#   {
#     print("Checking executable paths!")
#     path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
#        key = "path_type", data.table = TRUE)
#     if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
#       dir.create(path = path_layout[path_type == "output_dir", actual_path])
#     }
#     whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
#   }#preparing external paths for input
#   return(list2env(
#     x = list(read_layout = read_layout, 
#       path_layout = path_layout, 
#       adapters = adapters, 
#       whitelist = whitelist), envir = .GlobalEnv))
# }
#' @export
# stat_collector_old<-function(df, read_layout, mode = NULL){
#   forward_adapters<-read_layout[type %in% c("adapter", "flanking_adapter") & direction == "forward", id]
#   reverse_adapters<-read_layout[type %in% c("adapter", "flanking_adapter") & direction == "reverse", id]
#   df<-as_tidytable(df)
#   forward_zero_ids<-df %>%
#     tidytable::filter(query_id %in% forward_adapters) %>%
#     tidytable::group_by(id) %>%
#     tidytable::filter(all(best_edit_distance == 0)) %>%
#     tidytable::pull(id)
#   reverse_zero_ids<-df %>% 
#     tidytable::filter(query_id %in% reverse_adapters) %>%
#     tidytable::group_by(id) %>%
#     tidytable::filter(all(best_edit_distance == 0)) %>%
#     tidytable::pull(id)
#   if(mode == "stats"){
#     stats_reverse<-df %>%
#       tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
#       tidytable::group_by(query_id) %>%
#       tidytable::summarise(
#         null_distance = mean(best_edit_distance), 
#         sd_null = sd(best_edit_distance))
#     stats_forward<-df %>%
#       tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
#       tidytable::group_by(query_id) %>%
#       tidytable::summarise(
#         null_distance = mean(best_edit_distance), 
#         sd_null = sd(best_edit_distance))
#     null_distance<-rbind(stats_forward, stats_reverse)
#     return(null_distance)
#   }
#   if(mode == "graphics"){
#     stats_reverse<-df %>%
#       tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
#       tidytable::group_by(query_id)
#     stats_forward<-df %>%
#       tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
#       tidytable::group_by(query_id)
#     total_stats<-data.table::rbindlist(list(stats_forward, stats_reverse))
#     return(total_stats)
#   }
# }
#' @export
write_fastqas<-function(df, fn, type, append = FALSE){
  if(type == "fq"){
    fq_list<-df %>%
      dplyr::mutate(dummy = "+") %>%
      dplyr::select(id, seq, dummy, qual) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fq_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
  }
  if(type == "fa"){
    if(any(grepl(">", df$id)) == FALSE){
      df$id<-paste0(">",df$id)
      fa_list<-df %>%
        dplyr::select(id, seq) %>%
        as.matrix() %>%
        t() %>%
        as.character() %>%
        tibble::as_tibble()
      data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
    }
    fa_list<-df %>%
      dplyr::select(id, seq) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append)
  }
}
#' @export
read_fastqas<-function(fn, type, full_id = FALSE, ...){
  rfq_single<-function(fn, type, full_id = FALSE, ...){
    ext<-tools::file_ext(fn)
    if (ext == "gz"){
      cat_cmd<-"zcat"
    }
    else{
      cat_cmd<-"cat"
    }
    if(type == "fq"){
      res<-fread(cmd = glue::glue("{cat_cmd} {fn} | paste - - - - | cut -f1,2,4"), col.names = c("id", "seq","qual"),
        sep = "\t", header = FALSE, quote = "", ...)
      res<-res[,c(1,3,2)]
      if(full_id == TRUE){
        res$full_id<-res$id
        res$id<-str_split(res$id, pattern = " ") %>% sapply(., function(x){x[1]}, USE.NAMES = FALSE)
        return(res)
      } else {
        res$id<-str_split(res$id, pattern = " ") %>% sapply(., function(x){x[1]}, USE.NAMES = FALSE)
        return(res)
      }
    }
    if(type == "fa"){
      awk_cmd <- 'awk \'/^>/{if (NR>1) printf("\\n"); printf("%s\\t",$0); next} {printf("%s",$0);} END {printf("\\n");}\''
      if(ext == "gz"){
        res<-fread(cmd = glue::glue("gunzip -c {fn} | {awk_cmd}"), col.names = c("id", "seq"), sep = "\t", ...)
      } else {
        res<-fread(cmd = glue::glue("{awk_cmd} {fn}"), col.names = c("id", "seq"), sep = "\t", ...)
      }
      return(res)
    }
    
  }
    if(class(fn) == "list"){
      out<-pbapply::pblapply(X = fn, FUN = function(X){
        if(length(X) == 1){
            return(rfq_single(fn = X, full_id = FALSE, type = type, ...))
        } else {
            return(NA)
          }
      })
      return(out[-which(is.na(out))])
    } else {
      return(rfq_single(fn = fn, type = type, full_id = FALSE, ...))  
    }
}
#' @export
#legacy command--simplify and clean up
bajrun<-function(path_layout_form, read_layout_form, 
  test_mode = FALSE, nthreads_sigstrings = 1,
  nthreads_sigstract = 1, test_return_stage = NULL,
  external_sr_bc = FALSE, chunk_divisor = 50, 
  barcorrect = TRUE, jaccard_on = FALSE){
  prepare_to_anger(read_layout_form = read_layout_form, 
    external_path_form = path_layout_form)
  input_path<-path_layout["file","actual_path"]
  output_path<-path_layout["output_dir","actual_path"]
  if(external_sr_bc == TRUE){
    external_bc_path<-path_layout["external_sr_bc","actual_path"]
    external_bcs<-readLines(external_bc_path[[1]])
  }
  print(paste0("The files are coming in from ", input_path, 
    " and the files are going to be written to ", output_path))
  if(test_mode == TRUE && test_return_stage != "chunked_fastqs"){
    print("Testing that this actually works!")
    if(dir.exists(input_path[[1]])){
      print("Input path is a directory!")
      gz_files<-list.files(path = input_path[[1]], pattern = ".gz", full.names = TRUE)
      file_chunks<-split(gz_files, ceiling(seq_along(gz_files) / 50))
      lapply(seq_along(file_chunks), function(i){
        print("Testing to find all file chunks...")
        if(file.exists(file_chunks[[i]][1])){
          print(paste0("Chunk ",i," found!"))
        }
      })
      print("Now pulling in the first big chunk...")
      df<-lapply(X=file_chunks[[1]][1:20],FUN = function(X){jaeger::read_fastqas(fn =X, 
        type = "fq")}) %>% data.table::rbindlist(.)
    } else {
      if(file.exists(input_path[[1]])){
        print("Input path is a single fastq file!")
        df<-jaeger::read_fastqas(fn = input_path, type = "fq", full_id = TRUE)
      }
    }
    print("Now calculating misalignment distributions with first 100K sequences!")
    if(test_return_stage == "stats"){
      print("Returning stats mode! Generating  dataframe (1/3)...")
      total_nd<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
        sequences = df$fastq_files[1:100000],
        nthreads= 1)
      print("Returning stats mode! Generating null distances (2/3)...")
      null_distances<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
        sequences = df$fastq_files[1:100000],
        nthreads= 1) %>% jaeger::stat_collector(., read_layout, mode = "stats")
      print("Returning stats mode! Generating graphable df (3/3)...")
      total_stats<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
        sequences = df$fastq_files[1:100000],
        nthreads= 1) %>% jaeger::stat_collector(., read_layout, mode = "graphics")
      return(list2env(x = list(total_nd = total_nd, null_distances = null_distances, total_stats = total_stats)))
    }
    null_distance<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
      sequences = df$fastq_files[1:100000],
      nthreads= 1) %>% jaeger::stat_collector(., read_layout, mode = "stats")
    print(null_distance)
    if(test_return_stage == "import"){
      return(list2env(x = list(df = df, null_distance = null_distance),envir = .GlobalEnv))
    }
    print("Now giving all of the alignment code a test run...")
    print("Filtering the df for any sequences shorter than 150 bp...")
    df<-df[which(str_length(df$fastq_files) >= 150),]
    sigstrings<-vector(length = length(df$fastq_files))
    print(paste0("Generated a vector to hold your sigstrings that's ", 
      length(df$fastq_files)," spaces long!"))
    sigstrings<-bajalign_sigs(adapters, sequences = df$fastq_files,
       null_distance = null_distance, nthreads = nthreads_sigstrings)
    print("Done with sigstringing!")
    print(sigstrings[1:5])
    if(test_return_stage == "sigstrings"){
      return(list2env(x = list(df = df, sigstrings = sigstrings),envir = .GlobalEnv))
    }
    invisible(bajbatch(null_distance = null_distance, read_layout = read_layout, 
      sigstrings = sigstrings))
    print("Done with bajbatching!")
    print(sigstrings[1:5])
    if(test_return_stage == "bajbatch"){
      return(list2env(x = list(df = df, sigstrings = sigstrings),envir = .GlobalEnv))
    }
    print("Now doing baj_extract!")
    df_new<-baj_extract(sigstrings = sigstrings, whitelist_df = whitelist,
       df = df, verbose = FALSE, barcorrect = TRUE, nthreads = nthreads_sigstract)
    print("Done with extraction!")
    if(test_return_stage == "baj_extract"){
      return(list2env(x = list(df = df, df_new = df_new, sigstrings = sigstrings),envir = .GlobalEnv))
    }
  }
  print("Calculating misalignment distributions with first 100K sequences!")
  if(test_mode == TRUE){
    if(test_return_stage == "chunked_fastqs"){
      gz_files<-list.files(path = input_path[[1]], pattern = ".gz", full.names = TRUE)
      file_chunks<-split(gz_files, ceiling(seq_along(gz_files) / 50))
      file_chunks<-file_chunks[1:5]
      lapply(seq_along(file_chunks), function(i){
        print("Testing to find all file chunks...")
        if(file.exists(file_chunks[[i]][1])){
          print(paste0("Chunk ",i," found!"))
        }
      })
    }
  }
  if(dir.exists(input_path[[1]])){
    print("Input path is a directory!")
    gz_files<-list.files(path = input_path[[1]], pattern = ".gz", full.names = TRUE)
    file_chunks<-split(gz_files, ceiling(seq_along(gz_files) / chunk_divisor))
    lapply(seq_along(file_chunks), function(i){
      print("Testing to find all file chunks...")
      if(file.exists(file_chunks[[i]][1])){
        print(paste0("Chunk ",i," found!"))
      }
    })
  }
  #gotta fix this part--harcoded minimum chunk divisor of 20, and I think it can actually go faster
  df<-lapply(X=file_chunks[[1]][1:20],
    FUN = function(X){read_fastqas(fn =X, type = "fq")}) %>% 
    data.table::rbindlist(.)
  null_distance<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
    sequences = df$fastq_files[1:100000],
    nthreads= 1) %>% stat_collector(., read_layout, mode = "stats")
  print("Misalignment Thresholds as follows:")
  print(null_distance)
  write.csv(null_distance, file = paste0(output_path[[1]],"/misalignment_distances.csv"))
  pbapply::pblapply(seq_along(file_chunks), function(i){
    print(paste0("now on chunk ", i, " of processing!"))
    df<-lapply(X=file_chunks[[i]],FUN = function(X){jaeger::read_fastqas(fn =X, 
      type = "fq")}) %>% data.table::rbindlist(.)
    df<-df[which(str_length(df$fastq_files) >= 150),]
    sigstrings<-vector(length = length(df$fastq_files))
    print(paste0("Generated a vector to hold your sigstrings that's ", 
      length(df$fastq_files)," spaces long!"))
    sigstrings<-bajalign_sigs(adapters, sequences = df$fastq_files,
       null_distance = null_distance, nthreads_sigstrings) 
    print("Done with sigstringing!")
    invisible(bajbatch(null_distance = null_distance, 
      read_layout = read_layout, sigstrings = sigstrings))
    print("Done with bajbatching!")
    df_new<-baj_extract(sigstrings = sigstrings, whitelist_df = whitelist, 
      df = df, verbose = FALSE, barcorrect = barcorrect,
       nthreads = nthreads_sigstract, jaccard_on = jaccard_on)
    print("Done with baj_extracting!")
    barcodes<-df_new$barcode
    print(paste0(length(unique(barcodes)), " barcodes in total!"))
    if(external_sr_bc == TRUE){
      df_true<-df_new[which(df_new$barcode %in% external_bcs),]
      print(paste0(nrow(df_true), " barcodes in the external list!"))
    } else {
    df_true<-df_new$barcode %>%
      barcodes_to_bits %>% 
      {df_new[which(.%in%whitelist$whitelist_bcs),]}
    }
    
    
    df_true$id<-DescTools::StrExtractBetween(x = df_true$sig_id, left = "^", right = "\\*")
    df_true<-df_true %>% tidytable::unite("id", c("id", "umi", "barcode"), sep = "_", remove = FALSE) %>% 
      .[,c("id","filtered_read","filtered_qc")] %>% 
      data.table::setnames(c("id","seq","qual"))
    
    append_fastq<-file.exists(paste0(output_path, "/filtered_reads.fastq.gz"))
    fastqa_writer(df = df_true, fn = paste0(output_path, "/filtered_reads.fastq.gz"), 
      type = "fq", append = append_fastq)
    
    append_barcode<-file.exists(paste0(output_path,"/all_barcodes.txt"))
    write(x = barcodes, file = paste0(output_path,"/all_barcodes.txt"), 
      append = append_barcode)
    
    append_sigstrings<-file.exists(paste0(output_path, "sigsummary.txt"))
    write(x = sigstrings, file = paste0(output_path, "/sigsummary.txt"), 
      append = append_sigstrings)
    print("Done with this chunk!")
  })
}
#' @export
read_paf<-function(path, ...){
df<-data.table::fread(input = path, 
  fill = TRUE, sep = NULL, strip.white = FALSE, header = FALSE )
df<-df[[1]] %>% stringr::str_split(., pattern = "\t", simplify = TRUE) %>% {data.table::as.data.table(.)}
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
}
for (col in names(df)) {
  numeric_conversion <- suppressWarnings(as.numeric(df[[col]]))
  if (!any(is.na(numeric_conversion)) || all(is.na(numeric_conversion) == is.na(df[[col]]))) {
    df[[col]] <- numeric_conversion
  }
}
return(df)
}

setkey(df, "id")
str_extract(processed_sigstrings[1], pattern = "barcode:0:\\d+.?.:\\d+.?.?(?=\\|)") %>% str_split(., pattern = ":")
str_extract(processed_sigstrings[1], pattern = "<.+?.>") %>% str_split(., pattern = ":")

str_extract(processed_sigstrings[7], pattern = "<.+?.>") %>% str_split(., pattern = ":") %>% {.[[1]][2]}

barcode_extractor<-function(processed_sigstring){
  barcode_coords<-stringr::str_extract(processed_sigstring, 
    pattern = "barcode:0:\\d+.?.:\\d+.?.?(?=\\|)") %>% 
    stringr::str_split(., pattern = ":") %>%
    {.[[1]][3:4]}
  id<-stringr::str_extract(processed_sigstring, pattern = "<.+?.>") %>% 
    stringr::str_split(., pattern = ":") %>% {.[[1]][2]} %>% 
    {ifelse(test = grep(pattern = "+", x = .), 
    yes = stringr::str_remove(., pattern = "\\+.+.?$"), no = "")}
  barcode = substr(df[id, "seq"], start = barcode_coords[1], stop = barcode_coords[2]) 
  barcode<-ifelse(test = grepl(pattern = ":R>$", x = processed_sigstring), yes =revcomp(barcode), no = barcode)
  return(barcode)
}

extract_barcodes<-function(processed_sigstrings){
  out<-vector(mode = "list", length = length(processed_sigstrings))
  out<-lapply(processed_sigstrings, FUN = function(x){barcode_extractor(x)})
  return(out)
}

extract_components_dt <- function(dt, sig_col_name, components, id_col_name = "id") {
  pattern <- sprintf("(%s):0:(\\d+):(\\d+)(?=\\|)", paste(components, collapse = "|"))
  dt[, (components) := lapply(components, function(component) {
    matches <- stringi::stri_match_all_regex(get(sig_col_name), pattern)[[1]]
    sapply(1:nrow(matches), function(i) {
      comp_match <- matches[i, , drop = FALSE]
      if (comp_match[1, 1] == component) {
        start <- as.integer(comp_match[1, 3])
        end <- as.integer(comp_match[1, 4])
        seq <- dt[get(id_col_name) == id, seq]
        if (!is.null(seq) && start > 0 && end <= nchar(seq)) {
          return(substring(seq, start, end))
        }
      }
      return(NA_character_)
    })
  }), by = .(id_col_name)]
}

processed_sigstrings