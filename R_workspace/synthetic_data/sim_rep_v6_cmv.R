synth_gen<-function(number_of_cells, read_layout, whitelist){
  bases <- c("A", "C", "T", "G")
  data<-read_layout$seq
  names(data)<-read_layout$id
  if("poly_t" %in% names(data)){
    data[grep(pattern = "poly_t", x = names(data))]<-paste(rep("T", sample(15:30, size = 1)), collapse = "")
  } else if ("poly_a" %in% names(data)){
    data[grep(pattern = "poly_a", x = names(data))]<-paste(rep("A", sample(15:30, size = 1)), collapse = "")
  }
 reads<-grep(pattern = "read", x = names(data))
 for(i in 1:length(reads)){
   data[reads]<-replicate(1, paste(sample(c("A","C","T","G"), size = 300, replace = TRUE), collapse = ""))
 }
 barcodes<-grep(pattern = "barcode", x = names(data))
   subsample<-sample(whitelist$whitelist_bcs, size = 20000, replace = FALSE) %>% bits_to_sequence
   barcode_list<-sample(subsample, number_of_cells, replace = FALSE)
   cells<-pblapply(barcode_list, FUN = function(x){
     data[barcodes]<-x
     rng_umis<-replicate(10, paste(sample(c("A","C","T","G"), 
       size = read_layout$expected_length[which(read_layout$id == "umi")], 
       replace = TRUE), collapse = ""))
     umi_idx<-grep(pattern = "umi", x = names(data))
     forward_sequences<-list()
     for(i in 1:length(rng_umis)){
       data[umi_idx]<-rng_umis[i]
       forward_sequences[[i]]<-paste0(data, collapse = "")
       names(forward_sequences)[i]<-rng_umis[i]
     }
     reverse_sequences<-list()
     for(i in 1:length(rng_umis)){
       data[umi_idx]<-rng_umis[i]
       reverse_sequences[[i]]<-revcomp(paste0(data, collapse = ""))
       names(reverse_sequences)[i]<-rng_umis[i]
     }
     sequences<-append(forward_sequences, reverse_sequences)
     sequences<-data.table(id = paste0(x, "_", names(sequences),"_", 1:length(sequences)), seq = unlist(sequences))
   })
   cells<-rbindlist(cells)
   return(cells)
}