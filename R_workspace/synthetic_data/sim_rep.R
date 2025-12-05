total_simreps<-rbindlist(list(tcrb_simrep, sim_rep_tra))
total_simreps$revcomps<-revcomp(total_simreps$sequence)
total_sequences<-total_simreps$revcomps
random_barcodes<-sample(whitelist$V1, size = 2500, replace = FALSE)
sequences<-list()
bases<-c("A","C","T","G")

fasta<-rbindlist(sequences)
ec<-prop.table(table(anger_frame$error_code)) %>% data.frame(.)

sequences<-vector(length = 2500, mode = "list")
new_seqs<-pbmapply(x = random_barcodes, y = sequences, function(x,y){
  y<-paste0("TATGTTTATGTGTGTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGCATAGTTCTGCATGATGGGTTAGCAGCACCTCTACACGACGCTCTTCCGATCT",
         x,
         paste0(sample(bases, size = 12, replace = TRUE), collapse = ""),
         "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
         sample(total_sequences, size = sample(1:500, size = 1,replace = TRUE), replace = TRUE),
         "CCCATGTACTCTGCGTTGATACCACTGCTTAGGTGCTGCTAACCCATCATGCAGAACTATGCTTAACCTTAGCAA")
  return(y)
})
names(sequences)<-random_barcodes
sequences<-lapply(new_seqs, function(x){
  df<-as.data.table(x)
})


pbsim2<-"/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/cmv/miniconda3/envs/blazeseq/bin/pbsim"
for(i in 1:length(sequences)){
  dir.create(path = paste0("/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/synth_outs/",names(sequences)[[i]]))
  lord_fastquad(df = sequences[[i]], fn = paste0("/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/synth_outs/",
    names(sequences)[[i]],"/",names(sequences)[[i]],"_reference.fa"), type = "fa", append = FALSE)
  setwd(paste0("/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/synth_outs/",names(sequences)[[i]]))
  pbsim2_args<-c("--template-fasta",
                 paste0("/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/synth_outs/",names(sequences)[[i]],"/",
                   names(sequences)[[i]],"_reference.fa")
                 ,"--hmm_model","/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/cmv/miniconda3/envs/blazeseq/data/R103.model")
  processx::run(command = pbsim2, args = pbsim2_args, echo_cmd =TRUE, echo = FALSE, spinner = TRUE)
}

for(i in 1:length(sequences)){
  sequences[[i]]$id<-NA
  sequences[[i]]$id<-paste0(random_barcodes[i], "_",seq(nrow(sequences[[i]])))
  sequences[[i]]<-sequences[[i]][,c(2,1)]
  colnames(sequences[[i]])<-c("id","seq")
}
crosscheck<-table(sapply(seq_along(sequences), function(x){all(grepl(pattern = names(sequences)[x], x = sequences[[x]]$seq))}))





seq_length<-seq(1,16)
generate_mismatch<-function(barcode){
  sapply(X = seq_length, function(y){
    paste0(substr(x = barcode, 1, y[seq_along(y)]-1), bases, substr(x = barcode, y[seq_along(y)]+1, 16))
  }, simplify = FALSE) %>% unlist(.) %>% unique(.)
}
generate_mismatch<-Vectorize(generate_mismatch, "barcode")


la3_fastas<-list()
progress<-progress::progress_bar$new(total = length(la3), format = "  progress [:bar] :percent eta: :eta",)
for(i in 1:length(la3)){
  progress$tick(0)
  la3_fastas[[i]]<-read_fastqas(fn = la3[i], type = "fa")
  progress$tick()
}
names(la3_fastas)<-sapply(la3, function(x){basename(x)}) %>% gsub(pattern = ".correctedReads.fasta.gz", 
                                                                  replacement = "", x = .)
for(i in 1:length(la3_fastas)){
  la3_fastas[[i]]$id<-paste0(names(la3_fastas)[i],"_",seq(nrow(la3_fastas[[i]])))
}

anger_frame_iterants<-vector(length = 10, mode = "list")
for(i in 0:9){
  print(i)
  colnames(anger_frame)<-c("id","fastq_files","qc")
  anger_frame<-anger_frame[,c(1,3,2)]
  anger_frame_iterants[i]<-always_angry(anger_frame = anger_frame, primer_1 = primer_1, primer_2 = primer_2, 
    whitelist_path = whitelist_path, cores = 1, adapter_distance = i) %>%
    stripping_in_anger(.) %>%
    erred_in_anger(.)
}

anger_frame_iterants<-pblapply(seq(1:10), function(x){
  out<-always_angry(anger_frame = anger_frame, primer_1 = primer_1, primer_2 = primer_2, whitelist_path = whitelist_path, 
    cores = 1, adapter_distance = x) %>%
    stripping_in_anger(.) %>%
    erred_in_anger(.)
  return(out)
})

matched_files<-list.files(path = "/dartfs/rc/lab/A/AckermanLab/CMV/PostRot/RAGE-Seq_RDS/synth_outs", full.names = TRUE)
anger_frame<-pblapply(matched_files, function(x){
  files<-list.files(x, pattern = ".fastq", full.names = TRUE)
  files<-lapply(files, function(x){read_fastqas(fn = x, type = "fq")})
  files<-rbindlist(files)
  return(files)
})

anger_frame<-anger_frame_list[[1]]
table(which(names(anger_frame_list)[1] %in% anger_frame$putative_bcs))

total_outcomes<-rbindlist(anger_frame_list) %>% table(.$putative_bcs[which(.$cr_barcode == TRUE)] %in% names(.))

outcome<-list()
for(i in 1:length(anger_frame_list)){
  test<-anger_frame_list[[i]]
  outcome[[i]]<-table((test$putative_bcs[which(test$cr_barcode == TRUE)] %in% names(anger_frame_list)[i])) %>% data.frame(.)
}
