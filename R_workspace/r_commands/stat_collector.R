stat_collector<-function(df, read_layout, mode = "stats"){
  forward_adapters<-read_layout[type %in% c("static") & direction == "forward" & class != "poly_tail", .(class_id)]
  reverse_adapters<-read_layout[type %in% c("static") & direction == "reverse" & class != "poly_tail", .(class_id)]
  df<-tidytable::as_tidytable(df)
  forward_zero_ids<-df %>%
    tidytable::filter(query_id %in% forward_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance == 0)) %>%
    tidytable::pull(id)
  reverse_zero_ids<-df %>%
    tidytable::filter(query_id %in% reverse_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance == 0)) %>%
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