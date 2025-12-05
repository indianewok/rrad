#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include "/dartfs-hpc/rc/home/3/f005f43/rad_dev/src/edlib.h"

// Define the structure for storing alignment information
struct AlignmentInfo {
  int best_edit_distance;
  int best_start_pos;
  int best_stop_pos;
  int second_edit_distance;
  int second_start_pos;
  int second_stop_pos;
};

// Function to process a sequence and perform alignments
std::vector<AlignmentInfo> process_sequence(const std::string& sequence, const std::vector<std::string>& queries, int nthreads) {
  std::vector<AlignmentInfo> alignments;
  omp_set_num_threads(nthreads);
  
#pragma omp parallel for
  for (int j = 0; j < queries.size(); ++j) {
    const auto& query = queries[j];
    
    EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
    std::set<UniqueAlignment> uniqueAlignments;
    char* cquery = const_cast<char*>(query.c_str());
    char* csequence = const_cast<char*>(sequence.c_str());
    EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
    std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
    std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
    
    for (int& pos : start_positions) {
      ++pos;
    }
    for (int& pos : end_positions) {
      ++pos;
    }
    
    std::sort(start_positions.begin(), start_positions.end());
    std::sort(end_positions.begin(), end_positions.end());
    
    auto last = std::unique(start_positions.begin(), start_positions.end());
    start_positions.erase(last, start_positions.end());
    
    std::vector<int> unique_end_positions;
    for (const auto& start : start_positions) {
      auto pos = std::find(start_positions.begin(), start_positions.end(), start) - start_positions.begin();
      unique_end_positions.push_back(end_positions[pos]);
    }
    
    int edit_distance = cresult.editDistance;
    
    for (size_t k = 0; k < start_positions.size(); ++k) {
      UniqueAlignment ua = {edit_distance, start_positions[k], unique_end_positions[k]};
      uniqueAlignments.insert(ua);
    }
    
    auto it = uniqueAlignments.begin();
    UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
    ++it;
    UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
    
    alignments.push_back({
      best.edit_distance,
      best.start_position,
      best.end_position,
      secondBest.edit_distance,
      secondBest.start_position,
      secondBest.end_position
    });
    
    edlibFreeAlignResult(cresult);
  }
  
  return alignments;
}

// Function to stream FASTA file and perform alignments
// [[Rcpp::export]]
Rcpp::DataFrame stream_fasta_and_align(const std::string& fasta_path, Rcpp::CharacterVector adapters, int nthreads = 1) {
  std::ifstream fasta_file(fasta_path);
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  
  std::map<int, std::vector<AlignmentInfo>> all_alignments;
  std::string line, sequence;
  int sequence_counter = 1;
  
  while (std::getline(fasta_file, line)) {
    if (line[0] == '>') {
      if (!sequence.empty()) {
        auto alignments = process_sequence(sequence, queries, nthreads);
        all_alignments[sequence_counter] = alignments;
        sequence.clear();
        ++sequence_counter;
      }
    } else {
      sequence += line;
    }
  }
  
  if (!sequence.empty()) {
    auto alignments = process_sequence(sequence, queries, nthreads);
    all_alignments[sequence_counter] = alignments;
  }
  
  fasta_file.close();
  
  std::vector<int> ids;
  std::vector<std::string> query_ids;
  std::vector<int> best_edit_distances;
  std::vector<int> best_start_positions;
  std::vector<int> best_stop_positions;
  std::vector<int> second_edit_distances;
  std::vector<int> second_start_positions;
  std::vector<int> second_stop_positions;
  
  for (const auto& pair : all_alignments) {
    int sequence_id = pair.first;
    const auto& alignments = pair.second;
    for (size_t i = 0; i < alignments.size(); ++i) {
      ids.push_back(sequence_id);
      query_ids.push_back(Rcpp::as<std::string>(query_names[i]));
      best_edit_distances.push_back(alignments[i].best_edit_distance);
      best_start_positions.push_back(alignments[i].best_start_pos);
      best_stop_positions.push_back(alignments[i].best_stop_pos);
      second_edit_distances.push_back(alignments[i].second_edit_distance);
      second_start_positions.push_back(alignments[i].second_start_pos);
      second_stop_positions.push_back(alignments[i].second_stop_pos);
    }
  }
  
  Rcpp::DataFrame result = Rcpp::DataFrame::create(
    Rcpp::Named("id") = ids,
    Rcpp::Named("query_id") = query_ids,
    Rcpp::Named("best_edit_distance") = best_edit_distances,
    Rcpp::Named("best_start_pos") = best_start_positions,
    Rcpp::Named("best_stop_pos") = best_stop_positions,
    Rcpp::Named("second_edit_distance") = second_edit_distances,
    Rcpp::Named("second_start_pos") = second_start_positions,
    Rcpp::Named("second_stop_pos") = second_stop_positions
  );
  return result;
}
