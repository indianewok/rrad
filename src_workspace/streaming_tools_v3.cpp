#include "rad.h"
#include <fstream>
#include <zlib.h>
#include "rapidgzip/rapidgzip.hpp"
using namespace std;
using namespace Rcpp;

//WORK IN PROGRESS--STREAMING SEQUENCES

// [[Rcpp::export]]
Rcpp::DataFrame misalignment_stream(const std::string& fastq_path, Rcpp::CharacterVector adapters, 
  int nthreads = 1, int max_sequences = -1) {
  
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  // Filter out queries with "poly" in their names
  std::vector<std::string> filtered_queries;
  Rcpp::CharacterVector filtered_query_names;
  for (int i = 0; i < query_names.size(); ++i) {
    std::string name = Rcpp::as<std::string>(query_names[i]);
    if (name.find("poly") == std::string::npos) {
      filtered_queries.push_back(queries[i]);
      filtered_query_names.push_back(query_names[i]);
    }
  }
  
  std::map<int, std::vector<AlignmentInfo>> all_alignments;
  size_t initial_buffer_size = 98304;
  std::vector<char> buffer(initial_buffer_size);
  
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  auto read_line = [&](gzFile file, std::vector<char>& buffer) -> std::string {
    std::string result;
    while (gzgets(file, buffer.data(), buffer.size()) != Z_NULL) {
      result += std::string(buffer.data());
      if (result.back() == '\n') break;
      buffer.resize(buffer.size() * 2); // Double the buffer size if not complete
    }
    result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
    return result;
  };
  
  omp_set_num_threads(nthreads);
  
  while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL &&
    (max_sequences == -1 || processed_sequences < max_sequences)) {
    line = std::string(buffer.data());
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(1, line.find(' ') - 1); // Extract ID between '@' and first space
    sequence = read_line(fastq_file, buffer);
    if (sequence.empty()) break;
    
    plus = read_line(fastq_file, buffer);
    if (plus.empty() || plus[0] != '+') {
      Rcpp::Rcout << "Read too big or FASTQ formatting error. Skipping read.\n";
      // Skip the rest of the current sequence
      while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL && buffer[0] != '@');
      continue;
    }
    
    quality = read_line(fastq_file, buffer);
    if (quality.empty()) break;
    
    // Process the accumulated sequence
    std::vector<AlignmentInfo> alignments(filtered_queries.size());
    
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < filtered_queries.size(); ++j) {
      const auto& query = filtered_queries[j];
      
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
      
      alignments[j] = {
        best.edit_distance,
        best.start_position,
        best.end_position,
        secondBest.edit_distance
      };
      edlibFreeAlignResult(cresult);
    }
    
    all_alignments[sequence_counter] = alignments;
    ++sequence_counter;
    ++processed_sequences;
  }
  gzclose(fastq_file);
  std::vector<int> ids;
  std::vector<std::string> query_ids;  // to store the query ID for each alignment
  std::vector<int> best_edit_distances;
  std::vector<int> best_start_positions;
  std::vector<int> best_stop_positions;
  std::vector<int> second_edit_distances;
  
  for (const auto& pair : all_alignments) {
    int sequence_id = pair.first;
    const auto& alignments = pair.second;
    for (size_t i = 0; i < alignments.size(); ++i) {
      ids.push_back(sequence_id);
      query_ids.push_back(Rcpp::as<std::string>(filtered_query_names[i]));  // get the filtered query names correctly
      best_edit_distances.push_back(alignments[i].best_edit_distance);
      best_start_positions.push_back(alignments[i].best_start_pos);
      best_stop_positions.push_back(alignments[i].best_stop_pos);
      second_edit_distances.push_back(alignments[i].second_edit_distance);
    }
  }
  Rcpp::DataFrame result = Rcpp::DataFrame::create(
    Rcpp::Named("id") = ids,
    Rcpp::Named("query_id") = query_ids,
    Rcpp::Named("best_edit_distance") = best_edit_distances,
    Rcpp::Named("best_start_pos") = best_start_positions,
    Rcpp::Named("best_stop_pos") = best_stop_positions,
    Rcpp::Named("second_edit_distance") = second_edit_distances
  );
  return result;
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigalign_stream(const std::string& fastq_path, Rcpp::CharacterVector adapters, 
  const Rcpp::DataFrame& misalignment_threshold, int nthreads = 1, int max_sequences = -1) {
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double)misal_threshold[i], (double)misal_sd[i]);
  }
  
  char buffer[131072];
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  
  while (gzgets(fastq_file, buffer, sizeof(buffer)) != Z_NULL && (max_sequences == -1 || processed_sequences < max_sequences)) {
    line = std::string(buffer);
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(1, line.find(' ') - 1); // Extract ID between '@' and first space
    if (!gzgets(fastq_file, buffer, sizeof(buffer))) break;
    sequence = std::string(buffer);
    
    if (!gzgets(fastq_file, buffer, sizeof(buffer)) || buffer[0] != '+') {
      Rcpp::stop("Invalid FASTQ format: Missing '+' line.");
    }
    
    plus = std::string(buffer);
    if (!gzgets(fastq_file, buffer, sizeof(buffer))) break;
    quality = std::string(buffer);
    
    // Remove any trailing newline characters
    id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    quality.erase(std::remove(quality.begin(), quality.end(), '\n'), quality.end());
    
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      std::string query_name = Rcpp::as<std::string>(query_names[j]);
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTail(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
        std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
        std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
        for (int& pos : start_positions) { ++pos; }
        for (int& pos : end_positions) { ++pos; }
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
#pragma omp critical
          uniqueAlignments.insert(ua);
        }
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
        std::string best_signature_part = query_name + ":" + std::to_string(best.edit_distance) + ":" + std::to_string(best.start_position) + ":" + std::to_string(best.end_position);
#pragma omp critical
        temp_signature_parts.push_back(best_signature_part);
        if (!skip_second_best) {
          std::string second_best_signature_part = query_name + ":" + std::to_string(secondBest.edit_distance) + ":" + std::to_string(secondBest.start_position) + ":" + std::to_string(secondBest.end_position);
#pragma omp critical
          temp_signature_parts.push_back(second_best_signature_part);
        }
        edlibFreeAlignResult(cresult);
      }
    }
    // Add hardcoded signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string& a, const std::string& b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      });
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
    final_signature += "<" + std::to_string(length) + ":" + id + ":undecided>";
#pragma omp critical
    signature_map[sequence_counter] = final_signature;
    
    ++sequence_counter;
    ++processed_sequences;
  }
  gzclose(fastq_file);
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigstream(
    const std::string& fastq_path, 
    Rcpp::CharacterVector adapters, 
    const Rcpp::DataFrame& read_layout, 
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads = 1, 
    int max_sequences = -1,
    bool verbose = false, 
    bool import_only = false) {
  
  // Prepare for signature alignment
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double)misal_threshold[i], (double)misal_sd[i]);
  }
  
  size_t initial_buffer_size = 131072;
  std::vector<char> buffer(initial_buffer_size);
  
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  
  // Main processing loop
  while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL && (max_sequences == -1 || processed_sequences < max_sequences)) {
    line = std::string(buffer.data());
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(0, line.find(' ') - 1); // Extract ID between '@' and first space
    if (!gzgets(fastq_file, buffer.data(), buffer.size())) break;
    sequence = std::string(buffer.data());
    
    if (!gzgets(fastq_file, buffer.data(), buffer.size()) || buffer[0] != '+') {
      Rcpp::stop("Invalid FASTQ format: Missing '+' line.");
    }
    
    plus = std::string(buffer.data());
    if (!gzgets(fastq_file, buffer.data(), buffer.size())) break;
    quality = std::string(buffer.data());
    
    // Remove any trailing newline characters
    id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    quality.erase(std::remove(quality.begin(), quality.end(), '\n'), quality.end());
    
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      std::string query_name = Rcpp::as<std::string>(query_names[j]);
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTail(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
        std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
        std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
        for (int& pos : start_positions) { ++pos; }
        for (int& pos : end_positions) { ++pos; }
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
#pragma omp critical
          uniqueAlignments.insert(ua);
        }
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
        std::string best_signature_part = query_name + ":" + std::to_string(best.edit_distance) + ":" + std::to_string(best.start_position) + ":" + std::to_string(best.end_position);
#pragma omp critical
        temp_signature_parts.push_back(best_signature_part);
        if (!skip_second_best) {
          std::string second_best_signature_part = query_name + ":" + std::to_string(secondBest.edit_distance) + ":" + std::to_string(secondBest.start_position) + ":" + std::to_string(secondBest.end_position);
#pragma omp critical
          temp_signature_parts.push_back(second_best_signature_part);
        }
        edlibFreeAlignResult(cresult);
      }
    }
    // Add hardcoded signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string& a, const std::string& b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      });
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
    final_signature += "<" + std::to_string(length) + ":" + id + ":undecided>";
#pragma omp critical
    signature_map[sequence_counter] = final_signature;
    
    ++sequence_counter;
    ++processed_sequences;
  }
  gzclose(fastq_file);
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  if(!import_only){
    container = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
    VarScan(container, verbose);
    PositionFuncMap positionFuncMap;
    positionFuncMap = createPositionFunctionMap(container, verbose);
    std::vector<std::string> sigs = Rcpp::as<std::vector<std::string>>(signature_strings);
    return process_sigstrings(container, sigs, positionFuncMap, verbose, nthreads);
  } else {
    return signature_strings;
  }
}

// [[Rcpp::export]]
Rcpp::List sigstream_combined(
    const std::string& fastq_path, 
    Rcpp::CharacterVector adapters, 
    const Rcpp::DataFrame& read_layout, 
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads = 1, 
    int max_sequences = -1,
    bool verbose = false) {
  
  // Prepare for signature alignment
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double)misal_threshold[i], (double)misal_sd[i]);
  }
  size_t initial_buffer_size = 131072;
  std::vector<char> buffer(initial_buffer_size);
  
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  // Convert R DataFrame to C++ ReadLayout structure
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  // Process variable scanning and create position function map
  VarScan(readLayout, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  // Prepare the result list
  Rcpp::List result;
  std::vector<std::string> result_names;
  // Main processing loop
  while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL && (max_sequences == -1 || processed_sequences < max_sequences)) {
    line = std::string(buffer.data());
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(1, line.find(' ') - 1); // Extract ID between '@' and first space
    if (!gzgets(fastq_file, buffer.data(), buffer.size())) break;
    sequence = std::string(buffer.data());
    
    if (!gzgets(fastq_file, buffer.data(), buffer.size()) || buffer[0] != '+') {
      Rcpp::Rcout << "Read too big or FASTQ formatting error. Skipping read.\n";
      // Skip the rest of the current sequence
      while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL && buffer[0] != '@');
      continue;
    }
    
    plus = std::string(buffer.data());
    if (!gzgets(fastq_file, buffer.data(), buffer.size())) break;
    quality = std::string(buffer.data());
    
    // Remove any trailing newline characters
    id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    quality.erase(std::remove(quality.begin(), quality.end(), '\n'), quality.end());
    
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      std::string query_name = Rcpp::as<std::string>(query_names[j]);
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTail(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
        std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
        std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
        for (int& pos : start_positions) { ++pos; }
        for (int& pos : end_positions) { ++pos; }
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
#pragma omp critical
          uniqueAlignments.insert(ua);
        }
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
        std::string best_signature_part = query_name + ":" + std::to_string(best.edit_distance) + ":" + std::to_string(best.start_position) + ":" + std::to_string(best.end_position);
#pragma omp critical
        temp_signature_parts.push_back(best_signature_part);
        if (!skip_second_best) {
          std::string second_best_signature_part = query_name + ":" + std::to_string(secondBest.edit_distance) + ":" + std::to_string(secondBest.start_position) + ":" + std::to_string(secondBest.end_position);
#pragma omp critical
          temp_signature_parts.push_back(second_best_signature_part);
        }
        edlibFreeAlignResult(cresult);
      }
    }
    // Add hardcoded signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string& a, const std::string& b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      });
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
    final_signature += "<" + std::to_string(length) + ":" + id + ":undecided>";
    
    // Convert final_signature to a vector of strings
    std::vector<std::string> sig = {final_signature};
    // Process sigstrings
    Rcpp::CharacterVector final_sig = process_sigstrings(readLayout, sig, positionFuncMap, verbose, nthreads);
    
    // Get the processed signature string back
    std::string processed_signature = Rcpp::as<std::string>(final_sig[0]);
    // Parse the signature
    auto read_info_pos = processed_signature.find_last_of('<');
    std::string lengthTypeStr = processed_signature.substr(read_info_pos + 1, processed_signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    
    // Debugging statements
    if(verbose){
      Rcpp::Rcout << "Signature: " << processed_signature << std::endl;
      Rcpp::Rcout << "Length: " << lengthStr << ", Read ID: " << read_id << ", Type: " << type << std::endl;
    }
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    
    if(type == "undecided"){
      continue;
    }
    
    // Process the signature string
    std::stringstream ss(processed_signature);
    std::string token;
    Rcpp::List read_result;
    while (std::getline(ss, token, '|')) {
      std::string id;
      int editDistance, startPos, endPos;
      std::stringstream tokenStream(token);
      
      std::getline(tokenStream, id, ':');
      tokenStream >> editDistance;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> startPos;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> endPos;
      
      auto it = readLayout.get<id_tag>().find(id);
      if (it != readLayout.get<id_tag>().end() && it->type == "variable") {
        std::string extracted_seq;
        if (startPos > 0 && endPos >= startPos && (sequence.length() >= endPos)) {
          extracted_seq = sequence.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_seq = revcomp_cpp(extracted_seq);
          }
        } else {
          continue;
        }
        std::string classId = id;
        if (classId.substr(0, 3) == "rc_") {
          classId = classId.substr(3);
        }
        read_result[classId] = extracted_seq;
      }
    }
    result.push_back(read_result);
    result_names.push_back(id);
    
    ++sequence_counter;
    ++processed_sequences;
  }
  gzclose(fastq_file);
  result.attr("names") = result_names;
  return result;
}
