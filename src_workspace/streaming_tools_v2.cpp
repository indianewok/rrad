std::pair<std::vector<std::string>, std::vector<std::string>> stream_sequences_cpp(
    const std::string& fastq_path, 
    int chunk_size, 
    int max_sequences = -1) {
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    throw std::runtime_error("Failed to open the FASTQ file.");
  }
  
  char buffer[131072];
  std::vector<std::string> chunk_ids;
  std::vector<std::string> chunk_sequences;
  int processed_sequences = 0;
  
  while (gzgets(fastq_file, buffer, sizeof(buffer)) != Z_NULL && (max_sequences == -1 || processed_sequences < max_sequences)) {
    std::string line(buffer);
    if (line[0] != '@') {
      throw std::runtime_error("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    std::string id = line.substr(1, line.find(' ') - 1);
    if (!gzgets(fastq_file, buffer, sizeof(buffer))) break;
    std::string sequence(buffer);
    
    if (!gzgets(fastq_file, buffer, sizeof(buffer)) || buffer[0] != '+') {
      throw std::runtime_error("Invalid FASTQ format: Missing '+' line.");
    }
    
    if (!gzgets(fastq_file, buffer, sizeof(buffer))) break;
    
    // Remove trailing newline characters
    id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    
    chunk_ids.push_back(id);
    chunk_sequences.push_back(sequence);
    processed_sequences++;
    
    if (chunk_ids.size() == chunk_size) {
      break;
    }
  }
  
  gzclose(fastq_file);
  return {chunk_ids, chunk_sequences};
}

std::vector<std::string> sigalign_stream_cpp(
    const std::vector<std::string>& ids, 
    const std::vector<std::string>& sequences, 
    const std::vector<std::string>& queries,
    const std::vector<std::string>& query_names,
    const std::map<std::string, 
      std::pair<double, double>>& null_dist_map, int nthreads = 1) {
  
  std::map<int, std::string> signature_map;
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const std::string& id = ids[i];
    const std::string& sequence = sequences[i];
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      std::string query_name = query_names[j];
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
    signature_map[i + 1] = final_signature;
  }
  
  std::vector<std::string> signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

void process_and_write_cpp(
    const std::vector<std::string>& chunk_ids, 
    const std::vector<std::string>& chunk_sequences,
    const std::vector<std::string>& processed_signatures,
    const ReadLayout& readLayout,
    const PositionFuncMap& positionFuncMap,
    const std::string& output_fasta_path,
    const std::string& output_sig_path, bool verbose) {
  
  std::ofstream output_fasta(output_fasta_path, std::ios::app);
  std::ofstream output_sig(output_sig_path, std::ios::app);
  
  std::vector<std::string> chunk_new_ids(chunk_ids.size());
  std::vector<std::string> chunk_extracted_sequences(chunk_ids.size());
  
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < processed_signatures.size(); ++i) {
    const std::string& processed_signature = processed_signatures[i];
    const std::string& id = chunk_ids[i];
    const std::string& sequence = chunk_sequences[i];
    
    // Parse the signature
    auto read_info_pos = processed_signature.find_last_of('<');
    std::string lengthTypeStr = processed_signature.substr(read_info_pos + 1, processed_signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    
    if (verbose) {
      std::cout << "Signature: " << processed_signature << std::endl;
      std::cout << "Length: " << lengthStr << ", Read ID: " << read_id << ", Type: " << type << std::endl;
    }
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    
    if (type == "undecided") {
      continue;
    }
    
    // Process the signature string
    std::stringstream ss(processed_signature);
    std::string token;
    std::string new_id = id;
    std::string extracted_sequence;
    while (std::getline(ss, token, '|')) {
      std::string token_id;
      int editDistance, startPos, endPos;
      std::stringstream tokenStream(token);
      
      std::getline(tokenStream, token_id, ':');
      tokenStream >> editDistance;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> startPos;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> endPos;
      
      auto it = readLayout.get<id_tag>().find(token_id);
      if (it != readLayout.get<id_tag>().end() && it->type == "variable" && it->global_class != "read") {
        if (startPos > 0 && endPos >= startPos && (sequence.length() >= endPos)) {
          std::string extracted_seq = sequence.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_seq = revcomp_cpp(extracted_seq);
          }
          new_id += "|" + token_id + ":" + extracted_seq; // Append to new_id
        }
      } else if (it != readLayout.get<id_tag>().end() && it->global_class == "read") {
        if (startPos > 0 && endPos >= startPos && (sequence.length() >= endPos)) {
          extracted_sequence = sequence.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_sequence = revcomp_cpp(extracted_sequence);
          }
        }
      }
    }
    // Store results in vectors
    chunk_new_ids[i] = new_id;
    chunk_extracted_sequences[i] = extracted_sequence;
  }
  // Write chunk results to output files
  for (size_t i = 0; i < chunk_new_ids.size(); ++i) {
    output_fasta << ">" << chunk_new_ids[i] << "\n" << chunk_extracted_sequences[i] << "\n";
    output_sig << processed_signatures[i] << "\n";
  }
}

void sigstream_chunks_cpp(
    const std::string& fastq_path, 
    const std::vector<std::string>& adapters, 
    const ReadLayout& readLayout, 
    const std::map<std::string, 
      std::pair<double, double>>& null_dist_map, 
      const PositionFuncMap& positionFuncMap, 
      const std::string& output_fasta_path, 
      const std::string& output_sig_path, 
      int nthreads = 1, 
      int max_sequences = -1, 
      int chunk_size = 1000, 
      bool verbose = false) {
  
  int processed_sequences = 0;
  
  while (max_sequences == -1 || processed_sequences < max_sequences) {
    // Stream a chunk of sequences
    auto chunk = stream_sequences_cpp(fastq_path, chunk_size, max_sequences - processed_sequences);
    auto chunk_ids = chunk.first;
    auto chunk_sequences = chunk.second;
    
    if (chunk_ids.empty()) {
      break; // No more sequences to process
    }
    // Perform alignment
    auto final_signatures = sigalign_stream_cpp(chunk_ids, chunk_sequences, adapters, query_names, null_dist_map, nthreads);
    // Process and write results
    process_and_write_cpp(chunk_ids, chunk_sequences, final_signatures, readLayout, positionFuncMap, output_fasta_path, output_sig_path, verbose);
    processed_sequences += chunk_size;
  }
}

// [[Rcpp::export]]
Rcpp::List stream_sequences(const std::string& fastq_path, int chunk_size, int max_sequences = -1) {
  auto result = stream_sequences_cpp(fastq_path, chunk_size, max_sequences);
  return Rcpp::List::create(Rcpp::Named("ids") = result.first, Rcpp::Named("sequences") = result.second);
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigalign_stream_v2(const std::vector<std::string>& ids, const std::vector<std::string>& sequences, 
  Rcpp::CharacterVector adapters, const Rcpp::DataFrame& misalignment_threshold, int nthreads = 1) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double)misal_threshold[i], (double)misal_sd[i]);
  }
  
  auto result = sigalign_stream_cpp(ids, sequences, queries, query_names, null_dist_map, nthreads);
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
void sigstream_chunks(const std::string& fastq_path, Rcpp::CharacterVector adapters, 
  const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  const std::string& output_fasta_path, const std::string& output_sig_path, 
  int nthreads = 1, int max_sequences = -1, int chunk_size = 1000, bool verbose = false) {
  // Convert R DataFrame to C++ ReadLayout structure
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  // Process variable scanning and create position function map
  VarScan(readLayout, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double)misal_threshold[i], (double)misal_sd[i]);
  }
  sigstream_chunks_cpp(fastq_path, queries, readLayout, 
    null_dist_map, positionFuncMap, output_fasta_path, 
    output_sig_path, nthreads, max_sequences, chunk_size, verbose);
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigstream_parallel(
    const std::string& fastq_path,
    Rcpp::CharacterVector adapters,
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    int nthreads = 1,
    int max_sequences = -1,
    bool verbose = false, bool import_only = false) {
  
  Rcpp::Rcout << "Initializing rapidgzip reader..." << std::endl;
  UniqueFileReader file_reader = std::make_unique<StandardFileReader>(fastq_path);
  if (!file_reader) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  rapidgzip::ParallelGzipReader<> gzipReader(std::move(file_reader), nthreads);
  
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
  
  std::vector<char> buffer(131072);
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  
  auto read_line = [&gzipReader, &buffer]() -> std::string {
    std::string result;
    size_t bytesRead = 0;
    while ((bytesRead = gzipReader.read(buffer.data(), buffer.size())) > 0) {
      result.append(buffer.data(), bytesRead);
      if (result.back() == '\n') break;
      buffer.resize(buffer.size() * 2); // Double the buffer size if not complete
    }
    result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
    return result;
  };
  
  if (verbose) {
    Rcpp::Rcout << "Starting to read sequences..." << std::endl;
  }
  
  std::vector<std::string> all_signatures;
  std::vector<int> all_sequence_ids;
  
  while ((gzipReader.read(buffer.data(), buffer.size())) > 0 &&
    (max_sequences == -1 || processed_sequences < max_sequences)) {
    
    if (verbose) {
      Rcpp::Rcout << "Reading a new sequence..." << std::endl;
    }
    
    line = std::string(buffer.data(), buffer.size());
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(1, line.find(' ') - 1); // Extract ID between '@' and first space
    if (verbose) {
      Rcpp::Rcout << "Read ID: " << id << std::endl;
    }
    
    sequence = read_line();
    if (sequence.empty()) break;
    
    plus = read_line();
    if (plus.empty() || plus[0] != '+') {
      Rcpp::Rcout << "Read too big or FASTQ formatting error. Skipping read.\n";
      // Skip the rest of the current sequence
      while ((gzipReader.read(buffer.data(), buffer.size())) > 0 && buffer[0] != '@');
      continue;
    }
    
    quality = read_line();
    if (quality.empty()) break;
    
    Rcpp::Rcout << "Processing sequence of length " << sequence.size() << std::endl;
    
    int length = sequence.size();
    std::vector<std::string> temp_signature_parts;
    
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
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
{
  all_signatures.push_back(final_signature);
  all_sequence_ids.push_back(sequence_counter);
}

++sequence_counter;
++processed_sequences;
  }
  gzipReader.close();
  
  Rcpp::Rcout << "Generating output..." << std::endl;
  Rcpp::CharacterVector signature_strings(all_signatures.size());
  for (size_t i = 0; i < all_signatures.size(); ++i) {
    signature_strings[i] = all_signatures[i];
  }
  
  if (!import_only) {
    Rcpp::Rcout << "Processing read layout..." << std::endl;
    container = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
    VarScan(container, verbose);
    PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
    std::vector<std::string> sigs = Rcpp::as<std::vector<std::string>>(signature_strings);
    return process_sigstrings(container, sigs, positionFuncMap, verbose, nthreads);
  } else {
    return signature_strings;
  }
}
