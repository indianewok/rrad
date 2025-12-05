#include "rad.h"
using namespace std;
using namespace Rcpp;

KSEQ_INIT(gzFile, gzread)

std::string extract_corrected_barcode(const std::string& full_header) {
  size_t start_pos = full_header.find("CB:Z:");
  if (start_pos == std::string::npos) return "";
  
  size_t barcode_start = start_pos + 5; // Skip "CB:Z:"
  
  // Find the next space or tab after the barcode, stopping before any UB:Z: field
  size_t barcode_end = full_header.find_first_of("\t ", barcode_start);
  if (barcode_end == std::string::npos) {
    barcode_end = full_header.length(); // In case there is no space/tab after the barcode
  }
  
  return full_header.substr(barcode_start, barcode_end - barcode_start);
}

std::string extract_blaze_id(const std::string& header) {
  size_t start_pos = header.find("#");
  size_t end_pos = header.find("_", start_pos);
  
  if (start_pos != std::string::npos && end_pos != std::string::npos) {
    return header.substr(start_pos + 1, end_pos - start_pos - 1);
  }
  
  // Return the full header if the pattern is not found
  return header;
}

// [[Rcpp::export]]
Rcpp::DataFrame extract_blaze_barcode(std::string fastq_file, int print_every = 10000) {
  gzFile fp = gzopen(fastq_file.c_str(), "r");
  if (!fp) Rcpp::stop("Failed to open file.");
  
  kseq_t *seq = kseq_init(fp);
  int l;
  
  // Hash map to store the counts of corrected barcodes
  std::unordered_map<std::string, int> barcode_counts;
  
  int sequence_count = 0;
  
  // Read through the fastq file with kseq
  while ((l = kseq_read(seq)) >= 0) {
    sequence_count++;
    
    // Capture both the sequence name (header) and the comment (which might contain CB:Z)
    std::string header = "@" + std::string(seq->name.s);  // Add '@' symbol to the header
    std::string full_header = header;
    if (seq->comment.l > 0) {
      full_header += "\t" + std::string(seq->comment.s); // Append the comment section (CB:Z might be here)
    }
    
    // Extract the corrected barcode (CB:Z: field)
    std::string corrected_barcode = extract_corrected_barcode(full_header);
    
    // If the corrected barcode is missing, skip
    if (corrected_barcode.empty()) continue;
    
    // Increment the count for the corrected barcode
    barcode_counts[corrected_barcode]++;
    
    // Print progress every 'print_every' sequences
    if (sequence_count % print_every == 0) {
      std::cout << "Processed " << sequence_count << " sequences so far." << std::endl;
    }
  }
  
  kseq_destroy(seq);
  gzclose(fp);
  
  std::cout << "Total sequences processed: " << sequence_count << std::endl;
  
  // Prepare vectors for the DataFrame
  Rcpp::CharacterVector barcodes;
  Rcpp::IntegerVector total_counts;
  for (const auto& entry : barcode_counts) {
    barcodes.push_back(entry.first);              // Corrected barcode
    total_counts.push_back(entry.second);  // Total count
  }
  
  // Return the DataFrame with barcode and total counts
  return Rcpp::DataFrame::create(Rcpp::Named("barcode") = barcodes,
    Rcpp::Named("total") = total_counts);
}

// [[Rcpp::export]]
void extract_blaze_id_bc(std::string fastq_file, std::string output_file, int print_every = 10000) {
  gzFile fp = gzopen(fastq_file.c_str(), "r");
  if (!fp) Rcpp::stop("Failed to open file.");
  
  kseq_t *seq = kseq_init(fp);
  int l;
  
  // Vectors to store the extracted IDs and barcodes
  std::vector<std::string> ids;
  std::vector<std::string> barcodes;
  
  int sequence_count = 0;
  
  // Read through the fastq file with kseq
  while ((l = kseq_read(seq)) >= 0) {
    sequence_count++;
    
    // Capture both the sequence name (header) and the comment (which might contain CB:Z)
    std::string header = "@" + std::string(seq->name.s);  // Add '@' symbol to the header
    std::string full_header = header;
    if (seq->comment.l > 0) {
      full_header += "\t" + std::string(seq->comment.s); // Append the comment section (CB:Z might be here)
    }
    
    // Extract the corrected barcode (CB:Z: field)
    std::string corrected_barcode = extract_corrected_barcode(full_header);
    
    // If the corrected barcode is missing, skip
    if (corrected_barcode.empty()) continue;
    
    // Extract the read ID part between '#' and '_'
    std::string extracted_id = extract_blaze_id(header);
    
    // Add the extracted ID and corrected barcode to the vectors
    ids.push_back(extracted_id);
    barcodes.push_back(corrected_barcode);
    
    // Print progress every 'print_every' sequences
    if (sequence_count % print_every == 0) {
      std::cout << "Processed " << sequence_count << " sequences so far." << std::endl;
    }
  }
  
  kseq_destroy(seq);
  gzclose(fp);
  
  std::cout << "Total sequences processed: " << sequence_count << std::endl;
  
  // Write the results to a CSV file
  std::ofstream out_file(output_file);
  if (!out_file.is_open()) {
    Rcpp::stop("Failed to open output file.");
  }
  
  // Write the header
  out_file << "id,barcode\n";
  
  // Write the data
  for (size_t i = 0; i < ids.size(); ++i) {
    out_file << ids[i] << "," << barcodes[i] << "\n";
  }
  
  out_file.close();
  
  std::cout << "Data successfully written to " << output_file << std::endl;
}

// [[Rcpp::export]]
Rcpp::DataFrame read_seqs_cpp(std::string filename, std::string type, bool full_id = false, int max_reads = -1) {
  std::vector<std::string> ids;
  std::vector<std::string> full_ids;
  std::vector<std::string> seqs;
  std::vector<std::string> quals;
  bool is_fastq = (type == "fq" || type == "fastq");
  
  // Open the file
  gzFile fp = gzopen(filename.c_str(), "r");
  if (!fp) {
    Rcpp::stop("Could not open file: " + filename);
  }
  
  // Initialize kseq
  kseq_t *seq = kseq_init(fp);
  int read_counter = 0;
  int result;
  
  // Read sequences
  while ((result = kseq_read(seq)) >= 0) {
    // Handle full ID vs shortened ID
    std::string full_id_str = std::string(seq->name.s);
    if (seq->comment.l) {
      full_id_str += " " + std::string(seq->comment.s);
    }
    
    std::string short_id = full_id_str;
    size_t space_pos = short_id.find(' ');
    if (space_pos != std::string::npos) {
      short_id = short_id.substr(0, space_pos);
    }
    
    // Store sequence data
    ids.push_back(short_id);
    if (full_id) {
      full_ids.push_back(full_id_str);
    }
    seqs.push_back(std::string(seq->seq.s));
    
    // Handle quality scores for FASTQ
    if (is_fastq) {
      if (seq->qual.l == 0) {
        Rcpp::stop("Missing quality scores in FASTQ file");
      }
      quals.push_back(std::string(seq->qual.s));
    }
    
    read_counter++;
    if (max_reads > 0 && read_counter >= max_reads) {
      break;
    }
  }
  
  // Clean up
  kseq_destroy(seq);
  gzclose(fp);
  
  // Return appropriate DataFrame based on format and full_id option
  if (is_fastq) {
    if (full_id) {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("full_id") = full_ids,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("qual") = quals
      );
    } else {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("qual") = quals
      );
    }
  } else {
    if (full_id) {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("full_id") = full_ids,
        Rcpp::Named("seq") = seqs
      );
    } else {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("seq") = seqs
      );
    }
  }
}

// [[Rcpp::export]]
SEXP stream_fastqas(SEXP fn, std::string type = "fq", bool full_id = false, int max_reads = -1) {
  Rcpp::Function file_exists("file.exists");
  Rcpp::Function setDT("setDT", "data.table");
  Rcpp::Function pblapply("pblapply", "pbapply");
  Rcpp::Function rbindlist("rbindlist", "data.table");
  
  // Input validation
  if (type != "fa" && type != "fasta" && type != "fq" && type != "fastq") {
    Rcpp::stop("type must be one of: fa, fasta, fq, fastq");
  }
  
  // Normalize type to internal representation
  type = (type == "fasta") ? "fa" : (type == "fastq" ? "fq" : type);
  
  // Handle single file case
  auto read_single = [&](const std::string& filename) {
    if (!Rcpp::as<bool>(file_exists(filename))) {
      Rcpp::stop("File not found: " + filename);
    }
    
    Rcpp::DataFrame res = read_seqs_cpp(filename, type, full_id, max_reads);
    return setDT(res);
  };
  
  // Check if input is a vector/list
  if (TYPEOF(fn) == STRSXP && Rf_length(fn) > 1) {
    // Handle multiple files
    Rcpp::Function identity("identity");
    Rcpp::List out = pblapply(fn, identity);
    
    std::vector<SEXP> results;
    for (int i = 0; i < out.length(); i++) {
      results.push_back(read_single(Rcpp::as<std::string>(out[i])));
    }
    
    return rbindlist(Rcpp::wrap(results));
  } else {
    // Handle single file
    return read_single(Rcpp::as<std::string>(fn));
  }
}

// [[Rcpp::export]]
Rcpp::CharacterVector masked_sigalign(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  
  // Create masked versions of adapters
  std::vector<std::pair<std::string, int>> masked_adapters;
  for(const auto& query : queries) {
    // Store original length for position correction
    int orig_length = query.length();
    
    // Create masked version - keeping last 10 bases unmasked
    // Adjust these numbers based on your needs
    int unmask_length = std::min(10, (int)query.length());
    int mask_length = query.length() - unmask_length;
    
    std::string masked_query = std::string(mask_length, 'N') + 
      query.substr(query.length() - unmask_length);
    
    masked_adapters.push_back({masked_query, orig_length});
  }
  
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for(int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = 
      std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query_name = query_names[j];
      const auto& masked_adapter = masked_adapters[j];
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        // Configure edlib to treat N as equivalent to any base
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,                    // no max edit distance
          EDLIB_MODE_HW,        // semi-global alignment
          EDLIB_TASK_LOC,       // location tasks
          additionalEqualities,  // N matches any base
          4                      // number of equality pairs
        );
        
        char* cquery = const_cast<char*>(masked_adapter.first.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        
        EdlibAlignResult cresult = edlibAlign(
          cquery, 
          masked_adapter.first.size(), 
          csequence, 
          sequence.size(), 
          config
        );
        
        // Process alignment results with position correction
        std::vector<int> start_positions(
            cresult.startLocations, 
            cresult.startLocations + cresult.numLocations
        );
        std::vector<int> end_positions(
            cresult.endLocations, 
            cresult.endLocations + cresult.numLocations
        );
        
        // Adjust positions based on original adapter length
        for (int& pos : start_positions) {
          pos = pos + 1 - (masked_adapter.first.length() - masked_adapter.second);
        }
        for (int& pos : end_positions) { ++pos; }
        
        // Rest of your existing alignment processing code...
        std::sort(start_positions.begin(), start_positions.end());
        std::sort(end_positions.begin(), end_positions.end());
        
        auto last = std::unique(start_positions.begin(), start_positions.end());
        start_positions.erase(last, start_positions.end());
        
        std::vector<int> unique_end_positions;
        for (const auto& start : start_positions) {
          auto pos = std::find(start_positions.begin(), start_positions.end(), start) - 
            start_positions.begin();
          unique_end_positions.push_back(end_positions[pos]);
        }
        
        int edit_distance = cresult.editDistance;
        std::set<UniqueAlignment> uniqueAlignments;
        
        for (size_t k = 0; k < start_positions.size(); ++k) {
          UniqueAlignment ua = {
            edit_distance, 
            start_positions[k], 
            unique_end_positions[k]
          };
          uniqueAlignments.insert(ua);
        }
        
        // Extract best alignments and process results
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? 
        *it : UniqueAlignment{-1, -1, -1}; 
        ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? 
        *it : UniqueAlignment{-1, -1, -1};
        
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
        
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + 
    std::to_string(best.edit_distance) + ":" + 
    std::to_string(best.start_position) + ":" + 
    std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + 
      std::to_string(secondBest.edit_distance) + ":" + 
      std::to_string(secondBest.start_position) + ":" + 
      std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
edlibFreeAlignResult(cresult);
      }
    }
    
    // Add your hardcoded signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + 
      std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + 
      std::to_string(length));
    
    // Sort and join signature parts
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );
    
    std::string final_signature = std::accumulate(
      temp_signature_parts.begin(), 
      temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> string {
        return a + (a.length() > 0 ? "|" : "") + b;
      }
    );
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  // Calculate total length of filtered adapters
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  for(int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    // Silently skip if sequence is too short
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,
          EDLIB_MODE_HW,
          EDLIB_TASK_LOC,
          additionalEqualities,
          4
        );
        
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
        for (size_t i = 0; i < start_positions.size(); ++i) {
          UniqueAlignment ua = {edit_distance, start_positions[i], unique_end_positions[i]};
          uniqueAlignments.insert(ua);
        }
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1}; ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + std::to_string(best.edit_distance) + ":" + std::to_string(best.start_position) + ":" + std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + std::to_string(secondBest.edit_distance) + ":" + std::to_string(secondBest.start_position) + ":" + std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::DataFrame n_masked_misalignment_cigar_string(
    const std::string& input_path,
    Rcpp::CharacterVector adapters,
    int nthreads = 1,
    int max_sequences = 500000,
    size_t chunk_size = 50000) {
  std::vector<std::string> input_files = get_fastq_files(input_path);
  if (input_files.empty()) {
    Rcpp::warning("No FASTQ files found in the specified path.");
    return Rcpp::DataFrame::create();
  }
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  
  std::vector<std::string> filtered_queries;
  Rcpp::CharacterVector filtered_query_names;
  for (int i = 0; i < query_names.size(); ++i) {
    std::string name = Rcpp::as<std::string>(query_names[i]);
    if (name.find("poly") == std::string::npos) {
      filtered_queries.push_back(queries[i]);
      filtered_query_names.push_back(query_names[i]);
    }
  }
  
  std::vector<int> ids;
  std::vector<std::string> query_ids;
  std::vector<int> best_edit_distances;
  std::vector<int> best_start_positions;
  std::vector<int> best_stop_positions;
  std::vector<int> second_edit_distances;
  std::vector<std::string> aligned_sequences;
  std::vector<int> sequence_lengths;
  std::vector<std::string> cigar_strings;
  
  Rcpp::Rcout << "Generating misalignment threshold data!\n";
  
  omp_set_num_threads(nthreads);
  
  auto process_chunk = [&](const std::vector<std::string>& chunk_sequences, const std::vector<std::string>& chunk_ids, int sequence_counter, int remaining) {
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < std::min(static_cast<int>(chunk_sequences.size()), remaining); ++i) {
      const auto& sequence = chunk_sequences[i];
      std::vector<AlignmentInfo> alignments(filtered_queries.size());
      std::vector<std::string> local_aligned_seqs(filtered_queries.size());
      std::vector<std::string> local_cigars(filtered_queries.size());
      
      for (size_t j = 0; j < filtered_queries.size(); ++j) {
        const auto& query = filtered_queries[j];
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          4
        );
        
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), sequence.c_str(), sequence.size(), config);
        
        alignments[j] = {
          result.editDistance,
          result.startLocations[0] + 1,
          result.endLocations[0] + 1,
          (result.numLocations > 1) ? result.editDistance : -1
        };
        
        int start_pos = result.startLocations[0];
        int end_pos = result.endLocations[0] + 1;
        if (start_pos >= 0 && end_pos <= sequence.length()) {
          local_aligned_seqs[j] = sequence.substr(start_pos, end_pos - start_pos);
          
          // Generate CIGAR string using edlib's function
          if (result.alignment != nullptr) {
            char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
            if (cigar != nullptr) {
              local_cigars[j] = std::string(cigar);
              free(cigar);
            } else {
              local_cigars[j] = "";
            }
          } else {
            local_cigars[j] = "";
          }
        } else {
          local_aligned_seqs[j] = "";
          local_cigars[j] = "";
        }
        
        edlibFreeAlignResult(result);
      }
      
#pragma omp critical
{
  for (size_t j = 0; j < alignments.size(); ++j) {
    ids.push_back(sequence_counter + i);
    query_ids.push_back(Rcpp::as<std::string>(filtered_query_names[j]));
    best_edit_distances.push_back(alignments[j].best_edit_distance);
    best_start_positions.push_back(alignments[j].best_start_pos);
    best_stop_positions.push_back(alignments[j].best_stop_pos);
    second_edit_distances.push_back(alignments[j].second_edit_distance);
    aligned_sequences.push_back(local_aligned_seqs[j]);
    sequence_lengths.push_back(sequence.length());
    cigar_strings.push_back(local_cigars[j]);
  }
}
    }
  };
  
  int total_sequences = 0;
  for (const auto& file_path : input_files) {
    gzFile file = gzopen(file_path.c_str(), "rb");
    if (!file) {
      Rcpp::warning("Failed to open file: " + file_path);
      continue;
    }
    kseq_t* seq = kseq_init(file);
    std::vector<std::string> chunk_sequences;
    std::vector<std::string> chunk_ids;
    while (kseq_read(seq) >= 0 && (max_sequences == -1 || total_sequences < max_sequences)) {
      chunk_sequences.push_back(seq->seq.s);
      chunk_ids.push_back(seq->name.s);
      if (chunk_sequences.size() >= chunk_size || (max_sequences != -1 && total_sequences + chunk_sequences.size() >= max_sequences)) {
        int remaining = (max_sequences == -1) ? chunk_sequences.size() : std::min(static_cast<int>(chunk_sequences.size()), max_sequences - total_sequences);
        process_chunk(chunk_sequences, chunk_ids, total_sequences, remaining);
        total_sequences += chunk_sequences.size();
        chunk_sequences.clear();
        chunk_ids.clear();
        if (max_sequences != -1 && total_sequences >= max_sequences) break;
      }
    }
    
    if (!chunk_sequences.empty() && (max_sequences == -1 || total_sequences < max_sequences)) {
      int remaining = (max_sequences == -1) ? chunk_sequences.size() : std::min(static_cast<int>(chunk_sequences.size()), max_sequences - total_sequences);
      process_chunk(chunk_sequences, chunk_ids, total_sequences, remaining);
      total_sequences += chunk_sequences.size();
    }
    
    kseq_destroy(seq);
    gzclose(file);
    
    if (max_sequences != -1 && total_sequences >= max_sequences) break;
  }
  
  if (total_sequences == 0) {
    Rcpp::warning("No sequences were processed.");
    return Rcpp::DataFrame::create();
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("id") = ids,
    Rcpp::Named("query_id") = query_ids,
    Rcpp::Named("best_edit_distance") = best_edit_distances,
    Rcpp::Named("best_start_pos") = best_start_positions,
    Rcpp::Named("best_stop_pos") = best_stop_positions,
    Rcpp::Named("second_edit_distance") = second_edit_distances,
    Rcpp::Named("aligned_sequence") = aligned_sequences,
    Rcpp::Named("sequence_length") = sequence_lengths,
    Rcpp::Named("cigar") = cigar_strings
  );
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v2(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
      misal_sd[i],
      align_mean_start[i],
      align_mean_stop[i],
      align_sd_start[i],
      align_sd_stop[i]
    };
  }
  
  // Calculate total length of filtered adapters
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    // Silently skip if sequence is too short
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        
        // Calculate search boundaries based on alignment statistics
        int search_start = std::max(1, static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0)));
        int search_end = std::min(length, static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0)));

        // Set maximum edit distance k as floor of misalignment threshold + SD
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold + bounds.misal_sd));
        
        // Create substring for searching
        std::string search_region = sequence.substr(search_start - 1, search_end - search_start + 1);
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_LOC,
          additionalEqualities,
          4
        );
        
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csearch = const_cast<char*>(search_region.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csearch, search_region.length(), config);
        
        if (cresult.startLocations != nullptr) {
          std::vector<int> start_positions(cresult.startLocations, 
            cresult.startLocations + cresult.numLocations);
          std::vector<int> end_positions(cresult.endLocations, 
            cresult.endLocations + cresult.numLocations);
          
          // Adjust positions relative to original sequence
          for (int& pos : start_positions) { pos += search_start; }
          for (int& pos : end_positions) { pos += search_start; }
          
          std::sort(start_positions.begin(), start_positions.end());
          std::sort(end_positions.begin(), end_positions.end());
          
          auto last = std::unique(start_positions.begin(), start_positions.end());
          start_positions.erase(last, start_positions.end());
          
          std::vector<int> unique_end_positions;
          for (const auto& start : start_positions) {
            auto pos = std::find(start_positions.begin(), start_positions.end(), start) - 
              start_positions.begin();
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
          
          bool skip_second_best = false;
          double threshold = bounds.misal_threshold - bounds.misal_sd;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
          
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + 
    std::to_string(best.edit_distance) + ":" + 
    std::to_string(best.start_position) + ":" + 
    std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + 
      std::to_string(secondBest.edit_distance) + ":" + 
      std::to_string(secondBest.start_position) + ":" + 
      std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair :  signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v3(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    int min_continuous_matches = 8) {
  
  // Function to check if a position in a sequence is not an N
  auto is_non_n_match = [](const std::string& query, int pos) -> bool {
    return pos >= 0 && pos < query.length() && query[pos] != 'N';
  };
  
  // Function to check for continuous non-N matches in CIGAR string
  auto has_continuous_non_n_matches = [&is_non_n_match](
    const char* cigar, 
    const std::string& query, 
    int min_matches) -> bool {
      std::string number_buffer;
      int query_pos = 0;
      int continuous_matches = 0;
      
      for(int i = 0; cigar[i] != '\0'; i++) {
        if(std::isdigit(cigar[i])) {
          number_buffer += cigar[i];
        } else {
          int count = std::stoi(number_buffer);
          
          if(cigar[i] == '=') {
            for(int j = 0; j < count; j++) {
              if(is_non_n_match(query, query_pos + j)) {
                continuous_matches++;
                if(continuous_matches >= min_matches) {
                  return true;
                }
              } else {
                continuous_matches = 0;
              }
            }
          }
          
          // Update query position based on CIGAR operation
          if(cigar[i] == '=' || cigar[i] == 'X' || cigar[i] == 'D') {
            query_pos += count;
          }
          number_buffer.clear();
        }
      }
      return false;
    };
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  // Calculate total length of filtered adapters
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  for(int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = 
      std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    // Silently skip if sequence is too short
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          4
        );
        
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), 
          csequence, sequence.size(), config);
        
        if (cresult.startLocations != nullptr && cresult.alignment != nullptr) {
          char* cigar = edlibAlignmentToCigar(cresult.alignment, 
            cresult.alignmentLength,
            EDLIB_CIGAR_EXTENDED);
          
          if(cigar != nullptr) {
            bool has_enough_matches = has_continuous_non_n_matches(cigar, query, min_continuous_matches);
            free(cigar);
            
            if(has_enough_matches) {
              std::vector<int> start_positions(cresult.startLocations, 
                cresult.startLocations + cresult.numLocations);
              std::vector<int> end_positions(cresult.endLocations, 
                cresult.endLocations + cresult.numLocations);
              
              for (int& pos : start_positions) { ++pos; }
              for (int& pos : end_positions) { ++pos; }
              
              std::sort(start_positions.begin(), start_positions.end());
              std::sort(end_positions.begin(), end_positions.end());
              
              auto last = std::unique(start_positions.begin(), start_positions.end());
              start_positions.erase(last, start_positions.end());
              
              std::vector<int> unique_end_positions;
              for (const auto& start : start_positions) {
                auto pos = std::find(start_positions.begin(), start_positions.end(), start) - 
                  start_positions.begin();
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
              
              bool skip_second_best = false;
              auto null_data = null_dist_map.find(query_name);
              if (null_data != null_dist_map.end()) {
                double threshold = null_data->second.first - null_data->second.second;
                if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
                  skip_second_best = true;
                }
              }
              
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + 
    std::to_string(best.edit_distance) + ":" + 
    std::to_string(best.start_position) + ":" + 
    std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + 
      std::to_string(secondBest.edit_distance) + ":" + 
      std::to_string(secondBest.start_position) + ":" + 
      std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
            }
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    if(!temp_signature_parts.empty()) {
      temp_signature_parts.push_back("seq_start:0:1:1");
      temp_signature_parts.push_back("rc_seq_start:0:1:1");
      temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + 
        std::to_string(length));
      temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + 
        std::to_string(length));
      
      std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
        [](const std::string &a, const std::string &b) -> bool {
          int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
          int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
          return start_pos_a < start_pos_b;
        }
      );
      
      std::string final_signature = std::accumulate(
        temp_signature_parts.begin(), 
        temp_signature_parts.end(),
        std::string(),
        [](const std::string& a, const std::string& b) -> std::string {
          return a + (a.length() > 0 ? "|" : "") + b;
        }
      );
      
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
    }
  }
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v4(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
      misal_sd[i],
      align_mean_start[i],
      align_mean_stop[i],
      align_sd_start[i],
      align_sd_stop[i]
    };
  }
  
  // Calculate total length of filtered adapters
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    // Silently skip if sequence is too short
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        
        // Calculate search boundaries based on alignment statistics
        int search_start = std::max(1, static_cast<int>(length * (bounds.start_mean/100.0 - 1.5*bounds.start_sd/100.0)));
        int search_end = std::min(length, static_cast<int>(length * (bounds.stop_mean/100.0 + 1.5*bounds.stop_sd/100.0)));
        
        // Set maximum edit distance k as floor of misalignment threshold + SD
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold + bounds.misal_sd));
        
        // Create substring for searching
        std::string search_region = sequence.substr(search_start - 1, search_end - search_start + 1);
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_LOC,
          additionalEqualities,
          4
        );
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        // Search entire sequence
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          // First pass: collect all alignments meeting edit distance threshold
          for (int k = 0; k < cresult.numLocations; k++) {
            int start_pos = cresult.startLocations[k];
            int end_pos = cresult.endLocations[k];
            int edit_dist = cresult.editDistance;
            
            // Check if alignment is within expected region
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end - 1);
            
            if (edit_dist <= bounds.misal_threshold){ // + bounds.misal_sd) {
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, start_pos + 1, end_pos + 1});
            }
          }
          
          // Only proceed if we found at least one alignment in the expected region
          if (found_in_substring && !validAlignments.empty()) {
            // Merge overlapping alignments
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              // Check for overlap
              if (it->start_position <= current.end_position) {
                // Merge alignments
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
#pragma omp critical
{
  // Report all merged alignments
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
  }
}
          }
        }
        edlibFreeAlignResult(cresult);      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair :  signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v5(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
                     misal_sd[i],
                             align_mean_start[i],
                                             align_mean_stop[i],
                                                            align_sd_start[i],
                                                                          align_sd_stop[i]
    };
  }
  
  // Initialize adapter info with N counts
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      const auto& query = queries[j];
      int n_count = std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, static_cast<int>(query.length())};
    }
  }
  
  // Calculate total length of filtered adapters
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    // Silently skip if sequence is too short
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        auto adapter = adapter_info[query_name];
        
        // Calculate search boundaries based on alignment statistics
        int search_start = std::max(1, static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0)));
        int search_end = std::min(length, static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0)));
        
        // Set maximum edit distance k as floor of misalignment threshold + SD
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold + bounds.misal_sd));
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          4
        ); 
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        // Search entire sequence
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          // First pass: collect all alignments meeting criteria
          for (int k = 0; k < cresult.numLocations; k++) {
            int start_pos = cresult.startLocations[k];
            int end_pos = cresult.endLocations[k];
            int edit_dist = cresult.editDistance;
            
            // Get CIGAR string and calculate matches
            int total_matches = 0;
            if (cresult.alignment != nullptr) {
              char* cigar = edlibAlignmentToCigar(cresult.alignment, 
                                                  cresult.alignmentLength, 
                                                  EDLIB_CIGAR_EXTENDED);
              if (cigar != nullptr) {
                // Parse CIGAR for matches
                std::string number_buffer;
                for(int c = 0; cigar[c] != '\0'; c++) {
                  if(std::isdigit(cigar[c])) {
                    number_buffer += cigar[c];
                  } else if(cigar[c] == '=' || cigar[c] == 'M') {
                    total_matches += std::stoi(number_buffer);
                    number_buffer.clear();
                  } else {
                    number_buffer.clear();
                  }
                }
                free(cigar);
              }
            }
            
            // Adjust matches for N-masked positions
            int adjusted_matches = total_matches - adapter.n_count;
            
            // Calculate probability
            double p_match_at_position = std::pow(0.25, adjusted_matches);
            double p_no_match_at_position = 1 - p_match_at_position;
            int possible_positions = length - adjusted_matches + 1;
            double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
            double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
            bool is_significant = p_at_least_one_match <= 1.0;
            
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end - 1);
            
            if (edit_dist <= bounds.misal_threshold && 
                is_significant && 
                adjusted_matches > 0) {
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, start_pos + 1, end_pos + 1});
            }
          }
          
          // Only proceed if we found at least one alignment in the expected region
          if (found_in_substring && !validAlignments.empty()) {
            // Merge overlapping alignments
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              if (it->start_position <= current.end_position) {
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
#pragma omp critical
{
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
  }
}
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v6(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    bool verbose = false) {
  
  if(verbose) std::cout << "\n=== Starting Adapter Processing ===\n";
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  if(verbose) {
    std::cout << "Processing " << sequences.size() << " sequences with " 
              << queries.size() << " adapters\n\n";
  }
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  if(verbose) std::cout << "=== Initializing Alignment Bounds ===\n";
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
                     misal_sd[i],
                             align_mean_start[i],
                                             align_mean_stop[i],
                                                            align_sd_start[i],
                                                                          align_sd_stop[i]
    };
    if(verbose) {
      std::cout << "Adapter " << id << ":\n"
                << "  Misalignment threshold: " << misal_threshold[i] << " ± " << misal_sd[i] << "\n"
                << "  Alignment bounds: " << align_mean_start[i] << "% - " << align_mean_stop[i] << "%\n";
    }
  }
  
  if(verbose) std::cout << "\n=== Counting N's in Adapters ===\n";
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      const auto& query = queries[j];
      int n_count = std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, static_cast<int>(query.length())};
      if(verbose) {
        std::cout << "Adapter " << query_names[j] << ":\n"
                  << "  Sequence: " << query << "\n"
                  << "  N count: " << n_count << "\n"
                  << "  Total length: " << query.length() << "\n";
      }
    }
  }
  
  size_t total_adapter_length = 0;
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }
  
  if(verbose) std::cout << "\n=== Starting Sequence Processing ===\n";
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    int length = sequence.size();
    
    if(verbose) {
      std::cout << "\nProcessing sequence " << i + 1 << " (Length: " << length << ")\n";
      std::cout << "Sequence ID: " << ids[i] << "\n";
    }
    
    if(length < total_adapter_length) {
      if(verbose) std::cout << "  Sequence too short, skipping\n";
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if(verbose) {
        std::cout << "\n  Processing adapter: " << query_name << "\n";
        std::cout << "  Adapter sequence: " << query << "\n";
      }
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
          if(verbose) std::cout << "  Found poly-" << poly_base << " tail\n";
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        auto adapter = adapter_info[query_name];
        
        int search_start = std::max(1, static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0)));
        int search_end = std::min(length, static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0)));
        
        if(verbose) {
          std::cout << "  Search bounds: " << search_start << " - " << search_end << "\n";
          std::cout << "  N-count adjustment: " << adapter.n_count << "\n";
        }
        
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold + bounds.misal_sd));
        
        if(verbose) {
          std::cout << "  Max edit distance: " << max_edit_distance << "\n";
        }
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          4
        );
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          if(verbose) {
            std::cout << "  Found " << cresult.numLocations << " potential alignments\n";
          }
          
          for (int k = 0; k < cresult.numLocations; k++) {
            int start_pos = cresult.startLocations[k];
            int end_pos = cresult.endLocations[k];
            int edit_dist = cresult.editDistance;
            
            if(verbose) {
              std::cout << "\n    Alignment " << k + 1 << ":\n";
              std::cout << "    Position: " << start_pos + 1 << " - " << end_pos + 1 << "\n";
              std::cout << "    Edit distance: " << edit_dist << "\n";
            }
            
            int total_matches = 0;
            if (cresult.alignment != nullptr) {
              char* cigar = edlibAlignmentToCigar(cresult.alignment, 
                                                  cresult.alignmentLength, 
                                                  EDLIB_CIGAR_EXTENDED);
              if (cigar != nullptr) {
                if(verbose) std::cout << "    CIGAR: " << cigar << "\n";
                
                std::string number_buffer;
                for(int c = 0; cigar[c] != '\0'; c++) {
                  if(std::isdigit(cigar[c])) {
                    number_buffer += cigar[c];
                  } else if(cigar[c] == '=' || cigar[c] == 'M') {
                    total_matches += std::stoi(number_buffer);
                    number_buffer.clear();
                  } else {
                    number_buffer.clear();
                  }
                } 
                free(cigar);
              }
            }
            
            int adjusted_matches = total_matches - adapter.n_count;
            
            if(verbose) {
              std::cout << "    Total matches: " << total_matches << "\n";
              std::cout << "    Adjusted matches: " << adjusted_matches << "\n";
            }
            
            double p_match_at_position = std::pow(0.25, adjusted_matches);
            double p_no_match_at_position = 1 - p_match_at_position;
            int possible_positions = length - adjusted_matches + 1;
            double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
            double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
            bool is_significant = p_at_least_one_match < 1.0;
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end - 1);
            
            if(verbose) {
              std::cout << "      Detected within expected sequence range: " << (in_substring ? "YES" : "NO") << "\n";
              std::cout << "      Probability calculations:\n";
              std::cout << "      P(match at position): " << p_match_at_position << "\n";
              std::cout << "      P(at least one match): " << p_at_least_one_match << "%\n";
              std::cout << "      Significant: " << (is_significant ? "YES" : "NO") << "\n";
            }

            if ((edit_dist <= static_cast<int>(std::floor(bounds.misal_threshold)) && is_significant && in_substring) ||
              ((edit_dist < static_cast<int>(std::ceil(bounds.misal_threshold-bounds.misal_sd))) && is_significant && !in_substring)){
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, start_pos + 1, end_pos + 1});
              if(verbose) std::cout << "    Alignment ACCEPTED\n";
            } else {
              if(verbose) std::cout << "    Alignment REJECTED\n";
            }
          }
          
          if (found_in_substring && !validAlignments.empty()) {
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              if (it->start_position <= current.end_position) {
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
            if(verbose) {
              std::cout << "\n  Final merged alignments: " << mergedAlignments.size() << "\n";
            }
            
#pragma omp critical
{
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
    if(verbose) {
      std::cout << "  Adding signature: " << signature_part << "\n";
    }
  }
}
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  if(verbose) {
    std::cout << "\nFinal signature for sequence " << i + 1 << ":\n";
    std::cout << final_signature << "\n";
    std::cout << "----------------------------------------\n";
  }
  temp_signature_parts.clear();
}
  }
  
  if(verbose) std::cout << "\n=== Processing Complete ===\n";
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v7(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    bool verbose = false) {
  
  if(verbose) std::cout << "\n=== Starting Adapter Processing ===\n";
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  if(verbose) {
    std::cout << "Processing " << sequences.size() << " sequences with " 
              << queries.size() << " adapters\n\n";
  }
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector misal_mean_start = misalignment_threshold["misal_mean_start"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double misal_mean_start;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  if(verbose) std::cout << "=== Initializing Alignment Bounds ===\n";
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
      misal_sd[i],
      misal_mean_start[i],
      align_mean_start[i],
      align_mean_stop[i],
      align_sd_start[i],
      align_sd_stop[i]
    };
    if(verbose) {
      std::cout << "Adapter " << id << ":\n"
                << "  Misalignment threshold: " << misal_threshold[i] << " ± " << misal_sd[i] << "\n"
                << "  Alignment bounds: " << align_mean_start[i] << "% - " << align_mean_stop[i] << "%\n";
    }
  }
  
  if(verbose) std::cout << "\n=== Counting N's in Adapters ===\n";
  for(size_t j = 0; j < query_names.size(); ++j) {
    if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      const auto& query = queries[j];
      int n_count = std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, static_cast<int>(query.length())};
      if(verbose) {
        std::cout << "Adapter " << query_names[j] << ":\n"
                  << "  Sequence: " << query << "\n"
                  << "  N count: " << n_count << "\n"
                  << "  Total length: " << query.length() << "\n";
      }
    }
  }
  
  // size_t total_adapter_length = 0;
  // for(size_t j = 0; j < query_names.size(); ++j) {
  //   if(query_names[j] != "poly_a" && query_names[j] != "poly_t") {
  //     total_adapter_length += queries[j].length();
  //   }
  // }
  
  if (verbose) std::cout << "\n=== Counting N's in Adapters ===\n";
  int longest_adapter_length = 0;
  for (size_t j = 0; j < query_names.size(); ++j) {
    const auto& query = queries[j];
    int length = (int)query.size();
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      int n_count = (int)std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, length};
      if (length > longest_adapter_length) longest_adapter_length = length;
      if (verbose) {
        std::cout << "Adapter " << query_names[j] << ":\n"
                  << "  Sequence: " << query << "\n"
                  << "  N count: " << n_count << "\n"
                  << "  Total length: " << length << "\n";
      }
    } else {
      adapter_info[query_names[j]] = {0, length};
    }
  }
  
  size_t total_adapter_length = 0;
  for (size_t jj = 0; jj < query_names.size(); ++jj) {
    if (query_names[jj] != "poly_a" && query_names[jj] != "poly_t") {
      total_adapter_length += queries[jj].length();
    }
  }
  
  if (verbose) std::cout << "\n=== Determining Padding Size ===\n";
  // Heuristic: half of longest adapter length, at least 10
  int pad_size = std::max(10, longest_adapter_length / 2);
  std::string padding(pad_size, 'n');
  if (verbose) std::cout << "Using " << pad_size << " 'n' characters as padding on each end.\n";
  
  if (verbose) std::cout << "\n=== Starting Sequence Processing ===\n";
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& original_sequence = sequences[i];
    int length = original_sequence.size();
    
    if(verbose) {
      std::cout << "\nProcessing sequence " << i + 1 << " (Length: " << length << ")\n";
      std::cout << "Sequence ID: " << ids[i] << "\n";
    }
    
    if(length < total_adapter_length) {
      if(verbose) std::cout << "  Sequence too short, skipping\n";
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    std::string sequence = padding + original_sequence + padding;
    int padded_length = length + 2 * pad_size;
    
    if(verbose) {
      std::cout << "\n Padding sequence: original length is " << length << " and padded length is " << padded_length;
    }
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if(verbose) {
        std::cout << "\n  Processing adapter: " << query_name << "\n";
        std::cout << "  Adapter sequence: " << query << "\n";
      }
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(original_sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
          if(verbose) std::cout << "  Found poly-" << poly_base << " tail\n";
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        auto adapter = adapter_info[query_name];
        
        int thresholded_start = static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0));
        int thresholded_stop = static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0));
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold + bounds.misal_sd));
        
        int search_start = std::max(1, thresholded_start);
        int search_end = std::min(length, thresholded_stop);
        
        if(verbose) {
          std::cout << "  Search bounds: " << search_start << " - " << search_end << "\n";
          std::cout << "  N-count adjustment: " << adapter.n_count << "\n";
          std::cout << "  Max edit distance: " << max_edit_distance << "\n";
        }
        
        EdlibEqualityPair additionalEqualities[9] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, 
          {'n', 'A'}, {'n', 'C'}, {'n', 'G'}, {'n', 'T'}, 
          {'N', 'n'} //as far as i'm concerned this is the only critical one--otherwise,
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          9
        );
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          if(verbose) {
            std::cout << "  Found " << cresult.numLocations << " potential alignments\n";
          }
          
          for (int k = 0; k < cresult.numLocations; k++) {
            int padded_start_pos = cresult.startLocations[k];
            int padded_end_pos = cresult.endLocations[k];

            // Convert padded coordinates back to original
            int orig_start_zero_based = padded_start_pos - pad_size;
            int orig_end_zero_based = padded_end_pos - pad_size;
            
            int start_pos = orig_start_zero_based + 1;
            int end_pos = orig_end_zero_based + 1;
            
            int truncation_cover = std::max(0, -orig_start_zero_based) 
              + std::max(0, orig_end_zero_based - (length - 1));
            
            // Clamp coordinates to the original sequence boundaries
            if (start_pos < 1) start_pos = 1;
            if (end_pos > length) end_pos = length;
            int edit_dist = cresult.editDistance;
            
            if(verbose) {
              std::cout << "\n    Alignment " << k + 1 << ":\n";
              std::cout << "    Padded position: " << padded_start_pos + 1 << " - " << padded_end_pos + 1 << "\n";
              std::cout << "    Position: " << start_pos << " - " << end_pos << "\n";
              std::cout << "    Edit distance: " << edit_dist << "\n";
            }
            
            int total_matches = 0;
            if (cresult.alignment != nullptr) {
              char* cigar = edlibAlignmentToCigar(cresult.alignment, 
                                                  cresult.alignmentLength, 
                                                  EDLIB_CIGAR_EXTENDED);
              if (cigar != nullptr) {
                if(verbose) std::cout << "    CIGAR: " << cigar << "\n";
                
                std::string number_buffer;
                for(int c = 0; cigar[c] != '\0'; c++) {
                  if(std::isdigit(cigar[c])) {
                    number_buffer += cigar[c];
                  } else if(cigar[c] == '=' || cigar[c] == 'M') {
                    total_matches += std::stoi(number_buffer);
                    number_buffer.clear();
                  } else {
                    number_buffer.clear();
                  }
                } 
                free(cigar);
              }
            }
            
            int adjusted_matches = total_matches - adapter.n_count - truncation_cover;
            
            if(verbose) {
              std::cout << "    Total matches: " << total_matches << "\n";
              std::cout << "    Adapter N-adjustment: " << adapter.n_count << "\n";
              std::cout << "    Sequence n-padding adjustment: " << truncation_cover << "\n";
              std::cout << "    Adjusted matches: " << adjusted_matches << "\n";
            }
            
            double p_match_at_position = std::pow(0.25, adjusted_matches);
            double p_no_match_at_position = 1 - p_match_at_position;
            int possible_positions = length - adjusted_matches + 1;
            double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
            double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
            bool is_significant = p_at_least_one_match < 1.0;
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end);
            
            if(verbose) {
              std::cout << "      Detected within expected sequence range: " << (in_substring ? "YES" : "NO") << "\n";
              std::cout << "      Probability calculations:\n";
              std::cout << "      P(match at position): " << p_match_at_position << "\n";
              std::cout << "      P(at least one match): " << p_at_least_one_match << "%\n";
              std::cout << "      Significant: " << (is_significant ? "YES" : "NO") << "\n";
            }
            
            if ((edit_dist <= static_cast<int>(std::floor(bounds.misal_threshold)) && is_significant && in_substring) ||
                ((edit_dist < static_cast<int>(std::ceil(bounds.misal_threshold-bounds.misal_sd))) && is_significant && !in_substring)){
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, std::max(1, start_pos), std::min(length,end_pos)});
              if(verbose) std::cout << "    Alignment ACCEPTED\n";
            } else {
              if(verbose) std::cout << "    Alignment REJECTED\n";
            }
          }
          
          if (found_in_substring && !validAlignments.empty()) {
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              if (it->start_position <= current.end_position) {
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
            if(verbose) {
              std::cout << "\n  Final merged alignments: " << mergedAlignments.size() << "\n";
            }
            
#pragma omp critical
{
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
    if(verbose) {
      std::cout << "  Adding signature: " << signature_part << "\n";
    }
  }
}
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  if(verbose) {
    std::cout << "\nFinal signature for sequence " << i + 1 << ":\n";
    std::cout << final_signature << "\n";
    std::cout << "----------------------------------------\n";
  }
  temp_signature_parts.clear();
}
  }
  
  if(verbose) std::cout << "\n=== Processing Complete ===\n";
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}


//as of 12-11-24 this is the version that works (v8)

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v8(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    bool verbose = false) {
  
  if(verbose) std::cout << "\n=== Starting Adapter Processing ===\n";
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  if(verbose) {
    std::cout << "Processing " << sequences.size() << " sequences with " 
              << queries.size() << " adapters\n\n";
  }
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  struct AlignmentResult {
    int total_matches;
    int overlapping_n_count;
    int longest_clean_match;
    int current_clean_match;
  };
  
  auto parse_cigar = [&](const char* cigar_cstr,
                         const std::string& adapter_seq,
                         const std::string& target_seq,
                         int start_pos,
                         int end_pos) -> AlignmentResult {
                           AlignmentResult result = {0, 0, 0, 0};  // Initialize all fields to 0
                           
                           if (start_pos < 0 || end_pos > target_seq.size() || end_pos <= start_pos) {
                             if (verbose) {
                               std::cout << "Invalid alignment positions: " << start_pos << "-" << end_pos << "\n";
                             }
                             return result;
                           }
                           
                           std::string aligned_region = target_seq.substr(start_pos, end_pos - start_pos + 1);
                           
                           int adapter_idx = 0;
                           int target_idx = 0;
                           std::string cigar(cigar_cstr);
                           std::string aligned_adapter;
                           std::string aligned_target;
                           
                           if (verbose) {
                             std::cout << "      Processing CIGAR: " << cigar << "\n";
                             std::cout << "      Adapter sequence: " << adapter_seq << "\n";
                             std::cout << "      Target region: " << aligned_region << "\n";
                           }
                           
                           size_t i = 0;
                           while (i < cigar.length()) {
                             int count = 0;
                             while (i < cigar.length() && std::isdigit(cigar[i])) {
                               count = count * 10 + (cigar[i] - '0');
                               i++;
                             }
                             
                             if (i >= cigar.length()) break;
                             
                             char op = cigar[i++];
                             
                             switch (op) {
                             case 'M':
                             case '=':
                             case 'X':
                               for (int c = 0; c < count; ++c) {
                                 if (adapter_idx < adapter_seq.size() && target_idx < aligned_region.size()) {
                                   // Check if this position is an N/n overlap
                                   bool is_nn_overlap = (adapter_seq[adapter_idx] == 'N' && 
                                                         aligned_region[target_idx] == 'n');
                                   
                                   if (is_nn_overlap) {
                                     result.overlapping_n_count++;
                                     // Reset current clean match streak
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                     
                                     if (verbose) {
                                       std::cout << "      Found N/n overlap at adapter pos " << adapter_idx 
                                                 << " target pos " << target_idx 
                                                 << ", resetting clean match streak\n";
                                     }
                                   } else if (adapter_seq[adapter_idx] == aligned_region[target_idx]) {
                                     // If bases match and it's not an N/n overlap, increment clean match streak
                                     result.current_clean_match++;
                                     if (verbose && result.current_clean_match > result.longest_clean_match) {
                                       std::cout << "      New longest clean match: " << result.current_clean_match 
                                                 << " at adapter pos " << adapter_idx << "\n";
                                     }
                                   } else {
                                     // Mismatch - reset clean match streak
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                   }
                                   
                                   aligned_adapter += adapter_seq[adapter_idx];
                                   aligned_target += aligned_region[target_idx];
                                   result.total_matches++;
                                 }
                                 adapter_idx++;
                                 target_idx++;
                               }
                               break;
                               
                             case 'I':
                             case 'D':
                               // Any insertion or deletion breaks the clean match streak
                               result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                     result.current_clean_match);
                               result.current_clean_match = 0;
                               
                               if (op == 'I') {
                                 for (int c = 0; c < count; ++c) {
                                   if (adapter_idx < adapter_seq.size()) {
                                     aligned_adapter += adapter_seq[adapter_idx];
                                     aligned_target += '-';
                                   }
                                   adapter_idx++;
                                 }
                               } else { // 'D'
                                 for (int c = 0; c < count; ++c) {
                                   if (target_idx < aligned_region.size()) {
                                     aligned_adapter += '-';
                                     aligned_target += aligned_region[target_idx];
                                   }
                                   target_idx++;
                                 }
                               }
                               break;
                             }
                           }
                           
                           // Don't forget to check one last time at the end
                           result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                 result.current_clean_match);
                           
                           if (verbose) {
                             std::cout << "      Aligned Adapter: " << aligned_adapter << "\n";
                             std::cout << "      Aligned Target:  " << aligned_target << "\n";
                             std::cout << "      Total matches: " << result.total_matches << "\n";
                             std::cout << "      N/n overlaps: " << result.overlapping_n_count << "\n";
                             std::cout << "      Longest clean match: " << result.longest_clean_match << "\n";
                           }
                           
                           return result;
                         }; 
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  if(verbose) std::cout << "=== Initializing Alignment Bounds ===\n";
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
                     misal_sd[i],
                     align_mean_start[i],
                     align_mean_stop[i],
                     align_sd_start[i],
                     align_sd_stop[i]
    };
    if(verbose) {
      std::cout << "Adapter " << id << ":\n"
                << "  Misalignment threshold: " << misal_threshold[i] << " ± " << misal_sd[i] << "\n"
                << "  Alignment bounds: " << align_mean_start[i] << "% - " << align_mean_stop[i] << "%\n";
    }
  }
  
  if (verbose) std::cout << "\n=== Counting N's in Adapters ===\n";
  int longest_adapter_length = 0;
  for (size_t j = 0; j < query_names.size(); ++j) {
    const auto& query = queries[j];
    int length = (int)query.size();
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      int n_count = (int)std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, length};
      if (length > longest_adapter_length) longest_adapter_length = length;
      if (verbose) {
        std::cout << "Adapter " << query_names[j] << ":\n"
                  << "  Sequence: " << query << "\n"
                  << "  N count: " << n_count << "\n"
                  << "  Total length: " << length << "\n";
      }
    } else {
      adapter_info[query_names[j]] = {0, length};
    }
  }
  size_t total_adapter_length = 0;
  for (size_t jj = 0; jj < query_names.size(); ++jj) {
    if (query_names[jj] != "poly_a" && query_names[jj] != "poly_t") {
      total_adapter_length += queries[jj].length();
    }
  }
  
  if (verbose) std::cout << "\n=== Determining Padding Size ===\n";
  // Heuristic: half of longest adapter length, at least 10
  int pad_size = std::max(10, longest_adapter_length / 2);
  std::string padding(pad_size, 'n');
  if (verbose) std::cout << "Using " << pad_size << " 'n' characters as padding on each end.\n";
  
  if (verbose) std::cout << "\n=== Starting Sequence Processing ===\n";
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& original_sequence = sequences[i];
    int length = original_sequence.size();
    
    if(verbose) {
      std::cout << "\nProcessing sequence " << i + 1 << " (Length: " << length << ")\n";
      std::cout << "Sequence ID: " << ids[i] << "\n";
    }
    
    if(length < total_adapter_length) {
      if(verbose) std::cout << "  Sequence too short, skipping\n";
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if(verbose) {
        std::cout << "\n  Processing adapter: " << query_name << "\n";
        std::cout << "  Adapter sequence: " << query << "\n";
      }
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(original_sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
          if(verbose) std::cout << "  Found poly-" << poly_base << " tail\n";
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        auto adapter = adapter_info[query_name];
        
        int front_pad = (bounds.start_mean <= 25.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
        int end_pad   = (bounds.start_mean >= 75.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
        
        // Construct the modified sequence
        std::string front_padding_str = (front_pad > 0) ? padding : "";
        std::string end_padding_str   = (end_pad > 0)   ? padding : "";
        std::string sequence = front_padding_str + original_sequence + end_padding_str;
      
        int thresholded_start = static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0));
        int thresholded_stop = static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0));
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold));// + bounds.misal_sd));
        
        int search_start = std::max(1, thresholded_start);
        int search_end = std::min(length, thresholded_stop);
        
        if(verbose) {
          std::cout << "  Search bounds: " << search_start << " - " << search_end << "\n";
          std::cout << "  N-count adjustment: " << adapter.n_count << "\n";
          std::cout << "  Max edit distance: " << max_edit_distance << "\n";
        }
        
        EdlibEqualityPair additionalEqualities[9] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, 
          {'n', 'A'}, {'n', 'C'}, {'n', 'G'}, {'n', 'T'}, 
          {'N', 'n'} //as far as i'm concerned this is the only critical one
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          9
        );
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          if(verbose) {
            std::cout << "  Found " << cresult.numLocations << " potential alignments\n";
          }
          
          for (int k = 0; k < cresult.numLocations; k++) {
            int padded_start_pos = cresult.startLocations[k];
            int padded_end_pos = cresult.endLocations[k];
            // Convert padded coordinates back to original
            
            int orig_start_zero_based = (front_pad > 0) ? (padded_start_pos - front_pad) : padded_start_pos;
            int orig_end_zero_based   = (front_pad > 0) ? (padded_end_pos - front_pad) : padded_end_pos;
            
            int start_pos = orig_start_zero_based + 1;
            int end_pos = orig_end_zero_based + 1;
            
            int truncation_cover = (front_pad > 0) 
              ? std::max(0, -orig_start_zero_based)
                : std::max(0, orig_end_zero_based - (length - 1));
            
            // Clamp coordinates to the original sequence boundaries
            if (start_pos < 1) start_pos = 1;
            if (end_pos > length) end_pos = length;
            int edit_dist = cresult.editDistance;
            
            if(verbose) {
              std::cout << "\n    Alignment " << k + 1 << ":\n";
              std::cout << "    Padded position: " << padded_start_pos + 1 << " - " << padded_end_pos + 1 << "\n";
              std::cout << "    Position: " << start_pos << " - " << end_pos << "\n";
              std::cout << "    Edit distance: " << edit_dist << "\n";
            }
            
            int total_matches = 0;
            int longest_match = 0;
            int overlapping_n_count = 0;
            if (cresult.alignment != nullptr) {
              char* cigar = edlibAlignmentToCigar(cresult.alignment, 
                                                  cresult.alignmentLength, 
                                                  EDLIB_CIGAR_EXTENDED);
              if (cigar != nullptr) {
                if(verbose) std::cout << "    CIGAR: " << cigar << "\n";
                AlignmentResult parsed_result = parse_cigar(cigar, query, sequence, padded_start_pos, padded_end_pos);

                total_matches = parsed_result.total_matches;
                longest_match = parsed_result.longest_clean_match;
                overlapping_n_count = parsed_result.overlapping_n_count;
                free(cigar);
              }
            }
            
            int adjusted_matches = total_matches - adapter.n_count - truncation_cover + overlapping_n_count;
            
            if(verbose) {
              std::cout << "    Total matches: " << total_matches << "\n";
              std::cout << "    Adapter N-adjustment: " << adapter.n_count << "\n";
              std::cout << "    Sequence n-padding adjustment: " << truncation_cover << "\n";
              std::cout << "    Overlapping N/n count: " << overlapping_n_count << "\n";
              std::cout << "    Total adjusted matches: " << adjusted_matches << "\n";
            }
            
            double p_match_at_position = std::pow(0.25, adjusted_matches);
            double p_no_match_at_position = 1 - p_match_at_position;
            int possible_positions = length - adjusted_matches + 1;
            double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
            double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
            bool is_significant = p_at_least_one_match < 1.0;
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end);
            
            if(verbose) {
              std::cout << "      Detected within expected sequence range: " << (in_substring ? "YES" : "NO") << "\n";
              std::cout << "      Probability calculations:\n";
              std::cout << "      P(match at position): " << p_match_at_position << "\n";
              std::cout << "      P(at least one match): " << p_at_least_one_match << "%\n";
              std::cout << "      Significant: " << (is_significant ? "YES" : "NO") << "\n";
            }
            
            if ((edit_dist <= static_cast<int>(std::floor(bounds.misal_threshold)) && is_significant && in_substring) ||
                ((edit_dist < static_cast<int>(std::ceil(bounds.misal_threshold-bounds.misal_sd))) && is_significant && !in_substring)){
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, std::max(1, start_pos), std::min(length,end_pos)});
              if(verbose) std::cout << "    Alignment ACCEPTED\n";
            } else {
              if(verbose) std::cout << "    Alignment REJECTED\n";
            }
          }
          
          if (found_in_substring && !validAlignments.empty()) {
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              if (it->start_position <= current.end_position) {
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
            if(verbose) {
              std::cout << "\n  Final merged alignments: " << mergedAlignments.size() << "\n";
            }
            
#pragma omp critical
{
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
    if(verbose) {
      std::cout << "  Adding signature: " << signature_part << "\n";
    }
  }
}
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  if(verbose) {
    std::cout << "\nFinal signature for sequence " << i + 1 << ":\n";
    std::cout << final_signature << "\n";
    std::cout << "----------------------------------------\n";
  }
  temp_signature_parts.clear();
}
  }
  
  if(verbose) std::cout << "\n=== Processing Complete ===\n";
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v9(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    bool verbose = false) {
  
  // === Structures ===
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  struct AlignmentResult {
    int total_matches;
    int overlapping_n_count;
    int longest_clean_match;
    int current_clean_match;
  };
  
  struct ProbabilityMetrics {
    double p_match_at_position;
    double p_at_least_one_match;
    bool is_significant;
  };
  
  struct AlignmentDetails {
    int edit_distance;
    int start_position;
    int end_position;
    int longest_clean_match;
    bool is_valid;
    bool in_substring;
  };
  
  // === Helper Functions ===
  auto calculate_match_probability = [](int match_length, int sequence_length) -> ProbabilityMetrics {
    double p_match_at_position = std::pow(0.25, match_length);
    double p_no_match_at_position = 1 - p_match_at_position;
    int possible_positions = sequence_length - match_length + 1;
    double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
    double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
    
    return {
      p_match_at_position,
      p_at_least_one_match,
      p_at_least_one_match < 1.0
    };
  };
  
  auto parse_cigar_old = [&](const char* cigar_cstr,
                         const std::string& adapter_seq,
                         const std::string& target_seq,
                         int start_pos,
                         int end_pos) -> AlignmentResult {
                           AlignmentResult result = {0, 0, 0, 0};
                           
                           if (start_pos < 0 || end_pos > target_seq.size() || end_pos <= start_pos) {
                             if (verbose) std::cout << "Invalid alignment positions: " << start_pos << "-" << end_pos << "\n";
                             return result;
                           }
                           
                           std::string aligned_region = target_seq.substr(start_pos, end_pos - start_pos + 1);
                           int adapter_idx = 0;
                           int target_idx = 0;
                           std::string cigar(cigar_cstr);
                           std::string aligned_adapter, aligned_target;
                           
                           if (verbose) {
                             std::cout << "      Processing CIGAR: " << cigar << "\n";
                             std::cout << "      Adapter sequence: " << adapter_seq << "\n";
                             std::cout << "      Target region: " << aligned_region << "\n";
                           }
                           
                           size_t i = 0;
                           while (i < cigar.length()) {
                             int count = 0;
                             while (i < cigar.length() && std::isdigit(cigar[i])) {
                               count = count * 10 + (cigar[i] - '0');
                               i++;
                             }
                             
                             if (i >= cigar.length()) break;
                             
                             char op = cigar[i++];
                             
                             switch (op) {
                             case 'M':
                             case '=':
                             case 'X':
                               for (int c = 0; c < count; ++c) {
                                 if (adapter_idx < adapter_seq.size() && target_idx < aligned_region.size()) {
                                   bool is_nn_overlap = (adapter_seq[adapter_idx] == 'N' && 
                                                         aligned_region[target_idx] == 'n');
                                   
                                   if (is_nn_overlap) {
                                     result.overlapping_n_count++;
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                   } else if (adapter_seq[adapter_idx] == aligned_region[target_idx]) {
                                     result.current_clean_match++;
                                   } else {
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                   }
                                   
                                   aligned_adapter += adapter_seq[adapter_idx];
                                   aligned_target += aligned_region[target_idx];
                                   result.total_matches++;
                                 }
                                 adapter_idx++;
                                 target_idx++;
                               }
                               break;
                               
                             case 'I':
                             case 'D':
                               result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                     result.current_clean_match);
                               result.current_clean_match = 0;
                               
                               if (op == 'I') {
                                 for (int c = 0; c < count; ++c) {
                                   if (adapter_idx < adapter_seq.size()) {
                                     aligned_adapter += adapter_seq[adapter_idx];
                                     aligned_target += '-';
                                   }
                                   adapter_idx++;
                                 }
                               } else {
                                 for (int c = 0; c < count; ++c) {
                                   if (target_idx < aligned_region.size()) {
                                     aligned_adapter += '-';
                                     aligned_target += aligned_region[target_idx];
                                   }
                                   target_idx++;
                                 }
                               }
                               break;
                             }
                           }
                           
                           result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                 result.current_clean_match);
                           
                           if (verbose) {
                             std::cout << "      Aligned Adapter: " << aligned_adapter << "\n";
                             std::cout << "      Aligned Target:  " << aligned_target << "\n";
                             std::cout << "      Total matches: " << result.total_matches << "\n";
                             std::cout << "      N/n overlaps: " << result.overlapping_n_count << "\n";
                             std::cout << "      Longest clean match: " << result.longest_clean_match << "\n";
                           }
                           
                           return result;
                         };
  
  auto parse_cigar = [&](const char* cigar_cstr,
                         const std::string& adapter_seq,
                         const std::string& target_seq,
                         int start_pos,
                         int end_pos) -> AlignmentResult {
                           AlignmentResult result = {0, 0, 0, 0};
                           
                           if (start_pos < 0 || end_pos > target_seq.size() || end_pos <= start_pos) {
                             if (verbose) std::cout << "Invalid alignment positions: " << start_pos << "-" << end_pos << "\n";
                             return result;
                           }
                           
                           std::string aligned_region = target_seq.substr(start_pos, end_pos - start_pos + 1);
                           int adapter_idx = 0;
                           int target_idx = 0;
                           std::string cigar(cigar_cstr);
                           std::string aligned_adapter, aligned_target;
                           
                           if (verbose) {
                             std::cout << "      Processing CIGAR: " << cigar << "\n";
                             std::cout << "      Adapter sequence: " << adapter_seq << "\n";
                             std::cout << "      Target region: " << aligned_region << "\n";
                           }
                           
                           size_t i = 0;
                           while (i < cigar.length()) {
                             int count = 0;
                             while (i < cigar.length() && std::isdigit(cigar[i])) {
                               count = count * 10 + (cigar[i] - '0');
                               i++;
                             }
                             
                             if (i >= cigar.length()) break;
                             
                             char op = cigar[i++];
                             
                             switch (op) {
                             case 'M':
                             case '=':
                             case 'X':
                               for (int c = 0; c < count; ++c) {
                                 if (adapter_idx < adapter_seq.size() && target_idx < aligned_region.size()) {
                                   bool is_nn_overlap = (adapter_seq[adapter_idx] == 'N' && 
                                                         aligned_region[target_idx] == 'n');
                                   
                                   if (is_nn_overlap) {
                                     result.overlapping_n_count++;
                                     // Don't break the clean match streak for N/n overlaps
                                     result.current_clean_match++;
                                   } else if (adapter_seq[adapter_idx] == aligned_region[target_idx]) {
                                     // Matching bases increase the clean match streak
                                     result.current_clean_match++;
                                   } else {
                                     // Only mismatches break the clean match streak
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                   }
                                   
                                   aligned_adapter += adapter_seq[adapter_idx];
                                   aligned_target += aligned_region[target_idx];
                                   result.total_matches++;
                                 }
                                 adapter_idx++;
                                 target_idx++;
                               }
                               break;
                               
                             case 'I':
                               // Insertions don't break the clean match streak
                               for (int c = 0; c < count; ++c) {
                                 if (adapter_idx < adapter_seq.size()) {
                                   aligned_adapter += adapter_seq[adapter_idx];
                                   aligned_target += '-';
                                 }
                                 adapter_idx++;
                               }
                               break;
                               
                             case 'D':
                               // Deletions don't break the clean match streak
                               for (int c = 0; c < count; ++c) {
                                 if (target_idx < aligned_region.size()) {
                                   aligned_adapter += '-';
                                   aligned_target += aligned_region[target_idx];
                                 }
                                 target_idx++;
                               }
                               break;
                             }
                           }
                           
                           // Final check for the clean match streak
                           result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                 result.current_clean_match);
                           
                           if (verbose) {
                             std::cout << "      Aligned Adapter: " << aligned_adapter << "\n";
                             std::cout << "      Aligned Target:  " << aligned_target << "\n";
                             std::cout << "      Total matches: " << result.total_matches << "\n";
                             std::cout << "      N/n overlaps: " << result.overlapping_n_count << "\n";
                             std::cout << "      Longest clean match: " << result.longest_clean_match << "\n";
                           }
                           
                           return result;
                         };
  
  auto process_alignment = [&](const EdlibAlignResult& result, 
                               int k,
                               const std::string& query,
                               const std::string& sequence,
                               int front_pad,
                               int length,
                               int search_start,
                               int search_end,
                               bool verbose) -> AlignmentDetails {
                                 int padded_start_pos = result.startLocations[k];
                                 int padded_end_pos = result.endLocations[k];
                                 
                                 // Convert coordinates
                                 int orig_start_zero_based = (front_pad > 0) ? (padded_start_pos - front_pad) : padded_start_pos;
                                 int orig_end_zero_based = (front_pad > 0) ? (padded_end_pos - front_pad) : padded_end_pos;
                                
                                 int start_pos = std::max(1, orig_start_zero_based + 1);
                                 int end_pos = std::min(length, orig_end_zero_based + 1);
                                 
                                 AlignmentResult parsed_result = {0, 0, 0, 0};
                                 if (result.alignment != nullptr) {
                                   char* cigar = edlibAlignmentToCigar(result.alignment, 
                                                                       result.alignmentLength, 
                                                                       EDLIB_CIGAR_EXTENDED);
                                   if (cigar != nullptr) {
                                     parsed_result = parse_cigar(
                                       cigar,
                                       query,
                                       sequence,
                                       padded_start_pos,
                                       padded_end_pos);
                                     free(cigar);
                                   }
                                 }
                                 bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end);
                                 auto prob_metrics = calculate_match_probability(
                                   parsed_result.longest_clean_match,
                                   length);
                                 
                                 return {
                                   result.editDistance,
                                   start_pos,
                                   end_pos,
                                   parsed_result.longest_clean_match,
                                   prob_metrics.is_significant,
                                   in_substring
                                 };
                               };
  
  // === Initialize Variables ===
  if(verbose) std::cout << "\n=== Starting Adapter Processing ===\n";
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  // Extract threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  // Initialize alignment bounds and adapter info
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  // Set up alignment bounds
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
      misal_sd[i],
      align_mean_start[i], 
      align_mean_stop[i],
      align_sd_start[i],
      align_sd_stop[i]
    };
  }
  
  // Process adapters and count N's
  int longest_adapter_length = 0;
  size_t total_adapter_length = 0;
  
  for (size_t j = 0; j < query_names.size(); ++j) {
    const auto& query = queries[j];
    int length = query.size();
    
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      int n_count = std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, length};
      longest_adapter_length = std::max(longest_adapter_length, length);
      total_adapter_length += length;
    } else {
      adapter_info[query_names[j]] = {0, length};
    }
  }
  
  // Set up padding
  int pad_size = std::max(10, longest_adapter_length / 2);
  std::string padding(pad_size, 'n');
  
  // === Main Processing Loop ===
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& original_sequence = sequences[i];
    int length = original_sequence.size();
    
    if(length < total_adapter_length) {
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
      continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    // Process each query
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(original_sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
        continue;
      }
      
      // Regular adapter processing
      auto bounds = alignment_bounds[query_name];
      auto adapter = adapter_info[query_name];
      
      // Set up padding and search bounds
      int front_pad = (bounds.start_mean <= 25.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
      int end_pad = (bounds.start_mean >= 75.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
      std::string sequence = std::string(front_pad, 'n') + original_sequence + std::string(end_pad, 'n');
      
      int search_start = std::max(1, static_cast<int>(length * (bounds.start_mean/100.0 - bounds.start_sd/100.0)));
      int search_end = std::min(length, static_cast<int>(length * (bounds.stop_mean/100.0 + bounds.stop_sd/100.0)));
      
      // Perform alignment
      EdlibEqualityPair additionalEqualities[9] = {
        {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, 
        {'n', 'A'}, {'n', 'C'}, {'n', 'G'}, {'n', 'T'}, 
        {'N', 'n'}
      };
      
      EdlibAlignConfig config = edlibNewAlignConfig(
        static_cast<int>(
        std::floor(bounds.misal_threshold)),
        EDLIB_MODE_HW,
        EDLIB_TASK_PATH,
        additionalEqualities,
        9
      );
      
      EdlibAlignResult alignment_result = edlibAlign(
        query.c_str(),
        query.size(),
        sequence.c_str(),
        sequence.length(),
        config
      );
      
      if (alignment_result.startLocations != nullptr) {
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        for (int k = 0; k < alignment_result.numLocations; k++) {
          auto details = process_alignment(alignment_result, k, query, sequence,
                                           front_pad, length, search_start, search_end, verbose);
          
          if ((details.edit_distance <= static_cast<int>(std::floor(bounds.misal_threshold)) && 
              details.is_valid && details.in_substring) ||
              ((details.edit_distance < static_cast<int>(std::ceil(bounds.misal_threshold-bounds.misal_sd))) && 
              details.is_valid && !details.in_substring)) {
            
            if (details.in_substring) found_in_substring = true;
            validAlignments.insert({details.edit_distance, details.start_position, details.end_position});
          }
        }
        
        // Process valid alignments
        if (found_in_substring && !validAlignments.empty()) {
          // First, merge overlapping alignments
          std::vector<UniqueAlignment> mergedAlignments;
          auto it = validAlignments.begin();
          UniqueAlignment current = *it++;
          
          while (it != validAlignments.end()) {
            if (it->start_position <= current.end_position) {
              current.end_position = std::max(current.end_position, it->end_position);
              current.edit_distance = std::min(current.edit_distance, it->edit_distance);
            } else {
              mergedAlignments.push_back(current);
              current = *it;
            }
            ++it;
          }
          mergedAlignments.push_back(current);
          
#pragma omp critical
{
  // Now process the merged alignments instead of the raw valid alignments
  for (const auto& alignment : mergedAlignments) {
    temp_signature_parts.push_back(
      query_name + ":" + 
        std::to_string(alignment.edit_distance) + ":" + 
        std::to_string(alignment.start_position) + ":" + 
        std::to_string(alignment.end_position)
    );
    if(verbose) {
      std::cout << "\n  Final merged alignments: " << mergedAlignments.size() << "\n";
      std::cout << "  Adding signature: " << query_name << ":" 
                << alignment.edit_distance << ":" 
                << alignment.start_position << ":" 
                << alignment.end_position << "\n";
    }
  }
}
        }
      }
      edlibFreeAlignResult(alignment_result);
    }
    
    // Add standard signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    // Sort and join signature parts
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) {
                return std::stoi(a.substr(a.find_last_of(":") + 1)) < 
                  std::stoi(b.substr(b.find_last_of(":") + 1));
              });
    
#pragma omp critical
{
  std::string final_signature = std::accumulate(
    temp_signature_parts.begin(), 
    temp_signature_parts.end(),
    std::string(),
    [](const std::string& a, const std::string& b) {
      return a + (a.length() > 0 ? "|" : "") + b;
    }
  );
  
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
}
  }
  
  // Convert results to R
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign_v10(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads,
    bool verbose = false) {
  
  if(verbose) std::cout << "\n=== Starting Adapter Processing ===\n";
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  std::map<int, std::string> signature_map;
  
  if(verbose) {
    std::cout << "Processing " << sequences.size() << " sequences with " 
              << queries.size() << " adapters\n\n";
  }
  
  // Extract all threshold data
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };
  
  struct AdapterInfo {
    int n_count;
    int total_length;
  };
  
  struct AlignmentResult {
    int total_matches;
    int overlapping_n_count;
    int longest_clean_match;
    int current_clean_match;
  };
  
  auto parse_cigar = [&](const char* cigar_cstr,
                         const std::string& adapter_seq,
                         const std::string& target_seq,
                         int start_pos,
                         int end_pos) -> AlignmentResult {
                           AlignmentResult result = {0, 0, 0, 0};  // Initialize all fields to 0
                           
                           if (start_pos < 0 || end_pos > target_seq.size() || end_pos <= start_pos) {
                             if (verbose) {
                               std::cout << "Invalid alignment positions: " << start_pos << "-" << end_pos << "\n";
                             }
                             return result;
                           }
                           
                           std::string aligned_region = target_seq.substr(start_pos, end_pos - start_pos + 1);
                           
                           int adapter_idx = 0;
                           int target_idx = 0;
                           std::string cigar(cigar_cstr);
                           std::string aligned_adapter;
                           std::string aligned_target;
                           
                           if (verbose) {
                             std::cout << "      Processing CIGAR: " << cigar << "\n";
                             std::cout << "      Adapter sequence: " << adapter_seq << "\n";
                             std::cout << "      Target region: " << aligned_region << "\n";
                           }
                           
                           size_t i = 0;
                           while (i < cigar.length()) {
                             int count = 0;
                             while (i < cigar.length() && std::isdigit(cigar[i])) {
                               count = count * 10 + (cigar[i] - '0');
                               i++;
                             }
                             
                             if (i >= cigar.length()) break;
                             
                             char op = cigar[i++];
                             
                             switch (op) {
                             case 'M':
                             case '=':
                             case 'X':
                               for (int c = 0; c < count; ++c) {
                                 if (adapter_idx < adapter_seq.size() && target_idx < aligned_region.size()) {
                                   // Check if this position is an N/n overlap
                                   bool is_nn_overlap = (adapter_seq[adapter_idx] == 'N' && 
                                                         aligned_region[target_idx] == 'n');
                                   
                                   if (is_nn_overlap) {
                                     result.overlapping_n_count++;
                                     // Reset current clean match streak
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                     
                                     if (verbose) {
                                       std::cout << "      Found N/n overlap at adapter pos " << adapter_idx 
                                                 << " target pos " << target_idx 
                                                 << ", resetting clean match streak\n";
                                     }
                                   } else if (adapter_seq[adapter_idx] == aligned_region[target_idx]) {
                                     // If bases match and it's not an N/n overlap, increment clean match streak
                                     result.current_clean_match++;
                                     if (verbose && result.current_clean_match > result.longest_clean_match) {
                                       std::cout << "      New longest clean match: " << result.current_clean_match 
                                                 << " at adapter pos " << adapter_idx << "\n";
                                     }
                                   } else {
                                     // Mismatch - reset clean match streak
                                     result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                           result.current_clean_match);
                                     result.current_clean_match = 0;
                                   }
                                   
                                   aligned_adapter += adapter_seq[adapter_idx];
                                   aligned_target += aligned_region[target_idx];
                                   result.total_matches++;
                                 }
                                 adapter_idx++;
                                 target_idx++;
                               }
                               break;
                               
                             case 'I':
                             case 'D':
                               // Any insertion or deletion breaks the clean match streak
                               result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                     result.current_clean_match);
                               result.current_clean_match = 0;
                               
                               if (op == 'I') {
                                 for (int c = 0; c < count; ++c) {
                                   if (adapter_idx < adapter_seq.size()) {
                                     aligned_adapter += adapter_seq[adapter_idx];
                                     aligned_target += '-';
                                   }
                                   adapter_idx++;
                                 }
                               } else { // 'D'
                                 for (int c = 0; c < count; ++c) {
                                   if (target_idx < aligned_region.size()) {
                                     aligned_adapter += '-';
                                     aligned_target += aligned_region[target_idx];
                                   }
                                   target_idx++;
                                 }
                               }
                               break;
                             }
                           }
                           
                           // Don't forget to check one last time at the end
                           result.longest_clean_match = std::max(result.longest_clean_match, 
                                                                 result.current_clean_match);
                           
                           if (verbose) {
                             std::cout << "      Aligned Adapter: " << aligned_adapter << "\n";
                             std::cout << "      Aligned Target:  " << aligned_target << "\n";
                             std::cout << "      Total matches: " << result.total_matches << "\n";
                             std::cout << "      N/n overlaps: " << result.overlapping_n_count << "\n";
                             std::cout << "      Longest clean match: " << result.longest_clean_match << "\n";
                           }
                           
                           return result;
                         }; 
  
  std::map<std::string, AlignmentBounds> alignment_bounds;
  std::map<std::string, AdapterInfo> adapter_info;
  
  if(verbose) std::cout << "=== Initializing Alignment Bounds ===\n";
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
                     misal_sd[i],
                             align_mean_start[i],
                                             align_mean_stop[i],
                                                            align_sd_start[i],
                                                                          align_sd_stop[i]
    };
    if(verbose) {
      std::cout << "Adapter " << id << ":\n"
                << "  Misalignment threshold: " << misal_threshold[i] << " ± " << misal_sd[i] << "\n"
                << "  Alignment bounds: " << align_mean_start[i] << "% - " << align_mean_stop[i] << "%\n";
    }
  }
  
  if (verbose) std::cout << "\n=== Counting N's in Adapters ===\n";
  int longest_adapter_length = 0;
  for (size_t j = 0; j < query_names.size(); ++j) {
    const auto& query = queries[j];
    int length = (int)query.size();
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      int n_count = (int)std::count(query.begin(), query.end(), 'N');
      adapter_info[query_names[j]] = {n_count, length};
      if (length > longest_adapter_length) longest_adapter_length = length;
      if (verbose) {
        std::cout << "Adapter " << query_names[j] << ":\n"
                  << "  Sequence: " << query << "\n"
                  << "  N count: " << n_count << "\n"
                  << "  Total length: " << length << "\n";
      }
    } else {
      adapter_info[query_names[j]] = {0, length};
    }
  }
  size_t total_adapter_length = 0;
  for (size_t jj = 0; jj < query_names.size(); ++jj) {
    if (query_names[jj] != "poly_a" && query_names[jj] != "poly_t") {
      total_adapter_length += queries[jj].length();
    }
  }
  
  if (verbose) std::cout << "\n=== Determining Padding Size ===\n";
  // Heuristic: half of longest adapter length, at least 10
  int pad_size = std::max(10, longest_adapter_length / 2);
  std::string padding(pad_size, 'n');
  if (verbose) std::cout << "Using " << pad_size << " 'n' characters as padding on each end.\n";
  
  if (verbose) std::cout << "\n=== Starting Sequence Processing ===\n";
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& original_sequence = sequences[i];
    int length = original_sequence.size();
    
    if(verbose) {
      std::cout << "\nProcessing sequence " << i + 1 << " (Length: " << length << ")\n";
      std::cout << "Sequence ID: " << ids[i] << "\n";
    }
    
    if(length < total_adapter_length) {
      if(verbose) std::cout << "  Sequence too short, skipping\n";
#pragma omp critical
{
  signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
}
continue;
    }
    
    std::vector<std::string> temp_signature_parts;
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      
      if(verbose) {
        std::cout << "\n  Processing adapter: " << query_name << "\n";
        std::cout << "  Adapter sequence: " << query << "\n";
      }
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(original_sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
          if(verbose) std::cout << "  Found poly-" << poly_base << " tail\n";
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        auto bounds = alignment_bounds[query_name];
        auto adapter = adapter_info[query_name];
        
        int front_pad = (bounds.start_mean <= 25.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
        int end_pad   = (bounds.start_mean >= 75.0 && pad_size <= (adapter.total_length/2)) ? pad_size : 0;
        
        // Construct the modified sequence
        std::string front_padding_str = (front_pad > 0) ? padding : "";
        std::string end_padding_str   = (end_pad > 0)   ? padding : "";
        std::string sequence = front_padding_str + original_sequence + end_padding_str;
        
        int thresholded_start = static_cast<int>(length * (bounds.start_mean/100.0 - 1*bounds.start_sd/100.0));
        int thresholded_stop = static_cast<int>(length * (bounds.stop_mean/100.0 + 1*bounds.stop_sd/100.0));
        int max_edit_distance = static_cast<int>(std::floor(bounds.misal_threshold));// + bounds.misal_sd));
        
        int search_start = std::max(1, thresholded_start);
        int search_end = std::min(length, thresholded_stop);
        
        if(verbose) {
          std::cout << "  Search bounds: " << search_start << " - " << search_end << "\n";
          std::cout << "  N-count adjustment: " << adapter.n_count << "\n";
          std::cout << "  Max edit distance: " << max_edit_distance << "\n";
        }
        
        EdlibEqualityPair additionalEqualities[9] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, 
          {'n', 'A'}, {'n', 'C'}, {'n', 'G'}, {'n', 'T'}, 
          {'N', 'n'} //as far as i'm concerned this is the only critical one
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          max_edit_distance,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          9
        );
        
        std::set<UniqueAlignment> validAlignments;
        bool found_in_substring = false;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* cseq = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), cseq, sequence.length(), config);
        
        if (cresult.startLocations != nullptr) {
          if(verbose) {
            std::cout << "  Found " << cresult.numLocations << " potential alignments\n";
          }
          
          for (int k = 0; k < cresult.numLocations; k++) {
            int padded_start_pos = cresult.startLocations[k];
            int padded_end_pos = cresult.endLocations[k];
            // Convert padded coordinates back to original
            
            int orig_start_zero_based = (front_pad > 0) ? (padded_start_pos - front_pad) : padded_start_pos;
            int orig_end_zero_based   = (front_pad > 0) ? (padded_end_pos - front_pad) : padded_end_pos;
            
            int start_pos = orig_start_zero_based + 1;
            int end_pos = orig_end_zero_based + 1;
            
            int truncation_cover = (front_pad > 0) 
              ? std::max(0, -orig_start_zero_based)
                : std::max(0, orig_end_zero_based - (length - 1));
            
            // Clamp coordinates to the original sequence boundaries
            if (start_pos < 1) start_pos = 1;
            if (end_pos > length) end_pos = length;
            int edit_dist = cresult.editDistance;
            
            if(verbose) {
              std::cout << "\n    Alignment " << k + 1 << ":\n";
              std::cout << "    Padded position: " << padded_start_pos + 1 << " - " << padded_end_pos + 1 << "\n";
              std::cout << "    Position: " << start_pos << " - " << end_pos << "\n";
              std::cout << "    Edit distance: " << edit_dist << "\n";
            }
            
            int total_matches = 0;
            int longest_match = 0;
            int overlapping_n_count = 0;
            if (cresult.alignment != nullptr) {
              char* cigar = edlibAlignmentToCigar(cresult.alignment, 
                                                  cresult.alignmentLength, 
                                                  EDLIB_CIGAR_EXTENDED);
              if (cigar != nullptr) {
                if(verbose) std::cout << "    CIGAR: " << cigar << "\n";
                AlignmentResult parsed_result = parse_cigar(cigar, query, sequence, padded_start_pos, padded_end_pos);
                
                total_matches = parsed_result.total_matches;
                overlapping_n_count = parsed_result.overlapping_n_count;
                longest_match = parsed_result.longest_clean_match;
                free(cigar);
              }
            }
            
            // int adjusted_matches = (front_pad > 0) 
            //   ? (total_matches - truncation_cover + overlapping_n_count)
            //   : (total_matches - adapter.n_count - truncation_cover + overlapping_n_count);
            
            int adjusted_matches = total_matches - adapter.n_count - truncation_cover + overlapping_n_count;
            
            if(verbose) {
              std::cout << "    Total matches: " << total_matches << "\n";
              // if (front_pad > 0) {
              //   std::cout << "    Sequence n-padding adjustment: " << truncation_cover << "\n";
              // } else {
              std::cout << "    Adapter N-adjustment: " << adapter.n_count << "\n";
              std::cout << "    Sequence n-padding adjustment: " << truncation_cover << "\n";
              //}
              std::cout << "    Overlapping N/n count: " << overlapping_n_count << "\n";
              std::cout << "    Total adjusted matches: " << adjusted_matches << "\n";
              
            }
            
            double p_match_at_position = std::pow(0.25, adjusted_matches);
            double p_no_match_at_position = 1 - p_match_at_position;
            int possible_positions = length - adjusted_matches + 1;
            double p_no_match_anywhere = std::pow(p_no_match_at_position, possible_positions);
            double p_at_least_one_match = (1 - p_no_match_anywhere) * 100;
            bool is_significant = p_at_least_one_match < 1.0;
            bool in_substring = (start_pos >= search_start - 1 && end_pos <= search_end);
            
            if(verbose) {
              std::cout << "      Detected within expected sequence range: " << (in_substring ? "YES" : "NO") << "\n";
              std::cout << "      Probability calculations:\n";
              std::cout << "      P(match at position): " << p_match_at_position << "\n";
              std::cout << "      P(at least one match): " << p_at_least_one_match << "%\n";
              std::cout << "      Significant: " << (is_significant ? "YES" : "NO") << "\n";
            }
            
            if ((edit_dist <= static_cast<int>(std::floor(bounds.misal_threshold)) && is_significant && in_substring) ||
                ((edit_dist < static_cast<int>(std::ceil(bounds.misal_threshold-bounds.misal_sd))) && is_significant && !in_substring)){
              if (in_substring) found_in_substring = true;
              validAlignments.insert({edit_dist, std::max(1, start_pos), std::min(length,end_pos)});
              if(verbose) std::cout << "    Alignment ACCEPTED\n";
            } else {
              if(verbose) std::cout << "    Alignment REJECTED\n";
            }
          }
          
          if (found_in_substring && !validAlignments.empty()) {
            std::vector<UniqueAlignment> mergedAlignments;
            auto it = validAlignments.begin();
            UniqueAlignment current = *it++;
            
            while (it != validAlignments.end()) {
              if (it->start_position <= current.end_position) {
                current.end_position = std::max(current.end_position, it->end_position);
                current.edit_distance = std::min(current.edit_distance, it->edit_distance);
              } else {
                mergedAlignments.push_back(current);
                current = *it;
              }
              ++it;
            }
            mergedAlignments.push_back(current);
            
            if(verbose) {
              std::cout << "\n  Final merged alignments: " << mergedAlignments.size() << "\n";
            }
            
#pragma omp critical
{
  for (const auto& alignment : mergedAlignments) {
    std::string signature_part = query_name + ":" + 
      std::to_string(alignment.edit_distance) + ":" + 
      std::to_string(alignment.start_position) + ":" + 
      std::to_string(alignment.end_position);
    temp_signature_parts.push_back(signature_part);
    if(verbose) {
      std::cout << "  Adding signature: " << signature_part << "\n";
    }
  }
}
          }
        }
        edlibFreeAlignResult(cresult);
      }
    }
    
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
              [](const std::string &a, const std::string &b) -> bool {
                int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
                int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
                return start_pos_a < start_pos_b;
              }
    );
    
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
                                                  std::string(),
                                                  [](const std::string& a, const std::string& b) -> std::string {
                                                    return a + (a.length() > 0 ? "|" : "") + b;
                                                  });
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  if(verbose) {
    std::cout << "\nFinal signature for sequence " << i + 1 << ":\n";
    std::cout << final_signature << "\n";
    std::cout << "----------------------------------------\n";
  }
  temp_signature_parts.clear();
}
  }
  
  if(verbose) std::cout << "\n=== Processing Complete ===\n";
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

// [[Rcpp::export]]
NumericMatrix cmv_consensusMatrix(SEXP sequences) {
  CharacterVector seqs(sequences);
  int n_seq = seqs.length();
  if (n_seq == 0) return NumericMatrix(0);
  
  // Get length of first sequence and convert to std::string
  int seq_length = 0;
  for (int i = 0; i < n_seq; i++) {
    std::string curr_seq = std::string(seqs[i]);
    seq_length = std::max(seq_length, (int)curr_seq.length());
  }  
  // Initialize 4x seq_length matrix
  NumericMatrix consensus(4, seq_length);
  rownames(consensus) = CharacterVector::create("A", "C", "G", "T");
  // Count occurrences of each base at each position
  for (int i = 0; i < n_seq; i++) {
    std::string curr_seq = std::string(seqs[i]);
    
    for (int j = 0; j < seq_length; j++) {
      char base = curr_seq[j];
      switch(base) {
      case 'A':
      case 'a':
        consensus(0, j)++;
        break;
      case 'C':
      case 'c':
        consensus(1, j)++;
        break;
      case 'G':
      case 'g':
        consensus(2, j)++;
        break;
      case 'T':
      case 't':
        consensus(3, j)++;
        break;
      }
    }
  }
  return consensus;
}

// [[Rcpp::export]]
Rcpp::IntegerVector find_longest_matches(Rcpp::StringVector cigars) {
  std::vector<int> longest_matches;
  longest_matches.reserve(cigars.length());
  
  for(const auto& cigar : cigars) {
    std::string cigar_str = Rcpp::as<std::string>(cigar);
    int current_length = 0;
    int max_length = 0;
    std::string number_buffer;
    
    for(size_t i = 0; i < cigar_str.length(); i++) {
      char c = cigar_str[i];
      
      if(std::isdigit(c)) {
        number_buffer += c;
      } else {
        // Process operation
        if(c == '=') {
          // Found a match operation
          current_length = std::stoi(number_buffer);
          max_length = std::max(max_length, current_length);
        }
        number_buffer.clear();
      }
    }
    
    longest_matches.push_back(max_length);
  }
  
  return Rcpp::wrap(longest_matches);
}

// [[Rcpp::export]]
Rcpp::List analyze_cigar_patterns(Rcpp::StringVector cigars) {
  std::vector<int> longest_matches;
  std::vector<int> total_matches;
  std::vector<int> match_blocks;
  longest_matches.reserve(cigars.length());
  total_matches.reserve(cigars.length());
  match_blocks.reserve(cigars.length());
  
  for(const auto& cigar : cigars) {
    std::string cigar_str = Rcpp::as<std::string>(cigar);
    int current_length = 0;
    int max_length = 0;
    int total_match = 0;
    int blocks = 0;
    std::string number_buffer;
    
    for(size_t i = 0; i < cigar_str.length(); i++) {
      char c = cigar_str[i];
      
      if(std::isdigit(c)) {
        number_buffer += c;
      } else {
        int num = std::stoi(number_buffer);
        if(c == '=') {
          // Found a match operation
          current_length = num;
          max_length = std::max(max_length, current_length);
          total_match += num;
          blocks++;
        }
        number_buffer.clear();
      }
    }
    
    longest_matches.push_back(max_length);
    total_matches.push_back(total_match);
    match_blocks.push_back(blocks);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("longest_match") = Rcpp::wrap(longest_matches),
    Rcpp::Named("total_matches") = Rcpp::wrap(total_matches),
    Rcpp::Named("match_blocks") = Rcpp::wrap(match_blocks)
  );
}

// [[Rcpp::export]]
bool has_minimum_match_length(const std::string& cigar, int min_length) {
  std::string number_buffer;
  
  for(size_t i = 0; i < cigar.length(); i++) {
    char c = cigar[i];
    
    if(std::isdigit(c)) {
      number_buffer += c;
    } else {
      if(c == '=' && std::stoi(number_buffer) >= min_length) {
        return true;
      }
      number_buffer.clear();
    }
  }
  
  return false;
}

// [[Rcpp::export]]
List edlib_search(std::vector<std::string> queries, 
                                 std::vector<std::string> targets,
                                 int max_distance = 3,
                                 int num_threads = 0) {
  
  // Set number of threads (0 = auto-detect)
  if (num_threads <= 0) {
    num_threads = omp_get_max_threads();
  }
  omp_set_num_threads(num_threads);
  
  Rcout << "Starting multithreaded search with " << num_threads << " threads..." << std::endl;
  Rcout << "Processing " << queries.size() << " queries against " 
        << targets.size() << " targets (max distance: " << max_distance << ")" << std::endl;
  
  // Pre-allocate result list
  List results(queries.size());
  CharacterVector query_names(queries.size());
  
  // Progress tracking variables (need to be thread-safe)
  std::vector<int> progress_counter(num_threads, 0);
  int total_matches = 0;
  
  // Parallel loop over queries
#pragma omp parallel for schedule(dynamic)
  for (int q = 0; q < (int)queries.size(); q++) {
    
    int thread_id = omp_get_thread_num();
    progress_counter[thread_id]++;
    
    // Progress reporting (only from thread 0 to avoid race conditions)
    const std::string& query = queries[q];
    const char* query_cstr = query.c_str();
    int query_len = query.length();
    
    std::vector<std::string> matched_targets;
    std::vector<int> matched_distances;
    
    // Pre-create edlib config (thread-local)
    EdlibAlignConfig config = edlibNewAlignConfig(max_distance, 
                                                  EDLIB_MODE_NW, 
                                                  EDLIB_TASK_DISTANCE, 
                                                  NULL, 0);
    
    // Inner loop over targets (not parallelized to avoid nested parallelism)
    for (size_t t = 0; t < targets.size(); t++) {
      
      const std::string& target = targets[t];
      const char* target_cstr = target.c_str();
      int target_len = target.length();
      
      // Quick length-based filteringd
      if (abs(query_len - target_len) > max_distance) {
        continue;
      }
      
      EdlibAlignResult result = edlibAlign(query_cstr, query_len,
                                           target_cstr, target_len,
                                           config);
      if (result.editDistance != -1 && result.editDistance <= max_distance) {
        matched_targets.push_back(target);
        matched_distances.push_back(result.editDistance);
        total_matches++;
      }
      edlibFreeAlignResult(result);
    }
    // Create DataFrame for this query's matches (thread-safe since each thread writes to different index)
    if (matched_targets.size() > 0) {
      DataFrame query_matches = DataFrame::create(Named("target") = matched_targets,
                                                  Named("distance") = matched_distances);
      results[q] = query_matches;
    } else {
      // Empty DataFrame for queries with no matches
      results[q] = DataFrame::create(Named("target") = CharacterVector(0),
                                     Named("distance") = IntegerVector(0));
    }
    query_names[q] = query;
  }
  
  results.names() = query_names;
  Rcout << "Multithreaded search complete!" << std::endl;
  Rcout << "Found " << total_matches << " total matches across all queries." << std::endl;
  return results;
}
