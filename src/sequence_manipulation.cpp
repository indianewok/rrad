#include "rad.h"
#include "stringdist_api.h"

using namespace std;
using namespace Rcpp;

KSEQ_INIT(gzFile, gzread)

std::string revcomp_cpp(const std::string& sequence) {
  std::unordered_map<char, char> complement {
    {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
    {'N', 'N'}, {'n', 'n'}
  };
  std::string reversed(sequence.rbegin(), sequence.rend());
  for (char& nucleotide : reversed) {
    nucleotide = complement[nucleotide];
  }
  return reversed;
}

// [[Rcpp::export]]
Rcpp::StringVector revcomp(Rcpp::StringVector sequences) {
  int n = sequences.size();
  Rcpp::StringVector results(n);
  int j = 0; // Counter for results
  for (int i = 0; i < n; ++i) {
    if (sequences[i] == NA_STRING) {
      results[j] = NA_STRING;
      continue;
    }
    results[j] = revcomp_cpp(Rcpp::as<std::string>(sequences[i]));
    ++j;
  }
  return results;
}

std::vector<int64_t> sequence_to_bits_cpp(const std::string& sequence) {
  const int chunk_size = 32; // Number of bases per chunk
  std::vector<int64_t> results;
  for (size_t i = 0; i < sequence.size(); i += chunk_size) {
    int64_t result = 0;
    size_t end = std::min(i + chunk_size, sequence.size());
    for (size_t j = i; j < end; ++j) {
      char c = sequence[j];
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break;
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    }
    results.push_back(result);
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector sequence_to_bits(Rcpp::StringVector sequences) {
  int n = sequences.size();
  std::vector<int64_t> results;
  for(int i = 0; i < n; ++i) {
    std::string sequence = Rcpp::as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results.insert(results.end(), chunk_results.begin(), chunk_results.end());
  }
  return Rcpp::toInteger64(results);
}

// [[Rcpp::export]]
List sequence_to_bitlist(Rcpp::StringVector sequences) {
  int n = sequences.size();
  List results(n);
  for(int i = 0; i < n; ++i) {
    std::string sequence = as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results[i] = wrap(chunk_results);
  }
  return results;
}

std::string bits_to_sequence_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::string final_sequence;
  int chunk_size = 32;
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    std::string result;
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    final_sequence += result;
  }
  return final_sequence;
}

std::vector<unsigned int> bits_to_uint_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::vector<unsigned int> final_sequence;
  int chunk_size = 32; // Assuming 32 bases per int64_t value
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      unsigned int base = int64_code & 3; // Extract the last two bits
      final_sequence.push_back(base); // Add the base as an unsigned int to the sequence
      int64_code >>= 2; // Move to the next base
    }
  }
  return final_sequence;
}

// [[Rcpp::export]]
Rcpp::StringVector bits_to_sequence(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  // Convert the entire NumericVector to std::vector<int64_t>
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  Rcpp::StringVector final_sequences(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    // Process the single sequence
    std::string sequence = bits_to_sequence_cpp(bits_input, seqlength);
    final_sequences[i] = sequence;
  }
  return final_sequences;
}

std::string bits_to_hex_cpp(const std::vector<int64_t>& input, int sequence_length) {
  std::string final_sequence;
  const int chunk_size = 32;
  for (size_t i = 0; i < input.size(); ++i) {
    int64_t int64_code = input[i];
    std::string result;
    
    int k_mers_to_decode = (i < input.size() - 1) ? chunk_size : sequence_length % chunk_size;
    k_mers_to_decode = (k_mers_to_decode == 0) ? chunk_size : k_mers_to_decode;
    
    for (int j = 0; j < k_mers_to_decode; ++j) {
      int k_mer = int64_code & 0xF;
      if (k_mer < 10) {
        result.insert(result.begin(), '0' + k_mer);
      } else {
        result.insert(result.begin(), 'A' + (k_mer - 10));
      }
      int64_code >>= 4;
    }
    final_sequence += result;
  }
  return final_sequence;
}

// [[Rcpp::export]]
StringVector bits_to_hex(Rcpp::NumericVector input, int sequence_length = 16) {
  Rcpp::StringVector hex_sequences(input.size());
  std::vector<int64_t> cpp_input = Rcpp::as<std::vector<int64_t>>(input);
  int slen = sequence_length/2; 
  //for consistency and compression's sake, we convert from a rep of 2 bits->one char to 4 bits->one char
  for(size_t i = 0; i < cpp_input.size(); ++i) {
    std::vector<int64_t> single_input = {cpp_input[i]};
    std::string hex_sequence = bits_to_hex_cpp(single_input, slen);
    hex_sequences[i] = hex_sequence;
  }
  return hex_sequences;
}

std::vector<int64_t> revcomp_bits_cpp(const std::vector<int64_t>& input, int sequence_length = 16) {
  std::vector<int64_t> reversed(input.rbegin(), input.rend()); // Reverse the vector
  for (int64_t& chunk : reversed) {
    int64_t result = 0;
    for (int i = 0; i < sequence_length; i++) {
      result <<= 2; // Shift result to make room for the new bits
      int64_t base = chunk & 3; // Extract the two rightmost bits
      switch (base) {
      case 0: result |= 2; break; // A (00) becomes T (10)
      case 1: result |= 3; break; // C (01) becomes G (11)
      case 2: result |= 0; break; // T (10) becomes A (00)
      case 3: result |= 1; break; // G (11) becomes C (01)
      }
      chunk >>= 2; // Move to the next base in the chunk
    }
    // Assign the reversed and complemented chunk back
    chunk = result;
  }
  return reversed;
}

// [[Rcpp::export]]
Rcpp::NumericVector revcomp_bits(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  std::vector<int64_t> results;
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    std::vector<int64_t> subresult = revcomp_bits_cpp(bits_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::toInteger64(results);
}

// Function to calculate the number of differing bits for int64_t sequences
int hamming_distance(int64_t a, int64_t b) {
  int bit_diff_count = std::bitset<64>(a ^ b).count();
  return std::ceil(bit_diff_count / 2.0);  // Base Hamming distance
}

// [[Rcpp::export]]
int hamming_distance_rcpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  // Convert Rcpp::NumericVector (integer64) to std::vector<int64_t>
  std::vector<int64_t> a_int = Rcpp::fromInteger64(a);
  std::vector<int64_t> b_int = Rcpp::fromInteger64(b);
  
  if (a_int.size() != 1 || b_int.size() != 1) {
    Rcpp::stop("Expecting single integer64 values.");
  }
  
  // Return the Hamming distance between the first (and only) elements
  return hamming_distance(a_int[0], b_int[0]);
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp_heap(int64_t sequence, int sequence_length) {
  // Allocate mutations vector on the heap
  auto mutations = std::make_unique<std::vector<int64_t>>();
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations->push_back(mutated_sequence);
      }
    }
  }
  // Return the vector stored on the heap
  return *mutations;
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp(int64_t sequence, int sequence_length) {
  std::vector<int64_t> mutations;
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations.push_back(mutated_sequence);
      }
    }
  }
  return mutations;
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_quaternary_mutations(Rcpp::NumericVector sequences, Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  std::vector<int64_t> results;
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_quaternary_mutations_cpp(quaternary_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

// Function to generate recursive mutations
std::vector<int64_t> generate_recursive_quaternary_mutations_cpp(
    int64_t sequence, 
  int mutation_rounds, 
  int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::vector<int64_t> current_round{sequence};
  for (int round = 0; round < mutation_rounds; ++round) {
    std::vector<int64_t> next_round;
    for (const auto& seq : current_round) {
      auto mutations = generate_quaternary_mutations_cpp(seq, sequence_length);
      for (const auto& mutation : mutations) {
        if (all_mutations.insert(mutation).second) {
          next_round.push_back(mutation);
        }
      }
    }
    if (next_round.empty()) break;
    current_round = std::move(next_round);
  }
  all_mutations.erase(sequence);
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v2(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::queue<std::pair<int64_t, int>> sequences_to_mutate;
  sequences_to_mutate.push({sequence, 0});
  
  while (!sequences_to_mutate.empty()) {
    auto current = sequences_to_mutate.front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate.pop();
    
    if (current_round < mutation_rounds) {
      std::vector<int64_t> new_mutations = generate_quaternary_mutations_cpp(current_seq, sequence_length);
      for (const auto& mutation : new_mutations) {
        if (all_mutations.insert(mutation).second) {
          sequences_to_mutate.push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v3(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  // Use unique_ptr to allocate containers on the heap
  auto all_mutations = std::make_unique<std::unordered_set<int64_t>>();
  all_mutations->insert(sequence);
  
  auto sequences_to_mutate = std::make_unique<std::queue<std::pair<int64_t, int>>>();
  sequences_to_mutate->push({sequence, 0});
  
  while (!sequences_to_mutate->empty()) {
    auto current = sequences_to_mutate->front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate->pop();
    
    if (current_round < mutation_rounds) {
      // Store new mutations on the heap
      auto new_mutations = std::make_unique<std::vector<int64_t>>(generate_quaternary_mutations_cpp(current_seq, sequence_length));
      for (const auto& mutation : *new_mutations) {
        if (all_mutations->insert(mutation).second) {
          sequences_to_mutate->push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  
  // Return vector created from the heap-allocated unordered_set
  return std::vector<int64_t>(all_mutations->begin(), all_mutations->end());
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_recursive_quaternary_mutations(
    Rcpp::NumericVector sequences,
    int mutation_rounds = 2, 
  Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  std::vector<int64_t> results;
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_recursive_quaternary_mutations_cpp(quaternary_input, mutation_rounds, seqlength); 
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

// Helper function to print binary representation of int64_t
std::string int64_to_binary(int64_t value, int length) {
  std::string binary = std::bitset<64>(value).to_string();
  return binary.substr(64 - 2*length);
}

// Helper function to convert int64_t to quaternary string
std::string int64_to_quaternary(int64_t value, int length) {
  std::string result;
  for (int i = length - 1; i >= 0; --i) {
    int base = (value >> (2 * i)) & 0b11;
    switch(base) {
    case 0: result += 'A'; break;
    case 1: result += 'C'; break;
    case 2: result += 'T'; break;
    case 3: result += 'G'; break;
    }
  }
  return result;
}

std::vector<int64_t> generate_shifted_barcodes(int64_t original, int sequence_length, int max_shift = 2) {
  std::unordered_set<int64_t> shifted_barcodes_set;  // Use a set to avoid duplicates
  int64_t mask = (1LL << (2 * sequence_length)) - 1;  // Mask for the original sequence length
  
  // Generate right shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original >> (2 * shift)) & (mask >> (2 * shift))) | (padding << (2 * (sequence_length - shift)));
      shifted_barcodes_set.insert(shifted);  // Insert into set
    }
  }
  
  // Generate left shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original << (2 * shift)) & mask) | padding;
      shifted_barcodes_set.insert(shifted);  // Insert into set
    }
  }
  
  // Convert set to vector and return
  return std::vector<int64_t>(shifted_barcodes_set.begin(), shifted_barcodes_set.end());
}

// [[Rcpp::export]]
Rcpp::List shift_barcodes(
    SEXP sequences,
    int sequence_length,
    int max_shift = 2,
    bool verbose = false
  ) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  Rcpp::List result(cpp_input_list.size());
  for (size_t seq_index = 0; seq_index < cpp_input_list.size(); ++seq_index) {
    int64_t original = cpp_input_list[seq_index];
    if (verbose) {
      std::cout << "Processing sequence " << seq_index + 1 << ":\n";
      std::cout << "Original sequence: " << int64_to_quaternary(original, sequence_length) << "\n";
      std::cout << "Binary representation: " << int64_to_binary(original, sequence_length) << "\n";
    }
    std::vector<int64_t> shifted = generate_shifted_barcodes(original, sequence_length, max_shift);
    Rcpp::List seq_result(shifted.size());  // +1 to include the original sequence
    //seq_result[0] = Rcpp::wrap(original);  // Include the original sequence
    for (size_t i = 0; i < shifted.size(); ++i) {
      //seq_result[i + 1] = Rcpp::wrap(shifted[i]);
      seq_result[i] = Rcpp::wrap(shifted[i]);
      if (verbose) {
        std::cout << "Shifted " << (i < shifted.size()/2 ? "right" : "left") << ": " 
                    << int64_to_quaternary(shifted[i], sequence_length) << "\n";
        std::cout << "Binary: " << int64_to_binary(shifted[i], sequence_length) << "\n";
      }
    }
    result[seq_index] = seq_result;
  }
  return result;
}

// Helper function to convert int64_t to string representation
std::string int64_to_string(uint64_t value, int length) {
  std::string result;
  for (int i = length - 1; i >= 0; --i) {
    int base = (value >> (2 * i)) & 0b11;
    switch(base) {
    case 0: result += 'A'; break;
    case 1: result += 'C'; break;
    case 2: result += 'T'; break;
    case 3: result += 'G'; break;
    }
  }
  return result;
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::DataFrame barcode_correction_v3(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    int method = 2,
    bool verbose = false) {
  //auto start_time = std::chrono::high_resolution_clock::now();
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_map<int64_t, int> original_counts, corrected_counts;  // Track counts for original and corrected barcodes
  
  // Step 2: Deduplicate correct barcodes and consolidate counts
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    original_counts[correct_barcodes[i]] += Rcpp::as<std::vector<int>>(true_counts)[i];  // Consolidate counts
    corrected_counts[correct_barcodes[i]] = 0;  // Initialize corrected counts
  }
  
  // Deduplicate incorrect barcodes
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 3: Initialize incorrect barcode map (resolved, [putative_barcodes, mutation_type], dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 4: Generate mutations and shifts, including checking for both mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    std::unordered_set<int64_t> mutated_set;
    std::unordered_set<int64_t> shifted_set;
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(
      correct_barcode,
      mutation_rounds,
      sequence_length);
    
    mutated_set.insert(mutated_barcodes.begin(), mutated_barcodes.end());
    
    // Generate shifted barcodes
    if(max_shift > 0){
      std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
      shifted_set.insert(shifted_barcodes.begin(), shifted_barcodes.end());
    }
    // Combine both sets to track "mutated", "shifted", and "mutated_and_shifted"
    std::unordered_set<int64_t> unique_barcodes(mutated_set.begin(), mutated_set.end());
    unique_barcodes.insert(shifted_set.begin(), shifted_set.end());
    
    for (const auto& barcode : unique_barcodes) {
      if (incorrect_set.find(barcode) != incorrect_set.end()) {
        int mutation_type = 0;  // Default: mutated
        if (mutated_set.find(barcode) != mutated_set.end() && shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 2;  // mutated and shifted
        } else if (shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 1;  // shifted
        }
#pragma omp critical
        std::get<1>(incorrect_map[barcode]).emplace_back(correct_barcode, mutation_type);
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 5: Populate DL distances
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  std::cout << "Method is " << method << "\n";
  // 0 is OSA, 1 is LV, 2 is DL, 3 is hamming
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    if (putative_barcodes.size() > 1) {
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)  // Single thread
      );
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 6: Resolve collisions
  //size_t total_incorrect = incorrect_map.size();
  log_frequency = 50000;  // Log progress every 50,000 sequences
  std::atomic<size_t> processed_incorrect(0);
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    processed_incorrect++;
    
    if (putative_barcodes.empty()) {
      resolved = false;
      continue;
    }
    
    if (putative_barcodes.size() == 1) {
      resolved = true;
      continue;
    }
    
    // Find the minimum DL distance and handle uniqueness
    std::unordered_map<int64_t, double> barcode_to_min_dl_dist;
    for (size_t i = 0; i < putative_barcodes.size(); ++i) {
      int64_t barcode = putative_barcodes[i].first;
      double dist = dl_distances[i];
      barcode_to_min_dl_dist[barcode] = std::min(barcode_to_min_dl_dist[barcode], dist);
    }
    
    // Find the barcode with the smallest DL distance
    auto min_it = std::min_element(
      barcode_to_min_dl_dist.begin(),
      barcode_to_min_dl_dist.end(),
      [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
        return a.second < b.second;
      }
    );
    
    double min_dl_dist = min_it->second;
    std::vector<int64_t> min_barcodes;
    
    for (const auto& [barcode, dist] : barcode_to_min_dl_dist) {
      if (dist == min_dl_dist) {
        min_barcodes.push_back(barcode);
      }
    }
    
    if (min_barcodes.size() == 1) {
      for (const auto& pair : putative_barcodes) {
        if (pair.first == min_barcodes[0]) {
          putative_barcodes = {pair};
          break;
        }
      }
      resolved = true;
    } else {
      resolved = false;
    }
  }
  
  std::cout << "Collision resolution completed.\n";
  
  // Step 7: Prepare outputs
  std::vector<std::string> incorrect_barcode_vec(incorrect_map.size());
  std::vector<std::string> corrected_barcode_vec(incorrect_map.size());
  std::vector<bool> resolved_status_vec(incorrect_map.size());
  std::vector<std::string> mutation_type_vec(incorrect_map.size());
  
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    incorrect_barcode_vec[idx] = int64_to_string(incorrect_barcode, sequence_length);
    if (resolved) {
      corrected_barcode_vec[idx] = int64_to_string(putative_barcodes[0].first, sequence_length);
      mutation_type_vec[idx] = (putative_barcodes[0].second == 0) ? "mutated" :
        (putative_barcodes[0].second == 1) ? "shifted" : "mutated_and_shifted";
    } else {
      if (putative_barcodes.empty()) {
        mutation_type_vec[idx] = "unresolved_empty";
        corrected_barcode_vec[idx] = "";
      } else {
        mutation_type_vec[idx] = "unresolved_multiple";
        std::string unresolved_barcodes;
        for (size_t j = 0; j < putative_barcodes.size(); ++j) {
          unresolved_barcodes += int64_to_string(putative_barcodes[j].first, sequence_length);
          if (putative_barcodes.size() > 1 && j != putative_barcodes.size() - 1) {
            unresolved_barcodes += "|";
          }
        }
        corrected_barcode_vec[idx] = unresolved_barcodes;
      }
    }
    resolved_status_vec[idx] = resolved;
  }
  
  std::cout << "Resolving the outputs, now generating dataframe.\n";
  
  // Step 8: Return as DataFrame
  return Rcpp::DataFrame::create(
    Rcpp::Named("incorrect_barcodes") = Rcpp::wrap(incorrect_barcode_vec),
    Rcpp::Named("corrected_barcodes") = Rcpp::wrap(corrected_barcode_vec),
    Rcpp::Named("resolved_status") = Rcpp::wrap(resolved_status_vec),
    Rcpp::Named("mutation_type") = Rcpp::wrap(mutation_type_vec)
  );
}

// colon delim for makedb.py in airr is messy--consider changing the delim?