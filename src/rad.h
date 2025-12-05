#ifndef RAD
#define RAD
//here are the shared libraries for all of the rad tools. needs massive cleanup
//but it's mainly because I've never used a header file before so I don't even know which ones of these i can drop yet.

// [[Rcpp::plugins(openmp)]]  // Enable OpenMP in Rcpp
#include <Rcpp.h>
#include "edlib.h"
#include "RcppInt64.h"
#include "kseq.h"
#include <chrono>
#include <iomanip>
#include <queue>
#include <sys/resource.h>
#include <errno.h>
#include <boost/algorithm/string.hpp>
#include <mutex>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/bimap.hpp>
#include <boost/optional.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include <filesystem>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <functional>
#include <omp.h>
#include <set>
#include <unordered_set>
#include <cstdint>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <regex>
#include <numeric>
#include <cmath>
#include <iostream>
#include <zlib.h>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>

//structs for alignment tools
struct UniqueAlignment {
  int edit_distance;
  int start_position;
  int end_position;
  
  // Define how to compare two UniqueAlignments
  bool operator<(const UniqueAlignment& other) const {
    return std::tie(edit_distance, start_position, end_position) <
      std::tie(other.edit_distance, other.start_position, other.end_position);
  }
};

struct AlignmentInfo {
  int best_edit_distance;
  int best_start_pos;
  int best_stop_pos;
  int second_edit_distance;
  int second_start_pos;
  int second_stop_pos;
};

// Define tags for multi-indexing

struct id_tag {};
struct length_tag {};
struct type_tag {};
struct order_tag {};
struct direction_tag {};
struct global_class_tag {};

struct PositionInfo {
  std::tuple<std::string, std::string, std::string, std::string> position_data; // Tuple of primary and secondary start positions
};

struct ReadElement {
  std::string class_id;
  std::string seq; // Sequence information, there's nothing in here tbh
  boost::optional<int> expected_length; // Optional expected length
  std::string type; // Type of the element (e.g., "static", "variable")
  int order; // Order in the layout
  std::string direction; // Direction (e.g., "forward", "reverse")
  std::string global_class;
  boost::optional<std::tuple<int, int, int>> misalignment_threshold; // Optional misalignment threshold
  boost::optional<PositionInfo> position_data; // Optional position information for variable elements
  
  ReadElement(
    const std::string& class_id,
    const std::string& seq,
    const boost::optional<int>& expected_length,
    const std::string& type,
    int order,
    const std::string& direction,
    const std::string& global_class,
    const boost::optional<std::tuple<int, int, int>>& misalignment_threshold = boost::none,
    const boost::optional<PositionInfo>& position_data = boost::none
  ) : class_id(class_id),
  seq(seq),
  expected_length(expected_length),
  type(type),
  order(order),
  direction(direction),
  global_class(global_class),
  misalignment_threshold(misalignment_threshold),
  position_data(position_data) {}
};

typedef boost::multi_index::multi_index_container<
  ReadElement,
  boost::multi_index::indexed_by<
    boost::multi_index::ordered_unique<boost::multi_index::tag<id_tag>, boost::multi_index::member<ReadElement, std::string, &ReadElement::class_id>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<length_tag>, boost::multi_index::member<ReadElement, boost::optional<int>, &ReadElement::expected_length>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<type_tag>, boost::multi_index::member<ReadElement, std::string, &ReadElement::type>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<order_tag>, boost::multi_index::member<ReadElement, int, &ReadElement::order>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<direction_tag>, boost::multi_index::member<ReadElement, std::string, &ReadElement::direction>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<global_class_tag>, boost::multi_index::member<ReadElement, std::string, &ReadElement::global_class>>
  >
> ReadLayout;

struct SigElement {
  std::string class_id;
  std::string global_class;
  boost::optional<int> edit_distance; // Optional edit distance
  std::pair<int, int> position;
  std::string type;
  int order;
  std::string direction;
  boost::optional<bool> element_pass;
  boost::optional<std::string> seq;
  
  SigElement(
    std::string class_id,
    std::string global_class,
    boost::optional<int> edit_distance,
    std::pair<int, int> position,
    std::string type,
    int order,
    std::string direction,
    boost::optional<bool> element_pass = boost::none,
    boost::optional<std::string> seq = boost::none
  ) : class_id(std::move(class_id)),
  global_class(std::move(global_class)),
  edit_distance(edit_distance),
  position(position),
  type(std::move(type)),
  order(order),
  direction(std::move(direction)),
  element_pass(element_pass),
  seq(seq) {}
};

//manual memory release function for boost

struct sig_id_tag {};
struct sig_global_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};
struct sig_pass_tag {}; //added to be able to index by element_pass

// Define the multi-index container for SigString
typedef boost::multi_index::multi_index_container<
  SigElement,
  boost::multi_index::indexed_by<
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_id_tag>, boost::multi_index::member<SigElement, std::string, &SigElement::class_id>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_global_tag>, boost::multi_index::member<SigElement, std::string, &SigElement::global_class>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_ed_tag>, boost::multi_index::member<SigElement, boost::optional<int>, &SigElement::edit_distance>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_order_tag>, boost::multi_index::member<SigElement, int, &SigElement::order>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_dir_tag>, boost::multi_index::member<SigElement, std::string, &SigElement::direction>>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<sig_pass_tag>, boost::multi_index::member<SigElement, boost::optional<bool>, &SigElement::element_pass>>
  >
> SigString;

struct StringInfo {
  int str_length;
  std::string str_type;
  std::string str_id;
  boost::optional<std::pair<std::string, std::string>> split_types; // Stores the type of both split strings
  
  StringInfo(int length = 0, std::string id = "NA", std::string type = "undefined",
    boost::optional<std::pair<std::string, std::string>> splitTypes = boost::none)
    : str_length(length), str_type(std::move(type)), str_id(std::move(id)), split_types(splitTypes) {}
};

struct ReadData {
  std::unordered_map<std::string, StringInfo> string_info;
  SigString sigstring;
  boost::optional<SigString> secondary_sigstring; 
  
  void addStringInfo(const std::string& key, int length, const std::string& type, const std::string& id,
    boost::optional<std::pair<std::string, std::string>> splitTypes = boost::none) {
    string_info[key] = StringInfo(length, type, id, splitTypes);
  }
};

struct stat_elem {
  int forw_pass = 0;
  int rev_pass = 0;
  int unique_forw_pass = 0; 
  int unique_rev_pass = 0;
  int term_forw_elements = 0;
  int term_rev_elements = 0;
};

struct PositionCalcFunc {
  std::function<int(const ReadData&)> primaryStartFunc;
  std::function<int(const ReadData&)> secondaryStartFunc;
  std::function<int(const ReadData&)> primaryStopFunc;
  std::function<int(const ReadData&)> secondaryStopFunc;
};

using PositionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;

extern ReadLayout container;

struct dictionary{
  unsigned int key;
  unsigned int value;
  struct dictionary* next;
};

std::string findPolyTail(const std::string& sequence, char poly_base, int window_size, int min_count);

std::string findPolyTails(const std::string& sequence, char poly_base, int window_size, int min_count);

ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose);

void displayOrderedSigString(const ReadData& readData);

void pos_scan(ReadData& readData, const PositionFuncMap& positionFuncMap, const ReadLayout& readLayout, bool verbose);

void VarScan(ReadLayout& readLayout, bool verbose);

PositionFuncMap createPositionFunctionMap(const ReadLayout& readLayout, bool verbose);

void fillSigString(ReadData& readData, const ReadLayout& readLayout, 
  const std::string& signature, const PositionFuncMap& positionFuncMap, bool verbose);

Rcpp::CharacterVector process_sigstrings(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads);

std::vector<std::string> process_sigstrings_cpp(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads);

std::vector<std::string> generateSigString(ReadData& readData);

std::vector<std::string> processSigString(const std::string& sigstring);

std::string revcomp_cpp(const std::string& sequence);

std::vector<int64_t> sequence_to_bits_cpp(const std::string& sequence);

std::string bits_to_sequence_cpp(const std::vector<int64_t>& input, int sequence_length);

std::vector<int64_t> revcomp_bits_cpp(const std::vector<int64_t>& input, int sequence_length);

int hamming_bits_cpp(int64_t a, int64_t b, int sequence_length);

std::vector<std::string> get_fastq_files(const std::string& path);

std::vector<std::string> split(const std::string &s, char delimiter);

#endif