#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/variant.hpp>
#include <string>
#include <vector>
#include <tuple>
#include <Rcpp.h>
#include <cstdint>
#include <sstream>
#include <unordered_map>

using namespace boost::multi_index;
using namespace Rcpp;

// Declaring our sequence_to_bits_cpp function

std::vector<int64_t> sequence_to_bits_cpp(const std::string& sequence);

// Define a structure to hold combined information
struct ReadElement {
  std::string class_id;
  std::string seq;
  int expected_length;
  std::string type;
  int order;
  std::string direction;
  std::tuple<int, int, int> misalignment_threshold;

  ReadElement(
    std::string class_id,
    std::string seq,
    int expected_length,
    std::string type,
    int order,
    std::string direction,
    std::tuple<int, int, int> misalignment_threshold = std::make_tuple(0, 0, 0)
  ) : class_id(std::move(class_id)),
  seq(std::move(seq)),
  expected_length(expected_length),
  type(std::move(type)),
  order(order),
  direction(std::move(direction)),
  misalignment_threshold(misalignment_threshold) {}
};

struct id_tag {};
struct length_tag {};
struct type_tag {};
struct order_tag {};
struct direction_tag {};

// Define the multi-index container with appropriate indices for ReadLayout
typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<length_tag>, member<ReadElement, int, &ReadElement::expected_length>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>
  >
> ReadLayout;

// Define the SigElement structure
struct SigElement {
  std::string class_id;
  int edit_distance;
  std::pair<int, int> position;
  std::string type;
  int order;
  int expected_length;
  std::string direction;

  SigElement(
    std::string class_id,
    int edit_distance,
    std::pair<int, int> position,
    std::string type,
    int order,
    int expected_length,
    std::string direction
  ) : class_id(std::move(class_id)),
  edit_distance(edit_distance),
  position(position),
  type(std::move(type)),
  order(order),
  expected_length(expected_length),
  direction(std::move(direction)) {}
};

struct sig_id_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};

// Define the multi-index container for SigString
typedef multi_index_container<
  SigElement,
  indexed_by<
    ordered_unique<tag<sig_id_tag>, member<SigElement, std::string, &SigElement::class_id>>,
    ordered_non_unique<tag<sig_ed_tag>, member<SigElement, int, &SigElement::edit_distance>>,
    ordered_non_unique<tag<sig_order_tag>, member<SigElement, int, &SigElement::order>>,
    ordered_non_unique<tag<sig_dir_tag>, member<SigElement, std::string, &SigElement::direction>>
  >
> SigString;

struct StringInfo {
  int str_length;
  std::string str_type;

  StringInfo(int length = 0, std::string type = "undefined")
    : str_length(length), str_type(std::move(type)) {}
};

struct ReadData {
  std::unordered_map<std::string, StringInfo> string_info;
  SigString sigstring;

  void addStringInfo(const std::string& key, int length, const std::string& type) {
    string_info[key] = StringInfo(length, type);
  }
};

ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_thresholds) {
  ReadLayout container;

  // Assume 'query_id' matches 'class_id' from 'read_layout'
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_thresholds["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_thresholds["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_thresholds["misal_sd"];

  for (int i = 0; i < query_id.size(); ++i) {
    double threshold = misal_threshold[i];
    double sd = std::ceil(misal_sd[i]);

    int lower_bound = static_cast<int>(std::floor(threshold - sd));
    int upper_bound = static_cast<int>(std::ceil(threshold + sd));

    thresholdsMap[Rcpp::as<std::string>(query_id[i])] = std::make_tuple(
      lower_bound,
      static_cast<int>(threshold),
      upper_bound
    );
  }

  Rcpp::StringVector class_id = read_layout["class_id"];
  Rcpp::StringVector seq = read_layout["seq"];
  Rcpp::IntegerVector expected_length = read_layout["expected_length"];
  Rcpp::StringVector type = read_layout["type"];
  Rcpp::IntegerVector order = read_layout["order"];
  Rcpp::StringVector direction = read_layout["direction"];

  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);

    container.insert(ReadElement(
        cid,
        Rcpp::as<std::string>(seq[i]),
        expected_length[i],
        Rcpp::as<std::string>(type[i]),
        order[i],
        Rcpp::as<std::string>(direction[i]),
        thresholdsMap.count(cid) ? thresholdsMap[cid] : std::make_tuple(0, 0, 0)
    ));
  }
  return container;
}

void test_ref_layout(const ReadLayout& container) {
  // Test 1: Find an element by class_id 'read'
  {
    auto& class_id_index = container.get<id_tag>();
    auto range = class_id_index.equal_range("read");
    bool found_by_class_id = (range.first != class_id_index.end());
    Rcpp::Rcout << "Found by class_id (read): " << found_by_class_id << "\n";
  }
  // Test 2: Count elements of a specific type 'static'
  {
    auto& type_index = container.get<2>();
    std::string test_type = "static";
    int count_by_type = type_index.count(test_type);
    Rcpp::Rcout << "Count by type (static): " << count_by_type << "\n";
  }
  // Test 3: Retrieve misalignment data for type::static
  {
    auto& type_index = container.get<type_tag>();
    auto range = type_index.equal_range("static");
    for(auto it = range.first; it != range.second; ++it) {
      auto misalignment_data = it->misalignment_threshold;
      Rcpp::Rcout << "Misalignment Data for " << it->class_id << ": ("
                  << std::get<0>(misalignment_data) << ", "
                  << std::get<1>(misalignment_data) << ", "
                  << std::get<2>(misalignment_data) << ")\n";
    }
  }
  // Test 4: Print expected lengths of all elements (the real friends are the NAs we missed along the way)
  {
    auto& type_index = container.get<length_tag>(); //testing that the tag functioning works
    Rcpp::Rcout << "Expected Lengths:\n";
    for(auto it = type_index.begin(); it != type_index.end(); ++it) {
      Rcpp::Rcout << it->class_id << ": " << it->expected_length << "\n";
    }
  }
  // Test 5: See it all come together part one:(hopefully?)
  {
    auto printLayout = [&](const std::string& direction) {
      auto& dir_index = container.get<4>();
      auto range = dir_index.equal_range(direction);

      Rcpp::Rcout << direction << " sequence layout:\n";
      for(auto it = range.first; it != range.second; ++it) {
        if (it->type == "static") {
          Rcpp::Rcout << "[static:" << it->class_id << "]";
        } else {
          Rcpp::Rcout << "[variable:" << it->class_id << "]";
        }
      }
      Rcpp::Rcout << "\n";
    };

    printLayout("forward");
    printLayout("reverse");
  }
  // Test 5.5: See where the class ID order is being made--it's alphabetical, which is interesting
  {
    auto& id_index = container.get<id_tag>();
    Rcpp::Rcout << "Elements in Order by class_id:\n";
    for (auto it = id_index.begin(); it != id_index.end(); ++it) {
      Rcpp::Rcout << it->class_id << " ";
    }
    Rcpp::Rcout << "\n";
  }
  // Test 6: Given an element, return the class_id before and after it by order
  {
    auto& class_id_index = container.get<id_tag>();
    auto it = class_id_index.find("barcode");
    if (it == class_id_index.end()) {
      Rcpp::Rcout << "Class ID '" << "barcode" << "' not found.\n";
      return;
    }
    // Get the order of the found element
    int found_order = it->order;
    // Access the order index
    auto& order_index = container.get<order_tag>();
    // Find the element in the order index
    auto order_it = order_index.find(found_order);
    Rcpp::Rcout << "Found Element: " << order_it->class_id << " (Order: " << found_order << ")\n";
    // Print previous element if it exists
    if (order_it != order_index.begin()) {
      auto prev_it = std::prev(order_it);
      Rcpp::Rcout << "Previous Element: " << prev_it->class_id << " (Order: " << prev_it->order << ")\n";
    } else {
      Rcpp::Rcout << "No previous element.\n";
    }
    // Print next element if it exists
    auto next_it = std::next(order_it);
    if (next_it != order_index.end()) {
      Rcpp::Rcout << "Next Element: " << next_it->class_id << " (Order: " << next_it->order << ")\n";
    } else {
      Rcpp::Rcout << "No next element.\n";
    }
  }
}

// Function to fill the SigString in ReadData using ReadLayout as a reference
void fillSigString(ReadData& readData, const ReadLayout& readLayout, const std::string& signature) {
  std::stringstream ss(signature);
  std::string token;

  while (std::getline(ss, token, '|')) {
    std::stringstream tokenStream(token);
    std::string id;
    int editDistance, startPos, endPos;

    std::getline(tokenStream, id, ':');
    tokenStream >> editDistance;
    tokenStream.ignore(1); // Ignore the colon
    tokenStream >> startPos;
    tokenStream.ignore(1); // Ignore the colon
    tokenStream >> endPos;

    // Retrieve additional information from ReadLayout
    auto& class_id_index = readLayout.get<id_tag>();
    auto it = class_id_index.find(id);
    if (it != class_id_index.end()) {
      // Create SigElement with complete information
      SigElement element(
          id,
          editDistance,
          std::make_pair(startPos, endPos),
          it->type,
          it->order,
          it->expected_length,
          it->direction
      );
      readData.sigstring.insert(std::move(element));
    } else {
      // Handle case where id is not found in ReadLayout
      Rcpp::Rcerr << "Warning: id '" << id << "' not found in ReadLayout.\n";
    }
  }
  for (const auto& readElement : readLayout) {
    if (readData.sigstring.find(readElement.class_id) == readData.sigstring.end()) {
      SigElement element(
          readElement.class_id,
          0, // default edit distance
          std::make_pair(0, 0), // default positions
          readElement.type,
          readElement.order,
          readElement.expected_length,
          readElement.direction
      );
      readData.sigstring.insert(std::move(element));
    }
  }
}

// Display function for the ordered sigstring and multiple timepoints
void displayOrderedSigString(const SigString& sigstring) {
  auto printOrderedDirection = [&](const std::string& direction) {
    // Filter by direction and sort by order
    std::vector<SigElement> sorted_elements;
    auto& dir_index = sigstring.get<sig_dir_tag>();
    auto dir_range = dir_index.equal_range(direction);
    for (auto it = dir_range.first; it != dir_range.second; ++it) {
      sorted_elements.push_back(*it);
    }
    std::sort(sorted_elements.begin(), sorted_elements.end(),
      [](const SigElement& a, const SigElement& b) { return a.order < b.order; });

    // Print elements
    Rcpp::Rcout << direction << " sequence layout:\n";
    for (const auto& element : sorted_elements) {
      std::string element_type = (element.type == "static") ? "static" : "variable";
      std::string pass_fail = (element.edit_distance <= 0) ? ":pass" : ":fail";
      std::string position_info = "|start:" + std::to_string(element.position.first) +
        "|stop:" + std::to_string(element.position.second);
      Rcpp::Rcout << "[" << element_type << ":" << element.class_id << pass_fail
                  << position_info << "]";
    }
    Rcpp::Rcout << "\n";
  };

  printOrderedDirection("forward");
  printOrderedDirection("reverse");
}

// Scanning through for the variable sets and annotation positions
void VarScan(SigString& sigstring) {
  auto updatePositions = [&](std::vector<SigElement>& elements, const std::string& direction) {
    for (auto& element : elements) {
      if (element.type == "variable" && (element.class_id == "barcode" || element.class_id == "umi")) {
        // Initialize reference position
        int reference_pos = (direction == "forward") ? INT_MAX : INT_MIN;

        // Find the reference static element
        for (auto& static_element : elements) {
          if (static_element.type == "static" && static_element.class_id != "poly_a" && static_element.class_id != "poly_t") {
            if (direction == "forward" && static_element.order < element.order) {
              reference_pos = std::min(reference_pos, static_element.position.second);
            } else if (direction == "reverse" && static_element.order > element.order) {
              reference_pos = std::max(reference_pos, static_element.position.first);
            }
          }
        }

        // Update the position of the variable element
        if ((direction == "forward" && reference_pos != INT_MAX) ||
          (direction == "reverse" && reference_pos != INT_MIN)) {
          if (direction == "forward") {
            element.position.first = reference_pos + 1;
            element.position.second = element.position.first + element.expected_length - 1;
          } else {  // Reverse
            element.position.second = reference_pos - 1;
            element.position.first = element.position.second - element.expected_length + 1;
          }
        }
      }
    }
  };

  auto sortAndOperate = [&](const std::string& direction) {
    // ... Existing code to sort and update positions ...

    updatePositions(sorted_elements, direction); // Corrected call
  };

  sortAndOperate("forward");
  sortAndOperate("reverse");
}

// [[Rcpp::export]]
void prep_read_layout(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_thresholds, bool verbose = false) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_thresholds);
  if (verbose) {
    test_ref_layout(container);
  }
}

// [[Rcpp::export]]
void test_sigstring(const DataFrame& read_layout, const DataFrame& misalignment_thresholds, const std::string& sigstring, bool verbose = false) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_thresholds);
  ReadData readData;
  fillSigString(readData, container, sigstring);
  //setVariablePositions(readData.sigstring); // Setting variable positions
  VarScan(readData.sigstring);
  if(verbose){
    displayOrderedSigString(readData.sigstring);
  }
}
