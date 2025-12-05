#include <Rcpp.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace Rcpp;
using boost::multi_index_container;
using namespace boost::multi_index;

// Define tags for multi-indexing
struct id_tag {};
struct length_tag {};
struct type_tag {};
struct order_tag {};
struct direction_tag {};

// Define a structure to hold the combined information
struct ReadElement {
  std::string class_id;
  std::string seq; // Assuming l'seq' is relevant for both static and variable elements
  boost::optional<int> expected_length;
  std::string type;
  int order;
  std::string direction;
  boost::optional<std::tuple<int, int, int>> misalignment_threshold;

  ReadElement(
    std::string class_id,
    std::string seq,
    boost::optional<int> expected_length,
    std::string type,
    int order,
    std::string direction,
    boost::optional<std::tuple<int, int, int>> misalignment_threshold = boost::none
  ) : class_id(std::move(class_id)),
  seq(std::move(seq)),
  expected_length(expected_length),
  type(std::move(type)),
  order(order),
  direction(std::move(direction)),
  misalignment_threshold(misalignment_threshold) {}
};

// Define the multi-index container with appropriate indices
typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<length_tag>, member<ReadElement, boost::optional<int>, &ReadElement::expected_length>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>
  >
> ReadLayout;

ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold) {
  ReadLayout container;

  // Assume 'query_id' matches 'class_id' from 'read_layout'
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];

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
  Rcpp::StringVector class_column = read_layout["class"];
  Rcpp::StringVector seq = read_layout["seq"];
  Rcpp::IntegerVector expected_length = read_layout["expected_length"];
  Rcpp::StringVector type = read_layout["type"];
  Rcpp::IntegerVector order = read_layout["order"];
  Rcpp::StringVector direction = read_layout["direction"];

  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);
    std::string class_value = Rcpp::as<std::string>(class_column[i]); // Get the class value

    boost::optional<std::tuple<int, int, int>> misalignment_threshold_opt;
    if (thresholdsMap.count(cid)) {
        misalignment_threshold_opt = boost::optional<std::tuple<int, int, int>>(thresholdsMap[cid]);
      }else{
        misalignment_threshold_opt = boost::none;
    }

    boost::optional<int> expected_length_opt;
    if (class_value != "read") {
      expected_length_opt = expected_length[i];
    }

    container.insert(ReadElement(
        cid,
        Rcpp::as<std::string>(seq[i]),
        expected_length_opt,
        Rcpp::as<std::string>(type[i]),
        order[i],
        Rcpp::as<std::string>(direction[i]),
        misalignment_threshold_opt
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
      if(it->misalignment_threshold.is_initialized()){
      auto misalignment_data = it->misalignment_threshold.get();
      Rcpp::Rcout << "Misalignment Data for " << it->class_id << ": <"
                  << std::get<0>(misalignment_data) << ", "
                  << std::get<1>(misalignment_data) << ", "
                  << std::get<2>(misalignment_data) << ">\n";
      } else {
        Rcpp::Rcout << "No Misalignment Data for " << it->class_id << "\n";
      }
    }
  }
  // Test 4: Print expected lengths of all elements (the real friends are the NAs we missed along the way)
  {
    auto& type_index = container.get<length_tag>(); //testing that the tag functioning works
    Rcpp::Rcout << "Expected Lengths:\n";
    for(auto it = type_index.begin(); it != type_index.end(); ++it) {
      Rcpp::Rcout << it->class_id << ": ";
      if (it->expected_length.is_initialized()) {
        Rcpp::Rcout << it->expected_length.get();
      } else {
        Rcpp::Rcout << "None";
      }
      Rcpp::Rcout << "\n";
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

struct SigElement {
  std::string class_id;
  boost::optional<int> edit_distance; // Optional edit distance
  std::pair<int, int> position;
  std::string type;
  int order;
  boost::optional<int> expected_length;
  std::string direction;

  SigElement(
    std::string class_id,
    boost::optional<int> edit_distance,
    std::pair<int, int> position,
    std::string type,
    int order,
    boost::optional<int> expected_length,
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
    ordered_non_unique<tag<sig_ed_tag>, member<SigElement, boost::optional<int>, &SigElement::edit_distance>>,
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

void build_sigstring(ReadData& readData, const ReadLayout& readLayout, const std::string& signature) {
  std::stringstream ss(signature);
  std::string token;

  while (std::getline(ss, token, '|')) {
    std::stringstream tokenStream(token);
    std::string id;
    int startPos, endPos;

    std::getline(tokenStream, id, ':');

    // Retrieve additional information from ReadLayout
    auto& class_id_index = readLayout.get<id_tag>();
    auto it = class_id_index.find(id);
    if (it != class_id_index.end()) {
      boost::optional<int> editDistance;
      if (it->type != "variable") {
        int tempDistance;
        tokenStream >> tempDistance;
        editDistance = tempDistance;
      }
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> startPos;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> endPos;

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

  // Populate default elements for those not in the signature string
  for (const auto& readElement : readLayout) {
    if (readData.sigstring.find(readElement.class_id) == readData.sigstring.end()) {
      boost::optional<int> default_edit_distance = (readElement.type != "variable") ? boost::optional<int>(0) : boost::none;

      SigElement element(
          readElement.class_id,
          default_edit_distance,
          std::make_pair(0, 0), // Default positions
          readElement.type,
          readElement.order,
          readElement.expected_length,
          readElement.direction
      );
      readData.sigstring.insert(std::move(element));
    }
  }
}

void categorize_sigstring(ReadData& readData, const ReadLayout& readLayout) {
  std::string category = "undecided";
  auto& sigstring_index = readData.sigstring.get<sig_id_tag>();

  auto check_pass_fail = [&](const std::string& id) -> bool {
    auto it = sigstring_index.find(id);
    if (it != sigstring_index.end() && it->edit_distance) {
      auto threshold_it = readLayout.get<id_tag>().find(id);
      if (threshold_it != readLayout.get<id_tag>().end()) {
        int lower_bound = std::get<0>(threshold_it->misalignment_threshold.value_or(std::make_tuple(0, 0, 0)));
        return it->edit_distance.get() <= lower_bound;
      }
    }
    return false;
  };

  bool forward_primer_pass = check_pass_fail("forw_primer");
  bool rc_forward_primer_pass = check_pass_fail("rc_forw_primer");
  bool rev_primer_pass = check_pass_fail("rev_primer");
  bool rc_rev_primer_pass = check_pass_fail("rc_rev_primer");

  if (forward_primer_pass && !rc_forward_primer_pass) {
    category = "forward";
  } else if (!forward_primer_pass && rc_forward_primer_pass) {
    category = "reverse";
  } else if (forward_primer_pass && rc_forward_primer_pass && rev_primer_pass && rc_rev_primer_pass) {
    category = "true_resc_conc";
  } else if (forward_primer_pass && rc_forward_primer_pass && rev_primer_pass || rc_rev_primer_pass) {
    category = "one_resc_conc";
  }

  // Store the category in readData's str_type
  readData.string_info["str_type"] = StringInfo(0, category);
}

void print_sigstring(ReadData& readData, const ReadLayout& readLayout) {
  // Ensure 'str_type' exists in string_info
  if (readData.string_info.find("str_type") == readData.string_info.end()) {
    Rcpp::Rcerr << "Error: 'str_type' not found in string_info.\n";
    return;
  }

  const auto& category = readData.string_info["str_type"].str_type;
  const auto& sigstring = readData.sigstring;

  auto printOrderedDirection = [&](const std::string& direction) {
    auto& dir_index = sigstring.get<sig_dir_tag>();
    auto dir_range = dir_index.equal_range(direction);
    std::vector<SigElement> sorted_elements;

    for (auto it = dir_range.first; it != dir_range.second; ++it) {
      sorted_elements.push_back(*it);
    }

    std::sort(sorted_elements.begin(), sorted_elements.end(), [](const SigElement& a, const SigElement& b) { return a.order < b.order; });

    Rcpp::Rcout << category << " sequence layout:\n";
    for (const auto& element : sorted_elements) {
      std::string element_type = (element.type == "static") ? "static" : "variable";
      std::string pass_fail = element.edit_distance ? (element.edit_distance.get() <= 0 ? ":pass" : ":fail") : ":na";
      if (element.class_id == "poly_a" || element.class_id == "poly_t") {
        pass_fail = ":na";
      }
      std::string position_info = "|start:" + std::to_string(element.position.first) + "|stop:" + std::to_string(element.position.second);
      Rcpp::Rcout << "[" << element_type << ":" << element.class_id << pass_fail << position_info << "]";
    }
    Rcpp::Rcout << "\n";
  };

  printOrderedDirection(category);  // Print only the category that was determined
}

// [[Rcpp::export]]
void test_read_layout_container(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose = true) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  if (verbose) {
    test_ref_layout(container);
  }
}
// [[Rcpp::export]]
void test_sigstring_processing(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, const std::string& signature) {
  // Convert DataFrame to ReadLayout
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold);
  // Create ReadData
  ReadData readData;
  // Populate SigString
  build_sigstring(readData, readLayout, signature);
  // Categorize SigString
  categorize_sigstring(readData, readLayout);
  // Print SigString
  print_sigstring(readData, readLayout);
}
