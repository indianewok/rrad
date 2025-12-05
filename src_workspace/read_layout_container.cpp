#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/variant.hpp>
#include <string>
#include <tuple>
#include <Rcpp.h>

using namespace boost::multi_index;
using namespace Rcpp;

// Define a structure to hold the combined information
struct ReadElement {
  std::string class_id;
  std::string seq; // Assuming l'seq' is relevant for both static and variable elements
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
 struct direction_tag{};

// Define the multi-index container with appropriate indices
typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<length_tag>, member<ReadElement, int, &ReadElement::expected_length>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>
      // The index typing moves top-down from 0 to whatever
  >
> ReadLayout;

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

// [[Rcpp::export]]
void prep_read_layout(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_thresholds, bool verbose = false) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_thresholds);
  if (verbose) {
    test_ref_layout(container);
  }
}
