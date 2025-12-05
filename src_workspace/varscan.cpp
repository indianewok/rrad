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
#include <utility>
#include <functional>
#include <unordered_map>

using namespace Rcpp;
using boost::multi_index_container;
using namespace boost::multi_index;
using namespace std;

// Function to process variable elements for a given segment
void process_segment(const ReadLayout& readLayout, const std::vector<int>& segment_orders, const std::string& direction, bool verbose) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");

  std::unordered_map<std::string, std::tuple<std::string, int, int>> assignedPositions; // Tracks assigned positions

  auto parse_static_left = [&](const ReadElement& varElem, const std::vector<int>& orders) -> std::tuple<std::string, int, int> {
    if (verbose) {
      Rcpp::Rcout << "Processing varElem: " << varElem.class_id << ", Order: " << varElem.order << "\n";
    }
    if (!orders.empty() && varElem.order > orders.front()) {
      auto prevElem = readLayout.get<order_tag>().find(varElem.order - 1);
      if (prevElem != order_index.end()) {
        if (verbose) {
          Rcpp::Rcout << "Previous Element ID: " << prevElem->class_id << ", Type: " << prevElem->type << ", Order: " << prevElem->order << "\n";
        }
        if (prevElem->type == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
          auto assignedPos = assignedPositions[prevElem->class_id];
          int startOffset = 1;
          int stopOffset = std::get<2>(assignedPos) + (varElem.expected_length ? *varElem.expected_length - 1 : 0);
          return std::make_tuple(prevElem->class_id, startOffset, stopOffset);
        } else if (prevElem->type == "static") {
          int startOffset = 1;
          int stopOffset = varElem.expected_length ? 1 + (*varElem.expected_length - 1) : 1;
          return std::make_tuple(prevElem->class_id, startOffset, stopOffset);
        }
      }
    }
    if (verbose) {
      Rcpp::Rcout << "No previous element found for " << varElem.class_id << "\n";
    }
    return std::make_tuple("Start of Segment", 0, 0);
  };

  auto parse_static_right = [&](const ReadElement& varElem, const std::vector<int>& orders) -> std::tuple<std::string, int, int> {
    for (auto it = orders.rbegin(); it != orders.rend(); ++it) {
      if (*it > varElem.order) {
        auto staticElem = readLayout.get<order_tag>().find(*it);
        if (staticElem != order_index.end() && staticElem->type == "static") {
          int startOffset = varElem.expected_length ? -1 - (*varElem.expected_length - 1) : 1; // Default to 1 if no expected_length
          int stopOffset = -1;
          return std::make_tuple(staticElem->class_id, startOffset, stopOffset);
        }
      }
    }
    return std::make_tuple("End of Segment", 0, 0);
  };

  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
      continue;
    }
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
      std::tuple<std::string, int, int> proximal_primer;
      int read_order = (direction == "forward") ? *std::max_element(segment_orders.begin(), segment_orders.end()) :
        *std::min_element(segment_orders.begin(), segment_orders.end());
      if (it_var->order < read_order) {
        proximal_primer = parse_static_left(*it_var, segment_orders);
      } else {
        proximal_primer = parse_static_right(*it_var, segment_orders);
      }
      if (verbose) {
        Rcpp::Rcout << "Variable Element: " << it_var->class_id
                    << ", Direction: " << it_var->direction
                    << ", Proximal Primer ID: " << std::get<0>(proximal_primer)
                    << ", Start Offset: " << std::get<1>(proximal_primer)
                    << ", Stop Offset: " << std::get<2>(proximal_primer) << "\n";
      }
  }
}

// VarScan function
void VarScan(ReadLayout& readLayout, bool verbose) {

  std::unordered_map<std::string, PositionInfo> varPositions;
  auto& global_class_index = readLayout.get<global_class_tag>();
  auto& order_index = readLayout.get<order_tag>();
  auto read_range = global_class_index.equal_range("read");

  // Separate forward and reverse orders
  std::vector<int> forward_orders, reverse_orders;
  for (auto it_order = order_index.begin(); it_order != order_index.end(); ++it_order) {
    if (it_order->direction == "forward") {
      forward_orders.push_back(it_order->order);
    } else {
      reverse_orders.push_back(it_order->order);
    }
  }

  // Process each read separately
  for (auto it_read = read_range.first; it_read != read_range.second; ++it_read) {
    std::string read_direction = it_read->direction;
    int read_order = it_read->order;

    // Subdivide orders based on read position
    std::vector<int> left_vector, right_vector;
    if (read_direction == "forward") {
      auto read_order_here = std::find(forward_orders.begin(), forward_orders.end(), read_order);
      left_vector = std::vector<int>(forward_orders.begin(), read_order_here);
      right_vector = std::vector<int>(read_order_here, forward_orders.end());
    } else {
      auto read_order_here = std::find(reverse_orders.begin(), reverse_orders.end(), read_order);
      left_vector = std::vector<int>(reverse_orders.begin(), read_order_here);
      right_vector = std::vector<int>(read_order_here, reverse_orders.end());
    }

    // Process the segments to the left and right of the read
    process_segment(readLayout, left_vector, read_direction, verbose);
    process_segment(readLayout, right_vector, read_direction, verbose);
  }
}

// [[Rcpp::export]]
void test_read_layout_container_v2(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose = true) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  if (verbose) {
    VarScan(container, true);
    //test_ref_layout(container);
  }
}
