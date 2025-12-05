#include <Rcpp.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

// Define the ReadElement structure
struct ReadElement {
  int order;
  std::string class_id;
  std::string type;
  std::string direction;
  std::string global_class;
  std::optional<int> expected_length;
  std::string position_data;
};

// Define the index tags
struct order_tag {};
struct type_tag {};
struct direction_tag {};
struct global_class_tag {};
struct id_tag {};

// Define the ReadLayout container
using ReadLayout = boost::multi_index_container<
  ReadElement,
  boost::multi_index::indexed_by<
    boost::multi_index::ordered_unique<boost::multi_index::tag<order_tag>, BOOST_MULTI_INDEX_MEMBER(ReadElement, int, order)>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<type_tag>, BOOST_MULTI_INDEX_MEMBER(ReadElement, std::string, type)>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<direction_tag>, BOOST_MULTI_INDEX_MEMBER(ReadElement, std::string, direction)>,
    boost::multi_index::ordered_non_unique<boost::multi_index::tag<global_class_tag>, BOOST_MULTI_INDEX_MEMBER(ReadElement, std::string, global_class)>,
    boost::multi_index::ordered_unique<boost::multi_index::tag<id_tag>, BOOST_MULTI_INDEX_MEMBER(ReadElement, std::string, class_id)>
  >
>;

// Define the PositionInfo structure
struct PositionInfo {
  std::tuple<std::string, std::string, std::string, std::string> position_data;
};

// Initialization Functions

void initialize_indices(const ReadLayout& readLayout,
  const std::string& direction,
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto& direction_index = readLayout.get<direction_tag>();
  auto& global_class_index = readLayout.get<global_class_tag>();
}

void separate_orders(const ReadLayout& readLayout, 
  std::vector<int>& forward_orders, 
  std::vector<int>& reverse_orders) {
  auto& order_index = readLayout.get<order_tag>();
  for (auto it = order_index.begin(); it != order_index.end(); ++it) {
    if (it->direction == "forward") {
      forward_orders.push_back(it->order);
    } else {
      reverse_orders.push_back(it->order);
    }
  }
}

// Helper Functions

auto findNearestStatic = [&](auto begin, auto end, const std::string& direction) {
  return std::find_if(begin, end, [&](const ReadElement& e) {
    return e.type == "static" && e.direction == direction;
  });
};

void calculate_primary_positions(const ReadElement& elem, 
  const ReadElement& neighbor_elem, 
  bool parse_left, 
  std::string& primaryStart, 
  std::string& primaryStop, 
  const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, 
  bool static_mode) {
  primaryStart = neighbor_elem.class_id + (parse_left ? "|stop+1" : "|start-" + std::to_string(*elem.expected_length));
  primaryStop = neighbor_elem.class_id + (parse_left ? "|stop+" + std::to_string(*elem.expected_length) : "|start-1");
  
  if (static_mode && neighbor_elem.type == "variable" && assignedPositions.find(neighbor_elem.class_id) != assignedPositions.end()) {
    auto& assigned = assignedPositions.at(neighbor_elem.class_id);
    primaryStart += "|chained_start";
    primaryStop += "|chained_stop";
  }
}

void calculate_secondary_positions(const ReadElement& elem, 
  const ReadElement& neighbor_elem, 
  bool parse_left, 
  std::string& secondaryStart, 
  std::string& secondaryStop, 
  const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, 
  bool static_mode) {
  if (static_mode && neighbor_elem.type == "variable" && assignedPositions.find(neighbor_elem.class_id) != assignedPositions.end()) {
    auto& assigned = assignedPositions.at(neighbor_elem.class_id);
    auto next_elem = readLayout.get<order_tag>().find(elem.order + (parse_left ? 1 : -1));
    if (next_elem != readLayout.get<order_tag>().end() && next_elem->type == "static") {
      secondaryStart = next_elem->class_id + "|start-" + std::to_string(*next_elem->expected_length);
      secondaryStop = next_elem->class_id + "|start-1";
    } else {
      secondaryStart = std::get<0>(assigned) + "|chained_start";
      secondaryStop = std::get<2>(assigned) + "|chained_stop";
    }
  } else if (neighbor_elem.type == "static") {
    secondaryStart = neighbor_elem.class_id + "|start+1|terminal_linked";
    secondaryStop = neighbor_elem.class_id + "|stop+" + std::to_string(*elem.expected_length + 1) + "|terminal_linked";
  }
}

// Processing Functions

void process_static_elements(const ReadLayout& readLayout, 
  const std::vector<int>& segment_orders, 
  bool parse_left, 
  bool static_mode, 
  const std::string& direction, 
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  auto& order_index = readLayout.get<order_tag>();
  int read_order = parse_left ? *std::max_element(segment_orders.begin(), segment_orders.end()) : *std::min_element(segment_orders.begin(), segment_orders.end());
  
  for (const auto& elem : readLayout) {
    if (static_mode) {
      if (elem.type != "variable" || !elem.expected_length || elem.direction != direction || std::find(segment_orders.begin(), segment_orders.end(), elem.order) == segment_orders.end()) {
        continue;
      }
    } else {
      if (elem.global_class != "read") {
        continue;
      }
    }
    
    if ((parse_left && elem.order <= read_order) || (!parse_left && elem.order >= read_order)) {
      int neighbor_order = parse_left ? elem.order - 1 : elem.order + 1;
      auto neighbor_elem = order_index.find(neighbor_order);
      
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (neighbor_elem == order_index.end()) {
        primaryStart = parse_left ? "start_of_read" : "end_of_read";
        primaryStop = primaryStart + "+" + std::to_string(*elem.expected_length);
        secondaryStart = primaryStart;
        secondaryStop = primaryStop;
      } else {
        calculate_primary_positions(elem, *neighbor_elem, parse_left, primaryStart, primaryStop, assignedPositions, static_mode);
        calculate_secondary_positions(elem, *neighbor_elem, parse_left, secondaryStart, secondaryStop, assignedPositions, static_mode);
      }
      
      if (assignedPositions.find(elem.class_id) != assignedPositions.end()) {
        auto& existingEntry = assignedPositions.at(elem.class_id);
        if (std::get<0>(existingEntry) != primaryStart) {
          std::get<1>(existingEntry) = primaryStart;
        }
        if (std::get<2>(existingEntry) != primaryStop) {
          std::get<3>(existingEntry) = primaryStop;
        }
      } else {
        assignedPositions[elem.class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
      }
    }
  }
}

void process_variable_elements(const ReadLayout& readLayout, 
  const std::vector<int>& segment_orders, 
  bool parse_left, 
  const std::string& direction, 
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  auto& order_index = readLayout.get<order_tag>();
  int read_order = parse_left ? *std::max_element(segment_orders.begin(), segment_orders.end()) : *std::min_element(segment_orders.begin(), segment_orders.end());
  
  for (const auto& elem : readLayout) {
    if (elem.global_class != "read" || std::find(segment_orders.begin(), segment_orders.end(), elem.order) == segment_orders.end()) {
      continue;
    }
    
    if ((parse_left && elem.order <= read_order) || (!parse_left && elem.order >= read_order)) {
      int neighbor_order = parse_left ? elem.order - 1 : elem.order + 1;
      auto neighbor_elem = order_index.find(neighbor_order);
      
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (neighbor_elem == order_index.end()) {
        primaryStart = parse_left ? "start_of_read" : "end_of_read";
        primaryStop = primaryStart + "+" + std::to_string(*elem.expected_length);
        secondaryStart = primaryStart;
        secondaryStop = primaryStop;
      } else {
        calculate_primary_positions(elem, *neighbor_elem, parse_left, primaryStart, primaryStop, assignedPositions, false);
        calculate_secondary_positions(elem, *neighbor_elem, parse_left, secondaryStart, secondaryStop, assignedPositions, false);
      }
      
      if (assignedPositions.find(elem.class_id) != assignedPositions.end()) {
        auto& existingEntry = assignedPositions.at(elem.class_id);
        if (std::get<0>(existingEntry) != primaryStart) {
          std::get<1>(existingEntry) = primaryStart;
        }
        if (std::get<2>(existingEntry) != primaryStop) {
          std::get<3>(existingEntry) = primaryStop;
        }
      } else {
        assignedPositions[elem.class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
      }
    }
  }
}

void process_read_elements(const ReadLayout& readLayout, 
  const std::vector<int>& segment_orders, 
  const std::string& direction, 
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  process_static_elements(readLayout, segment_orders, true, true, direction, assignedPositions);
  process_static_elements(readLayout, segment_orders, false, true, direction, assignedPositions);
  process_variable_elements(readLayout, segment_orders, true, direction, assignedPositions);
  process_variable_elements(readLayout, segment_orders, false, direction, assignedPositions);
}

// Main Function

void CombinedVarScan_V2(ReadLayout& readLayout, bool verbose) {
  if (verbose) {
    Rcpp::Rcout << "Starting CombinedVarScan\n";
  }
  
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>> assignedPositions;
  std::vector<int> forward_orders, reverse_orders;
  separate_orders(readLayout, forward_orders, reverse_orders);
  
  auto& global_class_index = readLayout.get<global_class_tag>();
  auto read_range = global_class_index.equal_range("read");
  
  for (auto it_read = read_range.first; it_read != read_range.second; ++it_read) {
    std::string read_direction = it_read->direction;
    int read_order = it_read->order;
    
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
    
    process_read_elements(readLayout, left_vector, read_direction, assignedPositions);
    process_read_elements(readLayout, right_vector, read_direction, assignedPositions);
  }
  
  for (const auto& entry : assignedPositions) {
    const auto& class_id = entry.first;
    const auto& positions = entry.second;
    auto& id_index = readLayout.get<id_tag>();
    auto it = id_index.find(class_id);
    if (it != id_index.end()) {
      PositionInfo posInfo;
      posInfo.position_data = positions;
      id_index.modify(it, [&posInfo](ReadElement& elem) {
        elem.position_data = posInfo;
      });
      
      if (verbose) {
        Rcpp::Rcout << "Class ID: " << class_id << "\n"
                    << "  Primary Start: " << std::get<0>(positions) << "\n"
                    << "  Secondary Start: " << std::get<1>(positions) << "\n"
                    << "  Primary Stop: " << std::get<2>(positions) << "\n"
                    << "  Secondary Stop: " << std::get<3>(positions) << "\n";
      }
    }
  }
  if (verbose) {
    Rcpp::Rcout << "CombinedVarScan completed\n";
  }
}
