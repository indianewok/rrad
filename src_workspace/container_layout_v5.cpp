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
  std::string seq; // Sequence information
  std::string global_class;
  boost::optional<int> expected_length; // Optional expected length
  std::string type; // Type of the element (e.g., "static", "variable")
  int order; // Order in the layout
  std::string direction; // Direction (e.g., "forward", "reverse")
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

typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<length_tag>, member<ReadElement, boost::optional<int>, &ReadElement::expected_length>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>,
    ordered_non_unique<tag<global_class_tag>, member<ReadElement, std::string, &ReadElement::global_class>>
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
        Rcpp::as<std::string>(class_column[i]),
        misalignment_threshold_opt
    ));
  }

  //VarScan(container, true);
  return container;
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

void test_ref_layout(ReadLayout& container) {
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
  // Test 6: Count the number of 'read' elements in the read layout
  {
    auto& global_class_index = container.get<global_class_tag>();
    int count_reads = global_class_index.count("read");
    Rcpp::Rcout << "Number of 'read' elements: " << count_reads << "\n";
  }
  // Test 7: Count of "read" elements
  {
    auto& global_class_index = container.get<global_class_tag>();
    int read_count = global_class_index.count("read");
    Rcpp::Rcout << "Number of 'read' elements: " << read_count << "\n";
  }
  // Test 9: Variable Elements Relative to Reads
  {
    auto& global_class_index = container.get<global_class_tag>();
    auto read_elements = global_class_index.equal_range("read");

    for (auto it_read = read_elements.first; it_read != read_elements.second; ++it_read) {
      Rcpp::Rcout << "Read Element: " << it_read->class_id << " (Order: " << it_read->order << ")\n";

      // Variable Elements to the Left of the Read
      Rcpp::Rcout << "  Variable Elements to the Left:\n";
      auto& order_index = container.get<order_tag>();
      auto it_left = order_index.lower_bound(1);
      while (it_left != order_index.end() && it_left->order < it_read->order) {
        if (it_left->type == "variable" && it_left->direction == it_read->direction) {
          Rcpp::Rcout << "    " << it_left->class_id << "\n";
        }
        ++it_left;
      }

      // Variable Elements to the Right of the Read
      Rcpp::Rcout << "  Variable Elements to the Right:\n";
      auto it_right = std::next(order_index.find(it_read->order));
      while (it_right != order_index.end() && it_right->direction == it_read->direction) {
        if (it_right->type == "variable") {
          Rcpp::Rcout << "    " << it_right->class_id << "\n";
        }
        ++it_right;
      }
    }
  }
  // Test 10: MVP Primer for Each Variable
  {
    auto& type_index = container.get<type_tag>();
    auto var_elements = type_index.equal_range("variable");

    Rcpp::Rcout << "MVP Primer for Each Variable Element:\n";
    for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
      if (it_var->position_data) {
        auto& order_index = container.get<order_tag>();
        auto nearest_static_before = std::find_if(order_index.rbegin(), order_index.rend(),
          [&](const ReadElement& elem) {
            return elem.order < it_var->order && elem.type == "static" && elem.direction == it_var->direction;
          });
        if (nearest_static_before != order_index.rend()) {
          Rcpp::Rcout << "  Variable: " << it_var->class_id << ", MVP Primer: " << nearest_static_before->class_id << "\n";
        }
      }
    }
  }
  // Additional Test: Display Variable Elements and Their MVP Primers in Reverse Sequencing
  {
    auto& type_index = container.get<type_tag>();
    auto var_elements = type_index.equal_range("variable");

    Rcpp::Rcout << "Variable Elements MVP Primers (including reverse sequencing):\n";
    for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
      if (it_var->position_data) {
        auto& order_index = container.get<order_tag>();
        std::string mvp_id;

        if (it_var->direction == "forward") {
          // For forward sequencing, find nearest static element before the variable element
          auto nearest_static_before = std::find_if(order_index.rbegin(), order_index.rend(),
            [&](const ReadElement& elem) {
              return elem.order < it_var->order && elem.type == "static";
            });
          if (nearest_static_before != order_index.rend()) {
            mvp_id = nearest_static_before->class_id;
          }
        } else {
          // For reverse sequencing, find nearest static element after the variable element
          auto nearest_static_after = std::find_if(order_index.begin(), order_index.end(),
            [&](const ReadElement& elem) {
              return elem.order > it_var->order && elem.type == "static";
            });
          if (nearest_static_after != order_index.end()) {
            mvp_id = nearest_static_after->class_id;
          }
        }

        Rcpp::Rcout << "  Element: " << it_var->class_id << ", MVP Primer: " << mvp_id
                    << ", Direction: " << it_var->direction << "\n";
      }
    }
  }
  // Additional additional test: check VarScan position data
}

void parse_static_left(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, bool static_mode) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto& direction_index = readLayout.get<direction_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::max_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if(static_mode){
      if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
      continue;
    }
      if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    } else {
      if(it_var -> global_class != "read"){
        continue;
      }
    }
      if (it_var->order <= read_order) {
      // Find previous element
      if (verbose) {
        Rcpp::Rcout << "Processing varElem: " << it_var->class_id << "\n";
      }
      if (it_var->order > 1) {
        auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
        std::string prevType = prevElem->type;
        std::string primaryStart;
        std::string primaryStop;
        std::string secondaryStart;
        std::string secondaryStop;
        if (prevElem != order_index.end()) {
          // Calculate primary positions
          std::string primaryStart = prevElem->class_id + "|stop+1";
          if(static_mode){
            primaryStop = prevElem->class_id + "|stop+" + std::to_string(*it_var->expected_length - 1);
          } else {
              auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
              primaryStop = nextElem->class_id + "|start-1";
          }
          // Check for chained variables
          if(static_mode){
            if (prevType == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
            auto& prevAssigned = assignedPositions[prevElem->class_id];
            secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
            secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
            } else if (prevType == "static") {
            secondaryStart = prevElem->class_id + "|left_terminal_startlinked";
            secondaryStop = prevElem->class_id + "|left_terminal_stoplinked";
          }
          } else {
            if (prevType == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
              auto& prevAssigned = assignedPositions[prevElem->class_id];
              secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
              secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
            } else if (prevType == "static") {
              secondaryStart = prevElem->class_id + "|left_terminal_startlinked";
                auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
                  if(nextElem->global_class == "poly_tail"){
                    auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); // move one over from the poly tail if there is one
                    secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
                  } else {
                    secondaryStop = prevElem->class_id + "|left_terminal_stoplinked";
                }
            }
          }
          // Update assignedPositions
          if(assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
            auto& existingEntry = assignedPositions[it_var->class_id];
            if(std::get<0>(existingEntry) != primaryStart){
              std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
            }
            if(std::get<2>(existingEntry) != primaryStop){
              std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
            }
          } else {
            assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
          }
          if(verbose) {
            Rcpp::Rcout << "Parse Left: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
            Rcpp::Rcout << "Parse Left: Primary Stop "  << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
          }
        }
      }
      if (verbose) {
        Rcpp::Rcout << "Parse left varElem: " << it_var->class_id << "\n";
      }
    }
  }
}

void parse_static_right(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, bool static_mode) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::min_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = std::make_reverse_iterator(var_elements.second); it_var != std::make_reverse_iterator(var_elements.first); ++it_var) {
    if(static_mode){
      if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
      continue;
    }
      if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    } else {
      if(it_var -> global_class != "read"){
        continue;
      }
    }
    if (it_var->order >= read_order) {
      if (verbose) {
        Rcpp::Rcout << "Processing varElem: " << it_var->class_id << "\n";
      }
      if (it_var->order < segment_orders.back()) {
        auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
        if (nextElem != order_index.end()) {
          // Define nextType here
          std::string nextType = nextElem->type;

          std::string primaryStart;
          std::string primaryStop;
          std::string secondaryStart;
          std::string secondaryStop;

          // Calculate primary positions
          if(static_mode){
            primaryStart = nextElem->class_id + "|start-" + std::to_string(*it_var->expected_length - 1);
          } else {
              auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
              primaryStart = prevElem->class_id + "|stop+1";
          }
          primaryStop = nextElem->class_id + "|start-1";

          // Check for chained variables
          if(static_mode){
            if (nextType == "variable" && assignedPositions.find(nextElem->class_id) != assignedPositions.end()) {
              auto& nextAssigned = assignedPositions[nextElem->class_id];
              secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
              secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
            } else if (nextType == "static") {
              secondaryStart = nextElem->class_id + "|right_terminal_startlinked";
              secondaryStop = nextElem->class_id + "|right_terminal_stoplinked";
            }
          } else {
            if (nextType == "variable" && assignedPositions.find(nextElem->class_id) != assignedPositions.end()) {
              auto& nextAssigned = assignedPositions[nextElem->class_id];
              secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
              secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
            } else if (nextType == "static") {
                auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
              if(prevElem->global_class == "poly_tail"){
                auto skipElem = readLayout.get<order_tag>().find(it_var->order-2);
                secondaryStart = skipElem->class_id+"|stop+1|poly_skipped";
              } else {
                secondaryStart = prevElem->class_id + "|left_terminal_startlinked";
              }
                auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
                  if(nextElem->global_class == "poly_tail"){
                  auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); // move one over from the poly tail if there is one
                  secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
                } else {
                  secondaryStop = nextElem->class_id + "|right_terminal_stoplinked";
                }
              }
          }
          if(assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
            auto& existingEntry = assignedPositions[it_var->class_id];
            if(std::get<0>(existingEntry) != primaryStart){
              std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
            }
            if(std::get<2>(existingEntry) != primaryStop){
              std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
            }
          } else {
            assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
          }
          if(verbose) {
            Rcpp::Rcout << "Parse Right: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
            Rcpp::Rcout << "Parse Right: Primary Stop "  << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
          }
        }
      }
      if (verbose) {
        Rcpp::Rcout << "Parse right varElem: " << it_var->class_id << "\n";
      }
    }
  }
}

// VarScan function
void VarScan(ReadLayout& readLayout, bool verbose) {
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>> assignedPositions;
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

    //Parse static elements
    parse_static_left(readLayout, left_vector, read_direction, verbose, assignedPositions, true);
    parse_static_right(readLayout, right_vector, read_direction, verbose, assignedPositions, true);
    //Parse read elements
    parse_static_left(readLayout, left_vector, read_direction, verbose, assignedPositions, false);
    parse_static_right(readLayout, right_vector, read_direction, verbose, assignedPositions, false);

  }

  for (const auto& entry : assignedPositions) {
    const auto& class_id = entry.first;
    const auto& positions = entry.second;

    // Retrieve the ReadElement from readLayout using class_id
    auto& id_index = readLayout.get<id_tag>();
    auto it = id_index.find(class_id);
    if (it != id_index.end()) {
      // Update the position_data field
      PositionInfo posInfo;
      posInfo.position_data = positions;
      id_index.modify(it, [&posInfo](ReadElement& elem) {
        elem.position_data = posInfo;
      });
    }
  }
  if (verbose) {
    Rcpp::Rcout << "Final Assigned Positions:\n";
    for (const auto& entry : assignedPositions) {
      const auto& class_id = entry.first;
      const auto& positions = entry.second;
      Rcpp::Rcout << "Class ID: " << class_id << "\n"
                  << "  Primary Start: " << std::get<0>(positions) << "\n"
                  << "  Secondary Start: " << std::get<1>(positions) << "\n"
                  << "  Primary Stop: " << std::get<2>(positions) << "\n"
                  << "  Secondary Stop: " << std::get<3>(positions) << "\n";
    }
  }
}

ReadLayout vartest_cpp(const Rcpp::DataFrame& read_layout, bool verbose) {
  ReadLayout container;
  // Assume 'query_id' matches 'class_id' from 'read_layout'

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
        Rcpp::as<std::string>(class_column[i])
    ));
  }
  if(verbose){
    VarScan(container, true);
    auto& id_index = container.get<id_tag>();
    for (auto it = id_index.begin(); it != id_index.end(); ++it) {
      Rcpp::Rcout << "Class ID: " << it->class_id << "\n";
      Rcpp::Rcout << "Seq: " << it->seq << "\n";
      Rcpp::Rcout << "Type: " << it->type << "\n";
      Rcpp::Rcout << "Order: " << it->order << "\n";
      Rcpp::Rcout << "Direction: " << it->direction << "\n";

      if (it->position_data) {
        auto& posData = it->position_data.get().position_data;
        Rcpp::Rcout << "  Primary Start: " << std::get<0>(posData) << "\n";
        Rcpp::Rcout << "  Secondary Start: " << std::get<1>(posData) << "\n";
        Rcpp::Rcout << "  Primary Stop: " << std::get<2>(posData) << "\n";
        Rcpp::Rcout << "  Secondary Stop: " << std::get<3>(posData) << "\n";
      } else {
        Rcpp::Rcout << "  Position Data: Not available\n";
      }
      Rcpp::Rcout << "------------------------------------\n";
    }
  } else {
    VarScan(container, false);
  }
  return container;
}

// [[Rcpp::export]]
void test_read_layout_container_v2(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose = true) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  if (verbose) {
    test_ref_layout(container);
  }
}

// [[Rcpp::export]]
void VarTest (const Rcpp::DataFrame& read_layout, bool verbose = true){
   ReadLayout container = vartest_cpp(read_layout, verbose);
}

