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
  std::string seq; // Sequence information, there's nothing in here tbh
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
  std::string global_class;
  boost::optional<int> edit_distance; // Optional edit distance
  std::pair<int, int> position;
  std::string type;
  int order;
  std::string direction;
  boost::optional<bool> element_pass;

  SigElement(
    std::string class_id,
    std::string global_class,
    boost::optional<int> edit_distance,
    std::pair<int, int> position,
    std::string type,
    int order,
    std::string direction,
    boost::optional<bool> element_pass = boost::none
  ) : class_id(std::move(class_id)),
  global_class(std::move(global_class)),
  edit_distance(edit_distance),
  position(position),
  type(std::move(type)),
  order(order),
  direction(std::move(direction)),
    element_pass(element_pass) {}
};

struct sig_id_tag {};
struct sig_global_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};

// Define the multi-index container for SigString
typedef multi_index_container<
  SigElement,
  indexed_by<
    ordered_non_unique<tag<sig_id_tag>, member<SigElement, std::string, &SigElement::class_id>>,
    ordered_non_unique<tag<sig_global_tag>, member<SigElement, std::string, &SigElement::global_class>>,
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

struct PositionCalcFunc {
  std::function<int(const ReadData&)> primaryStartFunc;
  std::function<int(const ReadData&)> secondaryStartFunc;
  std::function<int(const ReadData&)> primaryStopFunc;
  std::function<int(const ReadData&)> secondaryStopFunc;
};

// diagnostic test function for the readlayout container
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

// Parse static elements to the left of a read
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
            primaryStop = prevElem->class_id + "|stop+" + std::to_string(*it_var->expected_length);
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
            secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
            secondaryStop = prevElem->class_id + "|none|left_terminal_linked";
          }
          } else {
            if (prevType == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
              auto& prevAssigned = assignedPositions[prevElem->class_id];
              secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
              secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
            } else if (prevType == "static") {
              secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
                auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
                  if(nextElem->global_class == "poly_tail"){
                    auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); // move one over from the poly tail if there is one
                    secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
                  } else {
                    secondaryStop = prevElem->class_id + "|none|left_terminal_linked";
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

//Parse static elements to the right of a read
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
            primaryStart = nextElem->class_id + "|start-" + std::to_string(*it_var->expected_length);
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
              secondaryStart = nextElem->class_id + "|none|right_terminal_linked";
              secondaryStop = nextElem->class_id + "|none|right_terminal_linked";
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
                secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
              }
                auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
                  if(nextElem->global_class == "poly_tail"){
                  auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); // move one over from the poly tail if there is one
                  secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
                } else {
                  secondaryStop = nextElem->class_id + "|none|right_terminal_linked";
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

// VarScan test cpp
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

// generating a function map
using PositionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;
PositionFuncMap createPositionFunctionMap(const ReadLayout& readLayout, bool verbose) {
  PositionFuncMap funcMap;
  for (const auto& element : readLayout) {
    if (element.type == "variable" && element.position_data) {
      const auto& positionData = element.position_data->position_data;
      
      auto primaryStart = std::get<0>(positionData);
      auto secondaryStart = std::get<1>(positionData);
      auto primaryStop = std::get<2>(positionData);
      auto secondaryStop = std::get<3>(positionData);

      auto parsePositionData = [verbose](const std::string& varElemClassId, const std::string& data, bool isPrimary) {
        std::stringstream ss(data);
        std::string refClassId, operation, additionalInfo;
        int offset = 0;

        std::getline(ss, refClassId, '|');
        std::getline(ss, operation, '|');
        if (!isPrimary) {
          std::getline(ss, additionalInfo, '|'); // Only for secondary positions
        }

        if(verbose) {
          Rcpp::Rcout << "Processing position data for variable element: " << varElemClassId << "\n";
          Rcpp::Rcout << "  Reference class ID: " << refClassId << "\n";
          Rcpp::Rcout << "  Operation: " << operation << "\n";
          Rcpp::Rcout << " Additional Information, if available:" << additionalInfo << "\n";
        }

        if (isPrimary || operation.find("+") != std::string::npos || operation.find("-") != std::string::npos) {
          std::string offsetStr = operation.substr(operation.find_first_of("+-"));
          offset = std::stoi(offsetStr);
        }
        return [varElemClassId, refClassId, offset, additionalInfo, &isPrimary, operation, verbose](const ReadData& readData) -> int {
          
          if(verbose){
            Rcpp::Rcout << "Made it into the lambda function!\n";
          }
          
          const auto& id_index = readData.sigstring.get<sig_id_tag>();
          
          auto refClass = id_index.find(refClassId);

          decltype(id_index.begin()) it;
          auto var_range = id_index.equal_range(varElemClassId);
          int var_order;
          
          for(auto it = var_range.first; it != var_range.second; ++it) {
             var_order = it->order;
          }
          
          auto range = id_index.equal_range(refClassId);
          size_t count = std::distance(range.first, range.second);
          if (count > 1) {
            if(verbose){
              Rcpp::Rcout << "Found more than one of this id: " << refClassId << ", and attempting to resolve for this variable position: " << varElemClassId << "\n";
            }
            auto closestIt = range.first; // Initialize with the first entry
            int closestDeltaOrder = std::numeric_limits<int>::max();
            for (auto it = range.first; it != range.second; ++it) {
              int delta_order = it->order - var_order;
              int absDeltaOrder = std::abs(delta_order);
              if ((delta_order < 0 && absDeltaOrder < closestDeltaOrder) || // Closer and to the left
                (delta_order > 0 && absDeltaOrder <= closestDeltaOrder && it->position.second > closestIt->position.second)) {
                closestIt = it;
                closestIt -> position.first;
                closestIt -> position.second;
                if(verbose){
                  Rcpp::Rcout << "Choosing " << closestIt->position.first << " as the closest start!\n";
                  Rcpp::Rcout << "Choosing " << closestIt->position.second << " as the closest stop!\n";
                }
              }
            }
          } else if (count == 1) {
            it = range.first;
          }
          if (it != id_index.end() && 
            ((refClass->type == "static" && refClass->element_pass && *refClass->element_pass) || 
            (refClass->type == "variable"))){
            if(verbose){
              Rcpp::Rcout << "Made it into the loop and past the if statements!";
            }
            int startPos = it->position.first;
            int stopPos = it->position.second;
            int basePos = (operation.find("start") != std::string::npos) ? it->position.first : it->position.second;
            int calculatedPos = basePos + offset;
            if (verbose) {
              Rcpp::Rcout << "  Variable Element: " << varElemClassId << "\n";
              Rcpp::Rcout << "  Reference Element: " << refClassId << ", Base Position: " << basePos << "\n";
              Rcpp::Rcout << " Positional Information from start to stop: "<< startPos << " and " << stopPos << "\n";
              Rcpp::Rcout << "  Calculated Position: " << calculatedPos << "\n";
            }
            return calculatedPos;
          } else if (!isPrimary && !additionalInfo.empty()) {
            if(verbose){
              string result = isPrimary ? "Primary" : "Secondary";
              Rcpp::Rcout << "Flagged for some failures, checking which one it is!\n";
              Rcpp::Rcout << result << "\n";
              Rcpp::Rcout << additionalInfo << " is the additional info!\n";
            }
            if (additionalInfo == "left_terminal_linked") {
              return 1; // Start of the sequence
            } else if (additionalInfo == "right_terminal_linked") {
              if(verbose){
                "Detected a right_terminal linkage!";
              }
              if(readData.string_info.find("length") != readData.string_info.end()){
                return readData.string_info.at("length").str_length; // End of the sequence
              }
            }
          }
          return 0; // Default return value if conditions are not met
        };
      };
      funcMap[element.class_id] = {
        parsePositionData(element.class_id, primaryStart, true),
        parsePositionData(element.class_id, secondaryStart, false),
        parsePositionData(element.class_id, primaryStop, true),
        parsePositionData(element.class_id, secondaryStop, false)
      };
    }
  }
  if (verbose) {
    Rcpp::Rcout << "Position function map created with " << funcMap.size() << " entries.\n";
  }
  return funcMap;
}

//parsing a sigstring and adding stuff to it
void fillSigString(ReadData& readData, const ReadLayout& readLayout, const PositionFuncMap& positionFuncMap, const std::string& signature, bool verbose) {
  // Process the signature and populate sigstring
  std::stringstream ss(signature);
  std::string token;
  auto& sig_id_index = readData.sigstring.get<sig_id_tag>();

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
    
  auto& class_id_index = readLayout.get<id_tag>();
  auto it = class_id_index.find(id);
  if (it != class_id_index.end()) {
    // Calculate misalignment threshold if available
    boost::optional<bool> element_pass;
    if (it->misalignment_threshold && it->type == "static") {
      int threshold_floor = std::get<0>(*it->misalignment_threshold);
      element_pass = (editDistance <= threshold_floor);
    }
      // Create SigElement with complete information
      SigElement element(
          id,
          it->global_class,
          editDistance,
          std::make_pair(startPos, endPos),
          it->type,
          it->order,
          it->direction,
          element_pass
      );
      readData.sigstring.insert(std::move(element));
      }
    }

  // Populate remaining elements from readLayout
  for (const auto& readElement : readLayout) {
    if (sig_id_index.find(readElement.class_id) == sig_id_index.end()) {
      SigElement element(
          readElement.class_id,
          readElement.global_class,
          0, // default edit distance
          std::make_pair(0, 0), // default positions
          readElement.type,
          readElement.order,
          readElement.direction
      );
      readData.sigstring.insert(std::move(element));
    }
  }

  // Update positions for variable elements using the function map
  for (auto& element : readData.sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
      if(verbose) {
        Rcpp::Rcout << "Dealing with everything that's not a read!\n";
        Rcpp::Rcout << "Working on " << element.class_id << "\n";
      }
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if(verbose){
        Rcpp::Rcout << newStartPos << "\n";
      }
       if(newStartPos <= 0){
         if(verbose){
           Rcpp::Rcout << "Triggered secondary func\n";
         }
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
       if(verbose){
         Rcpp::Rcout << newStopPos << "\n";
       }
        if(newStopPos <= 0){ 
          if(verbose){
            Rcpp::Rcout << "Triggered secondary func\n";
          }
          newStopPos = positionFuncs.secondaryStopFunc(readData);
      }
      sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
        elem.position.first = newStartPos;
        elem.position.second = newStopPos;
      });
      if (verbose) {
        Rcpp::Rcout << "Updated positions for variable class_id: " << element.class_id
                    << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
      }
    } else {
      continue;
    }
  }
  // I am incredibly lazy and wanted to defer evaluation of the reads until the end
  //and I really couldn't figure out a faster way to do it
  //even though i know there probably is a smarter implementation
  for (auto& element : readData.sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class == "read") {
      if(verbose) {
        Rcpp::Rcout << "Dealing with everything that's IS a read!\n";
      } //repetitive code--please fix this if you know what's good for you, you goddamn idiot
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if(newStartPos <= 0){
        //if(verbose){
        //  Rcpp::Rcout << "Triggered secondary func\n";
        //}
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
      if(newStopPos <= 0){ 
        newStopPos = positionFuncs.secondaryStopFunc(readData);
        }
      sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
        elem.position.first = newStartPos;
        elem.position.second = newStopPos;
      });
      if (verbose) {
        Rcpp::Rcout << "Updated positions for the read: " << element.class_id
                    << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
      }
    }
  }
}

//displaying a sigstring
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
      std::string pass_fail;
      if(element.type == "static" && element.global_class != "poly_tail"){
        pass_fail = (element.edit_distance <= 5) ? ":pass" : ":fail";
      } else {
        pass_fail = "";
      }
      if(element.position.first == 0 || element.position.second == 0){
        continue;
      }
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

// [[Rcpp::export]]
void SigTest (const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, const Rcpp::StringVector& sigstrings, bool verbose = false){
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  VarScan(container, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
  for(const auto& sigstring: sigstrings){
    ReadData readData;
    std::string signature = Rcpp::as<std::string>(sigstring);
    fillSigString(readData, container, positionFuncMap, signature, verbose);
    if(verbose){
      displayOrderedSigString(readData.sigstring);
    }
  }
}
