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
struct sig_pass_tag {}; //added to be able to index by element_pass

// Define the multi-index container for SigString
typedef multi_index_container<
  SigElement,
  indexed_by<
    ordered_non_unique<tag<sig_id_tag>, member<SigElement, std::string, &SigElement::class_id>>,
    ordered_non_unique<tag<sig_global_tag>, member<SigElement, std::string, &SigElement::global_class>>,
    ordered_non_unique<tag<sig_ed_tag>, member<SigElement, boost::optional<int>, &SigElement::edit_distance>>,
    ordered_non_unique<tag<sig_order_tag>, member<SigElement, int, &SigElement::order>>,
    ordered_non_unique<tag<sig_dir_tag>, member<SigElement, std::string, &SigElement::direction>>,
    ordered_non_unique<tag<sig_pass_tag>, member<SigElement, boost::optional<bool>, &SigElement::element_pass>>
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
              primaryStart = primaryStart + "|chained_start";
              primaryStop = primaryStop + "|chained_stop";
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
              primaryStart = primaryStart + "|chained_start";
              primaryStop = primaryStop + "|chained_stop";
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

//StatCounts from a sigstring, tabulate forward and reverse elements
stat_elem StatCounts(const SigString& sigstring) {
  stat_elem counts;
  std::set<std::string> seen_forw_pass_ids;
  std::set<std::string> seen_rev_pass_ids;
  for (const auto& element : sigstring) {
    if (element.global_class == "poly_tail" || element.type != "static") continue;
    if (element.direction == "forward" && element.element_pass && *element.element_pass) {
      ++counts.forw_pass;
      // Increment unique counter if it's a new ID
      if (seen_forw_pass_ids.insert(element.class_id).second) {
        ++counts.unique_forw_pass;
      }
    } else if (element.direction == "reverse" && element.element_pass && *element.element_pass) {
      ++counts.rev_pass;
      // Increment unique counter if it's a new ID
      if (seen_rev_pass_ids.insert(element.class_id).second) {
        ++counts.unique_rev_pass;
      }
    }
  }
  for(const auto& element : sigstring){
    if ((element.global_class == "forw_primer" || element.global_class == "rev_primer") && element.element_pass && *element.element_pass) {
      if (element.direction == "forward") {
        ++counts.term_forw_elements;
      } else if (element.direction == "reverse") {
        ++counts.term_rev_elements;
      }
    }
  }
  return counts;
}

//concatenation scanning
void concat_scan(ReadData& readData, bool verbose) {
  stat_elem counts = StatCounts(readData.sigstring);
  if(verbose){
    Rcpp::Rcout << "Forward pass: " << counts.forw_pass << "\nReverse pass: " << counts.rev_pass << "\n";
    Rcpp::Rcout << "Forward UNIQUE pass: " << counts.unique_forw_pass << "\nReverse UNIQUE pass: " << counts.unique_rev_pass << "\n";
    Rcpp::Rcout << "Forward terminal pass: " << counts.term_forw_elements << "\nReverse terminal pass: " << counts.term_rev_elements << "\n";
  }
  if(counts.forw_pass == 0 && counts.rev_pass > 0){
    readData.string_info["type"].str_type = "R";
  } else if (counts.rev_pass == 0 && counts.forw_pass > 0){
    readData.string_info["type"].str_type = "F";
  }
  if(counts.forw_pass != counts.unique_forw_pass && counts.rev_pass == 0){
    readData.string_info["type"].str_type = "F_F";
  } else if (counts.rev_pass != counts.unique_rev_pass && counts.forw_pass == 0){
    readData.string_info["type"].str_type = "R_R";
  }
  if(counts.forw_pass > 0 && counts.rev_pass > 0){
    readData.string_info["type"].str_type = "FR_RF";
  }
  if(verbose){
    Rcpp::Rcout << "This sigstring is a type " << readData.string_info["type"].str_type << "!\n";
    Rcpp::Rcout << "This read is " << readData.string_info["type"].str_length << " bases long!\n";
  }
}

void concat_solve(ReadData& readData, bool verbose) {
  std::vector<SigElement> sorted_elements;
  auto& order_index = readData.sigstring.get<sig_order_tag>();
  auto& id_index = readData.sigstring.get<sig_id_tag>();
  
  std::copy(order_index.begin(), order_index.end(), std::back_inserter(sorted_elements));
  // Assuming reads are identified and need position adjustment based on concatenation logic
  int readStart = -1;
  int readEnd = -1;
  for (const auto& element : sorted_elements) {
    // Identify the start and end positions of the read based on element order
    if (element.type == "read") {
      if (readStart == -1 || element.position.first < readStart) {
        readStart = element.position.first;
      }
      if (readEnd == -1 || element.position.second > readEnd) {
        readEnd = element.position.second;
      }
    }
  }
  for (auto& element : readData.sigstring) {
    if (element.type == "read") {
      Rcpp::Rcout << "Read Start and read end:" <<  readStart << " and " << readEnd << "\n";
    }
  }

  if (verbose) {
    Rcpp::Rcout << "Updated read positions based on concatenation logic.\n";
    Rcpp::Rcout << "New Start: " << readStart << ", New Stop: " << readEnd << "\n";
  }
}

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
        if (verbose) {
          Rcpp::Rcout << "Processing position data for variable element: " << varElemClassId << "\n";
          Rcpp::Rcout << "  Reference class ID: " << refClassId << "\n";
          Rcpp::Rcout << "  Operation: " << operation << "\n";
          Rcpp::Rcout << "  Additional Information, if available: " << additionalInfo << "\n";
        }
        if (operation.find("+") != std::string::npos || operation.find("-") != std::string::npos) {
          std::string offsetStr = operation.substr(operation.find_first_of("+-"));
          offset = std::stoi(offsetStr);
        }
        return [varElemClassId, refClassId, offset, additionalInfo, &isPrimary, 
          operation, verbose](const ReadData& readData) -> int {
            if (verbose) {
              Rcpp::Rcout << "Made it into the lambda function!\n";
            }
            
            const auto& id_index = readData.sigstring.get<sig_id_tag>();
            const auto& order_index = readData.sigstring.get<sig_order_tag>();
            
            auto refClass = id_index.find(refClassId);
            auto it = id_index.equal_range(refClassId).first;
            
            if (it != id_index.end() && 
              ((refClass->type == "static" && refClass->element_pass && *refClass->element_pass) || 
              (refClass->type == "variable") || (refClass->global_class == "poly_tail"))) {
              if (verbose) {
                Rcpp::Rcout << "Made it into the loop and past the if statements!\n";
              }
              int startPos = it->position.first;
              int stopPos = it->position.second;
              int basePos = (operation.find("start") != std::string::npos) ? it->position.first : it->position.second;
              int calculatedPos = basePos + offset;
              if (verbose) {
                Rcpp::Rcout << "  Variable Element: " << varElemClassId << "\n";
                Rcpp::Rcout << "  Reference Element: " << refClassId << ", Base Position: " << basePos << "\n";
                  Rcpp::Rcout << "  Positional Information from start to stop: " << startPos << " and " << stopPos << "\n";
                  Rcpp::Rcout << "  Calculated Position: " << calculatedPos << "\n";
              }
              return calculatedPos;
            } else if (!isPrimary && !additionalInfo.empty()) {
              
              std::string readType = readData.string_info.at("type").str_type;
              auto varClass = id_index.find(varElemClassId);
              
              if(verbose){
                string result = isPrimary ? "Primary" : "Secondary";
                Rcpp::Rcout << "Flagged for some failures, checking which one it is!\n";
                Rcpp::Rcout << result << "\n";
                Rcpp::Rcout << additionalInfo << " is the additional info!\n";
              }
              
              if (additionalInfo == "left_terminal_linked") {
                if(verbose){
                  "Detected a left_terminal linkage!";
                }
                if(varClass->global_class == "read"){
                   if(readType == "FR_RF"){
                  return 1;
                     } 
                else {
                  return 0;
                  }
               }
              } else if (additionalInfo == "right_terminal_linked") {
                if(verbose){
                  "Detected a right_terminal linkage!";
                  Rcpp::Rcout << "Global class is " << varClass->global_class << "!\n";
                }
                if(varClass->global_class == "read"){
                  if(readData.string_info.find("type") != readData.string_info.end()){
                    return readData.string_info.at("type").str_length; // End of the sequence
                  }
                } else {
                  return 0;
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

void pos_scan(ReadData& readData, const PositionFuncMap& positionFuncMap, bool verbose) {
  std::string readType = readData.string_info["type"].str_type;
  auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
  if(readType == "F" || readType == "R"){
    for (auto it = readData.sigstring.begin(); it != readData.sigstring.end(); ) {
      if ((readType == "F" && it->direction != "forward") || 
        (readType == "R" && it->direction != "reverse")) {
        it = readData.sigstring.erase(it); // Erase and move to the next element
      } else {
        ++it; // Move to the next element
      }
    }
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
  } else if(readType == "FR_RF"){
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
    //concat_solve(readData, verbose);
    // for (auto it = readData.sigstring.begin(); it != readData.sigstring.end(); ) {
    //   // Check if element_pass has a value and if the value is true
    //   if((it->element_pass && *(it->element_pass)) || it->type == "variable") {
    //     ++it; // Keep this element, move to the next
    //   }
    //   else {
    //     it = readData.sigstring.erase(it); // Erase this element
    //   }
    // }
  }
}

//parsing a sigstring and adding stuff to it
void fillSigString(ReadData& readData, const ReadLayout& readLayout, 
  const std::string& signature, const PositionFuncMap& positionFuncMap, bool verbose) {
  //adding information to "type"
  auto read_info_pos = signature.find_last_of('<');
  std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2); // Remove '<' and '>'
  std::stringstream lengthTypeStream(lengthTypeStr);
  std::string lengthStr, type;
  std::getline(lengthTypeStream, lengthStr, ':');
  std::getline(lengthTypeStream, type);
  
  int length = std::stoi(lengthStr);
  if(verbose){
    Rcpp::Rcout << "Length is " << length << "\n";
  }
  // Populate StringInfo with length and type
  readData.addStringInfo("type", length, type);
  
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
  pos_scan(readData, positionFuncMap, verbose);
}

//displaying a sigstring
void displayOrderedSigString(const ReadData& readData) {
  auto printOrderedDirection = [&](const std::string& direction) {
    // Filter by direction and sort by order
    std::vector<SigElement> sorted_elements;
    auto& dir_index = readData.sigstring.get<sig_dir_tag>();
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
        pass_fail = (element.element_pass == true) ? ":pass" : ":fail";
        if(pass_fail == ":fail"){
          continue;
        }
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
  if(readData.string_info.at("type").str_type == "F"){
    printOrderedDirection("forward");
  } else if(readData.string_info.at("type").str_type == "R"){
    printOrderedDirection("reverse");
  } else {
    printOrderedDirection("forward");
    printOrderedDirection("reverse");
  }
}


// [[Rcpp::export]]
void SigTest (const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  const Rcpp::StringVector& sigstrings, bool verbose = false){
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  VarScan(container, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
  for(const auto& sigstring: sigstrings){
    ReadData readData;
    std::string signature = Rcpp::as<std::string>(sigstring);
    fillSigString(readData, container, signature, positionFuncMap, verbose);
    if(verbose){
      displayOrderedSigString(readData);
    }
  }
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigalign_deprecated(
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
  for(int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  omp_set_num_threads(nthreads);
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence =sequences[i];
    std::ostringstream signature;
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    // Parallelize the inner loop through queries
    // Loop through each sequence
#pragma omp parallel for
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      int query_counter = j + 1;  // Local query_counter, initialized to the loop index + 1
      if (query_name == "poly_a" || query_name == "poly_t") {
        // Handle the regex query here.
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTail(sequence, poly_base, 14, 12);  // assuming window_size = 16 and min_count = 12
        if (!match_str.empty()){
#pragma omp critical
{
  temp_signature_parts.push_back(match_str);
}
        }
      } else {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
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
        // Remove duplicates from start_positions
        auto last = std::unique(start_positions.begin(), start_positions.end());
        start_positions.erase(last, start_positions.end());
        // Find minimal 'end' for each unique 'start'
        std::vector<int> unique_end_positions;
        for (const auto& start : start_positions) {
          auto pos = std::find(start_positions.begin(), start_positions.end(), start) - start_positions.begin();
          unique_end_positions.push_back(end_positions[pos]);
        }
        int edit_distance = cresult.editDistance;
        for (size_t i = 0; i < start_positions.size(); ++i) {
          UniqueAlignment ua = {edit_distance, start_positions[i], end_positions[i]};
          uniqueAlignments.insert(ua);
        }
        // Extract the two best unique alignments
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1}; ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;  // Variable to control whether to skip appending the second-best alignment
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
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        // Assumes the format is "query_name:edit_distance:start_position:end_position"
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
}
temp_signature_parts.clear();
  }
  // Convert signature_map to something R-friendly
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

