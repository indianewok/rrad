#include "rad.h"
using namespace std;
using namespace Rcpp;
//preparing the read layout container

ReadLayout container;
PositionFuncMap positionFuncMap;

ReadLayout prep_read_layout_cpp(
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    bool verbose) {
  ReadLayout container;
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  // Populate the thresholds map
  for (int i = 0; i < query_id.size(); ++i) {
    double threshold = misal_threshold[i];
    double sd = misal_sd[i];
    
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
  
  // Iterate over the read_layout and populate the ReadLayout container
  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);
    std::string class_value = Rcpp::as<std::string>(class_column[i]);
    
    boost::optional<std::tuple<int, int, int>> misalignment_threshold_opt;
    if (thresholdsMap.count(cid)) {
      misalignment_threshold_opt = boost::optional<std::tuple<int, int, int>>(thresholdsMap[cid]);
    } else {
      misalignment_threshold_opt = boost::none;
    }
    
    boost::optional<int> expected_length_opt;
    if (class_value != "read") {
      expected_length_opt = expected_length[i];
    }
    
    // Create a new ReadElement and insert it into the container
    ReadElement elem(
        cid,
        Rcpp::as<std::string>(seq[i]),
        expected_length_opt,
        Rcpp::as<std::string>(type[i]),
        order[i],
        Rcpp::as<std::string>(direction[i]),
        class_value,
        misalignment_threshold_opt
    );
    
    container.insert(elem);
    
    if(verbose){
      Rcpp::Rcout << "Inserted ReadElement:\n"
                  << "  class_id: " << elem.class_id << "\n"
                  << "  seq: " << elem.seq << "\n"
                  << "  expected_length: " << (elem.expected_length ? std::to_string(*elem.expected_length) : "none") << "\n"
                  << "  type: " << elem.type << "\n"
                  << "  order: " << elem.order << "\n"
                  << "  direction: " << elem.direction << "\n"
                  << "  global_class: " << elem.global_class << "\n"
                  << "  misalignment_threshold: "
                  << (elem.misalignment_threshold ? 
      std::to_string(std::get<0>(*elem.misalignment_threshold)) + ", " +
      std::to_string(std::get<1>(*elem.misalignment_threshold)) + ", " +
      std::to_string(std::get<2>(*elem.misalignment_threshold)) 
        : "none") << "\n";
    }
  }
  
  return container;
}

// [[Rcpp::export]]
void prep_read_layout_rcpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose) {
  prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
}

void parse_static_left(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::max_element(segment_orders.begin(), segment_orders.end());
  
  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if (it_var->type != "variable" || it_var->direction != direction) {
      if (verbose) {
        Rcpp::Rcout << "Skipping " << it_var->class_id << "\n";
      }
      continue;
    }
    
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    
    if (it_var->order <= read_order && it_var->order > 1) {
      auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
      auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
      auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
      auto next_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
      
      bool prevElemExists = (prevElem != order_index.end());
      bool nextElemExists = (nextElem != order_index.end());
      bool prePrevElemExists = (pre_prevElem != order_index.end());
      bool nextNextElemExists = (next_nextElem != order_index.end());
      
      //int nextElem_val = (nextElemExists && (nextElem->global_class == "stop")) ? 0 : 1;
      int prevElem_val = (prevElemExists && (prevElem->global_class == "start")) ? 0 : 1;
      int prevElem_stop_val = (prevElemExists && (prevElem->global_class == "start")) ? 1 : 0;
      
      
      if (verbose) {
        Rcpp::Rcout << "Processing varElem, global class: " << it_var->global_class << "\n";
        Rcpp::Rcout << "Processing varElem, class ID: " << it_var->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element: " << prevElem->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element, global_class: " << prevElem->global_class << "\n";
        if (prePrevElemExists) Rcpp::Rcout << "Pre-previous element: " << pre_prevElem->class_id << "\n";
        if (nextElemExists) Rcpp::Rcout << "Next element: " << nextElem->class_id << "\n";
        if (nextNextElemExists) Rcpp::Rcout << "Next-next element: " << next_nextElem->class_id << "\n";
      }
      
      int offset_result = 0;
      bool offset_available = true;
      bool is_read = false;
      if (it_var->expected_length && it_var->global_class != "read") {
        offset_result = (*it_var->expected_length);
      } else {
        offset_result = 0;
        offset_available = false;
        is_read = true;
      }
      if (verbose) {
        Rcpp::Rcout << "Offset result is " << std::to_string(offset_result) << "\n";
        Rcpp::Rcout << "Offset available: " << (offset_available ? "true" : "false") << "\n";
      }
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (prevElemExists) {
        primaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
        if (is_read) {
          if (verbose) {
            Rcpp::Rcout << "Found the read.\n";
            if (prevElemExists) Rcpp::Rcout << "prevElem exists: true\nprevElem class_id: " << prevElem->class_id << "\n";
            if (nextElemExists) Rcpp::Rcout << "nextElem exists: true\nnextElem class_id: " << nextElem->class_id << "\n";
          }
          bool prevPolyTail = prevElemExists && prevElem->global_class == "poly_tail";
          bool nextPolyTail = nextElemExists && nextElem->global_class == "poly_tail";
          bool prevShortStatic = prevElemExists && prevElem->type == "static" && *prevElem->expected_length <= 13 && prevElem->global_class != "poly_tail";
          bool nextShortStatic = nextElemExists && nextElem->type == "static" && *nextElem->expected_length <= 13 && nextElem->global_class != "poly_tail";
          if (prevPolyTail || prevShortStatic) {
            if (prePrevElemExists && pre_prevElem->type == "static") {
              if (prevPolyTail) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|left_terminal_linked";
              } else {
                int prevElem_offset = *prevElem->expected_length;
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset) + "|skipped";
              }
            }
          }
          if (nextPolyTail || nextShortStatic) {
            if (nextNextElemExists && next_nextElem->type == "static") {
              secondaryStop = next_nextElem->class_id + "|start-1";
              if (nextPolyTail) {
                secondaryStop += "|skipped";
              }
            }
          }
        } else {
          primaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result - prevElem_stop_val);
        }
        
        if (prevElemExists && prevElem->type == "variable" &! is_read) {
          auto& prevAssigned = assignedPositions[prevElem->class_id];
          primaryStart += "|chained_start";
          primaryStop += "|chained_stop";
          if (nextElemExists && nextElem->type == "static") {
            secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
            secondaryStop = nextElem->class_id + "|start-1";
          } else {
            secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
            secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
          }
        } else if (prevElemExists && prevElem->type == "static" &! is_read) {
          if (nextElemExists && nextElem->type == "static" && (nextElem->expected_length >= 13) && nextElem->global_class != "poly_tail") {
            secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
            secondaryStop = nextElem->class_id + "|start-1";
          } else {
            if (nextNextElemExists && next_nextElem->type == "static" && nextElem->global_class != "read" && next_nextElem->global_class != "poly_tail") {
              int nextElem_offset = *nextElem->expected_length;
              secondaryStart = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset+offset_result) + "|skipped";
              secondaryStop = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset) + "|skipped";
            } else {
              if (prevElem->global_class == "poly_tail" && prePrevElemExists && pre_prevElem->global_class != "read") {
                secondaryStart = pre_prevElem->class_id + "|stop+1|skipped";
              } else {
                secondaryStart = prevElem->class_id + "|stop+1|left_terminal_linked";
                secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result) + "|left_terminal_linked";
              }
            }
          }
        }
        
        // Update assignedPositions
        if (assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
          auto& existingEntry = assignedPositions[it_var->class_id];
          if (std::get<0>(existingEntry) != primaryStart) {
            std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
          }
          if (std::get<2>(existingEntry) != primaryStop) {
            std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
          }
        } else {
          assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
        }
        if (verbose) {
          Rcpp::Rcout << "Parse Left: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
          Rcpp::Rcout << "Parse Left: Primary Stop " << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
        }
      }
    }
    if (verbose) {
      Rcpp::Rcout << "Completed parsing left varElem: " << it_var->class_id << "\n";
    }
  }
}

void parse_static_right(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::min_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = std::make_reverse_iterator(var_elements.second); it_var != std::make_reverse_iterator(var_elements.first); ++it_var) {
    if (it_var->type != "variable" || it_var->direction != direction) {
      if (verbose) {
        Rcpp::Rcout << "Skipping " << it_var->class_id << "\n";
      }
      continue;
    }
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    
    if (it_var->order >= read_order) {
      auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
      auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
      auto next_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
      auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
      
      bool nextElemExists = (nextElem != order_index.end());
      bool prevElemExists = (prevElem != order_index.end());
      bool nextNextElemExists = (next_nextElem != order_index.end());
      bool prePrevElemExists = (pre_prevElem != order_index.end());
      
      int nextElem_val = (nextElemExists && (nextElem->global_class == "stop")) ? 0 : 1;
      int prevElem_val = (prevElemExists && (prevElem->global_class == "start")) ? 0 : 1;
      int nextElem_stop_val = (nextElemExists && (nextElem->global_class == "stop")) ? 1 : 0;
      
      if (verbose) {
        Rcpp::Rcout << "Processing varElem, global class: " << it_var->global_class << "\n";
        Rcpp::Rcout << "Processing varElem, class ID: " << it_var->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element: " << prevElem->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element, global class: " << prevElem->global_class << "\n";
        if (prePrevElemExists) Rcpp::Rcout << "Pre-previous element: " << pre_prevElem->class_id << "\n";
        if (nextElemExists) Rcpp::Rcout << "Next element: " << nextElem->class_id << "\n";
        if (nextNextElemExists) Rcpp::Rcout << "Next-next element: " << next_nextElem->class_id << "\n";
      }
      int offset_result = 0;
      bool offset_available = true;
      bool is_read = false;
      if (it_var->expected_length && it_var->global_class != "read") {
        offset_result = (*it_var->expected_length);
      } else {
        offset_available = false;
        is_read = true;
      }
      if (verbose) {
        Rcpp::Rcout << "Offset result is " << std::to_string(offset_result) << "\n";
        Rcpp::Rcout << "Offset available: " << (offset_available ? "true" : "false") << "\n";
      }
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (prevElemExists && nextElemExists) {
        primaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
        primaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val+nextElem_stop_val);
        secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val) + "|left_terminal_linked";
        secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
        // Handle read element separately
        if (is_read) {
          if (verbose) {
            Rcpp::Rcout << "Found the read.\n";
            if (prevElemExists) Rcpp::Rcout << "prevElem exists: true\nprevElem class_id: " << prevElem->class_id << "\n";
            if (nextElemExists) Rcpp::Rcout << "nextElem exists: true\nnextElem class_id: " << nextElem->class_id << "\n";
          }
          // Check for poly_tail or short static element in previous or next elements
          bool prevPolyTail = prevElemExists && prevElem->global_class == "poly_tail";
          bool nextPolyTail = nextElemExists && nextElem->global_class == "poly_tail";
          bool prevShortStatic = prevElemExists && prevElem->type == "static" && *prevElem->expected_length <= 13 && prevElem->global_class != "poly_tail";
          bool nextShortStatic = nextElemExists && nextElem->type == "static" && *nextElem->expected_length <= 13 && nextElem->global_class != "poly_tail";
          if(verbose){
            if(prevPolyTail) Rcpp::Rcout << "Previous element is a poly tail!\n";
            if(nextPolyTail) Rcpp::Rcout << "Next element is a poly tail!\n";
          }
          if (prevPolyTail || prevShortStatic) {
            if (prePrevElemExists) {
              if (prevPolyTail) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|skipped";
              } else {
                int prevElem_offset = *prevElem->expected_length;
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset) + "|skipped";
              }
            }
          }
          if (nextPolyTail || nextShortStatic) {
            if (nextNextElemExists) {
              if (nextPolyTail) {
                secondaryStop = next_nextElem->class_id + "|start-1|skipped";
              } else {
                int nextElem_offset = *nextElem->expected_length;
                secondaryStop = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset) + "|skipped";
              }
            }
          }
          else {
            secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
          }
        } else {
          primaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
        }
        
        if (nextElemExists && nextElem->type == "variable" && !is_read) {
          auto& nextAssigned = assignedPositions[nextElem->class_id];
          primaryStart += "|chained_start";
          primaryStop += "|chained_stop";
          if (prevElemExists && prevElem->type == "static") {
            secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
            secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result);
          } else {
            secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
            secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
          }
        } else if (nextElemExists && nextElem->type == "static" && !is_read) {
          if (prevElemExists && prevElem->type == "static") {
            secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
            secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result);
          } else {
            if (prePrevElemExists && pre_prevElem->type == "static" && prevElem->global_class != "read") {
              int prevElem_offset = *prevElem->expected_length;
              secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset + 1);
              secondaryStop = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset + offset_result);
            } else if (prePrevElemExists && pre_prevElem->type == "variable" && prevElem->global_class != "read") {
              if (prevElem->global_class == "poly_tail"||*prevElem->expected_length <= 13) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|skipped";
                secondaryStop = pre_prevElem->class_id + "|stop+" + std::to_string(offset_result) + "|skipped";
                } //added here
              } //else {
                else if(nextNextElemExists && next_nextElem->type == "variable" && nextElem->global_class != "read"){
                  int nextElem_offset = *nextElem->expected_length;
                  secondaryStart = next_nextElem->class_id + "|start-" + std::to_string(offset_result + nextElem_offset) + "|chained_start";
                  secondaryStop = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset) + "|chained_stop";
                } else {
                  secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result) + "|right_terminal_linked";
                  secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
                }
              //}
            }
          }
        //}removed
        
        // Update assignedPositions
        if (assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
          auto& existingEntry = assignedPositions[it_var->class_id];
          if (std::get<0>(existingEntry) != primaryStart) {
            std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
          }
          if (std::get<2>(existingEntry) != primaryStop) {
            std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
          }
        } else {
          assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
        }
        if (verbose) {
          Rcpp::Rcout << "Parse Right: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
          Rcpp::Rcout << "Parse Right: Primary Stop " << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
        }
      }
    }
    if (verbose) {
      Rcpp::Rcout << "Completed parsing right varElem: " << it_var->class_id << "\n";
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
    
    if (verbose) {
      Rcpp::Rcout << "Read direction: " << read_direction << ", Read order: " << read_order << "\n";
      Rcpp::Rcout << "Left vector: ";
      for (const auto& order : left_vector) {
        Rcpp::Rcout << order << " ";
      }
      Rcpp::Rcout << "\nRight vector: ";
      for (const auto& order : right_vector) {
        Rcpp::Rcout << order << " ";
      }
      Rcpp::Rcout << "\n";
    }
    
    //Parse static elements
    parse_static_left(readLayout, left_vector, read_direction, verbose, assignedPositions);
    parse_static_right(readLayout, right_vector, read_direction, verbose, assignedPositions);
    //Parse read elements
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

// [[Rcpp::export]]
void generate_position_setup_test(Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold) {
  container = prep_read_layout_cpp(read_layout, misalignment_threshold, true);
  VarScan(container, true);
  positionFuncMap = createPositionFunctionMap(container, true);
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

void filter_sigstrings(ReadData& readData, 
                         const PositionFuncMap& positionFuncMap, 
                         const ReadLayout& readLayout, 
                         bool verbose = false) {
    // Helper to find ReadElement by class_id
    auto findReadElementById = [&](const std::string& cid) -> const ReadElement* {
      auto& id_index = readLayout.get<id_tag>();
      auto it = id_index.find(cid);
      return (it != id_index.end()) ? &(*it) : nullptr;
    };
    
    // Check approximate length match (±20%)
    auto lengthMatchesExpected = [&](int actual_length, int expected_length) {
      double tolerance = 0.2;
      double lower = expected_length * (1.0 - tolerance);
      double upper = expected_length * (1.0 + tolerance);
      return (actual_length >= lower && actual_length <= upper);
    };
    
    // This lambda filters a single SigString and also checks for read/barcode presence
    auto filterLambda = [&](SigString& sig, bool primary = true) {
      std::vector<SigString::iterator> toRemove;
      bool read_or_barcode = false;
      for (auto it = sig.begin(); it != sig.end(); ++it) {
        if (it->type != "variable") {
          continue;
        }
        // Check if this is a read-like element
        if (it->global_class == "read"||it->global_class == "barcode") {
          read_or_barcode = true;
        }
        
        int startPos = it->position.first;
        int endPos = it->position.second;
        
        if (verbose) {
          std::cout << "Checking: " << it->class_id << " [" << startPos << "," << endPos << "]\n";
        }
        
        // Position validity
        if (startPos > endPos||(startPos < 0 && endPos < 0)||(it->global_class == "read" && endPos-startPos <= 100)) {
          if (verbose) std::cout << "Removing due to invalid positions.\n";
          toRemove.push_back(it);
          continue;
        }
        
        int actual_length = endPos - startPos + 1;
        // Check expected length if present
        if (const ReadElement* re = findReadElementById(it->class_id)) {
          if (re->expected_length && !lengthMatchesExpected(actual_length, *re->expected_length)) {
            if (verbose) {
              std::cout << "Removing " << it->class_id << " due to length mismatch. "
                        << "Actual: " << actual_length << ", Expected: " << *re->expected_length << "\n";
            }
            toRemove.push_back(it);
            continue;
          }
        }
      }
      
      // Remove invalid elements
      for (auto rit = toRemove.rbegin(); rit != toRemove.rend(); ++rit) {
        sig.erase(*rit);
      }
      if (readData.string_info.find("type") != readData.string_info.end()) {
        if(verbose) std::cout << readData.string_info["type"].str_type << "\n";
        if (!read_or_barcode) {
          if(primary){
            std::string readtype = readData.string_info["type"].str_type;
            readData.string_info["type"].str_type = readtype + "read_or_barcode_undecided";
          } else {
            std::string readtype = readData.string_info["type"].str_type;
            readData.string_info["type"].str_type = readtype + "read_or_barcode_undecided";
          }
        } 
      }
    };
    if (verbose) {
      std::cout << "Filtering the first string!\n";
    }
    filterLambda(readData.sigstring);
    // If there's a secondary sigstring, filter that as well
    if (readData.secondary_sigstring) {
      if (verbose) {
        std::cout << "Filtering the second string!\n";
      }
      filterLambda(*readData.secondary_sigstring, false);
    }
    
    if (verbose) {
      std::cout << "Filtering complete.\n";
    }
}

std::vector<std::string> processSigString(const std::string& sigstring) {
  std::vector<std::string> parts;
  size_t splitPos = sigstring.find("+");
  
  // Extracting the "<length:id:type>" part correctly
  size_t infoStartPos = sigstring.rfind("<") + 1;
  size_t infoEndPos = sigstring.find(">", infoStartPos);
  std::string infoPart = sigstring.substr(infoStartPos, infoEndPos - infoStartPos);
  // Parsing the length, ID, and type from the infoPart
  std::istringstream infoStream(infoPart);
  std::string lengthStr, id, type;
  std::getline(infoStream, lengthStr, ':');
  std::getline(infoStream, id, ':');
  std::getline(infoStream, type);
  if (splitPos != std::string::npos) {
    std::string firstPart = sigstring.substr(0, splitPos);
    std::string secondPart = sigstring.substr(splitPos + 1);
    // Finding the last ":" before "<" in each part to capture the correct length
    size_t lastColonInFirst = firstPart.rfind(":", firstPart.rfind("<"));
    std::string lengthFirst = firstPart.substr(lastColonInFirst + 1, firstPart.rfind("<") - lastColonInFirst - 1);
    size_t lastColonInSecond = secondPart.rfind(":", secondPart.rfind("<"));
    std::string lengthSecond = secondPart.substr(lastColonInSecond + 1, secondPart.rfind("<") - lastColonInSecond - 1);
    // Reconstructing the first and second parts with updated length and type
    size_t endFirstPart = firstPart.rfind("<");
    firstPart = firstPart.substr(0, endFirstPart) + "<" + lengthFirst + ":" + id + "+FR_RF:F>";
    size_t endSecondPart = secondPart.rfind("<");
    secondPart = secondPart.substr(0, endSecondPart) + "<" + lengthSecond + ":" + id + "+FR_RF:R>";
    parts.push_back(firstPart);
    parts.push_back(secondPart);
  } else {
    parts.push_back(sigstring); // Handling cases without "+"
  }
  return parts;
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

PositionFuncMap createPositionFunctionMap(const ReadLayout& readLayout, bool verbose) {
  PositionFuncMap funcMap;
  for (const auto& element : readLayout) {
    if (element.type == "variable" && element.position_data) {
      const auto& positionData = element.position_data->position_data;
      auto primaryStart = std::get<0>(positionData);
      auto secondaryStart = std::get<1>(positionData);
      auto primaryStop = std::get<2>(positionData);
      auto secondaryStop = std::get<3>(positionData);
      
      auto parsePositionData = [verbose](const std::string& varElemClassId, const std::string& data, bool isPrimary, bool isStart) {
        
        std::stringstream ss(data);
        std::string refClassId, operation, additionalInfo;
        int offset = 0;
        std::getline(ss, refClassId, '|');
        std::getline(ss, operation, '|');
        if (!isPrimary) {
          std::getline(ss, additionalInfo);
        }
        if (operation.find("+") != std::string::npos || operation.find("-") != std::string::npos) {
          std::string offsetStr = operation.substr(operation.find_first_of("+-"));
          offset = std::stoi(offsetStr);
        }
        if (verbose) {
          Rcpp::Rcout << "Processing position data for variable element: " << varElemClassId << "\n"
                      << "  Reference class ID: " << refClassId << "\n"
                      << "  Operation: " << operation << "\n"
                      << "  Additional Information: " << additionalInfo << "\n"
                      << "  Offset: " << offset << "\n";
        }
        
        return [varElemClassId, refClassId, offset, additionalInfo, isPrimary, isStart, operation, verbose]
        (const ReadData& readData) -> int {
          
          const auto& id_index = readData.sigstring.get<sig_id_tag>();
          const auto& order_index = readData.sigstring.get<sig_order_tag>();

          auto range = id_index.equal_range(refClassId);
          auto refClass = id_index.find(refClassId);
          auto varClass = id_index.find(varElemClassId);
          
         if (std::distance(range.first, range.second) > 1) {
            if(verbose){
              Rcpp::Rcout << "Multiple elements found with ID: " << refClassId << "\n";
              for (auto it = range.first; it != range.second; ++it) {
                Rcpp::Rcout << "Class ID: " << it->class_id << "\n";
                Rcpp::Rcout << "Type: " << it->type << "\n";
                Rcpp::Rcout << "Global Class: " << it->global_class << "\n";
                Rcpp::Rcout << "Start Position: " << it->position.first << "\n";
                Rcpp::Rcout << "Stop Position: " << it->position.second << "\n";
              }
            }
            
            if(varClass->direction == "reverse"){
              refClass = range.first;
              if(verbose){
                Rcpp::Rcout << "Solving for a reverse element--defaulting to the first element in range, assuming it'll be to our right.\n";
                Rcpp::Rcout << "Class ID: " << refClass->class_id << "\n";
                Rcpp::Rcout << "Type: " << refClass->type << "\n";
                Rcpp::Rcout << "Global Class: " << refClass->global_class << "\n";
                Rcpp::Rcout << "Start Position: " << refClass->position.first << "\n";
                Rcpp::Rcout << "Stop Position: " << refClass->position.second << "\n";
              }
            } else {
              refClass = std::prev(range.second);
              if(verbose){
                Rcpp::Rcout << "Not a reverse element--defaulting to the last element in range, assuming it'll be to our left.\n";
                Rcpp::Rcout << "Class ID: " << refClass->class_id << "\n";
                Rcpp::Rcout << "Type: " << refClass->type << "\n";
                Rcpp::Rcout << "Global Class: " << refClass->global_class << "\n";
                Rcpp::Rcout << "Start Position: " << refClass->position.first << "\n";
                Rcpp::Rcout << "Stop Position: " << refClass->position.second << "\n";
              }
            }
          }
          
          std::string readType = readData.string_info.at("type").str_type;
          bool isValidElement = true;
          bool isStatic = false;
          bool hasElementPass = false;
          bool elementPassIsTrue = false;
          
          if (verbose) {
            Rcpp::Rcout << "Processing " << varElemClassId << ": " << (isPrimary ? "Primary" : "Secondary") << 
              " mapping of the "<< (isStart ? "start" : "stop") << " position.\n";
            if (refClass != id_index.end()) {
              Rcpp::Rcout << refClass->class_id << " is the class of the referenced element.\n";
              Rcpp::Rcout << refClass->type << " is the type of the referenced element.\n";
            }
          }
          
          if (refClass == id_index.end()) {
            if (verbose) Rcpp::Rcout << "Invalid reference element: " << refClassId << " not found in id_index.\n";
            isValidElement = false;
          } else {
            if (refClass->type == "static") {
              isStatic = true;
              if (refClass->element_pass) {
                hasElementPass = true;
                if (*refClass->element_pass) {
                  elementPassIsTrue = true;
                }
              }
            }
            bool isStaticAndPass = isStatic && hasElementPass && elementPassIsTrue;
            if (!(isStaticAndPass || 
              refClass->type == "variable" || 
              refClass->global_class == "poly_tail" ||
              refClass->global_class == "start" ||
              refClass->global_class == "stop"
              )) {
              if (verbose) {
                if (!isStatic) 
                  Rcpp::Rcout << "Reference element: " << refClassId << " is not static.\n";
                if (!hasElementPass) 
                  Rcpp::Rcout << "Reference element: " << refClassId << " does not have element_pass.\n";
                if (!elementPassIsTrue) 
                  Rcpp::Rcout << "Reference element: " << refClassId << " has element_pass, but it is false.\n";
                if (refClass->type != "variable") 
                  Rcpp::Rcout << "Reference element: " << refClassId << " is not a variable element.\n";
                if (refClass->global_class != "poly_tail") 
                  Rcpp::Rcout << "Reference element: " << refClassId << " is not a poly_tail element.\n";
                 if (refClass->global_class != "start") 
                   Rcpp::Rcout << "Reference element: " << refClassId << " is not a start element.\n";
                 if (refClass->global_class != "stop") 
                  Rcpp::Rcout << "Reference element: " << refClassId << " is not a stop element.\n";
              }
              isValidElement = false;
            }
          }
          if (!isValidElement) {
            if (verbose) {
              Rcpp::Rcout << 
                "Invalid reference element: " << 
                refClassId << 
                " did not meet any valid criteria." <<
                " Now checking additional information:\n";
              if(!additionalInfo.empty()){
                Rcpp::Rcout << additionalInfo << "\n";
              }
            }
            if (additionalInfo.empty()||isPrimary) {
              return -1;
            }
          }
          
          int basePos = (operation.find("start") != std::string::npos) ? refClass->position.first : refClass->position.second;
          int calculatedPos = basePos + offset;
          
          if (verbose) {
            Rcpp::Rcout << " Putative base position: " << 
              basePos << 
              ", putative calculated position: " << 
              calculatedPos << "\n";
          }
          
          if (refClass->position.first < 0 && refClass->position.second < 0 && isPrimary) {
            if (verbose) Rcpp::Rcout << "Unannotated element, returning -1\n";
            return -1;
          }
          
          if (!isPrimary && !additionalInfo.empty()) {
    
            std::string readType = readData.string_info.at("type").str_type;
            if (verbose) {
              Rcpp::Rcout << "Processing additional info: " << additionalInfo << "\n"
                          << "Read type: " << readType << "\n"
                          << "Variable class: " << varClass->global_class << "\n";
            }
            if (
                (additionalInfo == "left_terminal_linked" || 
                 additionalInfo == "right_terminal_linked" || 
                 additionalInfo == "skipped")
              && readType != "FR_RF"
              && !isValidElement
            ) {
              auto startClassIt = (readType == "F") ? id_index.find("seq_start") : id_index.find("rc_seq_start");
              auto stopClassIt = (readType == "F") ? id_index.find("seq_stop") : id_index.find("rc_seq_stop");

              bool left_terminal = (additionalInfo == "left_terminal_linked");
              bool right_terminal = (additionalInfo == "right_terminal_linked");
              
              int start_order = startClassIt->order;
              int stop_order = stopClassIt->order;
              
              if (verbose) {
                Rcpp::Rcout << "Secondary position is invalid, and we're " <<
                  "trying to fix a secondary position.\nAttempting to overlook the crappy adapter.\n";
                if ((startClassIt == id_index.end() || stopClassIt == id_index.end())) {
                  Rcpp::Rcout << "Start or stop index position not found.\n";
                } else {
                  Rcpp::Rcout << "Start index position: " << start_order << "\n";
                  Rcpp::Rcout << "Stop index position: " << stop_order << "\n";
                }
              }
              bool allStatic = true;
              if (left_terminal||(additionalInfo == "skipped" && readType == "F")){
                if(allStatic && varClass->global_class == "read"){
                  if(isStart){
                  if(verbose){
                    Rcpp::Rcout << "This is a read that truncates at the start, it's left terminal-linked!\n";
                    Rcpp::Rcout << "Returning the start of string, position 1.\n";
                  }
                  return 1;
                  } else {
                    if(verbose){
                      Rcpp::Rcout << "This is a read that truncates at the end, it's right terminal-linked!\n";
                      Rcpp::Rcout << "Returning the end of string, position " <<  std::to_string(readData.string_info.at("type").str_length) << ".\n";
                    }
                    return readData.string_info.at("type").str_length;
                    
                  }
                }
                if (varClass != id_index.end() && startClassIt != id_index.end()) {
                  int varOrder = varClass->order;
                  int startOrder = startClassIt->order;
                  if (varOrder < startOrder) {
                    for (int order = varOrder + 1; order < startOrder; ++order) {
                      auto it = order_index.find(order);
                      if (it == order_index.end() || it->type != "static") {
                        allStatic = false;
                        break;
                      }
                    }
                  } else {
                    for (int order = startOrder + 1; order < varOrder; ++order) {
                      auto it = order_index.find(order);
                      if (it == order_index.end() || it->type != "static") {
                        allStatic = false;
                        break;
                      }
                    }
                  }
                }
                if (startClassIt != id_index.end() && std::abs(varClass->order - startClassIt->order) == 2){
                  auto nearestStaticIt = order_index.find(varClass->order - 1);
                  
                  if (nearestStaticIt != order_index.end() && 
                      nearestStaticIt->type == "static" && 
                      nearestStaticIt->element_pass && 
                      !*nearestStaticIt->element_pass) {
                      
                      if (verbose) {
                        Rcpp::Rcout << "This is a left-linked" << (isStart ? " start ":" stop ") << "position.\n";
                        Rcpp::Rcout << "Checking to see if the element to the left is in the direction of the read.\n";
                      }
                      
                      //built in here is the assumption that there's one adapter that's good enough to call it
                      if((nearestStaticIt->direction == "forward" && readType == "F") ||
                         (nearestStaticIt->direction == "reverse" && readType == "R")){
                        if (isStart) {
                          if(verbose){
                            Rcpp::Rcout << "Processed and returned the start position!\n";
                          }
                          return nearestStaticIt->position.second + 1;
                        } else {
                          if(verbose){
                            Rcpp::Rcout << "Processed and returned the stop position!\n";
                          }
                          return calculatedPos;
                        }
                      }
                  }
                }
              } else if(right_terminal||(additionalInfo == "skipped" && readType == "R")){
                if(allStatic && varClass->global_class == "read"){
                  if(verbose){
                    Rcpp::Rcout << "This is a read that truncates at the end, it's right terminal-linked!\n";
                    Rcpp::Rcout << "Returning the end of string, position " <<  std::to_string(readData.string_info.at("type").str_length) << ".\n";
                  }
                  return readData.string_info.at("type").str_length;
                }
                if (varClass != id_index.end() && startClassIt != id_index.end()) {
                  int varOrder = varClass->order;
                  int stopOrder = stopClassIt->order;
                  if (varOrder < stopOrder) {
                    for (int order = varOrder + 1; order < stopOrder; ++order) {
                      auto it = order_index.find(order);
                      if (it == order_index.end() || it->type != "static") {
                        allStatic = false;
                        break;
                      }
                    }
                  } else {
                    for (int order = stopOrder + 1; order < varOrder; ++order) {
                      auto it = order_index.find(order);
                      if (it == order_index.end() || it->type != "static") {
                        allStatic = false;
                        break;
                      }
                    }
                  }
                }
                if (stopClassIt != id_index.end() && std::abs(varClass->order - stopClassIt->order) == 2){
                  auto nearestStaticIt = order_index.find(varClass->order + 1);
                  
                  if (nearestStaticIt != order_index.end() && 
                      nearestStaticIt->type == "static" && 
                      nearestStaticIt->element_pass && 
                      !*nearestStaticIt->element_pass) {
                      
                      if (verbose) {
                        Rcpp::Rcout << "This is a right-linked" << (isStart ? " start ":" stop ") << "position.\n";
                        Rcpp::Rcout << "Checking to see if the element to the left is in the direction of the read.\n";
                      }
                      
                      if((nearestStaticIt->direction == "forward" && readType == "F") ||
                         (nearestStaticIt->direction == "reverse" && readType == "R")){
                        if(isStart) {
                          if(verbose){
                            Rcpp::Rcout << "Processed and returned the start position!\n";
                          }
                          return calculatedPos;
                        } else {
                          if(verbose){
                            Rcpp::Rcout << "Processed and returned the stop position!\n";
                          }
                          return nearestStaticIt->position.first - 1;
                        }
                      }
                  }
                }
              }
            }
          }
          return calculatedPos;
        };
      };
      
      funcMap[element.class_id] = {
        parsePositionData(element.class_id, primaryStart, true, true),
        parsePositionData(element.class_id, secondaryStart, false, true),
        parsePositionData(element.class_id, primaryStop, true, false),
        parsePositionData(element.class_id, secondaryStop, false, false)
      };
    }
  }
  if (verbose) {
    Rcpp::Rcout << "Position function map created with " << funcMap.size() << " entries.\n";
  }
  return funcMap;
}

void solve_single_sigstring(ReadData& readData, const PositionFuncMap& positionFuncMap, 
                            const ReadLayout& readLayout, bool verbose, bool isSecondary = false) {
  auto& sigstring = isSecondary ? *readData.secondary_sigstring : readData.sigstring;
  auto& sig_id_index = sigstring.get<sig_id_tag>();
  
  std::string readType = readData.string_info["type"].str_type;
  
  // Filter sigstring based on direction
  for (auto it = sigstring.begin(); it != sigstring.end(); ) {
    if ((readType == "F" && it->direction != "forward") || 
        (readType == "R" && it->direction != "reverse") ||
        (it->type == "static" && it->position.first == -1 && it->position.second == -1)) {
      //this removes elements that are R in an F, F in an R, and any leftover static elements that are unannotated. 
      it = sigstring.erase(it); // Erase and move to the next element
    } else {
      ++it; // Move to the next element
    }
  }
  
  // Process non-read variables
  for (auto& element : sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
      if (verbose) {
        Rcpp::Rcout << "Dealing with all non-read variables!\n";
        Rcpp::Rcout << "Working on " << element.class_id << "\n";
      }
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if (verbose) {
        Rcpp::Rcout << newStartPos << "\n";
      }
      if (newStartPos <= 0) {
        if (verbose) {
          Rcpp::Rcout << "Triggered the secondary function for start position:\n";
        }
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
      if (verbose) {
        Rcpp::Rcout << newStopPos << "\n";
      }
      if (newStopPos <= 0) {
        if (verbose) {
          Rcpp::Rcout << "Triggered the secondary function for stop position:\n";
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
    }
  }
  
  // Process read variables
  for (auto& element : sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class == "read") {
      if (verbose) {
        Rcpp::Rcout << "Dealing with everything that IS a read!\n";
      }
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if (newStartPos <= 0) {
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
      if (newStopPos <= 0) {
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
  if (verbose) {
    displayOrderedSigString(readData);
  }
}

void concat_solve_hybrid(ReadData& readData, const ReadLayout& readLayout, 
  const PositionFuncMap& positionFuncMap, bool verbose) {
  if(verbose){
    Rcpp::Rcout << "Moving into concatenate resolution!\n";
  }
  auto processSigElement = [&](const SigElement& element, const std::string& direction) {
    if (element.type != "variable" || positionFuncMap.count(element.class_id) == 0) {
      return;
    }
    const auto& positionFuncs = positionFuncMap.at(element.class_id);
    
    int newStartPos = positionFuncs.primaryStartFunc(readData);
    if (newStartPos == -1){
      newStartPos = positionFuncs.secondaryStartFunc(readData);
    }
    
    int newStopPos = positionFuncs.primaryStopFunc(readData);
    if (newStopPos == -1){
      newStopPos = positionFuncs.secondaryStopFunc(readData);
    }
    
    auto findAdjacentElement = [&](
      int currentOrder, 
      bool findStart, 
      int maxIterations = 3)->std::pair<boost::optional<SigElement>, int> {
        auto& order_index = readData.sigstring.get<sig_order_tag>();
        int step = (direction == "forward") ? 1 : -1;
        if (!findStart) step *= 1;
        int suitableElementsFound = 0;
        while (true) {
          currentOrder += step;
          auto adjacentElement = order_index.find(currentOrder);
          if (adjacentElement == order_index.end()) {
            break;
          }
          if(verbose){
            Rcpp::Rcout << "Looking at adjacent element: " << adjacentElement->class_id 
                        << " for " << (findStart ? "start" : "stop") << " position.\n";
          }
          bool isUnsuitable = 
            (adjacentElement->global_class == "poly_tail" &&
            (adjacentElement->position.second - adjacentElement->position.first <= 12)) ||
            //why do i have this like this, why am i a dumbass
            adjacentElement->global_class == "start" ||
            adjacentElement->global_class == "stop" ||
            (adjacentElement->type == "static" && adjacentElement->global_class != "poly_tail" &&
            ((adjacentElement->position.second - adjacentElement->position.first <= 12) || 
             (adjacentElement->element_pass && !*adjacentElement->element_pass)));
          
          if (!isUnsuitable) {
            suitableElementsFound++;
            if (suitableElementsFound > maxIterations) {
              break;
            }
            if (verbose) {
              Rcpp::Rcout << "Found suitable element: " << adjacentElement->class_id << "\n";
            }
            return {*adjacentElement, currentOrder};
          }
          if (verbose && isUnsuitable) {
            Rcpp::Rcout << "Skipping unsuitable element: " << adjacentElement->class_id << "\n";
          }
        }
        return {boost::none, currentOrder};
      };
    
    if(verbose){
      Rcpp::Rcout << "Current start position is " << newStartPos << " & current stop position is " << newStopPos << ".\n";
    }
    
    if (newStartPos <= 0) {
      auto [adjacentElement, newOrder] = findAdjacentElement(element.order, true);
      if (adjacentElement) {
        SigElement tempElement = *adjacentElement;
        tempElement.element_pass = true;
        auto& order_index = readData.sigstring.get<sig_order_tag>();
        order_index.replace(order_index.find(newOrder), tempElement);
        newStartPos = tempElement.position.second + 1;
        order_index.replace(order_index.find(newOrder), *adjacentElement);
      }
    }
    if (newStopPos <= 0) {
      auto [adjacentElement, newOrder] = findAdjacentElement(element.order, false);
      if (adjacentElement) {
        SigElement tempElement = *adjacentElement;
        tempElement.element_pass = true;
        auto& order_index = readData.sigstring.get<sig_order_tag>();
        order_index.replace(order_index.find(newOrder), tempElement);
        newStopPos = tempElement.position.first - 1;
        order_index.replace(order_index.find(newOrder), *adjacentElement);
      }
    }
        
        //added to basically prevent against whack stuff happening during concatenate resolution, where
    //the element might be in the class order of FR_RF, but isn't actually *there*.

    if(newStartPos >= newStopPos){
      auto [adjacentElement, newOrder] = findAdjacentElement(element.order, false);
      if (adjacentElement) {
        SigElement tempElement = *adjacentElement;
        tempElement.element_pass = true;
        auto& order_index = readData.sigstring.get<sig_order_tag>();
        order_index.replace(order_index.find(newOrder), tempElement);
        newStopPos = tempElement.position.first - 1;
        order_index.replace(order_index.find(newOrder), *adjacentElement);
      }
    }
    if (newStartPos > 0 && newStopPos > 0) {
      auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
      sig_id_index.modify(sig_id_index.find(element.class_id), 
        [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
      if (verbose) {
        Rcpp::Rcout << "Updated positions for " << element.class_id 
                    << ": Start: " << newStartPos << ", Stop: " << newStopPos << "\n";
      }
    } else if (verbose) {
      Rcpp::Rcout << "Failed to find suitable positions for " << element.class_id << "\n";
    }
  };
  
  // Process non-read elements
  for (const auto& element : readData.sigstring) {
    if (element.global_class != "read") {
      processSigElement(element, element.direction);
    }
  }
  
  // Process read elements
  for (const auto& element : readData.sigstring) {
    if (element.global_class == "read") {
      processSigElement(element, element.direction);
    }
  }
  
  // Generate and process new sigstrings
  std::vector<std::string> newSigStrings = generateSigString(readData);
  std::vector<std::string> deconcatenated = processSigString(newSigStrings[0]);
  
  if (verbose) {
    Rcpp::Rcout << "Original concatenated solution:\n" << newSigStrings[0] << "\n";
    Rcpp::Rcout << "Deconcatenated versions:\n";
    for (const auto& part : deconcatenated) {
      Rcpp::Rcout << part << "\n";
    }
  }
  
  // Clear existing data and fill with new sigstrings
  readData.sigstring.clear();
  readData.secondary_sigstring = boost::none;
  
  if (!deconcatenated.empty()) {
    fillSigString(readData, container, deconcatenated[0], positionFuncMap, verbose);
    if (deconcatenated.size() > 1) {
      fillSigString(readData, container, deconcatenated[1], positionFuncMap, verbose);
    }
  }
  
  readData.string_info["type"].str_type = "F-R";
  readData.string_info["type"].split_types = std::make_pair("F", "R");
}

void concat_solve_parallel(ReadData& readData, const ReadLayout& readLayout, 
  const PositionFuncMap& positionFuncMap, bool verbose) {
  std::string readType = readData.string_info["type"].str_type;
  std::string read_id = readData.string_info["type"].str_id;
  if (verbose) {
    std::cout << "Resolving parallel forward or reverse concatenation.\n";
    std::cout << read_id << " is the read ID.\n";
  }
  
  // Helper lambdas
  auto filterAndSortTerminalElements = [&](const SigString& sigstring) {
    std::vector<SigElement> terminal_elements;
    for (const auto& elem : sigstring) {
      if ((elem.global_class == "forw_primer" || elem.global_class == "rev_primer") 
        && elem.direction == "forward") {
        terminal_elements.push_back(elem);
      }
    }
    std::sort(terminal_elements.begin(), terminal_elements.end(),
      [](const SigElement& a, const SigElement& b) {
        return a.position.first < b.position.first;
      });
    return terminal_elements;
  };
  
  auto findSplitPosition = [&](const std::vector<SigElement>& elements) {
    int split_position = -1;
    bool found_first_pass = false;
    for (const auto& elem : elements) {
      if (elem.element_pass && *elem.element_pass) {
        if (!found_first_pass) {
          found_first_pass = true;
        } else {
          split_position = elem.position.first - 1; // Split before the second passing primer
          break;
        }
      }
    }
    return split_position;
  };
  
  auto distributeElements = [&](const SigString& sigstring, int split_position, bool verbose = false) {
    std::vector<SigElement> first_read_elements, second_read_elements;
    std::set<std::string> unique_class_ids;
    for (const auto& elem : sigstring) {
      if(elem.type == "variable"){
        continue;
      }
      unique_class_ids.insert(elem.class_id);
    }
    for (const auto& class_id : unique_class_ids) {
      // Gather all SigElements with the current class_id
      std::vector<SigElement> matching_elements;
      for (const auto& elem : sigstring) {
        if (elem.class_id == class_id) {
          matching_elements.push_back(elem);
        }
      }
      // Determine whether to assign elements to first or second read
      for (const auto& elem : matching_elements) {
        if (elem.position.first <= split_position) {
          if(verbose) std::cout << "Element with class_id " << elem.class_id << " sorted into the first read.\n";
          first_read_elements.push_back(elem);
        } else {
          if(verbose) std::cout << "Element with class_id " << elem.class_id << " sorted into the second read.\n";
          second_read_elements.push_back(elem);
        }
      }
    }
    
    // Sort the elements in each read by position
    std::sort(first_read_elements.begin(), first_read_elements.end(),
              [](const SigElement& a, const SigElement& b) {
                return a.position.first < b.position.first;
              });
    std::sort(second_read_elements.begin(), second_read_elements.end(),
              [](const SigElement& a, const SigElement& b) {
                return a.position.first < b.position.first;
              });
    return std::make_pair(first_read_elements, second_read_elements);
  };

  auto generate_parallel_sigstrings = [&](const std::vector<SigElement>& elements, int seq_stop_position, 
    int start_offset, const std::string& type, const std::string& str_id,
    const PositionFuncMap& positionFuncMap, bool verbose = false) {
    std::stringstream ss;
    ss << "seq_start:0:1:1";
    for (const auto& elem : elements) {
      if (elem.global_class != "start" && elem.global_class != "stop") {
        if(verbose){
          std::cout << "Element class ID is " << elem.class_id << ", element global class is " << elem.global_class <<
            ", start position is " <<
              elem.position.first << ", and stop position is " << elem.position.second << ".\n";
        }
        int adjusted_start = elem.position.first - start_offset + 1;
        int adjusted_end = elem.position.second - start_offset + 1;
        ss << "|" << elem.class_id << ":" 
           << elem.edit_distance.value_or(0) << ":" 
           << adjusted_start << ":" 
           << adjusted_end;
      }
    }
    // Add the seq_stop position after processing all elements
    int adjusted_stop = seq_stop_position - start_offset + 1;
    ss << "|seq_stop:0:" << adjusted_stop << ":" << adjusted_stop;
    ss << "<" << adjusted_stop << ":" << str_id << "+" << type << ":F>";
    return ss.str();
  };
  
  // Filter and sort terminal elements
  std::vector<SigElement> terminal_elements = filterAndSortTerminalElements(readData.sigstring);
  
  // Find the split position
  int split_position = findSplitPosition(terminal_elements);
  if (split_position == -1) {
    if (verbose) {
      std::cout << "No valid split point found.\n";
    }
    return;
  }
  
  // Distribute elements into two groups
  auto [first_read_elements, second_read_elements] = distributeElements(readData.sigstring, split_position);
  
  // Get `seq_stop` position for primary and secondary sigstrings
  auto seq_stop_primary = split_position;
  auto seq_stop_secondary = readData.sigstring.get<sig_id_tag>().find("seq_stop")->position.second;
  
  // Generate new sigstrings
  std::string first_sigstring = generate_parallel_sigstrings(first_read_elements, seq_stop_primary, 0, readType, read_id, positionFuncMap, verbose);
  std::string second_sigstring = generate_parallel_sigstrings(second_read_elements, seq_stop_secondary, split_position + 1, readType, read_id, positionFuncMap, verbose);
  
  if (verbose) {
    std::cout << "First Signature: " << first_sigstring << "\n";
    std::cout << "Second Signature: " << second_sigstring << "\n";
  }
  
  // Update ReadData with new sigstrings
  readData.sigstring.clear();
  fillSigString(readData, readLayout, first_sigstring, positionFuncMap, verbose);
  solve_single_sigstring(readData, positionFuncMap, readLayout, verbose);
  if (!readData.secondary_sigstring) {
    readData.secondary_sigstring = SigString();
  }
  readData.secondary_sigstring->clear();
  fillSigString(readData, readLayout, second_sigstring, positionFuncMap, verbose);
  solve_single_sigstring(readData, positionFuncMap, readLayout, verbose, true);
  // Update type information
  readData.string_info["type"].str_type = "F-F";
  readData.string_info["type"].split_types = std::make_pair("F", "F");
  
  if (verbose) {
    std::cout << "Completed parallel forward concatenation resolution.\n";
  }
}

//parsing a sigstring and adding stuff to it
void fillSigString(ReadData& readData, const ReadLayout& readLayout, 
      const std::string& signature, const PositionFuncMap& positionFuncMap, bool verbose) {
  if (signature.empty()) {
    throw std::runtime_error("sigstring is empty");
  }
  auto read_info_pos = signature.find_last_of('<');
  std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2); // Remove '<' and '>'
  std::stringstream lengthTypeStream(lengthTypeStr);
  std::string lengthStr, id, type;
  std::getline(lengthTypeStream, lengthStr, ':');
  std::getline(lengthTypeStream, id, ':');
  std::getline(lengthTypeStream, type);
  
  int length = std::stoi(lengthStr);
  if(verbose) {
    std::cout << "Length is " << length << ", ID is " << id << ", and type is " << type << ".\n";
  }  
  // Populate StringInfo with length and type
  readData.addStringInfo("type", length, id, type);
  SigString* targetSigString = &readData.sigstring; // Default to primary sigstring
  if (!readData.sigstring.empty() && !readData.secondary_sigstring) {
    // Condition 2: sigstring is populated and secondary_sigstring is not, populate secondary_sigstring.
    readData.secondary_sigstring = SigString(); // Initialize secondary_sigstring
    targetSigString = &(*readData.secondary_sigstring); // Target secondary sigstring
  } else if (!readData.sigstring.empty() && readData.secondary_sigstring) {
    // Condition 3: Both sigstring and secondary_sigstring are populated, empty sigstring and refill it with the new sigstring data.
    readData.sigstring.clear(); // Empty primary sigstring
    targetSigString = &readData.sigstring; // Target primary sigstring
  }
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
    if(verbose){
      std::cout << id << " is the id.\n";
      std::string exists = (it != class_id_index.end()) ? "Exists in the class id index.\n" : "Does not exist in the class id index.\n";
      std::cout << exists;
    }
    std::string readType = readData.string_info["type"].str_type;
    if (it != class_id_index.end()) {
      // Calculate misalignment threshold if available
      boost::optional<bool> element_pass;
      if (it->misalignment_threshold && it->type == "static" && (endPos-(startPos-1) >= 12)) {
        int misalignment_threshold = std::get<1>(*it->misalignment_threshold);
        if(readType != "FR_RF" && misalignment_threshold > 4){
          element_pass = (editDistance < misalignment_threshold);
          if(verbose){
            std::string pass = (element_pass)? "This element passed the edit distance being less than the threshold floor:\n" : "didn't pass threshold floor:\n";
            std::cout << pass << "\n";
            std::cout << misalignment_threshold << "\n";
          }
        } else {
          if(readType == "FR_RF" || misalignment_threshold <= 4){
            element_pass = (editDistance <= misalignment_threshold);
            if(verbose){
              std::string pass = (element_pass)? "This element passed the edit distance being less than the misalignment_threshold\n" : "didn't pass threshold floor\n";
              std::cout << pass << "\n";
              std::cout << misalignment_threshold << "\n";
            }
          }
        }
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
      targetSigString->insert(std::move(element));
    }
  }
  for (const auto& readElement : readLayout) {
    if (sig_id_index.find(readElement.class_id) == sig_id_index.end()) {
      SigElement element(
          readElement.class_id,
          readElement.global_class,
          0, // default edit distance
          std::make_pair(-1, -1), // default positions
          readElement.type,
          readElement.order,
          readElement.direction
      );
      targetSigString->insert(std::move(element));
    }
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame sig_extraction_updated(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold,
  Rcpp::DataFrame& df, const Rcpp::CharacterVector& processed_sigstrings, bool verbose) {
  // Initial setup
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  VarScan(readLayout, verbose);
  positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  
  // Convert processed_sigstrings to std::vector<std::string>
  std::vector<std::string> sigstrings = Rcpp::as<std::vector<std::string>>(processed_sigstrings);
  
  // Extract all unique class IDs
  std::set<std::string> uniqueClassIds;
  auto& class_id_index = readLayout.get<id_tag>();
  for (const auto& elem : class_id_index) {
    if (elem.type == "variable") {
      if (elem.class_id.substr(0, 3) == "rc_") {
        uniqueClassIds.insert(elem.class_id.substr(3));
      } else {
        uniqueClassIds.insert(elem.class_id);
      }
    }
  }
  
  // Modified sequence map to handle multiple sequences per ID
  std::unordered_map<std::string, std::vector<std::string>> seq_map;
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector seqs = df["seq"];
  for (int i = 0; i < ids.size(); ++i) {
    seq_map[Rcpp::as<std::string>(ids[i])].push_back(Rcpp::as<std::string>(seqs[i]));
  }
  
  // Store results with original indices to maintain order
  std::vector<std::string> result_ids;
  std::vector<int> sig_indices;
  std::vector<std::map<std::string, std::string>> result_values;
  
  int valid_sig_index = 0;
  for (size_t i = 0; i < sigstrings.size(); ++i) {  // Changed to size_t
    std::string signature = sigstrings[i];
    auto read_info_pos = signature.find_last_of('<');
    std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    
    if(verbose){
      Rcpp::Rcout << "Signature: " << signature << std::endl;
      Rcpp::Rcout << "Length: " << lengthStr << ", Read ID: " << read_id << ", Type: " << type << std::endl;
    }
    
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    
    if(type == "undecided"){
      continue;
    }
    
    result_ids.push_back(read_id);
    sig_indices.push_back(++valid_sig_index);
    
    std::map<std::string, std::string> current_values;
    std::stringstream ss(signature.substr(0, read_info_pos));
    std::string token;
    
    while (std::getline(ss, token, '|')) {
      std::string id;
      int editDistance, startPos, endPos;
      std::stringstream tokenStream(token);
      
      std::getline(tokenStream, id, ':');
      tokenStream >> editDistance;
      tokenStream.ignore(1);
      tokenStream >> startPos;
      tokenStream.ignore(1);
      tokenStream >> endPos;
      
      auto it = class_id_index.find(id);
      if (it != class_id_index.end() && it->type == "variable") {
        auto seq_it = seq_map.find(read_id);
        if (seq_it == seq_map.end()) {
          continue;
        }
        for (const auto& seq : seq_it->second) {
          if (startPos > 0 && endPos >= startPos && (seq.length() >= endPos)) {
            std::string extracted_seq = seq.substr(startPos - 1, endPos - startPos + 1);
            if (it->direction != "forward") {
              extracted_seq = revcomp_cpp(extracted_seq);
            }
            
            std::string classId = id;
            if (classId.substr(0, 3) == "rc_") {
              classId = classId.substr(3);
            }
            
            current_values[classId] = extracted_seq;
          }
        }
      }
    }
    result_values.push_back(current_values);
  }
  
  Rcpp::List data_table_list;
  data_table_list["id"] = result_ids;
  data_table_list["sig_index"] = sig_indices;
  
  for (const auto& classId : uniqueClassIds) {
    Rcpp::CharacterVector column(result_ids.size(), NA_STRING);
    for (size_t i = 0; i < result_values.size(); ++i) {
      const auto& current_values = result_values[i];
      auto it = current_values.find(classId);
      if (it != current_values.end()) {
        column[i] = it->second;
      }
    }
    data_table_list[classId] = column;
  }
  
  data_table_list.attr("class") = "data.table";
  data_table_list.attr("row.names") = Rcpp::seq(1, result_ids.size());
  return data_table_list;
}

std::vector<std::string> generateSigString(ReadData& readData) {  // Note: passing readData by non-const reference
  std::vector<std::string> sigStrings; // Vector to hold one or two sigstrings
  
  auto& typeInfo = readData.string_info["type"];  // Getting a modifiable reference to typeInfo
  auto generateFromSigString = [&](const SigString& sigstring, const std::string& typeOverride) -> std::string {
    std::stringstream sigStringStream;
    auto appendElements = [&](const std::string& direction) {
      std::vector<SigElement> sorted_elements;
      auto& dir_index = sigstring.get<sig_dir_tag>();
      auto dir_range = dir_index.equal_range(direction);
      for (auto it = dir_range.first; it != dir_range.second; ++it) {
        if (it->position.first <= 0 || it->position.second <= 0) continue;
        sorted_elements.push_back(*it);
      }
      std::sort(sorted_elements.begin(), sorted_elements.end(),
        [](const SigElement& a, const SigElement& b) { return a.order < b.order; });
      for (size_t i = 0; i < sorted_elements.size(); ++i) {
        const auto& element = sorted_elements[i];
        if (i > 0) sigStringStream << "|";
        sigStringStream << element.class_id << ":" << element.edit_distance.value_or(0)
                        << ":" << element.position.first << ":" << element.position.second;
      }
    };
    appendElements("forward");
    if (typeOverride == "R" || typeOverride == "FR_RF") {
      if (!sigStringStream.str().empty()) sigStringStream << "+";
      appendElements("reverse");
    }
    sigStringStream << "<" << typeInfo.str_length << ":" << typeInfo.str_id << ":" << typeOverride << ">";
    return sigStringStream.str();
  };
  // Logic to handle F-R type
  if(typeInfo.str_type == "F-R" || typeInfo.str_type == "F-F") {
    sigStrings.push_back(generateFromSigString(readData.sigstring, "F"));
    if(readData.secondary_sigstring) {
      if(typeInfo.str_type == "F-R"){
        sigStrings.push_back(generateFromSigString(*readData.secondary_sigstring, "R"));
      } else if (typeInfo.str_type == "F-F"){
        sigStrings.push_back(generateFromSigString(*readData.secondary_sigstring, "F"));
      }
    }
  } else {
    // Default behavior for other types
    sigStrings.push_back(generateFromSigString(readData.sigstring, typeInfo.str_type));
    if(readData.secondary_sigstring) {
      sigStrings.push_back(generateFromSigString(*readData.secondary_sigstring, typeInfo.str_type));
    }
  }
  return sigStrings;
}

void pos_scan(ReadData& readData, const PositionFuncMap& positionFuncMap, const ReadLayout& readLayout, 
  bool verbose) {
  std::string readType = readData.string_info["type"].str_type;
  if(readType == "F" || readType == "R"){
    solve_single_sigstring(readData, positionFuncMap, readLayout, verbose);
  } else if(readType == "FR_RF"){
    concat_solve_hybrid(readData, readLayout, positionFuncMap, verbose);
  } else if(readType == "F_F"||readType == "R_R"){
    concat_solve_parallel(readData, readLayout, positionFuncMap, verbose);
  }
}

Rcpp::CharacterVector process_sigstrings(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads = 1) {
  
  omp_set_num_threads(nthreads);
  std::vector<std::string> processed_sigstrings;
  
#pragma omp parallel
{
  std::vector<std::string> local_processed_sigstrings;
#pragma omp for nowait
  for (int i = 0; i < sigs.size(); ++i) {
    try {
      const auto& sigstring = sigs[i];
      ReadData readData;
      fillSigString(readData, readLayout, sigstring, positionFuncMap, verbose);
      concat_scan(readData, verbose);
      pos_scan(readData, positionFuncMap, readLayout, verbose);
      filter_sigstrings(readData, positionFuncMap, readLayout, verbose);
      std::vector<std::string> generated_sigstrings = generateSigString(readData);
      //filter_sigstring(generated_sigstrings);
      local_processed_sigstrings.insert(local_processed_sigstrings.end(), 
       generated_sigstrings.begin(), 
        generated_sigstrings.end());
    } catch (const std::exception& e) {
#pragma omp critical
{
  Rcpp::Rcerr << "Error processing sigstring " << i << ": " << e.what() << std::endl;
}
    } catch (...) {
#pragma omp critical
{
  Rcpp::Rcerr << "Unknown error processing sigstring " << i << std::endl;
}
    }
  }
#pragma omp critical
  processed_sigstrings.insert(processed_sigstrings.end(), local_processed_sigstrings.begin(), local_processed_sigstrings.end());
}
return wrap(processed_sigstrings);
}

std::vector<std::string> process_sigstrings_cpp(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads = 1) {
  omp_set_num_threads(nthreads);
  std::vector<std::string> processed_sigstrings;
#pragma omp parallel
{
  std::vector<std::string> local_processed_sigstrings;
#pragma omp for nowait
  for (int i = 0; i < sigs.size(); ++i) {
    const auto& sigstring = sigs[i];
    ReadData readData;
    fillSigString(readData, readLayout, sigstring, positionFuncMap, verbose);
    concat_scan(readData, verbose);
    pos_scan(readData, positionFuncMap, readLayout, verbose);
    filter_sigstrings(readData, positionFuncMap, readLayout, verbose);
    std::vector<std::string> generated_sigstrings = generateSigString(readData);
    local_processed_sigstrings.insert(local_processed_sigstrings.end(), 
      generated_sigstrings.begin(), 
      generated_sigstrings.end());
  }
#pragma omp critical
  processed_sigstrings.insert(processed_sigstrings.end(), local_processed_sigstrings.begin(), local_processed_sigstrings.end());
}
return processed_sigstrings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigrun(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  const Rcpp::StringVector& sigstrings, int nthreads = 1, bool verbose = false) {
  
  Rcpp::Rcout << "Starting sigrun..." << std::endl;
  
  try {
    Rcpp::Rcout << "Calling prep_read_layout_cpp..." << std::endl;
    container = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
    Rcpp::Rcout << "Completed prep_read_layout_cpp" << std::endl;
  } catch (const std::exception& e) {
    Rcpp::Rcerr << "Error in prep_read_layout_cpp: " << e.what() << std::endl;
    return Rcpp::CharacterVector();
  } catch (...) {
    Rcpp::Rcerr << "Unknown error in prep_read_layout_cpp" << std::endl;
    return Rcpp::CharacterVector();
  }
  
  try {
    Rcpp::Rcout << "Calling VarScan..." << std::endl;
    VarScan(container, verbose);
    Rcpp::Rcout << "Completed VarScan" << std::endl;
  } catch (const std::exception& e) {
    Rcpp::Rcerr << "Error in VarScan: " << e.what() << std::endl;
    return Rcpp::CharacterVector();
  } catch (...) {
    Rcpp::Rcerr << "Unknown error in VarScan" << std::endl;
    return Rcpp::CharacterVector();
  }
  
  PositionFuncMap positionFuncMap;
  try {
    Rcpp::Rcout << "Calling createPositionFunctionMap..." << std::endl;
    positionFuncMap = createPositionFunctionMap(container, verbose);
    Rcpp::Rcout << "Completed createPositionFunctionMap" << std::endl;
  } catch (const std::exception& e) {
    Rcpp::Rcerr << "Error in createPositionFunctionMap: " << e.what() << std::endl;
    return Rcpp::CharacterVector();
  } catch (...) {
    Rcpp::Rcerr << "Unknown error in createPositionFunctionMap" << std::endl;
    return Rcpp::CharacterVector();
  }
  
  std::vector<std::string> sigs;
  try {
    Rcpp::Rcout << "Converting sigstrings..." << std::endl;
    sigs = Rcpp::as<std::vector<std::string>>(sigstrings);
    Rcpp::Rcout << "Completed converting sigstrings" << std::endl;
  } catch (const std::exception& e) {
    Rcpp::Rcerr << "Error converting sigstrings: " << e.what() << std::endl;
    return Rcpp::CharacterVector();
  } catch (...) {
    Rcpp::Rcerr << "Unknown error converting sigstrings" << std::endl;
    return Rcpp::CharacterVector();
  }
  
  try {
    Rcpp::Rcout << "Calling process_sigstrings..." << std::endl;
    auto result = process_sigstrings(container, sigs, positionFuncMap, verbose, nthreads);
    Rcpp::Rcout << "Completed process_sigstrings" << std::endl;
    return result;
  } catch (const std::exception& e) {
    Rcpp::Rcerr << "Error in process_sigstrings: " << e.what() << std::endl;
    return Rcpp::CharacterVector();
  } catch (...) {
    Rcpp::Rcerr << "Unknown error in process_sigstrings" << std::endl;
    return Rcpp::CharacterVector();
  }
}
