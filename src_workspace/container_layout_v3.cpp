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
#include <utility> // for std::pair


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
  int start;
  int stop;
};


struct PositionInfo {
  std::function<std::pair<int, int>()> position_calculator;
};

// Define a structure to hold the combined information
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

// Define the multi-index container with appropriate indices
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

// Function to process variable elements for a given segment
void process_segment(const ReadLayout& readLayout, const std::vector<int>& segment_orders, const std::string& direction, bool verbose) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");

  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
      // Logic to find nearest static elements and calculate positions...
      continue;
    }
    // Check if the variable element is within the segment order range
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    auto parse_static_left = [&](const ReadElement& varElem, const std::vector<int>& orders) -> std::string {
      for (auto it = orders.rbegin(); it != orders.rend(); ++it) {
        if (*it < varElem.order) {
          auto staticElem = readLayout.get<order_tag>().find(*it);
          if (staticElem->type == "static") {
            if (verbose) {
              Rcpp::Rcout << "Static Left Found: " << staticElem->class_id
                          << ", Order: " << staticElem->order << "\n";
            }
            return staticElem->class_id;
          }
        }
      }
      return "Start of Segment";
    };
    auto parse_static_right = [&](const ReadElement& varElem, const std::vector<int>& orders) -> std::string {
      for (auto it = orders.begin(); it != orders.end(); ++it) {
        if (*it > varElem.order) {
          auto staticElem = readLayout.get<order_tag>().find(*it);
          if (staticElem->type == "static") {
            if (verbose) {
              Rcpp::Rcout << "Static Right Found: " << staticElem->class_id
                          << ", Order: " << staticElem->order << "\n";
            }
            return staticElem->class_id;
          }
        }
      }
      return "End of Segment";
    };

    std::string position_desc;
    int read_order = (direction == "forward") ? *std::max_element(segment_orders.begin(), segment_orders.end()) :
      *std::min_element(segment_orders.begin(), segment_orders.end());

    if (it_var->order < read_order) {
      position_desc = parse_static_left(*it_var, segment_orders);
    } else {
      position_desc = parse_static_right(*it_var, segment_orders);
    }

    if (verbose) {
      Rcpp::Rcout << "Variable Element: " << it_var->class_id
                  << ", Direction: " << it_var->direction
                  << ", Position Desc: " << position_desc << "\n";
    }
  }
}

std::string getPositionEquation(const ReadElement& varElement, const ReadLayout& readLayout, const std::string& direction) {
  auto& order_index = readLayout.get<order_tag>();

  // Find the nearest static elements to the left and right
  auto nearest_static_before = std::find_if(order_index.rbegin(), order_index.rend(),
    [&](const ReadElement& elem) {
      return elem.order < varElement.order && elem.type == "static" && elem.direction == direction;
    });
  auto nearest_static_after = std::find_if(order_index.begin(), order_index.end(),
    [&](const ReadElement& elem) {
      return elem.order > varElement.order && elem.type == "static" && elem.direction == direction;
    });

  // Construct the position equation
  std::string start_equation, stop_equation;
  if (nearest_static_before != order_index.rend()) {
    start_equation = nearest_static_before->class_id + " stop + 1";
  } else {
    start_equation = "Start of Segment";
  }
  if (nearest_static_after != order_index.end()) {
    stop_equation = nearest_static_after->class_id + " start - 1";
  } else {
    stop_equation = "End of Segment";
  }

  // Adjust for expected length (if available)
  if (varElement.expected_length) {
    stop_equation += " + " + std::to_string(*varElement.expected_length - 1);
  }

  return "Start: " + start_equation + ", Stop: " + stop_equation;
}

std::unordered_map<std::string, std::string> VarScan_V2(const ReadLayout& readLayout, bool verbose) {
  std::unordered_map<std::string, std::string> positionEquations;
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");

  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    std::string equation = getPositionEquation(*it_var, readLayout, it_var->direction);
    positionEquations[it_var->class_id] = equation;
  }
  if (verbose){
    for (const auto& [var_id, equation] : positionEquations) {
      Rcpp::Rcout << "Element: " << var_id << ", Position Equation: " << equation << "\n";
    }
  }
  return positionEquations;
}

// VarScan function
void VarScan_V1(ReadLayout& readLayout, bool verbose) {
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
  if (verbose) {
    Rcpp::Rcout << "Checking Variable Elements Positional Information:\n";
    for (const auto& element : readLayout) {
      if (element.type == "variable") {
        Rcpp::Rcout << "  Element: " << element.class_id;
        if (element.position_data) {
          auto pos_data = element.position_data.get();
          Rcpp::Rcout << ", Start Position: " << pos_data.start
                      << ", Stop Position: " << pos_data.stop << "\n";
        } else {
          Rcpp::Rcout << ", Position Data: Not Assigned\n";
        }
      }
    }
  }
}

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
}

// [[Rcpp::export]]
void VarTest(DataFrame read_layout, DataFrame misalignment_threshold) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  // Modify static elements with predefined positions
  auto& id_index = container.get<id_tag>();
  auto modify_static_position = [&](const std::string& id, int start, int stop) {
    auto it = id_index.find(id);
    if (it != id_index.end()) {
      container.modify(it, [&](ReadElement& elem) {
        elem.position_data = PositionInfo{start, stop};
      });
    }
  };

  // Assign dummy positions to static elements
  modify_static_position("forw_primer", 1, 23);
  modify_static_position("rev_primer", 60, 83);
  modify_static_position("rc_rev_primer", 1, 24);
  modify_static_position("rc_forw_primer", 60, 82);

  // Run VarScan function with verbose output
  VarScan_V1(container, true);
  VarScan_V2(container, true);

}

// [[Rcpp::export]]
void test_read_layout_container(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose = true) {
  if (verbose) {
    VarTest(read_layout, misalignment_threshold);
  }
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  if (verbose) {
    test_ref_layout(container);
  }
}

