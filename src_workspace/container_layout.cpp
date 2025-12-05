// [[Rcpp::plugins(cpp17)]]
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
#include <optional>

using namespace boost::multi_index;
using namespace Rcpp;

struct ReadElement {
  std::string class_id;
  std::string type;
  int order;
  std::string direction;
  std::optional<int> expected_length;
  std::optional<std::tuple<int, int, int>> misalignment_threshold;
  ReadElement(std::string class_id, std::string type, int order, std::string direction,
    const std::optional<int>& exp_length,
    const std::optional<std::tuple<int, int, int>>& misalign_threshold)
    : class_id(std::move(class_id)),
      type(std::move(type)),
      order(order),
      direction(std::move(direction)),
      expected_length(exp_length),
      misalignment_threshold(misalign_threshold) {}
};

struct id_tag {};
struct type_tag {};
struct order_tag {};
struct direction_tag {};

typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>
  >
> ReadLayout;
ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_thresholds) {
  ReadLayout container;
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_thresholds["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_thresholds["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_thresholds["misal_sd"];
  for (int i = 0; i < query_id.size(); ++i) {
    double threshold = misal_threshold[i];
    double sd = std::ceil(misal_sd[i]);
    int lower_bound = static_cast<int>(std::floor(threshold - sd));
    int upper_bound = static_cast<int>(std::ceil(threshold + sd));
    thresholdsMap[Rcpp::as<std::string>(query_id[i])] = std::make_tuple(lower_bound, static_cast<int>(threshold), upper_bound);
  }
  Rcpp::StringVector class_id = read_layout["class_id"];
  Rcpp::StringVector type = read_layout["type"];
  Rcpp::IntegerVector expected_length = read_layout["expected_length"];
  Rcpp::IntegerVector order = read_layout["order"];
  Rcpp::StringVector direction = read_layout["direction"];
  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);
    std::string t = Rcpp::as<std::string>(type[i]);
    std::optional<int> exp_length;
    std::optional<std::tuple<int, int, int>> misalign_threshold;
    if (t != "variable" && t != "read" && cid != "poly_a" && cid != "poly_t") {
      exp_length = expected_length[i];
      misalign_threshold = thresholdsMap.count(cid) ? thresholdsMap[cid] : std::make_tuple(0, 0, 0);
    }
    container.insert(ReadElement(
        cid, t, order[i], Rcpp::as<std::string>(direction[i]), exp_length, misalign_threshold
    ));
  }
  return container;
}
