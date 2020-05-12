#pragma once
#include "EliasFano.h"
#include <Rcpp.h>
#include <vector>
#include <set>

class QueryScore
{
public:
  typedef struct 
  {
    double tfidf;
    size_t index;
    int cartesian_product_sets;
  } GeneScore;
  friend class EliasFanoDB;
  int cells_in_query;
  int cell_types_in_query;
  double query_score;
  std::map<std::string, GeneScore> genes;
  // CellID (cell type , cell number)
  std::unordered_map<CellID , std::pair<std::vector<double>, int> > tfidf;

  
  

  
  QueryScore();
  void reset();
  void cell_type_relevance(const EliasFanoDB&, const Rcpp::List&, const std::set<std::string>&);
  void cell_tfidf(const EliasFanoDB&, const std::set<std::string>&);
  void estimateExpression(const Rcpp::List& gene_results, const EliasFanoDB& db, const Rcpp::CharacterVector& datasets, bool concsole_message);
  int calculate_cell_types(const std::set<std::string>&gene_set);
};
