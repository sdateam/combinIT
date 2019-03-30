#ifndef utils_h
#define utils_h
#include <Rcpp.h>
using namespace Rcpp;
List mlkmeans(const arma::mat& data, const int& clusters); 
List fastLm(const arma::vec & y, const arma::mat & X);

#endif