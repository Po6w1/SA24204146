#include <RcppArmadillo.h>

//' @title Proportion of Diagonally Dominant Rows
//' @description Computes the proportion of rows in a square matrix that satisfy diagonal dominance.
//' @name diagonal_dominance_cpp
//' @param mat Numeric matrix (assumed to be square).
//' @return A double representing the proportion of rows that satisfy diagonal dominance.
//' @useDynLib SA24204146
//' @import Rcpp
//' @export
// [[Rcpp::export]]
double diagonal_dominance_cpp(const Rcpp::NumericMatrix& mat) {

  // Convert Rcpp matrix to Armadillo matrix
  arma::mat A = Rcpp::as<arma::mat>(mat);
  
  // Check if the matrix is square
  if (A.n_rows != A.n_cols) {
    Rcpp::stop("Input matrix must be square.");
  }
  
  int n = A.n_rows;
  int dominant_rows = 0; // 统计满足对角占优的行数

  for (int i = 0; i < n; ++i) {
    double diag_value = std::abs(A(i, i)); // 取对角元素的绝对值
    double row_sum = 0.0;

    for (int j = 0; j < n; ++j) {
      if (i != j) {
        row_sum += std::abs(A(i, j)); // 非对角元素绝对值之和
      }
    }

    // 判断是否满足对角占优条件
    if (diag_value > row_sum) {
      dominant_rows++;
    }
  }

  // 返回满足对角占优的行比例
  return static_cast<double>(dominant_rows) / n;
}
