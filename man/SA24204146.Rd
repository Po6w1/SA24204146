\name{SA24204146}
\alias{SA24204146}
\docType{package}
\title{
  Analysis and Simulation of High-Dimensional Random Matrices
}
\description{
  The \code{SA24204146} package provides tools for generating, analyzing, and simulating high-dimensional random matrices. It includes functionalities for computing eigenvalues, analyzing spectral densities, comparing empirical results with theoretical predictions (e.g., Marčenko-Pastur law), and efficiently computing diagonal dominance using \code{Rcpp} and \code{RcppEigen}.
}
\details{
  This package is designed for researchers and students interested in random matrix theory and its applications. It combines computational efficiency (via Rcpp and Eigen) with detailed statistical analysis of eigenvalue distributions and diagonal dominance in high-dimensional settings.

  Key features include:
  \itemize{
    \item Simulation of random matrices with configurable dimensions and distributions.
    \item Computation and visualization of empirical spectral densities (ESD) compared to theoretical densities like Marčenko-Pastur.
    \item Kernel-based density estimation for eigenvalue distributions.
    \item Efficient computation of diagonal dominance for covariance matrices using \code{Rcpp}.
  }
}
\author{
  Zheng Gong

  Maintainer: Zheng Gong <gongzheng@mail.ustc.edu.cn>
}
\references{
  Bai, Z. D., & Silverstein, J. W. (2010). Spectral Analysis of Large Dimensional Random Matrices. Springer.

  Marčenko, V. A., & Pastur, L. A. (1967). Distribution of eigenvalues for some sets of random matrices. *Mathematics of the USSR-Sbornik*, *1*(4), 457–483.
}
\keyword{ package }
\examples{
  analyze_random_matrix_spectrum(n = 500, p = 400)
}
