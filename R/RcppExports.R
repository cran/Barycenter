# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Greenkhorn <- function(r, c, costm, lambda = 1, maxIter = 10000L, tolerance = 10^(-8)) {
    .Call('_Barycenter_Greenkhorn', PACKAGE = 'Barycenter', r, c, costm, lambda, maxIter, tolerance)
}

Sinkhorn <- function(a, b, costm, lambda = 1, maxIter = 10000L, tolerance = 10^(-8)) {
    .Call('_Barycenter_Sinkhorn', PACKAGE = 'Barycenter', a, b, costm, lambda, maxIter, tolerance)
}

Subgradient <- function(a, b, M, lambda, maxIter = 1000L, tolerance = 0.0001) {
    .Call('_Barycenter_Subgradient', PACKAGE = 'Barycenter', a, b, M, lambda, maxIter, tolerance)
}

