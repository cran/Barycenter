\name{Subgradient}
\alias{Subgradient}
\title{
Smoothed Dual Optima
}
\description{
Based on Sinkhorn's matrix scaling algorithm, the function returns the smoothed dual optima vector \eqn{\alpha}. This vector serves as a subgradient to compute a \link{Barycenter}.
}
\usage{
Subgradient(a, b, M, lambda, maxIter = 1000L, tolerance = 1e-04)
}
\arguments{
  \item{a}{
Normed vector of weights in a non-empty closed subset.
}
  \item{b}{
Normed vector of weights in a non-empty closed subset.
}
  \item{M}{
Cost matrix; the cost of mapping \eqn{a} to \eqn{b}.
}
  \item{lambda}{
Non-negative smoothing parameter.
}
  \item{maxIter}{
Maximum number of iterations.
}
  \item{tolerance}{
A threshold for the integrated stopping criterion.
}
}

\references{
Cuturi, M.: \code{Sinkhorn Distances: Lightspeed Computation of Optimal Transport}, Advances in Neural Information Processing Systems 26, 2013

Cuturi, M.: \code{Fast Computation of Wasserstein Barycenters}, Proceedings of the International Conference on Machine Learning, Beijing, China, 2014
}
\author{
Marcel Klatt
}


