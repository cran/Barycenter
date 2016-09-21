\name{eight}
\alias{eight}
\title{A small dataset to run the function Barycenter}
\format{List}
\usage{
eight
}
\description{
A list of five images of the digit eight. Each image is a \code{28 x 28} matrix, representing the weights of the corresponding pixels scaled s.t. they sum up to one. The images are included in a list. This dataset is based on the MNIST dataset.
}
\examples{
#Take a look into the dataset, for example the first image of the digit eight.
eight[[1]]
image(eight[[1]])
}
\keyword{dataset}
