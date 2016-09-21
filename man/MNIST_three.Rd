\name{MNIST_three}
\alias{MNIST_three}
\docType{data}
\title{
The full MNIST dataset of the digit three.
}
\description{
A list of 4351 images of the digit three. Each image is a 28 x 28 matrix, representing the weights of the corresponding pixels scaled s.t. they sum up to one. The images are included in a large list. This dataset is based on the MNIST dataset.
}
\usage{MNIST_three}
\format{List}
\examples{
#Take a look into the dataset, for example the 25th image of the digit three.
MNIST_three[[25]]
image(MNIST_three[[25]])
}
\keyword{datasets}
