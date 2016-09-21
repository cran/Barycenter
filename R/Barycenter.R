#'Computation of a Wasserstein Barycenter
#'
#'The function \code{Barycenter} takes in a list of matrices and computes the
#'corresponding Barycenter.
#'The list has to consist of matrices having all the same dimensions and each matrix represents the weights of the corresponding pixels of an image. The pixels should be scaled s.t. they sum up to one.
#'@author Marcel Klatt
#'@param images A list of matrices satisfying the prerequisites described above.
#'@param maxIter Maximum gradient iterations.
#'@param lambda Non-negative smoothing parameter (c.f. \link{Subgradient}).
#'@return The Barycenter of the images, represented by a \eqn{n x m} matrix. The function returns the Barycenter represented by a \eqn{n x m} matrix and prints also the corresponding image.
#'
#'Given the MNIST dataset, a Barycenter of the digit three is shown below. The Barycenter is based on 4351 images each represented by
#'a 28 x 28 pixel grid, respectively. The values for \code{lambda} and \code{maxIter} were set by default. The dataset is also available in this package (c.f. \link{MNIST_three}).
#'
#'\figure{threeMNIST.png}{test}
#'@references Cuturi, M.: \code{Fast Computation of Wasserstein Barycenters}, Proceedings of the International Conference on Machine Learning, Beijing, China, 2014
#'@examples #Computation of a Barycenter based on five images representing the digit eight, respectively.
#'\dontrun{Barycenter(eight)}
#'@export
Barycenter <- function(images, maxIter=10, lambda=60/median(costMatrix)){

  time <- proc.time() #to analyze the computation time

  #get dimension of the image as a grid embedded in [0,1]²
  dimension <- dim(images[[1]])
  n <- dimension[1]*dimension[2]

  #create a grid of the same dimension on [0,1]² and create the cost matrix
  coord1 <- seq(0,1,length.out = dimension[2])
  coord2 <- rev(seq(0,1,length.out = dimension[1]))
  coordinates <- expand.grid(coord1,coord2)
  costMatrix <- as.matrix(dist(coordinates, diag=TRUE, upper=TRUE))

  #initialize a_tild and a_hat
  a_tild <- rep(1/n,n)
  a_hat <- rep(1/n,n)
  t_0 <- 2
  t <- t_0


  #iteration using Sinkhorn_Distance (Algorithm 1)
  for(i in 1:maxIter){
    beta <- (t+1)/2
    a <- (1-1/beta) * a_hat + (1/beta) * a_tild
    #Form subgradient with Sinkhorn's Algorithm
    ALPHA <- 0
    for(j in 1:length(images)){
      #This step is based on RcppArmadillo. (Algorithm 3)
      ALPHA <- Subgradient(a,t(images[[j]]),costMatrix,lambda) + ALPHA  #Note, that we need to transpose each matrix, to be compatible with the coordinates defined above.
    }
    ALPHA <- (1/length(images)) * ALPHA
    a_tild <- a_tild*exp(-(t_0)*beta*ALPHA)
    a_tild <- a_tild/sum(a_tild)
    a_hat <- (1-1/beta)*a_hat + (1/beta) * a_tild
    t <- t+1
  }


  #Transforming Barycenter s.t. function "image" returns the correct orientation.
  a <- matrix(a,dimension[1],dimension[2],byrow=TRUE)
  a.temp <- a[,nrow(a):1]

  #Output
  print(image(a.temp))
  print(proc.time()-time) #Computation time
  return(a)
}
