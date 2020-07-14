
#' A function to subsample a dataset
#'
#' subSample() subsamples nProp observation and pProp variables of a n x p data.frame
#'
#'
#'@param data A n x p data.frame() with n the number of observations and p the number of varaibles
#'@param nProp The proportion of subsampled observations
#'@param pProp The proportion of subsampled variables
#'
#'
#' @export
subSample <- function(data, nProp, pProp, ...){
  # param = dimension
  n <- nrow(data)
  p <- ncol(data)
  rownames(data) <- 1:n
  nIndex <- sort(sample(n, nProp*n))
  pIndex <- sort(sample(p, pProp*p))
  perturbedData <- data[nIndex, pIndex]
  perturbedData
}

#' @export
randProjData <- function(data, randProjDim, randProjMethod, ...){
  rownames(data) <- 1:nrow(data)
  mat <- RPGenerate(p = ncol(data), d = randProjDim, method = randProjMethod, B2 = 1)
  perturbedData <- crossprod(mat, t(data.matrix(data)))
  perturbedData <- as.data.frame(t(perturbedData))
  rownames(perturbedData) <- rownames(data)
  perturbedData
}

#' @export
noiseGaussian <- function(data, noiseGaussianMean, noiseGaussianSD, ...) {
  n <- nrow(data)
  p <- ncol(data)
  rownames(data) <- 1:n
  noise <- rnorm(n*p, noiseGaussianMean, noiseGaussianSD)
  noiseMat <- matrix(noise, nrow = n, ncol = p)
  perturbedData <-  as.data.frame(as.double(as.matrix(data)) + noiseMat)
  colnames(perturbedData) <- colnames(data)
  rownames(perturbedData) <- rownames(data)
  perturbedData
}

