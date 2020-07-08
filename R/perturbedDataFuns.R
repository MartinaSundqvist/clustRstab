
#'
#' @export
subSample <- function(data, nProp, pProp, ...){
  # param = dimension
  n <- nrow(data)
  p <- ncol(data)
  nIndex <- sort(sample(n, nProp*n))
  pIndex <- sort(sample(p, pProp*p))
  perturbedData <- data[nIndex, pIndex]
  perturbedData
}

#' @export
randProjData <- function(data, randProjDim, randProjMethod, ...){
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
  noise <- rnorm(n*p, noiseGaussianMean, noiseGaussianSD)
  noiseMat <- matrix(noise, nrow = n, ncol = p)
  perturbedData <-  as.data.frame(as.double(as.matrix(data)) + noiseMat)
  colnames(perturbedData) <- colnames(data)
  rownames(perturbedData) <- rownames(data)
  perturbedData
}

