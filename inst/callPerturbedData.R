
subSample <- function(data, param = c(nProp, pProp)){
  # param = dimension
  n <- nrow(data)
  p <- ncol(data)
  nIndex <- sort(sample(n, param[1]*n))
  pIndex <- sort(sample(p, param[2]*p))
  data[nIndex, pIndex]
}

randProj <- function(i, data, param){
  # param = dimension
}

noiseGaussian <- function(data, param = c(0, 1)) {
  n <- nrow(data)
  p <- ncol(data)
  noise <- rnorm(n*p, param[1], param[2])
  noiseMat <- matrix(noise, nrow = n, ncol = p)
  perturbedData <-  as.data.frame(as.double(as.matrix(data)) + noiseMat)
  colnames(perturbedData) <- colnames(data)
  rownames(perturbedData) <- rownames(data)
  perturbedData
}
