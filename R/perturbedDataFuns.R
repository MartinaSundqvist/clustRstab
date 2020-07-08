
#'
#' @export
subSample <- function(data, param = c(nProp, pProp)){
  # param = dimension
  n <- nrow(data)
  p <- ncol(data)
  nIndex <- sort(sample(n, param[1]*n))
  pIndex <- sort(sample(p, param[2]*p))
  perturbedData <- data[nIndex, pIndex]
  perturbedData
}

randProj <- function(data, param){
  mat <- RPGenerate(p = ncol(data), d = param, method = "Haar", B2 = 1)
  perturbedData <- crossprod(mat, as.matrix(data))
  #perturbedData <- crossprod(mat, matrix(as.double(as.matrix(data)), nrow = nrow(data), ncol = ncol(data), byrow =T))
  #perturbedData <- t(mat) %*% as.double(as.matrix(data))
  as.data.frame(perturbedData)
}

#' @export
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

