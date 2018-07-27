options_default_kmeans <- list(
  iter.max = 10,
  nstart = 1,
  algorithm =  c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
)

options_default_gmm <- list(
  modelNames = NULL
)

options_default_hclust <- list(
  distance = "euclidean",
  method = "ward.D2"
)

options_default_subVar <- list(
  param = 0.8
)
