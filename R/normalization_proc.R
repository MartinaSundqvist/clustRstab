
#' NID normalization procedure
#'
#' @importFrom stats sd
#' @export
#'

NID.normalizaton <- function(score.pairs, score.pairs.boot){
  mean.normalized <- colMeans(score.pairs - score.pairs.boot) + 1
  sd.normalized <- apply(t(score.pairs) - colMeans(score.pairs.boot), 1, sd)
  return(list(mean.normalized = mean.normalized,
              sd.normalized = sd.normalized))
}


#' ARI normalization procedure
#'
#' @importFrom stats sd
#' @export

ARI.normalizaton <- function(score.pairs, score.pairs.boot){
  mean.normalized <- colMeans(score.pairs - score.pairs.boot)
  sd.normalized <- apply(t(score.pairs) - colMeans(score.pairs.boot), 1, sd)
  return(list(mean.normalized = mean.normalized,
              sd.normalized = sd.normalized))
}
