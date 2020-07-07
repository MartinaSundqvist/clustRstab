# ------------------------------------------------------------------------
# 4. A function to execute the perturbation and clustering of the data, and organize it by K
#-------------------------------------------------------------------------

# getClsByK() perturbes the data and compute clustering, organize the results after corresponding K

# - Get a list where each item corresponds to the classification of perturbed dataset df.x. cols = kVec, rowns = 1:n
# - From this create a list of length(kVec) items, where each item is n*nsim x 2 df
# - These df.s are then split each n (corresponding to each perturbed df), and cocanated in a n x length(kVec) df
# - These tables are then concatenated by rbind to a n*nsim x length(kVec) + 1 data.frame

getClsByK <- function(clsList, kVec){
  lapply(lapply(kVec, function(k)
    do.call(what = rbind,
            args = clsList
    )[,c(as.character(k), "df")]), function(j)
      sapply(split(x = j, f = j$df), function(i)
        i[,-ncol(i)]
      )
  )
}
