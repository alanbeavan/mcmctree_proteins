library(BAMMtools)

setwd("~/OneDrive - The University of Nottingham/seals/clock")
dists = read.csv("tree_lengths.csv")
hist(dists$dist, breaks=100, main=NULL, xlab="Tree Height")
breaks = getJenksBreaks(dists$dist, 4, subset = NULL)
abline(v = 0.02028777, col = "red")
abline(v = 0.28105305, col = "red")
abline(v = 0.67620711, col = "red")
abline(v = 2.05288159, col = "red")
min(dists$dist)

cluster = c()
for(i in 1:nrow(dists)) {
  len = length(cluster)
  for(j in seq(3)) {
    if(dists$dist[i] < breaks[j+1]) {
      cluster = c(cluster, j)
      break
    }
  }
  if(length(cluster) == len) {
    cluster = c(cluster, 4)
  }
}
dists$cluster = cluster
write.csv(file = "clusters.csv", dists, quote= FALSE, row.names=FALSE)
