centrotype = function(X, Sim, cluster_index)  {
  if(length(cluster_index) > 1){
    temp = which.max(apply(Sim[cluster_index, cluster_index], 1, sum))
    return(X[,cluster_index[temp]])
  }
  return(X[,cluster_index])
}