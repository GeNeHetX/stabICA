\name{centrotype}
\alias{centrotype}
\title{Compute centrotype}
\description{
  This is an R implementation of the centrotype function of stabilized-ICA algorithm
  of ncaptier (\url{https://github.com/ncaptier/stabilized-ica/}) to
  Compute the centrotype of the cluster of ICA components defined by cluster_index.
}
\usage{
centrotype(X, Sim,cluster_labels)
}
\arguments{
  \item{X}{a data matrix with rows representing observations
    and  columns representing variables.}
  \item{Sim}{a data similarity matrix (n_runs*n_component x n_runs*n_component)}
  \item{cluster_labels}{ a vector representing the labels of cluster f each components}
  
}
\details{
  \bold{Centrotype}
  
  Compute the centrotype of the cluster of ICA components defined by cluster_labels.
      centrotype : component of the cluster which is the most similar to the other components
                   of the cluster
}
\value{
  \item{centrotype}{A submatrix of X, which contains the most stable components. dim(centrotype)= n_observation x n_components}
}
\author{
  Audrey Beaufils
}