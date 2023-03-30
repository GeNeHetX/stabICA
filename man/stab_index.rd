\name{stab_index}
\alias{stab_index}
\title{Compute stability index}
\description{
  This is an R implementation of the _stability_index function of stabilized-ICA algorithm
  of ncaptier (\url{https://github.com/ncaptier/stabilized-ica/}) to
  Compute the stability index for the cluster of ICA components defined by cluster_index.
}
\usage{
  stab_index(Sim, cluster_index)
}
\arguments{
  \item{Sim}{a data similarity matrix (n_runs*n_component x n_runs*n_component)}
  \item{cluster_lindex}{ a vector representing the labels of cluster of each component}
}
\details{
  \bold{Stability_index} 
  Compute the difference between  average intra-cluster similarities and  average extra-cluster similarities
}
\value{
  \item{stab_index}{Float between 0 and 1
        stability index for the cluster of ICA components defined by cluster_labels}
}