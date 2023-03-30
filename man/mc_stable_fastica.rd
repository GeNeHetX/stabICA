
\name{mc_stable_fastica}
\alias{mc_stable_fastica }
\title{Stabilized ICA parallel}
\description{
  This is an R implementation of the stabilized-ICA algorithm
  of ncaptier (\url{https://github.com/ncaptier/stabilized-ica/}) to
  Compute the centrotype, the root mean square error and stability index of a data matrix .
}
\usage{
  mc_stable_fastica(table, n_comp, n_runs,...)
}
\arguments{
  \item{table}{a numeric data matrix with rows representing observation and columns representing varaibles }
  \item{n_comp}{int, number of component to extract of table}
  \item{maxit}{int, representing the maximum number of iteration fastICA algorithm must be run}
  \item{method}{charcater, "R" or "C" representing what language use fastICA algorithm}
  \item{bootstrap}{boolean, representing if algorithm must be do a boostrap (False by default)}
  \item{jackknife}{boolean, representing if algorithm must be sub sampling the data matrix(False by default). If TRUE is provide, the argument n_jacknife must be provide too.}
  \item{n_jackknife}{int, representing the size of subsampling, the number of column must be kepp for the following analysis}
  \item{nb_threads}{int, representing the number of thread, algorithm can be use (total core -1 by default)}
  \item{verbose}{ boolean, representing if steps must be print ( False by default)}
}
\details{
  \bold{mc_stable_fastica} 
  The algorithm can be split into 3 step :
  - Compute independent components for each run
  - Compute similarity matrix
  - Clustering components
  - Compute the stability index and the centrotype for each component
  - Return the obejct of class stableFastICA
}
\value{
  An object of S3 class \code{"stableFastICA"} is returned, with the following components:
  \item{rmse}{a vector of float representing the rmse of each component}
  \item{stab_component}{a matrix of stabilized components ( S matrix) with row representing observation and column representing components}
  \item{stab_index}{a vector of float representing the stability of each component}
}
\author{
  Audrey Beaufils
}
\examples{

  #---------------------------------------------------
  # Example 1 : A simple example
  #---------------------------------------------------
  set.seed(123)
  S <- matrix(runif(10000), 2500, 4)
  A <- matrix(sample(1:1000,100), 4, 25, byrow = TRUE)
  X <- S \%*\% A

  res<-mc_stable_fastica(table = X, n_comp = 4, n_runs = 1)

  #---------------------------------------------------
  # Example 2 : An example using jackknife
  #---------------------------------------------------

  set.seed(123)
  S <- matrix(runif(10000), 2500, 4)
  A <- matrix(sample(1:10000,1000), 4, 250, byrow = TRUE)
  X <- S \%*\% A
  res <- mc_stable_fastica(table = X, n_comp = 4, n_runs = 5, jackknife = TRUE, n_jackknife = 50)

  #---------------------------------------------------
  # Example 3 : An example ploting results
  #---------------------------------------------------

  set.seed(123)
  S <- matrix(runif(10000), 2500, 4)
  A <- matrix(sample(1:10000,1000), 4, 250, byrow = TRUE)
  X <- S \%*\% A

  stab_ind = matrix(nrow=0, ncol=2)
  rmse = matrix(nrow=0, ncol=2)
  res_li = list()
    
  for( i in 2:5){
    res <- mc_stable_fastica(table=X,n_comp=i,n_runs=5, jackknife= TRUE, n_jackknife =50 )
    stab_ind=rbind(stab_ind, cbind(res$stab_index,rep(i,i)))
    rmse = rbind(rmse, cbind(res$rmse,rep(i,5)))
    res_li =append(res_li,res)
  }

  colnames(rmse)=c("rmse","k")
  colnames(stab_ind)=c("stab_ind","k")
  rmse=as.data.frame(rmse)
  stab_ind=as.data.frame(stab_ind)

  ggpubr::ggboxplot(rmse,'k','rmse',color="k",title="Boxplot of rmse depends on number of component")
  ggpubr::ggboxplot(stab_ind,'k','stab_ind',color="k",title="Boxplot of rmse depends on number of component")

}