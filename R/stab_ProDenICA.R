ProDenICA_modified <-
  function(x, k=p, W0=NULL, whiten=FALSE, maxit = 20, thresh = 1e-7, restarts = 0,
           trace = FALSE, Gfunc=GPois, eps.rank=1e-7, ...)
  {
    this.call=match.call()
    p <- ncol(x)
    n <- nrow(x)
    x <- scale(x, T, F)## x should have mean zero
    if(whiten){## First sphere the data
      sx <- svd(x)
      ##Get the effective rank
      condnum=sx$d;condnum=condnum/condnum[1]
      good=condnum >eps.rank
      rank=sum(good)
      if(k>rank){
        warning(paste("Rank of x is ",rank,"; k reduced from",k," to ",rank,sep=""))
        k=rank
      }
      x <- sqrt(n) * sx$u[,good]# no need to rotate via  %*% t(sx$v[,good])
      whitener=sqrt(n)*scale(sx$v[,good],FALSE,sx$d[good])
    }
    else whitener=NULL
    ### Get a random start if needed
    if(is.null(W0))	W0 <- matrix(rnorm(p * k), p, k) else   k=ncol(W0)
    W0 <- ICAorthW(W0)
    ###Initialization
    GS <- matrix(0., n, k)
    gS <- GS
    gpS <- GS
    s <- x %*% W0
    flist <- as.list(1.:k)
    for(j in 1.:k)
      flist[[j]] <- Gfunc(s[, j], ...)
    flist0 <- flist
    crit0 <- mean(sapply(flist0, "[[", "Gs"))
    ### can try some better starts; only evaluated at first iteration
    while(restarts) {
      W1 <- matrix(rnorm(p * k), p, k)

      W1 <- ICAorthW(W1)

      s <- x %*% W1
      for(j in 1.:k){

        flist[[j]] <- Gfunc(s[, j], ...)

      }
      crit <- mean(sapply(flist, "[[", "Gs"))
      if(trace)
        cat("old crit", crit0, "new crit", crit, "\n")
      if(crit > crit0) {
        crit0 <- crit
        W0 <- W1
        flist0 <- flist
      }
      restarts <- restarts - 1
    }
    ###Here is the loop
    nit <- 0
    nw <- 10

    repeat {
      nit <- nit + 1
      gS <- sapply(flist0, "[[", "gs")
      gpS <- sapply(flist0, "[[", "gps")
      t1 <- t(x) %*% gS/n
      t2 <- apply(gpS, 2, mean)
      W1 <- t1 - scale(W0, F, 1/t2)
      W1 <- ICAorthW(W1)

      if(trace)
        cat("Iter", nit, "G", crit0, "crit", nw, "\n")
      W0 <- W1
      if((nit > maxit) | (nw < thresh))
        break
      s <- x %*% W0
      for(j in 1:k){

        flist0[[j]] <- Gfunc(s[, j], ...)

      }
      crit0 <- mean(sapply(flist0, "[[", "Gs"))
    }

    rl=list(W = W0, negentropy = crit0,s= x %*% W0,whitener=whitener,call=this.call)
    rl$density=lapply(flist0,"[[","density")
    class(rl)="ProDenICA"
    rl
  }


centrotype = function(X, Sim, cluster_index)  {
  if(length(cluster_index) > 1){
    temp = which.max(apply(Sim[cluster_index, cluster_index], 1, sum))
    return(X[,cluster_index[temp]])
  }
  return(X[,cluster_index])
}

stab_index = function(Sim,cluster_index)  {
  temp = Sim[cluster_index,cluster_index]
  ex_cluster = setdiff(seq(1, ncol(Sim)), cluster_index)

  aics = (1 / length(cluster_index) ^ 2) * sum(temp)
  aecs = (1 / (length(ex_cluster) * length(cluster_index))) * sum(Sim[cluster_index, ex_cluster])

    return(aics-aecs)
}


mc_stab_ProDenICA =function(table,k,n_runs,nb_threads=12){

  if(n_runs==1) {
    ica=NULL
    while(is.null(ica)) {
      ica = NULL
      err = tryCatch( {
        ica = ProDenICA_modified(table, k=k)
      },error = function(e){})
    }
    return(ica)
  }
  
  ICAs = parallel::mclapply(1 : n_runs, \(i){
    ica = NULL
    while(is.null(ica)) {
      ica = NULL
      err = tryCatch( {
        ica = ProDenICA_modified(table, k=k)
      },error = function(e){})
    }
    ica
  }, mc.cores = nb_threads)

  components = matrix(nrow = nrow(table), ncol = 0)
  negentropies = c()

  for(ica in ICAs){
    err = tryCatch({
      negentropies = c(negentropies, ica$negentropy)
      components = cbind(components, ica$s)
    },error = function(e){ 
      negentropies = c(negentropies, max(negentropies))
    })
  }

  best_ica = ICAs[[which.min(negentropies)]]


  Sim = abs(cor(components))
  hc = hclust(as.dist(1-Sim), method = 'average')
  cluster_index = cutree(hc, k=k)

  Stabilized = parallel::mclapply( 1:k, \(i){
    cluster_i = which(cluster_index==i)
    c(centrotype(components, Sim, cluster_i),
    stab_index(Sim, cluster_i))


  },mc.cores=nb_threads)


  res_table = do.call("rbind",Stabilized)
  res_table = res_table[order(res_table[,nrow(table)+1],decreasing=TRUE),]
  centrotype = t(res_table[,1:nrow(table)])
  stab_index = as.vector(res_table[,nrow(table)+1])

  return(list(best_ica=best_ica, negentropy = negentropies, stab_component=centrotype, stab_index=stab_index))


}


mc_stable_fastica <-
  function(table, n_comp, n_runs, maxit=2000, method=c("R","C"), 
          bootstrap = FALSE, jackknife=FALSE, verbose=FALSE, n_jackknife= NULL, nb_threads = parallel::detectCores()-1) {

  r = nrow(table)
  p = ncol(table)
  X = table

  if(n_comp > p) {
    print(paste0("Please choose a number of component smaller than ", p," (ncol of table)"))
    return(NULL)
  }

  print(paste0("The mc_stable_fastica function use ", min(nb_threads,n_runs), " cores"))

  if(verbose)  
    print("Compute independent components for each run")


  ICAs = parallel::mclapply(1:n_runs, \(i)  {  
    if(verbose)
      print(paste0("runs ", i, "/",n_runs))

    if(bootstrap)
      X = table[, sample(1:p, p,replace=TRUE)]
    if(jackknife && is.numeric(n_jackknife))
       X = table[, sample(1:p, n_jackknife)]
    
    ica = fastICA::fastICA(X, n.comp = n_comp, maxit = maxit, method = method)
    ica$rmse =  mean(sqrt(((ica$S %*% ica$A) - ica$X) ^ 2)) 
    ica
  
  },mc.cores = nb_threads)

  components = do.call("cbind",lapply(ICAs,\(x){x$S} ))
  rmse = do.call("c",lapply(ICAs,\(x){x$rmse}))

  if(verbose)
    print(" Compute similarity matrix")
    
  Sim = abs(cor(components))
  hc = hclust(as.dist(1 - Sim), method = 'average')
  cluster_index = cutree(hc, k = n_comp)

  if(verbose)
    print("Compute the stability index and the centrotype for each component")

  STABILIZED = parallel::mclapply( 1:n_comp, \(i) {
    cluster_i = which(cluster_index == i)
    c(centrotype(components, Sim, cluster_i), stab_index(Sim, cluster_i))


  },mc.cores = nb_threads)


  res_table = do.call("rbind", STABILIZED)
  res_table = res_table[order(res_table[, r + 1], decreasing = TRUE), ]
  centrotype = t(res_table[, 1:r])
  stab_index = as.vector(res_table[, r + 1])

  result = list( rmse = rmse, stab_component = centrotype, stab_index = stab_index)
  class(result)="stableFastICA"
  result

}
