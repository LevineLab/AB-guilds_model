##' Calculates \eqn{r_fk}.
##'
##' @param obj An object from \code{em_fast()} (remember,
##'   \code{aspect_bernoulli()} contains a list of such objects, from random
##'   restarts).
get_r_fk <- function(obj){
  ## Format beta and Gamma
  beta = obj$beta
  Gamma =  obj$gamma
  FF = ncol(beta)
  K = nguild = nrow(beta)
  G = nrow(Gamma)

  ## Add row and column names
  rownames(Gamma) <- 1:G
  colnames(Gamma) <- 1:K
  rownames(beta) <- 1:K
  colnames(beta) <- colnames(dat)

  ## Make genome-specific guild probabilities
  Vmult = Gamma %*% beta
  guild_prob_all = array(NA, dim = c(G, K, FF))
  for(f in 1:FF){
    for(g in 1:G){
      numers_by_guild = Gamma[g,] * beta[,f]
      denom = Vmult[g,f]
      if(denom == 0){
        guild_prob_all[g,,f] = 1/K
      } else {
        guild_prob_all[g,,f] = numers_by_guild/denom
      }
    }
  }
  dimnames(guild_prob_all)[[1]] <- 1:G
  dimnames(guild_prob_all)[[2]] <- 1:K
  dimnames(guild_prob_all)[[3]] <- colnames(beta)
  guild_prob = guild_prob_all %>% apply(c(2,3), mean)
  dimnames(guild_prob) = dimnames(beta)
  return(guild_prob)
}




##' Gets scores from an object.
##'
##' @param obj1 An object from \code{em_fast()} (remember,
##'   \code{aspect_bernoulli()} contains a list of such objects, from random
##'   restarts).
##' @return Table with guild function and scores.
get_scores <- function(obj1){

  ## Compute mapback genomes
  all_r_fk <- get_r_fk(obj1)

  KK = ncol(obj1$gamma)
  guild_score_list <- lapply(1:KK, function(kk){

    ## 1. Reduce each to A_k
    A_k = which(obj1$gamma[,kk] > 2/KK)
    ysums = colSums(dat[A_k,])
    q_fk = ysums / mean(ysums)

    ## 2. Compute scores $r_{fk}$ times $q_{fk}$
    s_fk = q_fk * all_r_fk[kk,]
  })
  scoremat = do.call(cbind, guild_score_list)
  return(scoremat)
}



##' Gets guilds and their scores from an object \code{obj}.
##'
##' @param obj An object from \code{em_fast()} (remember,
##'   \code{aspect_bernoulli()} contains a list of such objects, from random
##'   restarts).
##' @return Table with guild function and scores.
mapback <- function(obj){

  ## Compute mapback genomes
  all_r_fk <- get_r_fk(obj)

  KK = ncol(obj$gamma)
  guild_funs_list <- lapply(1:KK, function(kk){

    ## 1. Reduce each to A_k
    A_k = which(obj$gamma[,kk] > 2/KK)
    ysums = colSums(dat[A_k,])
    q_fk = ysums / mean(ysums)

    ## 2. Compute scores $r_{fk}$ times $q_{fk}$
    s_fk = q_fk * all_r_fk[kk,]

    ## 3. Further reduce to B_k
    FF = ncol(obj$gamma)
    oo = order(s_fk, decreasing = TRUE)
    for(ii in 1:FF){

      ## Top ii functions are oo[1:ii]
      guild_funs = oo[1:ii]

      ## Calculate number of mapback genomes
      mapbacks = sum(rowSums(dat[,guild_funs, drop = FALSE]) == ii)
      if(mapbacks <= 100) break
    }

    ## Guilds
    ii = ii-1
    guild_funs = oo[1:ii]
    df = data.frame(score = round(s_fk[guild_funs],3), fun_names = names(s_fk)[guild_funs])
    return(df)
  })
  return(guild_funs_list)
}




##' Fit the Aspect Bernoulli model
##'
##' Basic implementation of EM algorithm for the Aspect Bernoulli model,
##'
##' E[x] = gamma * beta
##'
##' @param x n-by-p binary matrix
##' @param k number of latent dimensions
##' @param gamma optional initialization of gamma
##' @param beta optional initialization of beta
##'
##' @return
##' - gamma an n-by-k matrix
##' - beta a k-by-p matrix
##'
##' @references
##' Bingham, E., Kabán, A., & Fortelius, M. (2009). The aspect Bernoulli model:
##' multiple causes of presences and absences. Pattern Analysis and Applications,
##' 12(1), 55-78.
aspect_bernoulli <- function(x, k, num_init = 10, num_steps = 100,
                             em_alg = em_fast, mc.cores = 1) {
  n <- nrow(x)
  p <- ncol(x)
  fit <- list()
  ## for (i in seq(num_init)) {
  fit = parallel::mclapply(seq(num_init), function(irep){

    # initialize
    gamma <- matrix(runif(n * k), n, k) ## Dirichlet?
    gamma <- gamma / rowSums(gamma)
    beta <- matrix(runif(k * p), k, p) ## Beta distribution?

    # run EM
    one_fit <- em_alg(x, gamma = gamma, beta = beta, num_steps = num_steps)
    ## fit[[i]] <- one_fit
    return(one_fit)
  }, mc.cores = mc.cores)

  return(fit)
}

em <- function(x, gamma, beta, num_steps) {
  # This is a vanilla implementation of the EM, but will be very slow because
  # of all the for loops and the responsibilities.  The em_fast algorithm is
  # much faster.
  n <- nrow(x)
  p <- ncol(x)
  k <- ncol(gamma)
  stopifnot(nrow(beta) == k)
  resp <- array(NA, c(n, p, k))
  obj <- rep(NA, num_steps)
  for (i in seq(num_steps)) {
    # E-step: update responsibilities
    for (g in seq(n)) {
      for (f in seq(p)) {
        for (l in seq(k)) {
          resp[g, f, l] <- gamma[g, l] * ifelse(x[g, f] == 1,
                                                beta[l, f],
                                                1 - beta[l, f])
        }
        sum_resp_gf <- sum(resp[g, f, ])
        resp[g, f, ] <- resp[g, f, ] / sum_resp_gf
      }
    }
    # M-step: update gamma and beta
    for (f in seq(p)) {
      for (l in seq(k)) {
        beta[l, f] <- sum(resp[, f, l] * x[, f]) / sum(resp[, f, l])
      }
    }
    for (g in seq(n)) {
      for (l in seq(k)) {
        gamma[g, l] <- sum(resp[g, , l]) / p
      }
    }
    # compute objective
    obj[i] <- log_lik_fast(gamma, beta, x)
    ## cat(obj[i], fill = TRUE)
  }
  list(gamma = gamma, beta = beta, obj = obj)
}

em_fast <- function(x, gamma, beta, num_steps, tol = 1E-6){
  # this is an implementation of the algorithm described in Section 2.2.1 of the
  # Aspect Bernoulli paper.
  # Our notation relative to theirs has a few transposes.
  # See handwritten notes aspect_Bernoulli.pdf for more on this.
  n <- nrow(x)
  p <- ncol(x)
  k <- ncol(gamma)
  stopifnot(nrow(beta) == k)
  obj <- rep(NA, num_steps)
  x = as.matrix(x)
  betalist = gammalist = list() ## temporary
  stopsign = FALSE
  for (i in seq(num_steps)) {
    ## printprogress(i, num_steps)
    beta_bar <- 1 - beta
    ratio <- x / (gamma %*% beta); ratio[is.nan(ratio)] <- 0 # 0/0 is 0 here
    ratio_bar <- (1 - x) / (gamma %*% beta_bar);ratio_bar[is.nan(as.matrix(ratio_bar))] <- 0 # 0/0 is 0 here
    gamma <- gamma * (ratio %*% t(beta) + ratio_bar %*% t(beta_bar))
    gamma <- gamma / rowSums(gamma)


    ratio <- x / (gamma %*% beta); ratio[is.nan(as.matrix(ratio))] <- 0 # 0/0 is 0 here
    ratio_bar <- (1 - x) / (gamma %*% beta_bar); ratio_bar[is.nan(as.matrix(ratio_bar))] <- 0 # 0/0 is 0 here

    beta <- beta * (t(gamma) %*% ratio)
    beta_bar <- beta_bar * (t(gamma) %*% ratio_bar)
    beta <- beta / (beta + beta_bar)
    ## betalist[[i]] = beta ## temporary
    ## gammalist[[i]] = gamma ## temporary

    # compute objective
    obj[i] <- log_lik_faster(gamma, beta, x)

    ## Every 10 steps, calculate convergence
    if((i > 1) & (i %% 10 == 0)){
      if(abs(obj[i] / obj[i-1] - 1) < tol) stopsign = TRUE
    }
    if(stopsign) break
  }
  obj = obj %>% na.omit() %>% as.numeric()
  list(gamma = gamma, beta = beta, obj = obj)
  ##   betalist = betalist, gammalist = gammalist) ## Temporary
}

log_lik <- function(gamma, beta, x) {
  # sum_gf log( sum_k gamma_gk beta_kf^y_gf(1-beta_kf)^(1-y_gf) )
  n <- nrow(gamma)
  k <- ncol(gamma)
  p <- ncol(beta)

  ll <- 0
  for (g in seq(n)) {
    for (f in seq(p)) {
      arg <- 0
      for (l in seq(k)) {
        arg <- arg + gamma[g, l] * ifelse(x[g, f] == 1,
                                          beta[l, f],
                                          1 - beta[l, f])
      }
      ll <- ll + log(arg)
    }
  }
  ll
}

log_lik_fast <- function(gamma, beta, x) {
  # sum_gf log( sum_k gamma_gk beta_kf^y_gf(1-beta_kf)^(1-y_gf) )
  n <- nrow(gamma)
  k <- ncol(gamma)
  p <- ncol(beta)

  ll <- 0
  for (g in seq(n)) {
    for (f in seq(p)) {
      ll <- ll + log(ifelse(x[g, f] == 1,
                    sum(gamma[g, ] * beta[, f]),
                    sum(gamma[g, ] * (1 - beta[, f]))))
    }
  }
  ll
}

log_lik_faster <- function(gamma, beta, x) {
  # sum_gf log( sum_k gamma_gk beta_kf^y_gf(1-beta_kf)^(1-y_gf) )
  n <- nrow(gamma)
  k <- ncol(gamma)
  p <- ncol(beta)

  pi <- gamma %*% beta
  pi_bar <- gamma %*% (1 - beta)
  sum(log(x * pi + (1 - x) * pi_bar))
}




##' From two probability matrices, form a (K x K) distance matrix of the
##' (n)-vectors. The distance between the vectors is the symmetric KL
##' divergence.
##'
##' @param mat1 Matrix 1 of size (n x K).
##' @param mat2 Matrix 2 of size (n x K).
##'
##' @return Return
form_symmetric_kl_distmat <- function(mat1, mat2){

  ## Manually add some small, in case some columns are all zero
  mat1 = (mat1 + 1E-10) %>% pmin(1)
  mat2 = (mat2 + 1E-10) %>% pmin(1)

  KK1 = ncol(mat1)
  KK2 = ncol(mat2)
  distmat = matrix(NA, ncol=KK2, nrow=KK1)
  for(kk1 in 1:KK1){
    for(kk2 in 1:KK2){
      ## if(kk1 > kk2) next
      ## mydist = sqrt(sum((mat1[,kk1, drop=TRUE]- mat2[,kk2, drop=TRUE])^2))
      mydist = symmetric_kl(mat1[,kk1, drop=TRUE], mat2[,kk2, drop=TRUE])
      distmat[kk1, kk2] = mydist
    }
  }
  stopifnot(all(!is.na(distmat)))
  return(distmat)
}

##' Symmetric KL divergence.
symmetric_kl <- function(vec1, vec2){
  stopifnot(all(vec1 <= 1) & all(vec1 >= 0))
  stopifnot(all(vec2 <= 1) & all(vec2 >= 0))
  kl <- function(vec1, vec2){
    sum(vec1 * log(vec1 / vec2))
  }
  return((kl(vec1, vec2) + kl(vec2, vec1))/2)
}


##' Get ordering of columns of mat1 that match mat2.
match_guilds <- function(mat1, mat2){

  ## Make distance matrix
  distmat <- form_symmetric_kl_distmat(mat1, mat2)

  ## Use Hungarian Algorithm to solve matching problem.
  fit <- clue::solve_LSAP(distmat)
  o <- as.numeric(fit)

  ## Return the optimal ordering of mat1's columns
  return(o)
}
