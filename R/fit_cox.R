#' Fit a Cox proportional hazards model to survival data
#'
#' @param x Design matrix of data, assumes survival times sorted in INCREASING order
#' @param censor Censoring vector of data
#'
#' @keywords survival, Cox proportional hazards
#'
#' @importFrom spatstat.utils revcumsum
#' @importFrom survival coxph
#' @importFrom bench mark
#' @import tidyr
#' @import ggbeeswarm
#' @import ggplot2
#'
#' @return Vector of estimated parameters for each covariate
#'
#' @examples
#' p = 5
#' n = 100
#' c0 = 0.5
#' beta = rnorm(p, 0, 1)
#' x = matrix(rnorm(n * length(beta)), ncol = length(beta))
#' yvar= rnorm(n, 0, 1)
#' times=rexp(n, c0 * exp(rowSums(t(t(x) * beta)) + yvar))
#' time.censor=rexp(n, c0 * exp(beta[1] * x[, 1] + yvar))
#' censorv=ifelse(times < time.censor, 1, 0)
#' y = ifelse(times < time.censor, times, time.censor)
#' ord = order(y, decreasing =F)
#' dataset = list(x = x[ord,], y = y[ord], censor = censorv[ord])
#' fit_cox(dataset$x, dataset$censor)
#' all.equal(fit_cox(dataset$x, dataset$censor), as.numeric(coxph(Surv(dataset$y, dataset$censor) ~ dataset$x)$coefficients))
#'
#' @export


fit_cox <- function(x, censor){
  p = ncol(x)
  beta_hat<-list()
  beta_hat[[1]] = matrix(rep(0,p)) ##initialize initial beta vector of zeroes
  theta = as.numeric(exp(x %*% beta_hat[[1]]))
  P = outer(theta,revcumsum(theta),'/') ##failure probabilities
  P[upper.tri(P)] <- 0
  W <- -P %*% diag(censor) %*% t(P)
  diag(W) <- diag(P %*% diag(censor) %*% t(1-P))
  xtwx = t(x) %*% W %*% x ##calculate the Hessian
  xtp = t(x) %*% (censor - P %*% censor) ##calculate the score

  ##iterate the process until convergence
  beta_hat[[2]] = beta_hat[[1]] + (solve(xtwx) %*% xtp)
  j=2
  while (norm(beta_hat[[j]] - beta_hat[[j-1]],type="f")>1e-5) {
    j=j+1
    theta = as.numeric(exp(x %*% beta_hat[[j-1]]))
    P = outer(theta,revcumsum(theta),'/')
    P[upper.tri(P)] <- 0
    W <- -P %*% diag(censor) %*% t(P)
    diag(W) <- diag(P %*% diag(censor) %*% t(1-P))
    xtwx = t(x) %*% W %*% x
    xtp = t(x) %*% (censor - P %*% censor)
    beta_hat[[j]]<-beta_hat[[j-1]] + solve(xtwx) %*% xtp
  }

  return(as.numeric(beta_hat[[j]]))
}
