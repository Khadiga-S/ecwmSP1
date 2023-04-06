logistic <- function(x) exp(x) / ( 1 + exp(x))


## One likelihood function for both prevalence and logistic regression

logl <- function(p2,  b, D1 = NULL, D2 = NULL, only1 = NULL){

  if (is.null(only1)) {

    pr    <- logistic(D1 %*% b) # pi

    theta <- 0

  }else{
    beta_hat <- matrix(b, nrow = ncol(D1) + ncol(D2)) # 1st for pt 2nd for theta

    pr      <- logistic(D1 %*% beta_hat[1:ncol(D1), ]) # pi

    theta   <- logistic(D2 %*% beta_hat[ncol(D1) + 1:ncol(D2), ])
  }

  pi_02 <- ((p2[,"p02"] * pr) + (1 - p2[,"p02"]) * (1 - pr))*(1 - theta)
  pi_11 <- (p2[, "p11"] + (1 - p2[,"p11"]) * theta) * pr + ((1 - p2[,"p11"]) + (p2[,"p11"] * theta)) * (1 - pr)

  lik <- ifelse(is.na(p2[,"p02"]), pi_11,  pi_02)

  - sum(log(lik))
}


## delta method to get the estimates and its se on a probability scale

prev <- function(coef, hess, level){

  jaco  <- numDeriv::jacobian(logistic, x = coef)

  vpar <- t(jaco) %*% solve(hess) %*% jaco

  se <- sqrt(diag(vpar))

  estimate <- data.frame(est   = c(logistic(coef)[1],
                                   logistic(coef)[level + 1]),

                         se_pi    = c(se[1], se[level + 1]),

                         min95 = c(ifelse(logistic(coef)[1] - 1.96 * se[1] < 0, 0, logistic(coef)[1] - 1.96 * se[1]),ifelse(logistic(coef)[level + 1] - 1.96 * se[level + 1] < 0, 0, logistic(coef)[level + 1] - 1.96 * se[level + 1])),

                         max95 = c(ifelse(logistic(coef)[1] + 1.96 * se[1] > 1, 1, logistic(coef)[1] + 1.96 * se[1]),ifelse(logistic(coef)[level + 1] + 1.96 * se[level + 1] > 1, 1, logistic(coef)[level + 1] + 1.96 * se[level + 1])),
                         row.names = c("prevalence",
                                       "One-sayer"))

  return(estimate %>% round(3))
}

fit_freq <- function(nobs, transMat, pihat, theta , only1 = NULL){
  if (is.null(only1)) {
    theta <- 0
    Q <-  diag(rep(c(sum(nobs[1:2])/sum(nobs), sum(nobs[3:4])/sum(nobs)), each = 2), 4) %*% transMat
  }else{
    theta <-  theta
    Q <- transMat
    Q[1, ] <- transMat[1, ] * (1 - theta)
    Q[2, ] <- transMat[2, ] + transMat[1, ] * theta

    Q[3, ] <- transMat[3, ] * (1 - theta)
    Q[4, ] <- transMat[4, ] + transMat[3, ] * theta

    Q <- diag(rep(c(sum(nobs[1:2])/sum(nobs), sum(nobs[3:4])/sum(nobs)), each = 2), 4) %*% Q
  }

  fitted <- sum(nobs) * Q  %*% c(pihat, 1 - pihat)

  return(data.frame(fitted = fitted))
  }
