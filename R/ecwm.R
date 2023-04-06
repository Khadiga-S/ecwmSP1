#' A logistic regression model for the ECWM accounting for One-Sayers
#' @description Fits a logistic regression to the
#'   extended crosswise model that accounts for self-protective One-Sayers,
#'   and predicts the probabilities of the sensitive characteristic and one-saying
#'   conditional on the predictors.
#'
#' @param prevalence a \code{prevalence} in the form
#' \code{response ~ terms}, where \code{response} is the response
#' variable and \code{terms} describes the linear predictors. The
#' response variable should be numeric with the values 1 for
#' "One" and 2 for "Two".
#' @param pyes name of the variable with the individual probabilities
#' of answering "Yes" to the innocuous question.
#' @param data data frame with individual records containing the response
#' and \code{pyes} variables, and the predictors.
#' @param only1 right hand formula for one-saying. The default is NULL means ECWM without one-saying. For the ECWM
#' with one-saying, should be in the form \code{~terms}.
#' @param iter the maximum number of iterations for the derivatives using the
#' method "Nelder-Mead". The Default is 500.
#'
#' @return A list with the following objects:
#' \describe{
#' \item{covsPrev}{data frame for the predictors of \code{prevalence}.}
#' \item{covsONLY}{data frame for the predictors of \code{only1}.}
#' \item{coefs}{data frame with the parameter estimates.}
#' \item{gof}{fit statistics including the log-likelihood, AIC and a likelihood-ratio
#' test against the intercept-only model.}
#' \item{fitted}{vector with the estimated probabilities.}
#' \item{designMat_Prev}{design matrix of \code{prevalence}.}
#' \item{designMat_ONLY}{design matrix of \code{only1}.}
#' \item{freqs}{data frame with the observed and estimated response frequencies.}
#' \item{hessian}{Hessian matrix of the parameter estimates.}
#' }
#' @importFrom stats optim pchisq xtabs model.matrix pt  model.frame get_all_vars
#' @importFrom dplyr mutate case_when %>%
#' @importFrom numDeriv jacobian
#'
#' @export
ecwm <- function(prevalence,  pyes, data, only1 = NULL, iter = 500){

  y <- data[[prevalence[[2]]]]

  y <- 2 - y

  if(!is.numeric(y))         stop("response variable is not of class numeric.")

  if(!all(y %in% 0:1))       stop("response variable has values other than 1 and 2.")

  if(sum(is.na(y)) > 0)      stop("response variable has missing values.")

  p <- data[[substitute(pyes)]]

  if(!is.numeric(p))         stop("pyes is not of class numeric.")

  if(sum(is.na(p)) > 0)      stop("pyes variable has missing values.")

  if(length(unique(p)) > 2)  warning("pyes variable has more than two unique values")

  p2 <- mutate(data.frame(y, p),
               p02 = case_when(p == min(p) & y == 0 ~ min(p),
                               p == max(p) & y == 0 ~ max(p)),
               p11  = case_when(p == min(p) & y == 1 ~ max(p),
                                p == max(p) & y == 1 ~ min(p)),
               .keep = "unused") %>%
    as.matrix()
  ##################
  d <- data.frame(xtabs(~ y + p))

  p1 <- min(p)

  P <- matrix(c(p1,    1 - p1,
                1 - p1,    p1,
                1 - p1,    p1,
                p1,    1 - p1), 4, 2)

  ##################

  D1 <- model.matrix(prevalence, data)

  mx <-  model.matrix(~p)

  if (!is.null(only1)) {

    D2 <-  model.matrix(only1, data)

  }else{
    D2 <- NULL
  }

  sat <- optim(c(0, 0), logl, p2 = p2, D1 = mx, only1 = NULL, control = list(maxit = iter))

  if (is.null(only1) & ncol(D1) == 1) {

    reg0 <- optim(0, logl, p2 = p2, D1 = D1, only1 = NULL, method = "Brent", lower = -20, upper = 20, hessian = TRUE)

    jaco  <- numDeriv::jacobian(logistic, x = reg0$par) # derivatives w.r.t. parameters

    vpar <- t(jaco) %*% solve(reg0$hessian) %*% jaco

    estimate <- data.frame(est   = logistic(reg0$par),

                           se    = sqrt(diag(vpar)),

                           min95 = ifelse(logistic(reg0$par) - 1.96 * sqrt(diag(vpar)) < 0, 0, logistic(reg0$par) - 1.96 * sqrt(diag(vpar))),

                           max95 = ifelse(logistic(reg0$par) + 1.96 * sqrt(diag(vpar)) > 1, 1, logistic(reg0$par) + 1.96 * sqrt(diag(vpar))),

                           row.names = c("MLE: Prevalence"))

    gof <- data.frame(G  = -2 * (sat$value - reg0$value),
                      df = 1,
                      p  = pchisq(-2 * (sat$value - reg0$value), df = 1, lower.tail = FALSE),
                      row.names = "")

    colnames(gof)[1] <- "G-squared"

    fitted <-  round(fit_freq(nobs = d$Freq, transMat = P, theta = 0, pihat = logistic(reg0$par)), 3)

    freqs <- data.frame(d, fitted)

    freqs <- mutate(freqs, y = 3 - as.numeric(y))

    colnames(freqs)[1:3] = c("response", "pyes", "obs")


  } else if (is.null(only1) & ncol(D1) > 1) {

    reg <- optim(rep(0, ncol(D1)), logl, p2 = p2, D1 = D1, only1 = NULL, control = list(maxit = iter))
    reg <- optim(reg$par, logl, p2 = p2, D1 = D1, only1 = NULL, method = "BFGS", hessian = TRUE)

    estimate <- data.frame(Coef.    = c(reg$par),
                           SE.      = c(sqrt(diag(solve(reg$hessian)))),
                           t.value  = c(reg$par)/c(sqrt(diag(solve(reg$hessian)))),
                           p.value  = 2 * pt(-abs(c(reg$par)/c(sqrt(diag(solve(reg$hessian))))), df = nrow(p2) - ncol(D1)),
                           row.names = colnames(D1))

    gof <-  data.frame(loglike = -reg$value %>%  round(3),
                       AIC     = 2 * (reg$value + ncol(D1)) %>%  round(1),
                       LR      = 2 * (optim(0, logl, p2 = p2, D1 = model.matrix(~1, data), only1 = NULL, method = "Brent", lower = -20, upper = 20)$value - reg$value) %>%  round(4),
                       df      = ncol(D1) - 1 ,
                       p       = pchisq(2 * (optim(0, logl, p2 = p2, D1 = model.matrix(~1, data), only1 = NULL, method = "Brent", lower = -20, upper = 20)$value - reg$value), df = ncol(D1) - 1 , lower.tail = FALSE) %>%  round(4),
                       row.names = "")

    probs           <- cbind(apply(D1 %*% c(reg$par), 2, logistic))
    colnames(probs) <- c("Prevalence")

    hessian <- reg$hessian


  } else if (!is.null(only1) & ncol(D1) == 1 & ncol(D2) == 1) {

    reg0 <- optim(matrix(c(0,-2), nrow = 2), logl, p2 = p2, D1  = D1, D2 = D2, only1 = only1, hessian = TRUE)

    jaco  <- numDeriv::jacobian(logistic, x = reg0$par)

    vpar <- t(jaco) %*% solve(reg0$hessian) %*% jaco

    estimate <- data.frame(est   = logistic(reg0$par),

                           se    = sqrt(diag(vpar)),

                           min95 = ifelse(logistic(reg0$par) - 1.96 * sqrt(diag(vpar)) < 0, 0, logistic(reg0$par) - 1.96 * sqrt(diag(vpar))),

                           max95 = ifelse(logistic(reg0$par) + 1.96 * sqrt(diag(vpar)) > 1, 1, logistic(reg0$par) + 1.96 * sqrt(diag(vpar))),

                           row.names = paste("MLE: ", c("Prevalence", "ONE Sayers")))

    gof <- data.frame(G  =  -2 * (reg0$value - optim(0, logl, p2 = p2, D1 = model.matrix(~1, data), only1 = NULL, method = "Brent", lower = -20, upper = 20)$value),
                      df = 1,
                      p  = pchisq(-2 * (reg0$value - optim(0, logl, p2 = p2, D1 = model.matrix(~1, data), only1 = NULL, method = "Brent", lower = -20, upper = 20)$value), df = 1, lower.tail = FALSE),
                      row.names = "")
    colnames(gof)[1] <- "G-squared"

    fitted <-  round(fit_freq(nobs = d$Freq, transMat = P, theta = logistic(reg0$par)[2], pihat = logistic(reg0$par)[1], only1 = 1), 3)

    freqs <- data.frame(d, fitted)

    freqs <- mutate(freqs, y = 3 - as.numeric(y))

    colnames(freqs)[1:3] = c("response", "pyes", "obs")

  }else{
    reg <- optim(matrix(rep(0, ncol(D1) + ncol(D2)), nrow = ncol(D1) + ncol(D2)), logl, p2 = p2, D1  = D1, D2 = D2, only1 = only1, control = list(maxit = iter))
    reg <- optim(reg$par, logl, p2 = p2, D1  = D1, D2 = D2, only1 = only1, method = "BFGS", hessian = TRUE)

    estimate <- data.frame(Coef.    = c(reg$par),

                           SE.      = c(sqrt(diag(solve(reg$hessian)))),

                           t.value  = c(reg$par)/c(sqrt(diag(solve(reg$hessian)))),

                           p.value  = 2 * pt(-abs(c(reg$par)/c(sqrt(diag(solve(reg$hessian))))), df = nrow(p2) - ncol(D1) - ncol(D2)),

                           row.names = c(paste0("Prevalence: ", colnames(D1)), paste0("ONLY ONE: ", colnames(D2))))

    gof <-  data.frame(loglike = -reg$value %>%  round(3),
                       AIC     = 2 * (reg$value + ncol(D1) + ncol(D2)) %>%  round(1),
                       LR      = 2 * (optim(matrix(c(0,-2), nrow = 2), logl, p2 = p2, D1  = model.matrix(~1, data), D2 = model.matrix(~1, data), only1 = only1)$value - reg$value) %>%  round(4),
                       df      = ncol(D1) + ncol(D2) - 2,
                       p       = pchisq(2 * (optim(matrix(c(0,-2), nrow = 2), logl, p2 = p2, D1  = model.matrix(~1, data), D2 = model.matrix(~1, data), only1 = only1)$value - reg$value), df = ncol(D1) + ncol(D2) - 2, lower.tail = FALSE) %>%  round(4),
                       row.names = "")

    probs           <- cbind(apply(D1 %*% matrix(reg$par[1:ncol(D1)], ncol = 1), 2, logistic), apply(D2 %*% matrix(reg$par[ncol(D1) + 1:ncol(D2)], ncol = 1), 2, logistic))

    colnames(probs) <- c("Prevalence", "ONE Sayers")
    hessian <- reg$hessian


  }


  cat(" ", "Call: ", paste(call("ecwm: prevalence", prevalence)), "\n", "\n")
  cat(" ", " ", paste(call("only1 ~", only1)), "\n", "\n")
  cat("Parameter estimates:", "\n")
  print(estimate %>% round(3))
  cat("\n")
  cat("\n")
  cat("Fit statistics")
  cat("\n")
  print(gof %>% round(3))
  cat("\n")
  if (ncol(D1) == 1){
    cat("observed and estimated frequencies \n \n")
    print(freqs)
  }


  invisible(list(covsPrev     = get_all_vars(prevalence, data),
                 covsONLY     = data.frame(ifelse(is.null(only1), "NULL", list(get_all_vars(only1, data)))),
                 coefs         = estimate,
                 gof          = gof,
                 fitted       = data.frame(fitted=ifelse(ncol(D1) == 1, estimate[,1], list(probs))),
                 designMat_Prev  = D1,
                 designMat_ONLY  = D2,
                 freqs = data.frame(freqs = ifelse(ncol(D1) == 1, freqs, data.frame(d))),
                 hessian = data.frame((var = ifelse(ncol(D1) > 1, list(hessian),list(reg0$hessian))) )

  )
  )

}
