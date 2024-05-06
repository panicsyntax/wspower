
#' @title Compute power for a One Factor ANOVA with n levels.
#' @description Takes means, sds, and sample sizes for each group.
#' Alpha is .05 by default, alternative values may be entered by user
#' @param means sample means of the n samples
#' @param sds sample standard deviations of the n samples
#' @param ns sample sizes of each sample
#' @param cs contrasts
#' @return The power of the contrasts and the overall F Test

MY_anova1f_nc <-
  function (means = c(0,1,2,3),
            sds = c(1,1,1,1),
            ns = c(10,5,10,20),
            alpha = 0.05,
            cs = c(-1, -1, 1, 1)) {
    K <- length(means)
    df <- data.frame(y = rep(NA, sum(ns)), group = factor(rep(LETTERS[1:K],ns)))
    dfy <- NULL
    for (k in 1:K){
      x <- stats::rnorm(ns[k], means[k], sds[k])
      X <- x
      MEAN <- means[k]
      SD <- sds[k]
      Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) *
        SD
      y <- MEAN + Z
      dfy <- c(dfy, y)
    }
    df$y <- dfy
    anova <- stats::aov(y ~ group, data = df)
    anova <- car::Anova(anova, type = "III")
    SSA <- anova[2, 1]
    SSwin <- anova[3, 1]
    dfwin <- anova[3, 2]
    mswin <- SSwin/dfwin
    dfbg <- anova[2, 2]
    eta2 <- SSA/(SSA + SSwin)
    f2 <- eta2/(1 - eta2)
    lambda <- f2 * dfwin
    minusalpha <- 1 - alpha
    Ft <- stats::qf(minusalpha, dfbg, dfwin)
    power <- 1 - stats::pf(Ft, dfbg, dfwin, lambda)
    delta <- (means %*% cs)/sqrt(mswin * sum(cs^2/ns))
    lambda.c = delta^2
    Ft.c <- stats::qf(minusalpha, 1, dfwin)
    power.contrast <- round(1 - stats::pf(Ft.c, 1, dfwin, lambda.c),
                            3)
    result <- cbind(ns, means, sds, cs)
    colnames(result) <- c("n", "mean", "sd", "contrast.coef")
    cat("Scenario:\n")
    prmatrix(result)
    message("Power against this contrast = ", power.contrast)
    message("Power of overall F test = ", power)
    invisible(list(scenario = result, power = power.contrast))
  }

#
#' @title Compute power for a One Factor ANOVA with four levels between three groups.
#' @description Takes means and sds for each group
#' assumes equal sample sizes for each group.
#' Alpha is .05 by default, alternative values may be entered by user
#' Can specify cohen's distance
#' Experiment-wise Significance Level
#'
#' @param Cohen.d The desired Cohen's Distance. Can be a vector
#' @param n desired sample size for each group. Can be a vector
#' @param sig.level experiment-wise significance level
#' @return the power for each desired cohen's d and n

ThreeGroupsMeansPower <- function(Cohen.d = c(.2, .5, .8, 1),
                               n = c(125, 139, 150, 160, 161, 165, 170, 184),
                               sig.level = 0.0015){
  myPowFun <- function(d, n){
    ans <-anova1f_3(m1 = d, m2 = 0, m3 = 0,
                    s1 = 1, s2 = 1, s3 = 1,
                    n1 = n, n2 = n, n3 = n,
                    alpha = sig.level)
    ans$Power
  }
  myVecPowFun <- Vectorize(myPowFun, vectorize.args = c("d", "n"))
  ThreeGroups <- outer(Cohen.d, n, FUN = myVecPowFun)
  dimnames(ThreeGroups) <- list(Cohen.d, n)
  ThreeGroups
}

#
#' @title Compute power for a One Factor ANOVA with four levels between three groups.
#' @description Takes means and sds for each group
#' assumes equal sample sizes for each group.
#' Alpha is .05 by default, alternative values may be entered by user
#' Can specify cohen's distance
#' The cohen's d for group 2 is half that of group 1
#' Experiment-wise Significance Level
#' @inheritParams ThreeGroupsMeansPower
#' @return the power for each desired cohen's d and n

ThreeGroupsMeansPower.stepped <- function(Cohen.d = c(.2, .5, .8, 1),
                                     n = c(125, 139, 150, 160, 161, 165, 170, 184),
                                     sig.level = 0.0015){
  myPowFun <- function(d, n){
    ans <-anova1f_3(m1 = d, m2 = d/2, m3 = 0,
                    s1 = 1, s2 = 1, s3 = 1,
                    n1 = n, n2 = n, n3 = n,
                    alpha = sig.level)
    ans$Power
  }
  myVecPowFun <- Vectorize(myPowFun, vectorize.args = c("d", "n"))
  ThreeGroups <- outer(Cohen.d, n, FUN = myVecPowFun)
  dimnames(ThreeGroups) <- list(Cohen.d, n)
  ThreeGroups
}

#
#' @title Compute power for a One Factor ANOVA with four levels between three groups.
#' @description Takes means and sds for each group
#' assumes equal sample sizes for each group.
#' Alpha is .05 by default, alternative values may be entered by user
#' Can specify cohen's distance
#' Contrasts can be specified for each group
#' Experiment-wise Significance Level
#' @inheritParams ThreeGroupsMeansPower
#' @param c1-c3 contrast for groups 1-3
#' @return the power for each desired cohen's d and n
#'
ThreeGroupsMeansPower.contrast <- function(Cohen.d = c(.2, .5, .8, 1),
                               n = c(125, 139, 150, 160, 161, 165, 170, 184),
                               sig.level = 0.0015,
                               c1=2, c2=-1, c3=-1){
  myPowFun <- function(d, n){
    ans <-anova1f_3c(m1 = d, m2 = 0, m3 = 0,
                     s1 = 1, s2 = 1, s3 = 1,
                     n1 = n, n2 = n, n3 = n,
                     alpha = sig.level,
                     c1 = c1, c2 = c2, c3 = c3) # contrast
    ans$Power
  }
  myVecPowFun <- Vectorize(myPowFun, vectorize.args = c("d", "n"))
  ThreeGroups <- outer(Cohen.d, n, FUN = myVecPowFun)
  dimnames(ThreeGroups) <- list(Cohen.d, n)
  ThreeGroups
}

#' @title Compute power for correlations among n groups.
#' @description The following power.3cor method uses the Fisher Z transform of the
#' correlations and a SS-type statistic having approximately a chisqure
#' distribution under H0 and a non-central chisquare distribution under H1.
#' This tests for ANY difference among the correlations. H0: r1 = r2 = r3
#' @param r vector sample correlation coefficients
#' @param n vector of corresponding sample sizes
#' @param sig.level significance level. Default is .05
#' @param groupNames optional: can specify the names of the groups
#' @return The power to detect ANY difference among 3 correlation coefficients
power.3cor <- function(r, n,
                       sig.level = 0.05,
                       groupNames = NULL){
  K <- length(r)
  if (length(n) != K) stop("r and n must be the same length")
  z <- atanh(r) # Fisher Z transform of r
  varz <- 1/(n-3)
  # commonmean <- weighted.mean(z, n/sum(n))
  commonmean <- weighted.mean(z, (1/varz)/sum(1/varz))
  W <- sum((z-commonmean)^2/varz)
  # W ~ chisq on K-1 df under H0
  # W ~ chisq on K-1 df with ncp ncp under H1
  ncp <- W
  cv <- qchisq(1-sig.level/2, K-1)
  power <- 1-pchisq(cv, K-1, ncp = ncp)
  out <- cbind(r, n)
  if (!is.null(groupNames)) rownames(out) <- groupNames else
    rownames(out) <- paste("Group", 1:K)
  colnames(out) <- c("Correlation", "n")
  cat("Scenario:\n")
  prmatrix(out)
  cat("The power to detect ANY difference is:", power)
  invisible(list(power = power, scenario = out))
}

#' @title Power for correlations among n groups.
#' @description The following power.3cor.sim uses via regression to test for
#' ANY difference in correlations
#' @param Nsim Number of simulations. Default is 1000
#' @inheritParams power.3cor
#' @return The power to detect ANY difference among 3 correlation coefficients with linear regression and ANOVA
power.3cor.sim <- function(Nsim = 1000, r, n, sig.level = 0.05,
                           groupNames = NULL){
  require(mvtnorm)
  require(lmtest)
  K <- length(r)
  if (length(n) != K) stop("r and n must be the same length")
  reject.LR <- reject.F <- vector("numeric", Nsim)
  dat <- matrix(NA, nrow = sum(n), ncol = 3)
  index.start <- c(0, cumsum(n[1:(K-1)])) + 1
  index.end <- index.start + n - 1
  for (i in 1:Nsim){
    for (k in 1:K){
      samp <- rmvnorm(n[k], c(0,0), matrix(c(1, r[k], r[k], 1), nrow = 2, byrow = TRUE))
      dat[index.start[k]:index.end[k],1:2] <- samp
    }
    df <- as.data.frame(dat)
    names(df) <- c("Y", "X", "source")
    df$source <- as.factor(rep(1:K, n))
    lmfit.full <- lm(Y ~ X*source, data = df)
    lmfit.reduced <- lm(Y ~ X + source, data = df)
    # LR test
    testinfo <- lrtest(lmfit.full, lmfit.reduced)
    pval <- testinfo$`Pr(>Chisq)`[2]
    reject.LR[i] <- ifelse(pval < sig.level, 1, 0)
    # ANOVA F test
    testinfo <- anova(lmfit.full, lmfit.reduced)
    pval <- testinfo$`Pr(>F)`[2]
    reject.F[i] <- ifelse(pval < sig.level, 1, 0)
  }
  power.LR <- mean(reject.LR)
  power.F <- mean(reject.F)
  out <- cbind(r, n)
  if (!is.null(groupNames)) rownames(out) <- groupNames else
    rownames(out) <- paste("Group", 1:3)
  colnames(out) <- c("Correlation", "n")
  cat("Scenario:\n")
  prmatrix(out)
  cat("The power by LR to detect ANY difference is:", power.LR, "\n")
  cat("The power by ANOVA to detect ANY difference is:", power.F, "\n")
  invisible(list(power.LR = power.LR, power.F = power.F, scenario = out))
}

#' @title Test linear combinations of correlations
#' @description The following power.3core.sim.c uses simulation via regression to
#' test H0: r %*% ccoef = 0, i.e. to test null hypotheses that
#' specify specific linear combinations of correlations are zero.
#' @inheritParams power.3cor.sim
#' @param ccoef contrast coefficients. can be a vector
power.3cor.sim.c <- function(Nsim = 1000, r, n, ccoef, sig.level = 0.05,
                             groupNames = NULL){
  require(mvtnorm)
  require(gmodels)
  K <- length(r)
  ccoef <- if (is.matrix(ccoef)) ccoef else matrix(ccoef, ncol = K, byrow = TRUE)
  if ((length(n) != K) || (ncol(ccoef) != K)) stop("r, n and ccoef must be compatible")
  if (!all(rowSums(ccoef) == 0)) stop("contrast coefficients must some to 0")
  ncontrasts <- nrow(ccoef)
  coefmat <- matrix(c(rep(0, ncontrasts + ncontrasts + ncontrasts*(K-1)),
                      t(ccoef[,2:K])), nrow = ncontrasts, byrow = FALSE)
  reject.F <- vector("numeric", Nsim)
  dat <- matrix(NA, nrow = sum(n), ncol = 2)
  index.start <- c(0, cumsum(n[1:(K-1)])) + 1
  index.end <- index.start + n - 1
  for (i in 1:Nsim){
    for (k in 1:K){
      samp <- rmvnorm(n[k], c(0,0), matrix(c(1, r[k], r[k], 1), nrow = 2, byrow = TRUE))
      dat[index.start[k]:index.end[k],] <- samp
    }
    df <- as.data.frame(dat)
    names(df) <- c("Y", "X")
    if (is.null(groupNames))
      df$source <- as.factor(rep(1:K, n)) else
        df$source <- as.factor(rep(groupNames, n))
    lmfit.full <- lm(Y ~ X*source, data = df)
    # F test
    testinfo <- glh.test(lmfit.full, coefmat)
    pval <- testinfo$p.value
    reject.F[i] <- ifelse(pval < sig.level, 1, 0)
  }
  power.F <- mean(reject.F)
  out <- cbind(r, n, t(ccoef))
  if (!is.null(groupNames)) rownames(out) <- groupNames else
    rownames(out) <- paste("Group", 1:3)
  colnames(out) <- c("Correlation", "n", paste("Contrast", 1:ncontrasts))
  cat("\nScenario:\n")
  prmatrix(out)
  if (ncontrasts == 1)
    cat("The power by F test against this contrast is:", power.F, "\n")
  else
    cat("The power by F test against these contrasts is:", power.F, "\n")
  invisible(list(power.F = power.F, scenario = out))
}

#' @title Test Moderation effects
#' @description The following power.3core.sim.c.M uses simulation via regression to
#' test H0:
#' prM1 = Pr(M = 1) for each group
#' Y in the code is CVD risk; X is WS score.
#' H0: there is no moderation effect
#' H0: a contrast of the moderation effects is zero
#' @inheritParams power.3cor.sim.c
#' @param prM1 Probability that the Moderator is 1 for each group. Assumed equal among groups
#' @param bWS description
#' @return The Power by F Test

power.3cor.sim.c.M <- function(Nsim = 1000, r, n, ccoef,
                               prM1,
                               bWS = NULL,
                               sig.level = 0.05,
                               groupNames = NULL){
  require(mvtnorm)
  require(gmodels)
  K <- length(r)
  ccoef <- if (is.matrix(ccoef)) ccoef else matrix(ccoef, ncol = K, byrow = TRUE)
  if ((length(n) != K) || (ncol(ccoef) != K)) stop("r, n and ccoef must be compatible")
  if (!all(rowSums(ccoef) == 0)) stop("contrast coefficients must some to 0")
  ncontrasts <- nrow(ccoef)
  coefmat <- matrix(c(rep(0, ncontrasts + ncontrasts + ncontrasts*(K-1)),
                      t(ccoef[,2:K])), nrow = ncontrasts, byrow = FALSE)
  if (is.null(bWS)) bWS <- rep(0, K)
  reject.F <- vector("numeric", Nsim)
  dat <- matrix(NA, nrow = sum(n), ncol = 3)
  index.start <- c(0, cumsum(n[1:(K-1)])) + 1
  index.end <- index.start + n - 1
  for (i in 1:Nsim){
    for (k in 1:K){
      samp <- rmvnorm(n[k], c(0,0), matrix(c(1, r[k], r[k], 1), nrow = 2, byrow = TRUE))
      colnames(samp) <- c("Y", "X")
      logitprM1 <- prM1[k]/(1 - prM1[k]) + bWS[k] * samp[,"X"]
      M <- rbinom(n[k], size = 1, prob = plogis(logitprM1))
      dat[index.start[k]:index.end[k],] <- cbind(samp, M)
    }
    df <- as.data.frame(dat)
    names(df) <- c("Y", "X", "M")
    if (is.null(groupNames))
      df$source <- as.factor(rep(1:K, n)) else
        df$source <- as.factor(rep(groupNames, n))
    lmfit.full <- lm(Y ~ X*source*M, data = df)
    # coefs: Int, X, source2 - sourceK, M, X*source2 - X*sourceK, X*M,
    #        M*source2 - M*sourceK, X*M*source2 - X*M*sourceK
    # length is  2+(K-1) + 1 + (K-1) + 1 + (K-1) + (K-1) = 4+4(K-1) = 4K
    #
    print(summary(lmfit.full))
    # F test
    testinfo <- glh.test(lmfit.full, coefmat)
    pval <- testinfo$p.value
    reject.F[i] <- ifelse(pval < sig.level, 1, 0)
  }
  power.F <- mean(reject.F)
  out <- cbind(r, n, t(ccoef))
  if (!is.null(groupNames)) rownames(out) <- groupNames else
    rownames(out) <- paste("Group", 1:3)
  colnames(out) <- c("Correlation", "n", paste("Contrast", 1:ncontrasts))
  cat("\nScenario:\n")
  prmatrix(out)
  if (ncontrasts == 1)
    cat("The power by F test against this contrast is:", power.F, "\n")
  else
    cat("The power by F test against these contrasts is:", power.F, "\n")
  invisible(list(power.F = power.F, scenario = out))
}

#' @title Returns the power of a moderation effect given N total sample size
#' @param N Total Sample Size
#' @param pt0 proportion in group 0
#' @param pm0 proportion with binary moderator = 0, assumed the same in each of group 0 and group 1
#' @param m   standardized size of the moderation effect: For the normally distribution response with sd sigma, if the population mean
#' outcome for (T = t, Moderator = k) is mu[t,k], then the moderation
#' effect is M = (mu[1,1] - mu[0,1] - (mu[1,0] - mu[0,0])) in the units
#' of the response, and the standardized moderation effect
#' is m = M/sigma.
#' @return The power of a moderation effect given N total sample size
power.moderation <- function(N, pt0 = 1/2, pm0 = 1/2,
                             m = 0.5, sig.level = 0.05)
{
  N0 <- ceiling(N*pt0)
  N1 <- N - N0
  n00 <- ceiling(N0*pm0)
  n10 <- ceiling(N1*pm0)
  n01 <- N0 - n00
  n11 <- N1 - n10

  var.m <- sum(1/c(n00, n10, n01, n11))
  sd.m <- sqrt(var.m)

  cv <- qnorm(1-sig.level/2)
  1-pnorm(cv - m/sd.m) + pnorm(-cv - m/sd.m)

}

#' Two-group comparison, continuous outcome, binary moderator.
#' n = 160/experimental group; moderator balanced.
#' @param std.trt.effect standardized treatment effect. Default is .5
#' @param std.mod.effect standardized moderation effect. Default is 0
#' @param nt0m0 sample size of trt=0 moderation=0
#' @param nt0m1 sample size of trt=0 moderation=1
#' @param nt1m0 sample size of trt=1 moderation=0
#' @param nt1m1 sample size of trt=1 moderation=1
#' @param alpha experiment-wise significance level
#' @param ct0m0-ct1m1 contrasts for each group
#' @param st0m0-st1m1 standard deviations for each group
#' @return Returns the power


moderator.pow <- function(std.trt.effect = 0.5,
                          std.mod.effect = 0,
                          nt0m0 = 80, nt1m0 = 80, nt0m1 = 80, nt1m1 = 80, #first index is trt, second is moderator
                          alpha = 0.0017,
                          ct0m0 = -1, ct1m0 = 1, ct0m1 = 1, ct1m1 = -1, #contrast coefs
                          st0m0 = 1, st1m0 = 1, st0m1 = 1, st1m1 = 1 #standard devs
){
  mt0m0 <- 0; mt1m0 <- mt0m0 + std.trt.effect
  mt0m1 <- 0; mt1m1 <- mt0m1 + std.trt.effect + std.mod.effect
  anova1f_4c(m1 = mt0m0, m2 = mt1m0, m3 = mt0m1, m4 = mt1m1,
             s1 = st0m0, s2 = st1m0, s3 = st0m1, s4 = st1m1,
             n1 = nt0m0, n2 = nt1m0, n3 = nt0m1, n4 = nt1m1,
             alpha = alpha,
             c1 = ct0m0, c2 = ct1m0, c3 = ct0m1, c4 = ct1m1)
}

#' @title Generates moderation data
#' @param n sample size per arm
#' @return Matrix of Trt Effect, Moderator Effect, Trt by Moderator Effect
#' @return Variance of the moderator coefficient estimate
#' @return Correlation of Response and Moderator for each Treatment Group
#' @return Parameter Estimates for Each Variable

genModData <-
  function(n = 8# n per arm
  ){
    # two arms
    Trt <- rep(c(1, 0), c(n, n))
    # Mod <- ifelse(runif(2*n) > .5,
    #               rnorm(2*n, 0, 1),
    #               rnorm(2*n, 2, 1))
    Mod <- rnorm(2*n, 0, 2) # standard normal
    TbyM <- Trt * Mod
    X <- cbind(1, Trt, Mod, TbyM)
    dm0 <- density(Mod[Trt == 0])
    plot(dm0)
    dm1 <- density(Mod[Trt == 1])
    plot(dm1)
    plot(range(c(dm0$x, dm1$x)),
         range(c(dm0$y, dm1$y)),
         type = "n", xlab = "moderator",
         ylab = "Density")
    lines(dm0, lty = 2, col = "black")
    lines(dm1, lty = 2, col = "red")
    XtXinv <- solve(t(X) %*% X)
    print(XtXinv)
    cat("\nVariance of the moderator coefficient estimate: ",
        XtXinv[4,4], "\n")
    Y <- rnorm(2*n, mean = X %*% c(0, 0, 0, .5))
    print(cor(Y[Trt == 0], Mod[Trt == 0]))
    print(cor(Y[Trt == 1], Mod[Trt == 1]))
    fit <- lm(Y ~ X)
    fit
  }

#' @title Powers for WD
#' @param siglevel experiment-wise significance level
#' @param alphamoderation alpha of moderation effect
#' @param nWD number in trt
#' @param nNegCntrl number in negative control
#' @param nNeuCntrl number in neutral control
#' @param effectSizes effect sizes
#' @param m standardized moderation effect
#' @return Matrix of Power for different Scenarios
WDpowers <- function(
    siglevel = 0.05/10, alphamoderation = 0.05,
    nWD = 160, nNegCntrl = 160, nNeuCntrl = 160,
    # effectSizes contains Cohen effect for mean difference between WD and combined control groups
    effectSizes = c(.2, .5, .8, 1.0),
    m = .8){
  ans <- matrix(NA, nrow = length(effectSizes), ncol = 12)
  colnames(ans) <- c("nWD", "nNegCntrl", "nNeuCntrl", "d", "pow1", "pow2", "pow3", "pow4",
                     "pow5", "m", "pow6", "pow7")
  i <- 1
  for (E in effectSizes){
    # power for WD vs combined control groups
    pow1 <- anova1f_3c(m1 = E, m2 = 0, m3 = 0,
                       s1 = 1, s2 =  1, s3 = 1,
                       n1 = nWD, n2 = nNegCntrl, n3 = nNeuCntrl,
                       alpha = siglevel,
                       c1 = 1, c2 = -.5, c3 = -.5)
    # power for WD vs combined control groups when the control groups differ by E/2
    pow2 <- anova1f_3c(m1 = E, m2 = E/2, m3 = 0,
                       s1 = 1, s2 =  1, s3 = 1,
                       n1 = nWD, n2 = nNegCntrl, n3 = nNeuCntrl,
                       alpha = siglevel,
                       c1 = 1, c2 = -.5, c3 = -.5)
    # power for difference between the two control groups of size E
    pow3 <- anova1f_3c(m1 = E, m2 = E, m3 = 0,
                       s1 = 1, s2 =  1, s3 = 1,
                       n1 = nWD, n2 = nNegCntrl, n3 = nNeuCntrl,
                       alpha = siglevel,
                       c1 = 0, c2 = 1, c3 = -1)
    # power for difference between the two control groups of size E/2
    pow4 <- anova1f_3c(m1 = E, m2 = E/2, m3 = 0,
                       s1 = 1, s2 =  1, s3 = 1,
                       n1 = nWD, n2 = nNegCntrl, n3 = nNeuCntrl,
                       alpha = siglevel,
                       c1 = 0, c2 = 1, c3 = -1)
    # power for trend when there is a trend
    pow5 <- anova1f_3c(m1 = E, m2 = E/2, m3 = 0,
                       s1 = 1, s2 = 1, s3 = 1,
                       n1 = nWD, n2 = nNegCntrl, n3 = nNeuCntrl,
                       alpha = siglevel,
                       c1 = 1, c2 = 0, c3 = -1)
    # power for moderation effect of size m*E in 2-group study with balanced moderator
    pow6 <- moderator.pow.cont(std.trt.effect = E, std.mod.effect = m*E,
                               pt0 = 0.5, pm0 = 0.50,
                               Ntotal = nWD + nNegCntrl,
                               alpha = alphamoderation,
                               ct0m0 = 1, ct1m0 = -1, ct0m1 = -1, ct1m1 = 1)
    # power for moderation effect of size m*E in 2-group study with unbalanced moderator 40-60
    pow7 <- moderator.pow.cont(std.trt.effect = E, std.mod.effect = m*E,
                               pt0 = 0.5, pm0 = 0.40,
                               Ntotal = nWD + nNegCntrl,
                               alpha = alphamoderation,
                               ct0m0 = 1, ct1m0 = -1, ct0m1 = -1, ct1m1 = 1)
    ans[i,] <- c(nWD, nNegCntrl, nNeuCntrl,
                 E,
                 pow1$Power, pow2$Power, pow3$Power, pow4$Power, pow5$Power,
                 m, pow6$Power, pow7$Power)
    i <- i + 1
  }
  ans
}

#' @title Binary outcome of social withdrawal with 2 groups. Studies 2 and 3.
#' @param RR relative risk. Can be a vector
#' @param p0 probability of group 0
#' @param power Power. Can  be a vector
#' @param sig.level experiment-wise significance level
#' @return Sample Size per Group

TwoGroupsProps<- function(RR = c(1.4, 1.5, 2.0),
                          p0 = c(.5, .55, .6, .65, .7),
                          power = c(.90, .85, .80),
                          sig.level = 0.0015){
  myPowFun <- function(p0,p1,sig.level,pow){
    ans <- power.prop.test(p1 = p0, p2 = p1, sig.level = sig.level, power = pow)
    ceiling(ans$n)
  }
  myVecPowFun <- Vectorize(myPowFun, vectorize.args = c("p0", "p1", "pow"))

  ans <- expand.grid(RR, p0, power)
  names(ans) <- c("RR", "p0", "power")
  ans$p1 <- pmin(ans$RR*ans$p0, 1)
  n <- myVecPowFun(p0 = ans$p0, p1 = ans$p1, pow = ans$power, sig.level = sig.level)
  ans$npergroup <- n
  ans
}

#' @title The power for a continuous outcome with a balanced moderator and two treatment groups
#' @param std.trt.effect Standardized Treatment Effect. Default = .5
#' @param std.mod.effect Standardized Moderation Effect. Default = 0
#' @param pt0 Probability that treatment is 0. Default = .5
#' @param pm0 Probability that moderation is 0. Default = .5
#' @param Ntotal Total Sample Size. Default = 320
#' @param alpha experiment-wise significance level. Default = .0017
#' @param ct0m0-ct1m1 contrast coefficients for each group
#' @param st0m0-stm1m1 standard deviations for each group
#' @return Power

moderator.pow.cont <- function(std.trt.effect = 0.5,
                               std.mod.effect = 0,
                               pt0 = 0.5, #P(T = 0)
                               pm0 = 0.5, #P(M = 0) = P(M = 0 | T = 0) = P(M = 0 | T = 1)
                               Ntotal = 4*80,
                               alpha = 0.0017,
                               ct0m0 = 1, ct1m0 = -1, ct0m1 = 1, ct1m1 = -1, #contrast coefs
                               st0m0 = 1, st1m0 = 1, st0m1 = 1, st1m1 = 1 #standard devs
){
  pt0m0 <- pt0*pm0; pt1m0 <- pm0*(1-pt0); pt0m1 <- pt0*(1-pm0); pt1m1 <- (1-pt0)*(1-pm0)
  nt0m0 <- pt0m0*Ntotal; nt1m0 <- pt1m0*Ntotal; nt0m1 <- pt0m1*Ntotal; nt1m1 <- pt1m1*Ntotal
  mt0m0 <- 0; mt1m0 <- mt0m0 + std.trt.effect
  mt0m1 <- 0; mt1m1 <- mt0m1 + std.trt.effect + std.mod.effect
  anova1f_4c(m1 = mt0m0, m2 = mt1m0, m3 = mt0m1, m4 = mt1m1,
             s1 = st0m0, s2 = st1m0, s3 = st0m1, s4 = st1m1,
             n1 = nt0m0, n2 = nt1m0, n3 = nt0m1, n4 = nt1m1,
             alpha = alpha,
             c1 = ct0m0, c2 = ct1m0, c3 = ct0m1, c4 = ct1m1)
}


#' @title The power for a continuous moderator and two treatment groups allowing for unbalanced moderator; continuous outcome; 3 groups
#' @param std.trt.effect Standardized Treatment Effect. Default = .5
#' @param std.mod.effect Standardized Moderation Effect. Default = 0
#' @param pt0 Probability that treatment is 0. Default = 1/3
#' @param pt1 Probability that treatment is 1. Default = 1/3
#' @param pm0 Probability that moderation is 0. Default = .5
#' @param Ntotal Total Sample Size. Default = 480
#' @param alpha significance level. Default = .05
#' @param ct0m0-ct1m1 contrast coefficients for each group
#' @param st0m0-stm1m1 standard deviations for each group
#' @return Power

moderator.pow.cont.unbalanced <- function(std.trt.effect = 0,
                                 std.mod.effect = 0.8,
                                 pt0 = 1/3, #P(T = 0)
                                 pt1 = 1/3, #P(T = 1)
                                 #P(M = 0) = P(M = 0 | T = k), k = 0, 1, 2:
                                 pm0 = 0.5,
                                 Ntotal = 3*160,
                                 alpha = 0.05,
                                 #contrast coefs:
                                 ct0m0 = 1, ct1m0 = -1, ct2m0 = 0,
                                 ct0m1 = 1, ct1m1 = -1, ct2m1 = 0,
                                 #standard devs
                                 st0m0 = 1, st1m0 = 1, st2m0 = 1,
                                 st0m1 = 1, st1m1 = 1, st2m1 = 1
){
  require(pwr2ppl)
  pt0m0 <- pt0*pm0; pt1m0 <- pm0*pt1; pt2m0 <- pm0*(1-pt0 - pt1)
  pt0m1 <- pt0*(1-pm0); pt1m1 <- pt1*(1-pm0); pt2m1 <- (1-pt0-pt1)*(1-pm0)
  nt0m0 <- pt0m0*Ntotal; nt1m0 <- pt1m0*Ntotal; nt2m0 <- pt2m0*Ntotal;
  nt0m1 <- pt0m1*Ntotal; nt1m1 <- pt1m1*Ntotal; nt2m1 <- pt2m1*Ntotal
  mt0m0 <- 0; mt1m0 <- mt0m0 + std.trt.effect/2; mt2m0 <- mt0m0 + std.trt.effect
  mt0m1 <- 0; mt1m1 <- mt0m1 + std.trt.effect/2 + std.mod.effect/2
  mt2m1 <- mt0m1 + std.trt.effect + std.mod.effect

  ans <- MY_anova1f_c(means = c(mt0m0, mt1m0, mt2m0, mt0m1, mt1m1, mt2m1),
                      sds = c(st0m0, st1m0, st2m0, st0m1, st1m1, st2m1),
                      ns = c(nt0m0, nt1m0, nt2m0, nt0m1, nt1m1, nt2m1),
                      alpha, cs = c(ct0m0, ct1m0, ct2m0, ct0m1, ct1m1, ct2m1)
  )
}

#' Two-group comparison, binary outcome, binary moderator.
#' n = 160/experimental group; moderator balanced.
#' dichotomous moderator effect for 2 groups, binary outcome:
#' use a chisquare goodness of fit test
#' The power for a continuous moderator and two treatment groups allowing for unbalanced moderator; continuous outcome; 3 groups
#' @param trt.effect Additive Treatment Effect. Default = .2
#' @param mod.effect Additive Moderation Effect. Default = 0
#' @param mod.main Additive Main Effect of Moderator. Default = .01
#' @param pt0m0 Response Rate given Trt=0, Moderator = 0
#' @param pt0 Probability of assignment to treatment 0. Default = .5
#' @param pm0 Probability that moderation is 0 in both treatment groups. Default = .5
#' @param Ntotal Total Sample Size. Default = 480
#' @param alpha significance level. Default = .05
#' @return Power

moderator.pow.binary <- function(trt.effect = 0.20,
                                 mod.effect = 0,
                                 mod.main = 0.01,
                                 pt0m0 = .70, #response rate given trt = 0, moderator = 0
                                 pm0 = 0.50, #prob. moderator = 0 (in both trt groups)
                                 pt0 = 0.50, #prob. of assignment to trt = 0
                                 alpha = 0.05,
                                 Ntotal = 4*80
){
  # H0: no moderation effect
  pt1m0 <- pt0m0 + trt.effect
  pt0m1 <- pt0m0 + mod.main
  pt1m1 <- pt0m1 + trt.effect
  Ps0 <- c(pt0m0, pt0m1, pt1m0, pt1m1) #conditional probabilities of response given trt and moderator
  P0 <- rbind(Ps0, 1-Ps0) #matrix of probs given trt and moderator
  P0 <- P0 * matrix(c(pm0, 1-pm0), ncol = 4, nrow = 2, byrow = TRUE)#joint prob of response and moderator given trt
  P0 <- P0 * matrix(c(pt0, pt0, 1-pt0, 1-pt0), ncol = 4, nrow = 2, byrow = TRUE)#P(Y,T,M)
  E0 <- P0*Ntotal

  # H1 with moderation effect s.t. treatment effect when M=1 is augmented by mod.effect
  pt1m1 <- pt0m1 + trt.effect + mod.effect
  Ps1 <- c(pt0m0, pt0m1, pt1m0, pt1m1)
  P1 <- rbind(Ps1, (1-Ps1))#P(Y | T, M)
  P1 <- P1 * matrix(c(pm0, 1-pm0), ncol = 4, nrow = 2, byrow = TRUE)
  P1 <- P1 * matrix(c(pt0, 1-pt0), ncol = 4, nrow = 2, byrow = FALSE)
  E1 <- P1*Ntotal

  dimnames(P0) <- dimnames(P1) <- dimnames(E0) <- dimnames(E1) <-
    list(c("respond", "not respond"), c("T0M0", "T0M1", "T1M0", "T1M1"))
  names(Ps0) <- colnames(P0)
  cat("\nH0 is no moderation effect.\n")
  cat(" Conditional probabilities of response are:\n")
  print(Ps0)
  cat("\n Marginal probabilities are:\n")
  prmatrix(P0)
  cat("\n Expected counts are:\n")
  prmatrix(E0)

  names(Ps1) <- colnames(P1)
  cat("\n\nH1 is moderation effect of ", mod.effect, "\n")
  cat(" Conditional probabilities of response are:\n")
  print(Ps1)
  cat("\n Marginal probabilities are:\n")
  prmatrix(P1)
  cat("\n Expected counts are:\n")
  prmatrix(E1)

  w <- ES.w1(P0, P1)
  # cat("\nEffect size is w = ", w, "\n")
  #
  # myX2 <- sum((E1 - E0)^2/E0)
  # cat("\nChi-square stat is X2 = ", myX2, "\n")

  ans <- pwr.chisq.test(w = w, N = Ntotal, df = 3, sig.level = alpha, power = NULL)
  ans$conditionalH0 <- Ps0
  ans$marginalH0 <- P0
  ans$E0 <- E0
  ans$conditionalH1 <- Ps1
  ans$marginalH1 <- P1
  ans$E1 <- E1
  ans
}

#' @title Test for trend in proportions among k groups
#' @param pi proportions for k groups. Vector
#' @param alpha significance level
#' @param power Power
#' @return Sample size per group
CAsamplesize <- function(pi, alpha = 0.05, power = 0.80){
  k <- length(pi)
  di <- seq(0,k-1, by = 1)
  ci <- (0:(k-1)) - (k-1)/2
  D <- sum(ci*pi)
  za <- qnorm(1-alpha)
  zb <- qnorm(power)
  p <- mean(pi)
  q <- mean(1-pi)
  x <- sum(ci^2*pi*(1-pi))
  nstar <- (za *sqrt(k*(k^2-1)*p*q/12) + zb*sqrt(x))^2/D^2
  n <- (nstar/4)*(1 + sqrt(1 + 2/D/nstar))^2
  n
}

#' @title Computing Sample Size for a given power and effect size in Weight Stigma Study using Mean with WS
#' @description Main effect of WS on LE8 at baseline.
#' assuming LE8 is normally distributed with sd = 15 (from Lloyd-Jones et al.)
#' prWS = Pr(WS = 1), i.e. prob. a subject reports experiencing WS
#' meanNoWS = mean of LE8 when WS = 0
#' Cohen.d = [(mean LE8 with WS) - (mean LE8 without WS)]/sd
#' @param sd standard deviation. Assumed to be equal among groups
#' @param Cohen.d Cohen's Distance. Can be a vector
#' @param prWS Proportion with Weight Stigma Variable
#' @param power Power. Default = .8
#' @param alpha Significance Level. Default = .05
#' @param meanNoWS mean without Weight Stigma Variable. Default = 60
#' @param writeExcel Default=True. To create an Excel file
#' @param file Default="Baseline" name of output file
#' @return a dataframe that includes the sample size for the specified power. Also excel file if specified.

ss <- function(sd = 15, Cohen.d = c(.2, .3, .4, .5),
               prWS = c(.5, .6, .7, .8), power = .80,
               alpha = 0.05, meanNoWS = 60,
               writeExcel = TRUE, file = "Baseline"){
  require(pwr)
  require(pwr2ppl)
  require(Superpower)
  require(pwrss)
  ans <- matrix(NA, nrow = length(prWS)*length(Cohen.d),
                ncol = 8)
  colnames(ans) <- c("prWS", "meanNoWS", "meanWS", "sd", "d", "nWS", "nNoWS", "N")
  r <- 1
  for (prws in prWS){
    for (d in Cohen.d){
      meanWS <- meanNoWS - sd*d
      N <- pwrss.t.2means(mu1 = meanNoWS, mu2 = meanWS,
                          sd1 = sd, sd2 = sd,
                          paired = FALSE,
                          kappa = prws/(1-prws),
                          margin = 0, power = power, alpha = alpha,
                          alt = "not equal")$n
      nWS <- N[1]
      nNoWS <- N[2]

      ans[r,] <- c(prws, meanNoWS, meanWS, sd, d, nWS, nNoWS, nWS + nNoWS)
      r <- r + 1
    }
  }
  if (writeExcel) write.csv(ans, file = paste0(file, "pow", power, ".csv"))
  ans
}

#' The function p2 computes the proportion in the second group given the
#' proportion in the first group and the Cohen h effect size.
#' h = 0.2 is "small", h = 0.5 is "medium", h = 0.80 is "large"
#' @param p1 proportion in the first group
#' @param h Cohen h effect size
#' @return proportion in the second group
p2 <- function(p1, h = 0.3){
  half.h <- h/2
  asinsqrtp1 <- asin(sqrt(p1))
  ans <- (sin(asinsqrtp1 - half.h))^2
  ans
}

#' @title Computing Sample Size for a given power and effect size in Weight Stigma Study using Proportion with WS
#' @description Main effect of WS on LE8 at baseline.
#'assuming LE8 is normally distributed with sd = 15 (from Lloyd-Jones et al.)
#' prWS = Pr(WS = 1), i.e. prob. a subject reports experiencing WS
#' @param sd standard deviation. Assumed to be equal among groups
#' @param Cohen.h Cohen's h. Can be a vector
#' @param prWS Proportion with Weight Stigma Variable
#' @param power Power. Default = .8
#' @param alpha Significance Level. Default = .05
#' @param propNoWS proportion without Weight Stigma Variable. Default = .70
#' @param writeExcel Default=True. To create an Excel file
#' @param file Default="Baseline" name of output file
#' @return a dataframe that includes the sample size for the specified power. Also excel file if specified.

ss.prop <- function(Cohen.h = c(.2, .3, .4, .5),
                    prWS = c(.5, .6, .7, .8), power = .80,
                    alpha = 0.05, propNoWS = .70,
                    writeExcel = TRUE, file = "Baseline"){
  ans <- matrix(NA, nrow = length(prWS)*length(Cohen.h),
                ncol = 7)
  colnames(ans) <- c("prWS", "propNoWS", "propWS", "h", "nWS", "nNoWS", "N")
  r <- 1
  for (prws in prWS){
    for (h in Cohen.h){
      propWS <- p2(propNoWS, h)
      N <- pwrss.z.2props(p1 = propNoWS, p2 = propWS,
                          arcsin.trans = FALSE,
                          kappa = prws/(1-prws),
                          margin = 0, power = power, alpha = alpha,
                          alt = "not equal")$n
      nWS <- N[1]
      nNoWS <- N[2]

      ans[r,] <- c(prws, propNoWS, propWS, h, nWS, nNoWS, nWS + nNoWS)
      r <- r + 1
    }
  }
  if (writeExcel) write.csv(ans, file = paste0(file, "-Proppow", power, ".csv"))
  ans
}

#' @title Computing Sample Size for detecting Moderation Effect
#' @description Main effect of WS on LE8 at baseline.
#'assuming LE8 is normally distributed with sd = 15 (from Lloyd-Jones et al.)
#' prWS = Pr(WS = 1), i.e. prob. a subject reports experiencing WS
#' meanNoWS = mean of LE8 when WS = 0
#' Cohen.d = [(mean LE8 with WS) - (mean LE8 without WS)]/sd
#' @param sd standard deviation. Assumed to be equal among groups
#' @param Cohen.m Cohen's m standaridized moderation effect. Can be a vector
#' @param prWS Proportion with Weight Stigma Variable
#' @param prM0 Probability that the moderator is 0. Can be a vector
#' @param alpha Significance Level. Default = .05
#' @param N Total Sample Size
#' @param writeExcel Default=True. To create an Excel file
#' @param file Default="Baseline" name of output file
#' @return a dataframe that includes the sample size for the specified power. Also excel file if specified.

ss.mod <- function(Cohen.m = c(.3, .4, .5, .6, .7, .8),
                   prWS = c(.5, .6, .7, .8),
                   prM0 = c(.5, .75),
                   N = seq(400, 800, by = 100),
                   alpha = 0.05,
                   writeExcel = TRUE, file = "Baseline"){
  ans <- matrix(NA, nrow = length(prWS)*length(prM0)*
                  length(Cohen.m)*length(N),
                ncol = 7)
  colnames(ans) <- c("prWS", "prM0", "m", "nWS", "nNoWS", "N", "power")
  r <- 1
  for (prws in prWS){
    for (prm0 in prM0){
      for (m in Cohen.m){
        for (ntot in N){

          pow <- power.moderation(N = ntot,
                                  pt0 = prws, # prWS
                                  pm0 = prm0, # pr(M = 0)
                                  m = m,
                                  sig.level = alpha)

          nWS <- pow$n0
          nNoWS <- pow$n1
          power <- pow$power

          ans[r,] <- c(prws, prm0, m, nWS, nNoWS, ntot, power)
          r <- r + 1
        }
      }
      if (writeExcel) write.csv(ans, file = paste0(file, "Modpow", ".csv"))
    }
  }
  ans
}
