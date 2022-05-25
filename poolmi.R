MultipleImputationCoefficientTable <- function(models, large.sample.df = FALSE)
{
  m <- length(models)
  .coef <- function(object)
  {
    object$coef
  }
  coefs <- sapply(models, .coef)
  vars <- sapply(models, FUN = function(x) diag(vcov(x)))
  coef.mean <- apply(coefs, 1, mean, na.rm = FALSE)
  ses <- multipleImputationStandardErrors(coefs, vars)
  df.c <- df.residual(models[[1]])
  dfs <- multipleImputationDegreesOfFreedom(coefs, vars, df.c, large.sample.df)
  tvals <- coef.mean / ses
  correct <- models[[1]]$correction
  pvals <-  pvalAdjust(2 * pt(abs(tvals), dfs, lower.tail = FALSE), correct)
  results <- cbind(coef.mean, ses, tvals, dfs, pvals)
  row.rm <- which(apply(results, 1, function(x){all(is.na(x))}))
  if (length(row.rm) > 0)
    results <- results[-row.rm,]
  coef.names <- rownames(models[[1]]$summary$coef)
  if(models[[1]]$type == "Multinomial Logit")
  {
    coef.names <- models[[1]]$original$coefnames
    levs <- models[[1]]$original$lev[-1]
    n.levs <- length(levs)
    k <- length(pvals)
    n.coef.names <- length(coef.names)
    coef.names <- rep(coef.names, rep(n.levs, n.coef.names))
    #levs <- levs[rep(1:n.levs, rep(k / n.levs, n.levs))]
    coef.names <- paste(coef.names, levs)
  }
  dimnames(results) <- list(coef.names,
                            c("Estimate", "Std. Error", "t value", "Degrees of Freedom", "Pr(>|t|)"))
  results
}


multipleImputationDegreesOfFreedom <- function(coefs, vars, df.complete, large.sample.df = FALSE)
{
  #  https://www.stata.com/manuals13/mi.pdf    Page 71 of the pDF  Rubin 1987
  m <- ncol(coefs)
  B <- apply(coefs, 1, var)
  W <- apply(vars, 1, FUN = mean, na.rm = FALSE)
  r <- (1 + 1 / m) * B / W
  df.large <- (m - 1) * (1 + 1 / r) ^ 2
  if (large.sample.df)
    return(df.large)
  if (missing(df.complete))
    stop("'df.complete': Degrees of freedom of the complete model.")
  T.denom <- W + (1 + 1 / m) * B
  lambda <- (1 + 1 / m) * B / T.denom
  df.obs <- df.complete * (df.complete + 1) * ( 1 - lambda) / (df.complete + 3)
  df.small <- (1 / df.large + 1 / df.obs) ^ -1
  df.small
}

multipleImputationStandardErrors <- function(coefs, vars)
{
  B <- apply(coefs, 1, var)
  W <- apply(vars, 1, FUN = mean, na.rm = FALSE)
  T <- W + (1 + 1 / ncol(coefs)) * B
  sqrt(T)
}

multipleImputationRelativeImportance <- function(models)
{
  coefs <- sapply(models, function(object) object$coef)
  tmp.coefs <- unname(apply(coefs, 1, mean, na.rm = FALSE))
  signs <- if (models[[1]]$type == "Ordered Logit") sign(tmp.coefs[1:(models[[1]]$n.predictors)])
  else sign(tmp.coefs[-1])
  
  result <- list()
  correct <- models[[1]]$correction
  models.raw.importance <- sapply(models, function(m) m$relative.importance$raw.importance)
  result$raw.importance <- apply(models.raw.importance, 1, mean, na.rm = FALSE)
  result$importance <- signs * 100 * prop.table(result$raw.importance)
  vars <- sapply(models, function(m) m$relative.importance$standard.errors ^ 2)
  result$standard.errors <- multipleImputationStandardErrors(models.raw.importance, vars)
  df.c <- df.residual(models[[1]])
  result$df <- multipleImputationDegreesOfFreedom(models.raw.importance, vars, df.c, FALSE)
  result$statistics <- signs * result$raw.importance / result$standard.errors
  result$p.values.raw <-  2 * pt(abs(result$statistics), result$df, lower.tail = FALSE)
  result$p.values <- pvalAdjust(result$p.values.raw, correct)
  result$statistic.name <- "t"
  result
}

#' @importFrom stats pchisq
multipleImputationCrosstabInteraction <- function(models, relative.importance)
{
  n <- nrow(models[[1]]$interaction$bb)
  m <- ncol(models[[1]]$interaction$bb)
  split.size <- models[[1]]$interaction$split.size
  correction <- models[[1]]$correction
  res <- list(label = models[[1]]$interaction$label,
              split.size = split.size, pvalue = NA,
              relative.importance = relative.importance,
              original.r2 = mean(sapply(models, function(m){m$interaction$original.r2})),
              full.r2 = mean(sapply(models, function(m){m$interaction$full.r2})))
  
  bb.all <- sapply(models, function(m){m$interaction$bb})
  bc.all <- sapply(models, function(m){m$interaction$bc})
  ss.all <- sapply(models, function(m){m$interaction$ss^2})
  sc.all <- sapply(models, function(m){m$interaction$sc^2})
  bb <- apply(bb.all, 1, mean, na.rm=T)
  bc <- apply(bc.all, 1, mean, na.rm=T)
  ss <- multipleImputationStandardErrors(bb.all, ss.all)
  sc <- multipleImputationStandardErrors(bc.all, sc.all)
  coef.sign <- compareCoef(matrix(bb, nrow=n), matrix(bc, nrow=n),
                           matrix(ss, nrow=n), matrix(sc, nrow=n),
                           split.size[1:m], correction, relative.importance)
  res$coef.pvalues <- coef.sign$pvalues
  res$coef.tstat <- coef.sign$tstat
  if (relative.importance)
    res$coef.pFDR <- coef.sign$pFDR
  
  net.coef.all <- sapply(models, function(m){m$interaction$net.coef})
  net.coef <- apply(net.coef.all, 1, mean)
  bb <- matrix(bb, nrow=n)
  
  # Report normalised relative importance scores but use raw scores for p-values
  if (relative.importance)
    bb <- apply(bb, 2, function(x){x/sum(abs(x))*100})
  
  combined.coefs <- cbind(bb, net.coef)
  colnames(combined.coefs) <- names(split.size)
  res$coefficients <- combined.coefs
  
  
  if (relative.importance)
    return(res)
  
  if (models[[1]]$type %in% c("Linear", "Quasi-Poisson"))
  {
    fstat <- mean(sapply(models, function(m){m$interaction$anova.fstat}))
    df1 <- models[[1]]$interaction$anova.df1
    df2 <- models[[1]]$interaction$anova.df2
    res$anova.test <- "F-test"
    res$pvalue <-  pf(fstat, df1, df2, lower.tail=F)
  } else
  {
    cstat <- mean(sapply(models, function(m){m$interaction$anova.dev}))
    df1 <- models[[1]]$interaction$anova.df1
    res$anova.test <- "Chi-square test"
    res$pvalue <- pchisq(cstat, df1, lower.tail=F)
  }
  return(res)
}

