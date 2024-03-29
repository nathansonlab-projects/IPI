
# ---------------------------------------------------------------------- #
jointMeta <- function(model.list, geno.varname, env.varname, snpname, alpha = 0.05, W = NULL)
# perform meta analysis for data with snp * environment interaction term from n datasets
# implementation of methods described in Becker & Wu, 2007, Manning et al., 2011.
#
# input:
#   model.list (list): entries are the glm objects from fitting each data set
#     eg the output of glm
#   geno.varname (string): name of the variable containing genotype data
#   env.varname (string): name of variable containing the environmental effect that
#     interacts with genotype
#   snpname (string): snp of interest  (for labeling output)
#   alpha (numeric): alpha-level for confidence intervals. default to 95% 
#   W (matrix): optional; weight-matrix if you using anything but the standard identity matrix
#   e.g., if including terms other than interaction and genotype main effect
#
# output:
#   out (data.frame): df containing estimates of snp and snp * env effects with confidence 
#   intervals, test statistics and p-values for joint coefficient test and test of heterogeneity, 
#   and the number of studies included for each snp
{
  library(matrixcalc)
  library(Matrix)
  
  int.varname <- paste(geno.varname, env.varname, sep = ":")
  nv <- c(geno.varname, int.varname)
  b <- c()
 
  k <- length(model.list)
  
  # reduce model.list to only valid models
  # setup b, vector of snp and snp * env coefficient estimates from the individual models
  
  i = 1    # the current index
  ct = 0  # count number of iterations; stop when reach k
  
  repeat
  {
    coef <- summary(model.list[[i]])$coefficients
    coef <- coef[!is.na(coef$Estimate),]
    # if the interaction term is NA, remove this model from the list
    # occurs when there is only 1 unique combination of genotype * env
    if( !(int.varname %in% rownames(coef)) )
    {
      model.list[[i]] <- NULL
    } else
    {
      b.tmp <- coef[nv,1]
      b <- c(b, b.tmp)
      i <- i + 1
    }
    
    ct <- ct + 1
    if( ct  == k )
    {
      break
    }
  } 
  
  # recalculate k now that invalid models are removed
  k <- length(model.list)
  
  # I needs to be 2x2 at minimum, even with only one study included
  # if( k == 1)
  # {
  #   I <- diag(1, nrow=2, ncol=2)
  # } else 
  # {
  #   I <- diag(1, nrow=k, ncol=k)
  # }
  
  I <- diag(1, nrow=2, ncol=2)
  
  sigma.list <- list()
  
  if( is.null(W) )
  {
    W <- c()
    for( i in 1:k)
    {
      W <- rbind(W, I)
    }
    
  }
  
  
  for( i in 1:k)
  {
    
    # get covariance matrices
    # vcov is identical to:
    #   X <- cbind(1, genotype, prior, genotype * prior)
    #   V <- diag(fit$fit * (1 - fit$fit))
    #   sigma.hat <- solve( t(X) %*% V %*% X )
    sigma.list[[i]] <- vcov(model.list[[i]])[nv,nv]
  }
  
  sigma <- bdiag(sigma.list)
  sigma.inv <- solve(sigma)
  
  # if k = 1, this is identical to the lm/glm solution
  covB <- solve( t(W) %*% sigma.inv %*% W)
  B <- covB %*% t(W) %*% sigma.inv %*% b
  
  # T statistics for joint hypothesis test
  T <- as.numeric(t(B) %*% solve(covB) %*% B)
  p <- pchisq(T, df = 2, lower.tail = F)
  
  # test of heterogeneity
  Q = as.numeric(t(b - W %*% B) %*% sigma.inv %*% (b - W %*% B))
  q.p <- pchisq(Q, df= 1, lower.tail = F)
  
  z <- qnorm(1 - (alpha/2))
  out <- data.frame(snp = snpname,
                    bp = strsplit(snpname, ":")[[1]][2],
                    geno.b = B[1,1], 
                    geno.b.se = covB[1,1], 
                    geno.b.ci.lo = B[1] - z * sqrt(covB[1,1]),
                    geno.b.ci.hi = B[1] + z * sqrt(covB[1,1]),
                    int.b = B[2,1],
                    int.b.se = covB[2,2],
                    int.b.ci.lo = B[2] - z,
                    int.b.ci.hi = B[2] + z,
                    joint.p = p,
                    het.p = q.p,
                    n.studies = k
  )
  
  
  return(out)
}
# ---------------------------------------------------------------------- #
