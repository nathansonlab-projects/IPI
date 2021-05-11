
# pluta 5/11/21
# v0.1

# setwd("~/Documents/nathansonlab/IPI/Meta/ipinivo_priorint/")
library(BEDMatrix) 
library(parallel)
library(pbapply)
library(matrixcalc)
library(Matrix)

source("~/IPI/jointMeta.R")
source("~/IPI/iraeAssoc.R")
source("~/gwas_tools/alignSnps.R")


# ================================================================================= #
# ======================================= functions =============================== #
# ================================================================================= #

# --------------------------------------------------------------------------------- #
prepData <- function(BEDFILE, COVFILE, COVARS)
# function to set up genotypes with phenotype and covariate data
# input:
#   BEDFILE (BEDMatrix object), genotype data
#   COVFILE (string), name of the file containing covariate data
#   COVARS (string), comma separated string of covariates, corresponding to columns in
#       COVFILE
#
# output:
#   list of data.frames, containing genotype and covariate data, truncated to the covariates
#   of interest
{
  # the list of variables we will always consider, plus the additional specified covariates
  if( COVARS == "NA")
  {
    covar.names <- c("PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
  } else
  {
    covar.names <- c(strsplit(COVARS, ",")[[1]], "PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
  }
  
  if( !file.exists(BEDFILE))
  {
    stop(paste0(BEDFILE, " not found. Exiting"))
  }
  
  if( !file.exists(COVFILE))
  {
    stop(paste0(COVFILE, " not found. Exiting"))
  }
  
  # much more efficient  than trying to use the .raw file
  print("reading input...")
  geno.dat <- BEDMatrix(BEDFILE, simple_names = TRUE)
  cov.dat <- read.table(COVFILE, header = T, sep = ",")
  
  if(all(is.na(match(rownames(geno.dat), cov.dat$GWASID))))
  {
    stop("could not match geno.dat subject ids to cov.dat subjects ids")
  }
  
  if(!all(covar.names %in% colnames(cov.dat)))
  {
    print("the following covars were not found in the covariate data:")
    print(covars[covars %in% colnames(cov.dat)])
    stop()
  }
  
  cov.dat <- cov.dat[match(rownames(geno.dat), cov.dat$GWASID),]
  cov.dat <- cov.dat[ ,colnames(cov.dat) %in% covar.names,]
  print("done")
  
  
  return(list(geno.dat,cov.dat))
  
}
# ---------------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------------- #
runJointMetaSNP <- function(snp, dat1, dat2)
# run joint meta analysis for a given snp
# joint meta requires the full variance-covariance matrix between the two data sets;
# full model needs to be run
# 
# input: snp (string), snp of interest, should be a column of dat1 and dat2
# output: df (data.frame), summary statistics for the joint meta analysis of snp
{
  ind <- colnames(dat1[[1]]) == snp
  
  # this needs to be a model.obj, oops
  fit1 <- iraeAssoc.priorint( dat1[[1]][,ind], dat1[[2]], fit.only = TRUE)
  fit2 <- iraeAssoc.priorint( dat2[[1]][,ind], dat2[[2]], fit.only = TRUE)
  df <- jointMeta(list(fit1, fit2), "geno.dat", "Prior", colnames(dat1[[1]])[ind]) 
  return(df)
}
# ---------------------------------------------------------------------------------- #
# ================================================================================== #
# ================================================================================== #
# ================================================================================== #



# ================================================================================== #
# ================================= MAIN  ========================================== #
# ================================================================================== #

args = commandArgs(trailingOnly=TRUE)

#  perform joint meta analysis of 2 datasets
if( length(args) < 7 )
{
  print("need to provide 6 arguments: BEDFILE1 COVFILE1 COVARS1 BEDFILE2 COVFILE2 COVARS2 OUTNAME")
  print(" ")
  stop()
}

#  how can i generalize this to n studies?
BEDFILE1 = args[1]
COVFILE1 = args[2]
COVARS1 = args[3]

BEDFILE2 = args[4]
COVFILE2 = args[5]
COVARS2  = args[6]

OUTNAME = args[7]
#BEDFILE1 = "chr22-ipinivo-QC2.bed"
#BIMFILE1 = "chr22-ipinivo-QC.bim"
#COVFILE1 = "../../ipi.nivo.pheno.txt"
#COVARS1 = "studyarm,NDoseIpi_2L"
#BEDFILE2 = "nivo-chr22.qc.bed"
#BIMFILE2 = "nivo-chr22.qc.bim"
#COVFILE2 =  "../../nivo.pheno.txt"
#COVARS2 = "NDose.Nivo,Stage"

print("reading first set of files...")
print(paste0("BEDFILE1 = ", BEDFILE1))

BIMFILE1 <- paste0(substr(BEDFILE1, 1, nchar(BEDFILE1) - 3), "bim")
print(paste0("BIMFILE1 = ", BIMFILE1))
print(paste0("COVFILE1 = ", COVFILE1))
print(paste0("COVARS1 = ", COVARS1))

BIM1 <- read.table(BIMFILE1, header = F)
#COVFILE1 <- read.table(COVFILE1, header = T, sep = ",")
dat1 <- prepData(BEDFILE1, COVFILE1, COVARS1)
print('done')

print('reading second set of files...')
print(paste0("BEDFILE2 = ", BEDFILE2))
#COVFILE2 <- read.table(COVFILE2, header = T, sep = ",")
BIMFILE2 <- paste0(substr(BEDFILE2, 1, nchar(BEDFILE2) - 3), "bim")
print(paste0("BIMFILE2 = ", BIMFILE2))
print(paste0("COVFILE2 = ", COVFILE2))
print(paste0("COVARS2 = ", COVARS2))

BIM2 <- read.table(BIMFILE2, header = F)
dat2 <- prepData(BEDFILE2, COVFILE2, COVARS2)
print('done')

# need to remove duplicate snps; 
dat1[[1]] <- dat1[[1]][,!duplicated(colnames(dat1[[1]]))]
dat2[[1]] <- dat2[[1]][,!duplicated(colnames(dat2[[1]]))]

# for now, reduce to the common set of snps ****
snps <- intersect(colnames(dat1[[1]]), colnames(dat2[[1]]))
snps <- snps[1:10]
dat1[[1]] <- dat1[[1]][ ,colnames(dat1[[1]]) %in% snps]
dat2[[1]] <- dat2[[1]][ ,colnames(dat2[[1]]) %in% snps]

# do i need to align snps? wont hetereogenity test catch this?
print("aligning snps...")

# cant do this with apply because i need to specify both the column and the snp name
# maybe there is a better solution
#dat2[[1]] <- apply(dat2[[1]], 2, alignSnps, BIM1, BIM2, colnames(dat2[[1]])[j])

for( i in 1:length(snps))
{
  dat2[[1]][,i] <- alignSnps(dat2[[1]][,i], BIM1, BIM2, colnames(dat2[[1]])[i])
}
#dat2[[1]] <- unlist(dat2[[1]])
#dat2[[2]] <- unlist(dat2[[2]])
print("done")


print("setting up parallel processing..")
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")

# export necessary functions
clusterExport(cl, list("iraeAssoc.priorint", "jointMeta", "bdiag"))
print("done")

print("running meta-analysis...")
out <- pblapply(snps, runJointMetaSNP, dat1, dat2, cl = cl)

# convert list output to df
out <- do.call(rbind.data.frame, out)
print("done")

stopCluster(cl)

write.table(OUTNAME, out, col.names = TRUE, row.names = FALSE, append = FALSE)

# ================================================================================== #
