rm(list = ls())

# pluta 9/14/20
# general script  to run all 4 kinds of ipi models (noprior - irae, priorint - irae, 
# noprior - coxph, priorint - coxph)
# writes output  with everything required for meta/pathway analysis

source("~/IPI/iraeAssoc.R")

# ------------------------------------ user input  ------------------------------ #
args = commandArgs(trailingOnly = TRUE)

if( length(args) < 7)
{
  print("need to provide 7 arguments: BEDFILE COVARS COVFILE CHR MODEL MAF OUTNAME")
  print("BEDFILE: plink .bed file of genotype data, must have corresponding .bim and .fam files")
  print("COVARS: list of covariates to use in the model, in a comma seperated list")
  print("COVFILE: file of covariate data (ipi.nivo.pheno.txt)")
  print("CHR: numeric chromosome of interest")
  print("Model: choice of 1-5:")
  print("1: irae, no prior treatment")
  print("2: irae, prior interaction")
  print("3: coxph, no prior treatemt")
  print("4: coxph, prior interaction")
  print("5: irae, all subjects")
  print("OUTPREFIX: prefix attached to all output files")
  print("")
  print("If no covariates are used, set COVARS to NA")
  print("")
  print("")
  print(paste0("you provided ", length(args), " arguments:"))
  print(paste0("BEDFILE = ", args[1]))
  print(paste0("COVARS = ", args[2]))
  print(paste0("COVFILE = ", args[3]))
  print(paste0("CHR = ", args[4]))
  print(paste0("MDL = ", args[5]))
  print(paste0("MAF = ", args[6]))
  print(paste0("OUTPREFIX = ", args[7]))
  stop("missing necessary arguments")
}

BEDFILE = args[1]
COVARS = args[2]
COVFILE = args[3]
CHR = args[4]
MDL=args[5]
MAF=args[6]
OUTPREFIX=args[7]
# ------------------------------------------------------------------------------- #


# =============================================================================== #
# ================================== MAIN ======================================= #
# =============================================================================== #
print(paste0("BEGIN: ", date()))

print("running with the following parameters:")
print(paste0("BEDFILE = ", BEDFILE))
print(paste0("COVARS = ", COVARS))
print(paste0("COVFILE = ", COVFILE))
print(paste0("CHR = ", CHR))
print(paste0("MAF = ", MAF))
print(paste0("OUTPREFIX = ", OUTPREFIX))
print("")

# the list of variables we will always consider, plus the additional specified covariates
if( COVARS == "NA")
{
  covar.names <- c("Prior", "Surv_Months", "Vital_Status_2yrs")
} else
{
  covar.names <- c(strsplit(COVARS, ",")[[1]], "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
}

if( !file.exists(BEDFILE))
{
  stop(paste0(BEDFILE, " not found. Exiting"))
}

if( !file.exists(COVFILE))
{
  stop(paste0(COVFILE, " not found. Exiting"))
}


print("reading input...")
geno.dat <- BEDMatrix(BEDFILE, simple_names = TRUE)
cov.dat <- read.table(COVFILE, header = T, sep = ",")

if(all(is.na(match(rownames(geno.dat), cov.dat$GWASID))))
{
  stop("could not match geno.dat subject ids to cov.dat subjects ids")
}

cov.dat <- cov.dat[match(rownames(geno.dat), cov.dat$GWASID),]
cov.dat <- cov.dat[ ,colnames(cov.dat) %in% covar.names,]
print("done")

print("setting up cores for parallel processing....")

# read BIM file for A1 and A2
BIM <- read.table( paste0( substr(BEDFILE, 1, nchar(BEDFILE) - 4), ".bim"), header = FALSE)
colnames(BIM) <- c("chr", "MarkerName", "GD", "bp", "A1", "A2")

# macOS workaround
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")
print("done")



# ------------- irae, no prior treatment ------------- #
if( MDL == 1)
{
  # make sure you specified the right covar names
  if(any(!(covar.names %in% colnames(cov.dat))))
  {
    ind <- which( !(covar.names %in% colnames(cov.dat)))
    print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
    stop()
  }
  
  # if irae, drop survival data
  cov.dat <- cov.dat[,!(colnames(cov.dat) %in% c("Surv_Months", "Vital_Status_2yrs"))]
  
  print("run glms for irae noprior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-noprior.irae.txt")
  runSNP.glm.parallel( geno.dat[ cov.dat$Prior == 0,], 
                       cov.dat[ cov.dat$Prior == 0,], 
                       BIM,
                       "iraeAssoc",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ---------------- irae, prior treatment -------------- #
if( MDL == 2)
{
  print("run glms for irae prior-interaction model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-priorint.irae.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "iraeAssoc.priorint",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ------------------ coxph no prior ------------------- #
if( MDL == 3)
{
  print("run glms for coxph no prior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-noprior.coxph.txt")
  runSNP.glm.parallel( geno.dat[ cov.dat$Prior == 0,], 
                       cov.dat[ cov.dat$Prior == 0,], 
                       BIM,
                       "coxPH",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ----------------- coxph prior int ------------------- #
if( MDL == 4)
{
  print("run glms for coxph prior-interaction model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-priorint.coxph.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "coxPH.priorint",
                       MAF,
                       outname)
}
# ----------------------------------------------------- #


# ----------------- irae ------------------- #
if( MDL == 5)
{
  # make sure you specified the right covar names
  if(any(!(covar.names %in% colnames(cov.dat))))
  {
    ind <- which( !(covar.names %in% colnames(cov.dat)))
    print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
    stop()
  }
  
  # if irae, drop survival data
  cov.dat <- cov.dat[,!(colnames(cov.dat) %in% c("Surv_Months", "Vital_Status_2yrs"))]
  
  print("run glms for irae noprior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-all.irae.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "iraeAssoc",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


print("stopping clusters...")
stopCluster(cl)
print("done")

print(paste0("analysis complete: ", date()))
# =============================================================================== #
# =============================================================================== #
# =============================================================================== #
