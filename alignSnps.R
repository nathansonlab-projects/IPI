
# pluta 5/7/21

# set of functions for aligning snps between two datasets

# ---------------------------------------------------------------------- #
checkGenoFlip <- function(bim1, bim2)
# compare genotypes from two datasets to see if they match
# input: bim1, bim2 (data.frames)- bim files from plink
# output: boolean- TRUE if both alleles match
{
  # if the same alleles are in both studies but in reverse order, need to flip
  if( sum( c(bim1$V5, bim1$V6) %in% c(bim2$V5, bim2$V6)) == 2)
  {
    if( bim1$V5 == bim2$V6 )
    {
      return(TRUE)
    }
  } 
  return(FALSE)
}
# ---------------------------------------------------------------------- #


# ---------------------------------------------------------------------- #
flipGeno <- function( geno.vec )
# if change minor allele to major allele, also flip geno types
# 0 becomes 2, 2 becomes 0
# input: geno.vec (integer), vector of genotype values (0,1,2)
# output: geno.vec (integer), vector of genotypes with 0 and 2 reversed
{
  geno.vec[ geno.vec == 0] <- -9
  geno.vec[ geno.vec == 2] <- 0
  geno.vec[ geno.vec == -9] <- 2
  return(geno.vec)
}
# ---------------------------------------------------------------------- #

# --------------------------- flipAllele ------------------------------------------- #
flipAllele <- function( allele )
# input: allele (character), one of AGCT
# output: allele (character), one of AGCT, flipped from input
{
  print(paste0("allele = ", allele))
  allele <- as.character(allele)
  
  if( allele == "A")
  {
    return("T")
  } else
    if( allele == "G")
    {
      return("C")
    } else
      if( allele == "C" )
      {
        return("G")
      } else
        if( allele == "T")
        {
          return("A")
        } else
          return(allele)
}
# ---------------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------------- #
alignSnps <- function( snp, BIM1, BIM2, geno.dat )
# align snps in BIM1 to BIM2, flip corresponding genotypes in geno.dat
# BIM1, BIM2 (data.frame), bim files from plink
# geno.dat (data.frame), genotype data imported by BEDMatrix
{
  if( length(which(BIM1$V2 == snp)) == 0)
  {
    print(paste0("dim(BIM1) = ", dim(BIM1)))
    print(paste0("dim(BIM2) = ", dim(BIM1)))
    print(paste0("snp = ", snp))
  }
  
  BIM1$V2 <- as.character(BIM1$V2)
  BIM2$V2 <- as.character(BIM2$V2)
  
  bim1 <- BIM1[ which(BIM1$V2 == snp),]
  bim2 <- BIM2[ which(BIM2$V2 == snp),]
  geno <- geno.dat[,snp]
  
  # align alleles
  # case one, the ref and alt allele are flipped in direction
  # eg C T  and T C
  if( checkGenoFlip(bim1, bim2))
  {
    geno <- flipGeno(geno)
  }
  
  # case two, the ref and alt allele are flipped in both name and/or direction
  # eg C T  and G A
  
  print(paste0("bim2$V5 == ", bim2$V5))
  print(paste0("bim2$V6 == ", bim2$V6))
  if( sum( c(bim1$V5, bim1$V6) %in% c(bim2$V5, bim2$V6)) == 0)
  {
    
    bim2$V5 <- flipAllele(bim2$V5)
    bim2$V6 <- flipAllele(bim2$V6)
    
    if( checkGenoFlip(bim1, bim2) )
    {
      geno <- flipGeno(geno)
    }
  }
  
  return(geno)
}
# ---------------------------------------------------------------------------------- #

# sample call:
#test = lapply(colnames(dat1[[1]][,1:10]), alignSnps, BIM1, BIM2, dat1[[1]])
#out <- matrix(unlist(test), nrow=226,ncol=10)