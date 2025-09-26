# Check that rows in the survival
# model sum to one

check_row_sums <- function(){
  psiFYB <- runif(1)
  psiAB <- runif(1)
  psiBA <- runif(1)
  
  phiFY <- runif(1)
  phiA <- runif(1)
  phiB <- runif(1)
  
  pA <- runif(1)
  pB <- runif(1)
  
  po <- ps <- array(NA, dim=c(4,4))
  ps[ 1,1 ] <- 0
  ps[ 1,2 ] <- phiFY  * (1-psiFYB ) 
  ps[ 1,3 ] <- phiFY  * psiFYB 
  ps[ 1,4 ] <- (1-phiFY )
  
  ps[ 2,1 ] <- 0
  ps[ 2,2 ] <- phiA  * (1-psiAB )
  ps[ 2,3 ] <- phiA  * psiAB 
  ps[ 2,4 ] <- (1-phiA )
  
  ps[ 3,1 ] <- 0
  ps[ 3,2 ] <- phiB  * psiBA 
  ps[ 3,3 ] <- phiB  * (1-psiBA )
  ps[ 3,4 ] <- (1-phiB )
  
  ps[ 4,1 ] <- 0
  ps[ 4,2 ] <- 0
  ps[ 4,3 ] <- 0
  ps[ 4,4 ] <- 1
  
  # Define probabilities of O(t) given S(t)
  po[ 1,1 ] <- 1 
  po[ 1,2 ] <- 0
  po[ 1,3 ] <- 0
  po[ 1,4 ] <- 0
  
  po[ 2,1 ] <- 0
  po[ 2,2 ] <- pA 
  po[ 2,3 ] <- 0
  po[ 2,4 ] <- (1-pA )
  
  po[ 3,1 ] <- 0
  po[ 3,2 ] <- 0
  po[ 3,3 ] <- pB 
  po[ 3,4 ] <- (1-pB )
  
  po[ 4,1 ] <- 0
  po[ 4,2 ] <- 0
  po[ 4,3 ] <- 0
  po[ 4,4 ] <- 1
  
  return(
  list(ps_sums=rowSums(ps),
  po_sums=rowSums(po) )
  )
}

# Check row sums simulating different probs
# 10 times
for (i in 1:10){ print(paste0("iteration=",i)); check_row_sums() |> print() }

