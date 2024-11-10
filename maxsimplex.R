#' Title simplex for maximisation
#'
#' @param A A matrix containing the coefficients of the constraints.
#' @param b A vector containing the values of the constraints.
#' @param z A vector containing the coefficients of the objective function
#'
#' @return A matrix representing the linear programming problem
#' @export
#'
#' @examples
#'A <- matrix(c(2,1,1,3),nrow=2,byrow=TRUE)
#' b <- c(10,15)
#' z <- c(3,5)
#'
maxsimplex <- function(A,b,z){

  # dimensions of problem
  n <- ncol(A) #number of variables
  m <- nrow(A) #number of constraints

  #initialization
  matrice_simplexe <- cbind(A,diag(m),b)
  matrice_simplexe <- rbind(matrice_simplexe,c(-z,rep(0,m),0))
  colnames(matrice_simplexe) <- c(paste0("x", 1:n),paste0("s",1:m),"RHS")
  rownames(matrice_simplexe) <- c(paste0("contraintes",1:m),"z")

  # simplex algorithm
  iteration <- 0
  max_iteration <- 1000 #limiter le nombre d'iteration pour les boucles infinies
  repeat{
    iteration <- iteration + 1

    # pivot  column
    colonne_pivot <- which.min(matrice_simplexe[m+1, 1:(n+m)])
    if(matrice_simplexe[m+1,colonne_pivot]>=0){
      break # Solution optimale
    }

    # pivot row
    ratios <- matrice_simplexe[1:m,n+m+1]/matrice_simplexe[1:m,colonne_pivot]
    ratios[matrice_simplexe[1:m,colonne_pivot] <=0] <- Inf
    if(all(is.infinite(ratios))){
      stop(" solution non bornee")
    }
    ligne_pivot <- which.min(ratios)

    # pivot  Operation
    pivot <- matrice_simplexe[ligne_pivot,colonne_pivot]
    matrice_simplexe[ligne_pivot,] <- matrice_simplexe[ligne_pivot,]/pivot
    for (i in 1:(m+1)) {
      if(i!=ligne_pivot){
        matrice_simplexe[i,] <- matrice_simplexe[i,]-matrice_simplexe[i,colonne_pivot]*matrice_simplexe[ligne_pivot,]
      }
    }
    if(iteration>max_iteration){
      stop("Number of iteration reached")
    }
  }

  # extraction of  solution
  solution <- numeric(n)
  for(i in 1:n){
    colonne <- which(matrice_simplexe[, 1:(n+m)]==1, arr.ind = TRUE)
    ligne <- which(colonne[,2]==i)
    if(length(ligne)>0){
      solution[i] <- matrice_simplexe[colonne[ligne[1], 1],n+m+1]
    }
  }
  z_optimal <- matrice_simplexe[m+1,n+m+1]
  return(list(solution=solution,z_optimal=z_optimal,matrice_finale=matrice_simplexe))
}

