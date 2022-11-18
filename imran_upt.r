rb <- function(th,k=2) {
  k*(th[2]-(th[1])^2)^2 + (1-th[1])^2 
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}




newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  #Initializing the values before iterating
  f <- func(theta,...) #Defining the value of the function at the current theta.
  if (is.infinite(f)){stop("The objective value is infinite.")} #If the objective
  #is infinite, we stop and return a warning message.
  
  g <- grad(theta) #Defining the value of the gradient at the current theta.
  if (sum(is.infinite(g)) != 0) {stop("At least one derivative value is infinite.")}
  #If there is infinite element in gradient, we stop and return a warning message.
  iter <- 0 #Starting the number of iterations with zero.
  
  for (i in 1:maxit) {
    f <- func(theta,...) #Updating the objective in each iteration.
    g <- grad(theta) #Updating the gradient in each iteration.
    if (is.null(hess)){ #If hessian is not provided, then we will approximate it
      #by finite differencing as below. This method is taken from lecture notes.
      g0 <- grad(theta) ## gradient of the function at theta
      h <- matrix(0,length(theta),length(theta)) ## finite difference Hessian
      for (i in 1:length(theta)) { ## loop over parameters
        theta1 <- theta; theta1[i] <- theta1[i] + eps ## increase theta[i] by eps
        g1 <- grad(theta1) ## compute resulting gradient
        h[i,] <- (g1 - g0)/eps ## approximate second derivatives
      }
    }
    else {h <- hess(theta)} #If it is provided, define the current hessian as h.
    hf <- h #As we may need to perturb h to be positive definite, defining another
    #hessian matrix to check in the end whether it is positive definite or not.
    R <- try(chol(h), silent = TRUE) # try to decompose the matrix by Cholesky decomposition
    if(inherits(R, "try-error")) # if decomposition fails, then current hessian is
      #not positive definite and we have to perturb it to be so.
    {
      multiplier <- norm(h)*1e-6 #The first multiplier for perturbation
      while(inherits(R, "try-error")) # if it fails again, we keep on perturbing
        #until the hessian becomes positive definite.
      {
        h <- h + diag(length(theta))*multiplier #Perturbing by multiplying with multiplier.
        R <- try(chol(h)) #Checking positive definiteness again, which will be
        #the argument inside of the while loop.
        multiplier <- 10*multiplier #Updating multiplier if it's still not pos.def.
      }
    }
    
    delta <- backsolve(R, forwardsolve(t(R), -g)) #Finding the descent direction.
    if (sum(abs(g) < tol*(abs(f) + fscale)) == length(theta)) {break} #Judging
    #convergence. If it converges, then we break the loop.
    
    if (func(theta + delta,...) >= f) { #If our descent direction can't improve 
      #the objective, we have to halve the direction vector.
      for(j in 1:max.half) { #Halving max max.half times.
        delta <- delta/2 #Halving
        if(func(theta + delta,...) < f & is.finite(func(theta+delta)) == TRUE)
          {break} #If the direction starts improving, we break halving.
      }
    }
    if (func(theta + delta,...) >= f) {warning("max.half is reached without convergence.")
      break} #If the direction still can't improve the objective after max.half
    #times, we break the optimization iterations and show a warning.
    if(is.finite(func(theta+delta,...)) == FALSE){
      warning("max.half is reached with the trial objective function not finite", call. = FALSE)}
    #If the objective is not finite (-Inf, NaN, NA etc.), we show a warning.
    
    theta <- theta + delta #After finding a decreasing direction, updating theta.
    iter <- iter+1 #As we updated theta, we increase the total number of iterations.
  }
  
  if ((sum(abs(g) < tol*(abs(f) + fscale)) != length(theta)) & iter == maxit) {
    warning("maxit is reached without convergence.") #If all elements of the gradient
    #vector still do not converge to zero after maxit iteration, we show a warning.
  } 
  
  if (!inherits(try(chol(hf), silent = TRUE), "try-error")) { #We are checking
    #whether our last not-perturbed hessian matrix is positive definite by using
    #Cholesky decomposition.
    Rf <- try(chol(hf))
    Hi <- backsolve(Rf, forwardsolve(t(Rf), diag(length(theta)))) #If it is positive
    #definite, we are finding its inverse by Cholesky decomp. and return it
  }
  else{
    Hi <- NULL #If the last hessian is not pos. definite, we do not return it.
    warning("Hessian matrix is not positive definite", call. = FALSE)
  }
  
  return(list("f"=f,"theta"=theta, "g"=g, "iter"=iter, "Hi"=Hi)) #Returning what we have.
}

theta <- c(200,200)
n1 <- newt(theta,rb,gb,hess=hb,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
n2 <- newt(theta,rb,gb,hess=NULL,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
system.time(newt(theta,rb,gb,hess=hb,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6))




installed.packages(priority="base")

installed.packages(priority="recommended")
