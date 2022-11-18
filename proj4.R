# Group members:
# 1. Jie Xin Ng s1859154
# 2. Qiqi Xiao s2323497
# 3. Imran Toprak s2430699

# GitHub Repo: https://github.com/narmi1998/proj4

# Team member contributions: 
# the three of us did all the questions and helped each other out when stuck
# we compared and modified our codes into this final version
# each of us contributed equally towards this assignment, 33% each 

### Overview ####
# This project aims to create a function 'newt' 
# 'newt' is a function implementing Newton’s method for minimisation of functions
# Newton’s method minimises successive quadratic approximations (formula in 'newt')
# to the objective function, D 
# The parameters are initially guessed, which lead to finding the quadratic function 
# matching D, the gradient vector and The Hessian matrix at that guess. 
# The quadratic is minimised to find an improved guess and to update the hessian matrix
# until we find a guess at the parameters which leads to gradient = 0
# The condition that the gradient vector must be zero and the 
# Hessian matrix has to be positive definite must be satisfied to find the minimum.
# To avoid the possibility of divergence, D is checked on whether it is reduced
# at each step. 
# If D is not reduced, the algorithm backtracks to the previous value to find a 
# value which reduces D.
# Another important check is to check whether the updated Hessian matrix is 
# positive definite at each step.
# This is to ensure that the quadratic has a proper minimum.
# If the Hessian matrix is not positive definite, it is perturbed by adding a 
# multiple of the identity matrix until it is so that it could be used in the algorithm. 
# In 'newt', the convergence (minimum) is estimated by testing if all elements 
# of the gradient vector have absolute value less than a constant, 
# tol (convergence tolerance) times the absolute value of D
# plus fscale (estimate of magnitude of D near optimum)
# If the initial guess does not lead to non-finite D or gradient, 
# 'newt' will always return the estimates
# warnings are issued to inform user when the returned value is far 
# from supposed minimum, the warnings are disccussed further in 'newt' 
#################





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

f <- function(th, k = 2)
{
  (th[1])^2 + th[1]*k
}

g <- function(th, k = 2)
{
  2*th[1] + k
}

theta <- c(10)
diag(length(theta))

newt(c(10),f,g,hess=hb,k = 3,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)



## loop over parameters
for (i in 1:length(theta)) { 
  theta1 <- theta; theta1[i] <- theta1[i] + eps ## increase theta[i] by eps
  g1 <- grad(theta1,...) ## compute resulting gradient
  h[i,] <- (g1 - g0)/eps ## approximate second derivatives
}


newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  
  # newt is a function created to find the estimated minimisation of a function 
  # This is implemented using Newton's method of optimisation 
  # 5 types of warnings are issued when:
  # warning 1: func or grad are not finite at the initial theta (stops the algorithm)
  # warning 2: If step fails to reduce func after trying max.half step halvings
  # warning 3: If maxit is reached without convergence
  # warning 4: If hess is not positive definite at convergence
  # warning 5: If func is not finite at an iteration
  #
  # @param theta - vector of initial values for the optimization parameters
  # @param func - objective function to minimize with params theta and ...
  # @param grad - gradient vector of func w.r.t. parameter vector
  # @param hess - Hessian matrix function w.r.t. parameter vector, approximated if NULL
  # @param ... - arguments of func, grad and hess after parameter vector are passed
  # @param tol - convergence tolerance
  # @param fscale - rough estimate of func near the optimum (for convergence testing)
  # @param maxit - maximum number of iterations to try before giving up
  # @param max.half - maximum number of times a step is halved before concluding 
  # that the step has failed to improve the objective.
  # @param eps - the finite difference intervals to estimate hess if it is not provided.
  #
  # If warning 1 does not occur,
  # @return list containing:
  # f = estimated minimum of the objective function
  # theta = the value of the parameters at estimated minimum
  # iter = number of iterations to find estimated minimum
  # g = the gradient vector at estimated minimum 
  # Hi = inverse of Hessian matrix at estimated minimum 
  
  
  #Initializing the values before iterating
  f <- func(theta,...) #Defining the value of the function at the current theta.
  if (is.infinite(f)){stop("The objective value is infinite.")} #If the objective
  #is infinite, we stop and return a warning message.
  
  g <- grad(theta,...) #Defining the value of the gradient at the current theta.
  if (sum(is.infinite(g)) != 0) {stop("At least one derivative value is infinite.")}
  #If there is infinite element in gradient, we stop and return a warning message.
  iter <- 0 #Starting the number of iterations with zero.
  
  for (i in 1:maxit) {
    f <- func(theta,...) #Updating the objective in each iteration.
    g <- grad(theta,...) #Updating the gradient in each iteration.
    
    if (is.null(hess)){ 
      # If hess is not provided, approximate it by finite differencing. 
      # source: Statistical Programming lecture notes.
      g0 <- grad(theta,...) # gradient of the function at theta
      h <- matrix(0,length(theta),length(theta)) # initialise Hessian matrix
      
      ## loop over parameters
      for (i in 1:length(theta)) { 
        theta1 <- theta; theta1[i] <- theta1[i] + eps ## increase theta[i] by eps
        g1 <- grad(theta1,...) ## compute resulting gradient
        h[i,] <- (g1 - g0)/eps ## approximate second derivatives
      }
      
    }
    else {h <- hess(theta,...)} #If it is provided, define the current hessian as h.
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
        print(h)
        print(diag(length(theta)))
        h <- h + diag(length(theta))*multiplier #Perturbing by multiplying with multiplier.
        R <- try(chol(h), silent = TRUE) #Checking positive definiteness again, which will be
        #the argument inside of the while loop.
        multiplier <- 10*multiplier # Updating multiplier if it's still not pos.def.
      }
    }
    
    delta <- backsolve(R, forwardsolve(t(R), -g)) # Finding the descent direction.
    if (sum(abs(g) < tol*(abs(f) + fscale)) == length(theta)) {break} 
    # Judging convergence. If it converges, then we break the loop.
    
    if (func(theta + delta,...) >= f) { 
      #If our descent direction can't improve func, we halve the direction vector.
      for(j in 1:max.half) { #Halving max max.half times.
        delta <- delta/2 #Halving
        if(func(theta + delta,...) < f & is.finite(func(theta+delta,...)) == TRUE)
          {break} #If the direction starts improving, we break halving.
      }
    }
    
    
    if (func(theta + delta,...) >= f) {
      # warning 2 (step fails to reduce func after trying max.half step halvings)
      warning("max.half is reached without convergence.")
      break # we break the optimization iterations if warning 2 happens
      } 
    
    
    if(is.finite(func(theta+delta,...)) == FALSE){
      #warning 5 (If objective is not finite at iteration)
      warning("max.half is reached with the trial objective function not finite", 
              call. = FALSE)}
    
    theta <- theta + delta #After finding a decreasing direction, update theta.
    iter <- iter + 1 # As we updated theta, we increase the total number of iterations.
  }
  
  
  if ((sum(abs(g) < tol*(abs(f) + fscale)) != length(theta)) & iter == maxit) {
    # warning 3 (maxit reached without convergence)
    warning("maxit is reached without convergence.") 
    # all elements of the gradient vector still do not converge to zero after maxit
  } 

  
  if (!inherits(try(chol(hf), silent = TRUE), "try-error")) { 
    # Check whether last non-perturbed hessian matrix is positive definite by using
    # Cholesky decomposition.
    Rf <- try(chol(hf), silent = TRUE)
    # If it is pos. definite, find its inverse by Cholesky decomp and return it
    Hi <- backsolve(Rf, forwardsolve(t(Rf), diag(length(theta)))) 
  }
  
  
  else{
    # warning 4 (hess is not pos. definite)
    Hi <- NULL #If the last hessian is not pos. definite, do not return it.
    warning("Hessian matrix is not positive definite", call. = FALSE)
  }
  
  # Returning estimations
  return(list("f"=f,"theta"=theta, "iter"=iter, "g"=g, "Hi"=Hi)) 
}

theta <- c(1,100)
n1 <- newt(theta,rb,gb,hess=hb,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
n2 <- newt(theta,rb,gb,hess=NULL,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
system.time(newt(theta,rb,gb,hess=hb,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6))




installed.packages(priority="base")

installed.packages(priority="recommended")
