hess_fd <- function(theta)
{
  gv <- grad(theta)
  hess <- matrix(0, length(theta), length(theta))
  for (i in 1:length(theta)){
    theta1 <- theta
    theta1[i] <- theta1[i] + eps
    gv1 <- grad(theta1)
    hess[i, ] <- (gv1 - gv)/eps
  }
  hess
}

hess_pd <- function(hess)
{
  R <- try(chol(hess), silent = TRUE) # try to decompose the matrix by Cholesky decomposition
  if(inherits(R, "try-error")) # if decomposition fails
  {
    multiple <- 1e-6 + norm(hess)
    hess <- hess + multiple*diag(dim(hess)[1])
    R <- try(chol(hess))
    while(inherits(R, "try-error"))
    {
      multiple <- multiple*10
      hess <- hess + multiple*diag(dim(hess)[1])
      R <- try(chol(hess))
    }
  }
  list(hess_pd = hess_pd, R = R)
}


newt <- function(theta, func, grad, hess, k=2, tol=1e-8, fscale=1, maxit=100, max.half=20, eps=1e-6)
{
  if (!is.finite(func(theta)))
  { stop("The objective is not finite at the initial theta", call. = FALSE) }
  if (!all(is.finite(grad(theta))))
  { stop("The derivatives are not finite at the initial theta", call. = FALSE) }
  
  iter = 1
  while(iter <= maxit)
  {
    if (is.null(hess)){
      hess <- hess_fd(theta) # the finite difference approximation to hessian matrix  
      hess <- (t(hess) + hess)/2
    }
    else{
      hess <- hess(theta)
    }
    hess_pd <- hess_pd(hess) # check if hessian matrix is positive definite
  
    dd <- backsolve(hess_pd$R, forwardsolve(t(hess_pd$R), -grad(theta)))
    halve_times <- 1 # count for number of times halving the descent direction
    while(func(theta+dd) >= func(theta) & halve_times <= max.half)
    {
      dd <- dd/2
      halve_times <- halve_times + 1
    }
    if(func(theta+dd) >= func(theta)){
      warning("fail to reduce the objective function", call. = FALSE)
    }
    
    theta <- theta + dd
    criterion <- tol*func(theta)+fscale
    if(all(abs(grad(theta)) < criterion))
    {
      f <- func(theta)
      g <- grad(theta)
      
      if (is.null(hess)){
        hess <- hess_fd(theta) # the finite difference approximation to hessian matrix  
        hess <- (t(hess) + hess)/2
      }
      else{
        hess <- hess(theta)
      }
      hess_pd <- hess_pd(hess) # check if hessian matrix is positive definite
      if (hess_pd == hess){
        Hi <- backsolve(hess_pd$R, forwardsolve(t(hess_pd$R), diag(dim(hess_pd)[1])))
        newton <- list(f = f, theta = theta, iter = iter, g = g, Hi = Hi)
        return(newton)
      }
      else{
        newton <- list(f = f, theta = theta, iter = iter, g = g)
        return(newton)
        warning("Hessian matrix is not positive definite", call. = FALSE)
      }
      break
    }
    
    iter <- iter + 1
  }
  
  if (iter > maxit)
  { stop("fail to reach convergence", call. = FALSE) }
}

rb <- function(th, k=2)
{
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <- function(th, k=2)
{
  c(-2*(1-th[1]) - k*4*th[1]*(th[2]-th[1]^2), k*2*(th[2]-th[1]^2))
}

hb <- function(th, k=2){
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

print(newt(c(10, 0.1), func=rb, grad=gb, hess=hb))

