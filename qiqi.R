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
    multiplier <- 1e-6 + norm(hess)
    while(inherits(R, "try-error"))
    {
      hess <- hess + multiplier*diag(length(theta))
      R <- try(chol(hess))
      multiplier <- multiple*10
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
      hm <- hess_fd(theta) # the finite difference approximation to hessian matrix  
      hm <- (t(hm) + hm)/2
    }
    else {hm <- hess(theta)}
    
    hess_pd <- hess_pd(hm) # check if hessian matrix is positive definite
    f <- func(theta)
    g <- grad(theta)
    
    criterion <- tol*(abs(f)+fscale)
    if(all(abs(g) < criterion)) {break}
  
    dd <- backsolve(hess_pd$R, forwardsolve(t(hess_pd$R), -g))
    for (i in 1: max.half)
    {
      dd <- dd/2
      a = func(theta+dd)
      if(func(theta+dd) < f & all(is.finite(func(theta+dd))) == TRUE) {break}
    }
    if(func(theta+dd) >= f){
      warning("fail to reduce the objective function", call. = FALSE)}
    if(all(is.finite(func(theta+dd))) == FALSE){
      warning("fail to let the objective function be finite", call. = FALSE)}
    
    theta <- theta + dd
    iter <- iter + 1
  }
  
  if (iter > maxit)
  { stop("fail to reach convergence", call. = FALSE) }
  
  if (!inherits(try(chol(hm), silent = TRUE), "try-error")) {
    Hi <- backsolve(hess_pd$R, forwardsolve(t(hess_pd$R), diag(length(theta))))
    
  }
  else{
    Hi <- NULL
    warning("Hessian matrix is not positive definite", call. = FALSE)
  }
  
  newton <- list(f = f, theta = theta, iter = iter, g = g, Hi = Hi)
  return(newton)
  
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

print(newt(theta=c(35,35), func=rb, grad=gb, hess=hb))

