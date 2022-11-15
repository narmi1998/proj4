library(MASS)

rb <- function(th,k=2) {
  k*(th[2]-(th[1])^2)^2 + (1-th[1])^2 
  }
theta <- c(35, 35)

f <- rb(theta)
h<-hb(theta)
g <- gb(theta)

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
  f <- func(theta)
  if (is.infinite(f)){stop("The objective value is infinite.")}
  g <- grad(theta)
  if (sum(is.infinite(g)) != 0) {stop("At least one derivative value is infinite.")}

  iter <- 0
  
  for (i in 1:maxit) {
    f <- func(theta)
    g <- grad(theta)
    if (is.null(hess)){
      g0 <- gb(theta) ## grad of nll at th0
      h <- matrix(0,length(theta),length(theta)) ## finite diference Hessian
      for (i in 1:length(theta)) { ## loop over parameters
        theta1 <- theta; theta1[i] <- theta1[i] + eps ## increase th0[i] by eps
        g1 <- gb(theta1) ## compute resulting nll
        h[i,] <- (g1 - g0)/eps ## approximate second derivs
      }
    }
    else {h <- hess(theta)}

    if (sum(eigen(h)$values>=0) != length(theta)) {
      multiplier <- norm(h)*1e-6
      while(sum(eigen(h)$values>=0) != length(theta)) {
        h <- h + diag(length(theta))*multiplier
        multiplier <- 10*multiplier 
        }
    }
    
    delta <- -ginv(h)%*%g
    if (sum(abs(g) < tol*(abs(f) + fscale)) == length(theta)) {break} 
        
    if (func(theta + delta) >= f) {
      for(j in 1:max.half) {
        delta <- delta/2
        if(func(theta + delta) < f) {break}
      }
    }
    if (func(theta + delta) >= f) {warning("max.half is reached without convergence.")}
    theta <- theta + delta
    iter <- iter+1
  }
  if (sum(eigen(h)$values>=0) != length(theta)) {
    warning("The hessian at convergence is not positive definite.")
    Hi <- NULL}
    
  else {Hi <- ginv(h)}
  return(list("f"=f,"theta"=theta, "g"=g, "iter"=iter, "Hi"=Hi))
  
}

n1 <- newt(theta,rb,gb,hess=hb,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
n2 <- newt(theta,rb,gb,hess=NULL,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)

  
installed.packages(priority="base")

installed.packages(priority="recommended")
