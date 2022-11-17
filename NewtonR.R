newt<-function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
{
  
}


# return
f <- 1 #the value of the objective function at the minimum.
t <- 1 #theta the value of the parameters at the minimum.
iter<- 1 #iter the number of iterations taken to reach the minimum.
g <- 1 #g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
Hi <- 1 #Hi the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).

?deriv
install(Ryacas0)
library(Ryacas0)
library(calculus)
# test functions
subst <- function(e, ...) do.call(substitute, list(e[[1]], list(...)))
rb <- function(th,k=2) {
  alpha <- th[1]
  beta <- th[2]
  #expression(k*(th[2]-th[1]^2)^2 + (1-th[1])^2)
  f <- (expression(k*(beta-alpha^2)^2 + (1-alpha)^2))
  substitute(k*(beta-alpha^2)^2 + (1-alpha)^2, list(k = 2)) # or k = k
}
# testing
ff<-force(expression(x^2))
subst(ff, x = 2)
D <- funcx <- force(rb(c('x','y'),k=2))
grad1 <- force(D(force(funcx), 'alpha')) # force it ***
grad1 <- substitute((-(2 * (1 - alpha) + 2 * (2 * (2 * alpha * (beta - alpha^2))))), list(alpha = 4, beta = 0.1)) ###
grad1 <- substitute(expression(force(grad1)), alpha = 4, beta = 0.1)
grad2 <- D(funcx, 'beta')
grad2 <- substitute((2 * (2 * (beta - alpha^2))), list(alpha = 4, beta = 0.1)) ##!!!
G <- c(grad1, grad2)
nll <- deriv(funcx, c("alpha","beta"), function.arg=c("alpha","beta"), hessian=TRUE)
H <- attr(nll, "hessian")
v <- -H^(-1) %*% G
## dist(V) = .......
quad <- D + t(V) %*% G + 0.5 * t(V) %*%H%*%V  ## + dist(V)
GV <- c(-(2 * (1 - alpha) + 2 * (2 * (2 * alpha * (beta -alpha^2)))), 2 * (2 * (beta - alpha^2)))
GV <- force(c(force(D(force(funcx), 'alpha')), D(funcx, 'beta')))
class(GV[[1]])
D(GV[[1]],'alpha')
symbols
abc <- deparse(GV[[1]])
a <- expression(10 + x + y)
D(GV[[1]],'alpha')

########################
##### HESSIAN
th <- C(4,0.1)
n <- length(th)
symbols <- c('alpha','beta')
D <- funcx <- force(rb(c('x','y'),k=2))
H <- matrix(0,2,2)
GV <- force(c(force(D(force(funcx), 'alpha')), D(funcx, 'beta')))
for (i in 1:n)
{
  a <- (D(GV[[i]], symbols[1]))
  b <- (D(GV[[i]], symbols[2]))
  c<- eval(do.call("substitute", list(a, list(alpha = th[1],beta = th[2]))))
  d <- eval(do.call("substitute", list(b, list(alpha = th[1],beta = th[2]))))
  H[i,] <- c(c,d)
}

H
#######################

D(funcx,'alpha')
H <- array(0,dim = c(2,2)) # 2 = length of th
H[1,1] <- D(D(funcx,'alpha'),'alpha')
H[1,2] <- H[2,1] <- D(D(funcx,'alpha'),'beta')
H[2,2] <- D(D(funcx,'beta'),'beta')


Hessian <- function(th, gv, symbols)
{
  n <- length(th)
  mat <- t(matrix(GV,n,n))
  
  for (i in 1:n)
  {
    mat[,i] <- dev(mat[,i], symbols[i])
  }
  # sub loop ***
  mat <- substitute(mat, list(symbols[1] = th[1], symbols[2] = th[2]))
  return(mat)
  
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

D(abc,'x')
abc <- hb(x,k=2)
f=expression(x^2 + 2)
D(f,'x')
theta <- c(1,2)

ex <- expression(hb(x,k=2))
D(ex,'x')
funcx <- rb(c(4,0.1),k=2)



nlli <- deriv( ## what to differentiate...
  expression(-y*(log(alpha)+beta*t)+alpha*exp(beta*t)+lgamma(y+1)),
  c("alpha","beta"), ## differentiate w.r.t. these
  function.arg=c("alpha","beta","t","y"),## return function - with these args
  hessian=TRUE) ## compute the Hessian as well as the gradient



nll4 <- function(th,t,y) {
  ## negative log lik function for optimization by 'nlm', where computations are done
  ## using nlli, produced by 'deriv'. th=c(alpha,beta).
  nli <- nlli(th[1],th[2],t,y) ## call nlli to get -ve log lik and derivative terms
  nll <- sum(nli) ## sum negative log likelihood terms
  attr(nll,"gradient") <- colSums(attr(nli,"gradient")) ## sum term grad vectors
  attr(nll,"hessian") <- apply(attr(nli,"hessian"),c(2,3),sum) ## same for Hessians
  nll
} 
th0 <- c(10,.1)
t80 <- 1:13 
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases
nnn <- nll4(th0,t80,y)
attr(nlli(th0[1],th0[2],t80,y),"gradient")
abc <-nlli(th0[1],th0[2],t80,y)

