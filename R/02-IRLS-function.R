
rats <- read_csv(here("data", "rats.csv"))

# p^x (1-p)^1-x ; 0 <= p <= 1, x = 0,1,2,3,...,n
# E[x] = np = 5 * (1-exp{-b0 - b1x - b2 x^2})
# var(x) = np(1-p)

# write the model as a linear function of log(E[x])
# # log(mu) = B^TX
# mu = exp(B^TX) > 0
# log(mu) = b0 + b1x1 + b2x2 + b3x3 + error
# E(yi) = mu_i = exp(B^Txi)
# var(yi) = sig_i^2 = exp(B^TXi)
# Xi = design matrix with column of 1's if intercept
# yi = exp(B^TXi) + ei
# wi = 1/exp(B^Txi)

JacfW_binom = function(beta,X){
    # Required elements for Poisson regression
    # beta is a vector of regression parameters (including beta_0,...,beta_k)
    # X is a matrix of observed covariates (n by k), including n observations on k covariates
    # Note, if there is an intercept,X must contain a column of 1's in its first column
    # On output, E(y)=f, the Jacobian, and a weight matrix (inverse of variances) is output
    xb = -(X %*% beta) # n by 1
    f = 5 * (1 - exp(xb)) # n by 1
    J = diag(as.vector(f)) %*% X
    var = exp(xb)*f
    W = diag(1/as.vector(var))
    list(f=f, J=J, W=W)
}

X <- rats %>% select(`dose-ppm`, `animals-tested`) %>%  mutate(inter = 1, .before=1)
X1 <- as.matrix(X)
X1[,3] <- X1[,3]^2
y <- rats %>% select(`tumor-incidence`)
y1 <-as.matrix(y)


b0 <- c(0.009960264,0.003557714,-5.19416e-05)
b0 <- c(5.82531,-0.01123,-0.04826)

# xb = -(X1 %*% b0); xb
# f = 5 * (1 - exp(xb)); f
# var = exp(xb)*f; var
# diag(1/as.vector(var))




# JacfW_binom(b0,X1)

#b0 = b1 = b2 = 1
# p = 1 - exp(-0.009960264 - 0.003557714*X1 - -5.19416e-05 * X1^2)
# n = nrow(X1)
# 
# 
# n %*% p

# xb = -(1 + X1 + X1^2) %*% b0 # n by 1
# f = 5 * (1 - exp(xb)) # n by 1
# J = diag(as.vector(f)) %*% X1
# var = exp(xb)*f
# W = diag(1/as.vector(var))



GN = function(y,X, beta0, Jac, Wt = 1, maxit, IRLS = TRUE){
    # y is the response variable
    # X is a matrix of covariates
    # beta0, is initial value for model parameters
    # Note, if there is an intercept,X must contain a column of 1's in its first column
    # Jac is name of a function that computes the expectation (Jac$f), the n by p 
    # Jacobian matrix (Jac$J), where n is the
    # number of cases and p is the number of parameters, It should have two input variables
    # beta and X, where beta is a vector of parameters. 
    # If IRLS is TRUE, Jac should also compute the weigh 
    # matrix, usually inverse of the variance, (Jac$Wt) as output (an n by n matrix).
    # Note, output these in the Jac function as a lst with names f, J, abd W
    # If IRLS is FALSE, Weight is an n by n matrix on entry. The default value of 
    # of W in this case is the identity matrix (no weight)
    # ... allows entry of the parameters needed for Jac and Wt
    
    fJ = function(beta,X) Jac(beta,X)
    for (it in 1: maxit){
        a = fJ(beta0,X)
        J = a$J
        f = a$f
        if (IRLS==TRUE){
            Wt = a$W
        }
        JW = t(J) %*% Wt
        print(sprintf('it = %3.0f   b0 = %6.6f b1 = %6.6f b2 = %6.6f, grad=%2.2e',
                      it,beta0[1],beta0[2],beta0[3],norm(JW %*% (y-f))))
        dir = solve(JW%*%J)%*%JW %*% (y-f)
        # if IRLS = TRUE, no step-having
        beta0 = beta0 + dir
    }
    beta0
}

maxit = 20
b0 <- c(0.009960264,0.003557714,-5.19416e-05)
a <- GN(y1,X1, b0, JacfW_binom, Wt = 1, maxit, IRLS = TRUE)

