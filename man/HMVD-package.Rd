\name{HMVD-package}
\alias{HMVD-package}
\docType{package}
\title{
Group Association Test using a Hidden Markov Model
}
\description{
HMVD performs an association test between a group of variables and the outcome. 
Posterior probabilities are provided for each variable indicating how likely each variable is associated with the outcome.
}
\details{
\tabular{ll}{
Package: \tab HMVD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-05-12\cr
License: \tab GPL-3 \cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Maintainer: Yichen Cheng<ycheng@fredhutch.org>
}
\references{
Cheng, Y., Dai, J. and Kooperberg, C. (2015). Group association test using hidden Markov model. 
Biostatistics, in pres.
}
\keyword{ package }
\examples{
#############################################################################
#### compute the p-value and do parameter estimation for continuous outcome 
n = 4000; m = 20
X = matrix(rnorm(n*m),n)
Y = rowMeans(X[,1:4])*.2 + rnorm(n)
HMVD(Y,X)$p.value #### approximate p-value

HMVD(Y,X,nperm.max = 20)$p.value.perm #### p-value based on permutations
#### in practice we would use more permutations

out = HMVD(Y,X,method='estimation')
round(out$prob,2) ###posterior probability
out$theta ### common effect size

#### compute the p-value and do parameter estimation for binary outcome 
n = 4000; m = 20
X = matrix(rnorm(n*m),n)
p = rowMeans(X[,1:4])*.4
Y = rbinom(n,1,p = exp(p)/(1+exp(p)))
HMVD(Y,X,model.type='D')$p.value #### approximate p-value

HMVD(Y,X,nperm.max = 20)$p.value.perm #### p-value based on permutations
#### in practice we would use more permutations

out = HMVD(Y,X,model.type='D',method='estimation')
round(out$prob,2) ###posterior probability
out$theta ### common effect size
}
