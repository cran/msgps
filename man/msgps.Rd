\name{msgps}
\alias{msgps}
\alias{dfgps}
\alias{aicc.dfgps}
\alias{bic.dfgps}
\alias{cp.dfgps}
\alias{gcv.dfgps}
\alias{print.msgps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
msgps (Degrees of Freedom of Elastic Net, Adaptive Lasso and Generalized Elastic Net)
}
\description{
This package computes the degrees of freedom of the lasso, elastic net, generalized elastic net and adaptive lasso based on the generalized path seeking algorithm.  The optimal model can be selected by model selection criteria including Mallows' Cp, bias-corrected AIC (AICc), generalized cross validation (GCV) and BIC.
}
\usage{
msgps(X,y,penalty="enet", alpha=0, gamma=1, lambda=0.001, tau2, STEP=20000, 
STEP.max=200000,  DFtype="MODIFIED",  p.max=300, intercept=TRUE, stand.coef=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
predictor matrix
}
  \item{y}{
%%     ~~Describe \code{y} here~~
response vector
}
  \item{penalty}{
%%     ~~Describe \code{penalty_type} here~~
The penalty term.  The \code{"enet"} indicates the elastic net:
 \deqn{\alpha/2||\beta||_2^2+(1-\alpha)||\beta||_1.} 
Note that \code{alpha=0} is the lasso penalty. The \code{"genet"} is the generalized elastic net: 
 \deqn{log(\alpha+(1-\alpha)||\beta||_1).} 
The \code{"alasso"} is the adaptive lasso, which is a weighted version of the lasso given by
 \deqn{w_i||\beta||_1,}
 where  \eqn{w_i} is \eqn{1/(\hat{\beta}_i)^{\gamma}}. Here \eqn{\gamma>0} is a tuning parameter, and \eqn{\hat{\beta}_i} is the ridge estimate with regularization parameter being \eqn{\lambda \ge 0}.  
}
  \item{alpha}{
%%     ~~Describe \code{ex_para} here~~
The value of \eqn{\alpha} on \code{"enet"} and \code{"genet"} penalty. 
}
  \item{gamma}{
%%     ~~Describe \code{ex_para} here~~
The value of \eqn{\gamma} on \code{"alasso"}. 
}
  \item{lambda}{
%%     ~~Describe \code{ex_para} here~~
The value of  regularization parameter \eqn{\lambda \ge 0} for ridge regression, which is used to calculate the weight vector of  \code{"alasso"} penalty.  Note that the ridge estimates can be ordinary least squared estimates when \code{lambda=0}.
}
  \item{tau2}{
%%     ~~Describe \code{tau2} here~~
Estimator of error variance for Mallows' Cp.  The default is the unbiased estimator of error vairance of the most complex model.  When the unbiased estimator of error vairance of the most complex model is not available (e.g., the number of variables exceeds the number of samples), \code{tau2} is the variance of response vector.
}
  \item{STEP}{
%%     ~~Describe \code{STEP} here~~
The approximate number of steps. 
}
  \item{STEP.max}{
%%     ~~Describe \code{STEP.max} here~~
The number of steps in this algorithm can often exceed \code{STEP}.  When the number of steps exceeds \code{STEP.max}, this algorithm stops.  
}
  \item{DFtype}{
%%     ~~Describe \code{DFtype} here~~
\code{"MODIFIED"} or \code{"NAIVE"}.  The \code{"MODIFIED"} update is much more efficient thatn \code{"NAIVE"} update.
}
  \item{p.max}{
%%     ~~Describe \code{p.max} here~~
If the number of selected variables exceeds \code{p.max}, the algorithm stops.
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
When intercept is \code{TRUE}, the result of intercept is included.
}
  \item{stand.coef}{
%%     ~~Describe \code{stand.coef} here~~
When stand.coef is \code{TRUE}, the standardized coefficient is displayed.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
%% ~put references to the literature/web site here ~
Friedman, J. (2008).  Fast sparse regression and classification. \eqn{Technical report}, Standford University.\cr
Hirose, K., Tateishi, S. and Konishi, S.. (2011).  Efficient algorithm to select tuning parameters in sparse regression modeling with regularization. arXiv:1109.2411 (arXiv).
}
\author{
Kei Hirose\cr
\email{mail@keihirose.com}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{coef.msgps}, \code{plot.msgps}, \code{predict.msgps} and \code{summary.msgos} objects.}
\examples{
#data
X <- matrix(rnorm(100*8),100,8)
beta0 <- c(3,1.5,0,0,2,0,0,0)
epsilon <- rnorm(100,sd=3)
y <- X \%*\% beta0 + epsilon
y <- c(y)

#lasso
fit <- msgps(X,y)
summary(fit) 
coef(fit) #extract coefficients at t selected by model selection criteria
coef(fit,c(0, 0.5, 2.5)) #extract coefficients at some values of t
predict(fit,X[1:10,]) #predict values at t selected by model selection criteria
predict(fit,X[1:10,],c(0, 0.5, 2.5)) #predict values at some values of t
plot(fit,criterion="cp") #plot the solution path with a model selected by Cp criterion

#elastic net
fit2 <- msgps(X,y,penalty="enet",alpha=0.5)
summary(fit2) 

#generalized elastic net
fit3 <- msgps(X,y,penalty="genet",alpha=0.5)
summary(fit3)

#adaptive lasso
fit4 <- msgps(X,y,penalty="alasso",gamma=1,lambda=0)
summary(fit4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
