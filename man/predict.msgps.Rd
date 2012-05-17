\name{predict.msgps}
\alias{predict.msgps}
\alias{coef.msgps}
\alias{coef.dfgps}
\alias{coef.step.dfgps}
\alias{coefmat.dfgps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make predictions from a "msgps" object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions predicts fitted values via msgps function.  
}
\usage{
\method{predict}{msgps}(object, X, tuning,...)
\method{coef}{msgps}(object, tuning,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
Fitted \code{"msgps"} model object.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
Matrix of vector of new input \code{x}.
}
  \item{tuning}{
%%     ~~Describe \code{tuning} here~~
Tuning  parameter vector \code{t} where predictions are required.  If \code{tuning} is missing, solutions selected by Cp, bias-corrected AIC (AICC), generalized cross validation (GCV) and BIC are displayed.}
  \item{\dots}{Other parameters}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The object returned depends on type.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kei Hirose\cr
\email{mail@keihirose.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
#data
X <- matrix(rnorm(100*8),100,8)
beta0 <- c(3,1.5,0,0,2,0,0,0)
epsilon <- rnorm(100,sd=3)
y <- X \%*\% beta0 + epsilon
y <- c(y)

#fit
fit <- msgps(X,y)
coef(fit) #extract coefficients at t selected by model selection criteria
coef(fit,c(0, 0.5, 2.5)) #extract coefficients at some values of t
predict(fit,X[1:10,]) #predict values at t selected by model selection criteria
predict(fit,X[1:10,],c(0, 0.5, 2.5)) #predict values at some values of t
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
