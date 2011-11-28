\name{summary.msgps}
\alias{summary.msgps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A summary of  "msgps" object..
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions summarizes the "msgps" object.  
}
\usage{
\method{summary}{msgps}(object, digits=max(3, getOption("digits") - 3), num.result = 20, 
coef.result=100,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
Fitted \code{"msgps"} model object.
}
   \item{digits}{
%%     ~~Describe \code{object} here~~
The digits of the output.
}
 \item{num.result}{
%%     ~~Describe \code{num.result} here~~
The number of tuning parameter and the corresponding degrees of freedom displayed in this code. 
}
  \item{coef.result}{
%%     ~~Describe \code{coef.result} here~~
If the coef.result exceeds the number of variables, the result of coefficient is not described in this code. 
}
  \item{\dots}{Other parameters on summary}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{df}{The degrees of freedom for each tuning parameter.}
\item{tuning.max}{Maximum value of tuning parameter.}
\item{ms.coef}{The coefficient selected by each model selection criterion.}
\item{ms.tuning}{The values of tuning parameter of models selected by each model selection criterion.}
\item{ms.df}{The degerees of freedom selected of models each model selection criterion.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
summary(fit) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
