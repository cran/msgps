\name{plot.msgps}
\alias{plot.msgps}
\alias{plot.df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot the solution path from a "msgps" object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions predicts fitted values from a "msgps" object.  
}
\usage{
\method{plot}{msgps}(x, criterion="cp", xvar="norm", yvar="coef", yvar.dflasso=TRUE, 
stand.coef=TRUE, plot.step = 1000, col=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{object} here~~
Fitted \code{"msgps"} model object.
}
%%  \item{type}{
%%     ~~Describe \code{tuning} here~~
%%The type of \code{"solution"} depicts the solution path along with the tuning parameter selected by model selection criterion of \code{criteion}. The type of \code{"df"} shows the degrees of freedom.
%%}
  \item{criterion}{
%%     ~~Describe \code{tuning} here~~
The code \code{criterion} plots the value of tuning parameter of each criterion (\code{"cp"}, \code{"aicc"}, \code{"gcv"}, \code{"bic"}).  The code \code{"none"} does not depict the tuning parameter.
}
  \item{xvar}{
%%     ~~Describe \code{tuning} here~~
The type of x variable.   \code{"xvar=norm"} is max|beta|/|beta|, \code{"xvar=sum"} is max|beta|, \code{"xvar=step"} is the number of steps, and  \code{"xvar=t"} is tuning parameter.
}
  \item{yvar}{
%%     ~~Describe \code{tuning} here~~
The type of y variable.   \code{"yvar=coef"} is the standardized coefficients, and  \code{"tvar=df"} is the degrees of freedom.  
}
  \item{yvar.dflasso}{
%%     ~~Describe \code{tuning} here~~
For lasso penalty, the degrees of freedom of the lasso (the number of non-zero parameters) is given when  \code{"yvar=df"} and \code{"yvar.dflasso=TRUE"}.
}
  \item{stand.coef}{
%%     ~~Describe \code{tuning} here~~
The standardized coefficients and tuning parameters are dipicted if "stand.coef=TRUE".
}
  \item{plot.step}{
%%     ~~Describe \code{tuning} here~~
The number of steps to plot the solution of df.  As \code{plot.step} increases, the picture will be well-looking whereas the file size of the picture will increase.
}
  \item{col}{
%%     ~~Describe \code{tuning} here~~
The color option.}
  \item{\dots}{Other graphical parameters to plot}
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

\seealso{\code{coef.msgps}, \code{predict.msgps} and \code{summary.msgps} objects.}

\examples{
#data
X <- matrix(rnorm(100*8),100,8)
beta0 <- c(3,1.5,0,0,2,0,0,0)
epsilon <- rnorm(100,sd=3)
y <- X \%*\% beta0 + epsilon
y <- c(y)

#fit
fit <- msgps(X,y)
plot(fit,criterion="cp") #plot the solution path with a model selected by Cp criterion
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
