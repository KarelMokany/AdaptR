\name{AdaptR}
\alias{AdaptR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Species distribution modelling with genetic adaptation
}
\description{
Incorporate genetic adaptation into the projected distribution of a species.
}
\usage{
AdaptR(Pfile = Parameter.File.Path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pfile}{
%%     ~~Describe \code{Pfile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns output files in the specified output directory
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Karel Mokany & Alex Bush
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Pfile = Parameter.File.Path) 
{
    .C("AdaptR", argv = as.character(c(Pfile)), arg_i_catch = as.integer(c(0, 
        0)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
