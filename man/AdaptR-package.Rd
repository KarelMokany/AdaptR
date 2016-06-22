\name{AdaptR-package}

\alias{AdaptR-package}

\docType{package}

\title{ 
Overview of the functions in the AdaptR package
}

\description{
AdaptR is a dynamic modelling technique for incorporating the capacity for genetic adaptation in projecting changes in the distributions of species over time. The \pkg{AdaptR} package currently provides basic functions to prepare the imput files, run the model, and generate some simple summary plots of the outputs. Future updates will improve the flexibility of data inputs and results processing. 
}
\details{
The functions in the \pkg{AdaptR} package provide the tools necessary for running the AdaptR model, including functions to prepare the environmental data. Major functionality includes:
  \itemize{
    \item Compressing temporal series of environmental grids for multiple variables to \pkg{AdaptR}'s compressed grid format, used for efficient model implementation  
    \item running AdaptR
    \item Plotting key outputs
  }

To see the preferable citation of the package, type \code{citation("AdaptR")}.  
}

\author{The \pkg{AdaptR} development team is Karel Mokany and Alex Bush. Where others have contributed to
individual functions, credits are provided in function help pages.

The maintainers of the R version of \pkg{AdaptR} are Karel Mokany <Karel.Mokany@csiro.au> and Alex Bush <Alex.Bush@csiro.au>.
}
