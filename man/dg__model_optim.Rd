% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dg__model_optim.R
\name{dg__model_optim}
\alias{dg__model_optim}
\title{calls the optimzation routine to estimate dG models}
\usage{
dg__model_optim(start_par, parlist, varlist, maxit, trace_optim = FALSE)
}
\arguments{
\item{start_par}{vector of starting parameters for modeling}

\item{parlist}{list with all model parameters}

\item{varlist}{list with all variables}

\item{maxit}{maximum number of iterations for optim}

\item{trace_optim}{logical, if TRUE, shows trace = 3, default: FALSE}
}
\value{
returns a model in data.table format
}
\description{
calls the optimzation routine to estimate dG models
}
