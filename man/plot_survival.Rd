% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_survival.R
\name{plot_survival}
\alias{plot_survival}
\title{Plot a Kaplan-Meier Curve}
\usage{
plot_survival(
  dbms = c("oracle", "postgresql", "redshift", "sql server", "pdw", "netezza",
    "bigquery", "sqlite"),
  user = NULL,
  password = NULL,
  server,
  port,
  schema = NULL
)
}
\arguments{
\item{dbms}{one of the following supported Database Management Systems: c("oracle","postgresql","redshift","sql server","pdw", "netezza","bigquery","sqlite")}

\item{user}{username for OMOP v5.0 instance}

\item{password}{password for OMOP v5.0 instance}

\item{server}{OMOP v5.0 server}

\item{port}{NULL by default.}

\item{schema}{NULL by default.}
}
\description{
Plot a Kaplan-Meier Curve
}
