% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classification.R
\name{load_brcaness_classifier}
\alias{load_brcaness_classifier}
\title{Load pre-trained BRCAness random forest model}
\usage{
load_brcaness_classifier()
}
\value{
A trained random forest model object
}
\description{
This function loads a pre-trained random forest model for BRCAness classification from a saved .rda file.
}
\examples{
brcaness_model <- load_brcaness_classifier()
}
