% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miic-sdg.R
\name{miicsdg}
\alias{miicsdg}
\title{Create synthetic data trough the MIIC algorithm}
\usage{
miicsdg(
  original_data,
  n_synthetic_samples = NA,
  method_DAG = "MIIC_to_DAG",
  root_index = 1,
  ids_columnsDiscrete = NA,
  state_order = NA,
  keep_continuous_variable_boundaries = T,
  confidence_threshold = 1,
  consistent = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{original_data}{A data frame containing the original data we want to mimic}

\item{n_synthetic_samples}{An integer representing the number of synthetic samples to be created. Default value: equal to the number of samples in the input data}

\item{method_DAG}{A string ('DFS', 'CPDAG_to_DAG', 'MIIC_to_DAG') corresponding to the method used to transform the MIIC network into a DAG}

\item{root_index}{An integer corresponding to the position of the feature in the original data frame that should be used as root in the DFS visit if DFS is chosen as method to build the DAG}

\item{keep_continuous_variable_boundaries}{A boolean, if TRUE, original boundaries are kept for continuous variables}

\item{confidence_threshold}{A floating point number corresponding to the threshold used to filter the less probable edges following the skeleton step. See the MIIC algorithm for more details}

\item{consistent}{A Boolean. If TRUE the method iterates over skeleton and orientation steps to ensure consistency of the network. See the MIIC algorithm for more details}

\item{verbose}{A boolean indicating if the execution should be documented with more informations}

\item{ids_columns_discrete}{A vector representing the index of features that should be used as non-numerical}
}
\value{
A list containig 5 objects: \itemize{
\item{\emph{synthetic_data}:}{ a data frame containing the synthetic data}
\item{\emph{adjacency_matrix_DAG}:}{ an adjacency matrix representing the dag used for creating the data (to read in the sense row -> column)}
\item{\emph{data_types}:}{ a data frame containing the information on the natture of variables (0 = discrete, 1 = continuous) useful for the visualization of the DAG on the miic web server}
\item{\emph{edges_miic_server}:}{ the edges file that can be used for the visualization of the DAG on the miic web server. You can use the function miicsdg::writeToFile to save this file plus the data_types file}
\item{\emph{miic}:}{ the object returned from running miic}
}
}
\description{
Create synthetic data trough the MIIC algorithm
}
\examples{
library(miicsdg)
library(datasets)
data(iris)

miicsdg_list = miicsdg::miicsdg(iris)

miicsdg_synthetic_data = miicsdg_list[['synthetic_data']]
miicsdg_adjacencyMatrixDAG = miicsdg_list[['adjacency_matrix_DAG']]
miicsdg_dataTypes = miicsdg_list[['data_types']]
miicsdg_edgesMiicServer = miicsdg_list[['edges_miic_server']]
miic_sdg_miic = miicsdg_list[['miic']]

miicsdg::writeToFile(miicsdg_list , '~/testMIICSDG')


}
