# MIIC-SDG
This repository contains the R source code for MIIC-SDG, an algorithm for the generation of synthetic data that mimic the original data. MIIC-SDG have been shown (submitted paper) to generate data with a promising balance between data re-usability (corresponding to data quality) and privacy preservation. The MIIC-SDG package contains 2 functions:
- miicsdg: for the generation of the synthetic data. The input parameters are described using ?miicsdg::miicsdg
- writeToFile: to write to a folder the two files representing the DAG network, that can be visualized at the following address: https://miic.curie.fr/vis_NL.php


## Installation
From GitHub:
```R
# install.packages("remotes")
remotes::install_github("miicTeam/miic-sdg")
```

## Quick start
```R
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
```
## Output
The object returned by MIIC-SDG contains 5 objects:
- the synthetic data as a data frame
- the DAG used to sample the synthetic data as an adjacency matrix
- a data frame containing the data type (discrete or continuous) of the input variables as a data frame
- a data frame representing the DAG object with a specific format; it can be used to visualize the network using the MIIC web server
- an object corresponding to the output of MIIC.

The returned network is very important in real scenarios when the investigator wants to inspect the set of associations between variables to better understand the structure of the data frame in terms of correlations. The graph can be uploaded in the https://miic.curie.fr/vis_NL.php web page. For this step it is possible to use the second function of the package, which saves two files to a folder defined by the user. The two outputs can be uploaded on the web-server.

## Authors
- Sella Nadir
- Florent Guinot
- Nikita Lagrange
- Laurent-Philippe Albou
- Jonathan Desponds
- Hervé Isambert
