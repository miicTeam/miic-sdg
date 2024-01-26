# MIIC-SDG
This repository contains the R source code for MIIC-SDG, an algorithm for the generation of synthetic data that mimic the original data. MIIC-SDG have been shown (submitted paper) to generate data with a promising balance between data re-usability (corresponding to data quality) and privacy preservation. The MIIC-SDG package contains 2 functions:
- miicsdg: for the generation of the synthetic data
- writeToFile: to write to a folder the two files representing the DAG network, that can be visualized at the following address: https://miic.curie.fr/vis_NL.php


## Installation
From GitHub:
```R
# install.packages("remotes")
remotes::install_github("miicTeam/miicsdg")
```

## Quick start
```R
library(miicsdg)
library(datasets)
data(iris)

miicsdg_list = miicsdg(iris)

miicsdg_synthetic_data = miicsdg_list[['synthetic_data']]
miicsdg_adjacencyMatrixDAG = miicsdg_list[['adjacency_matrix_DAG']]
miicsdg_dataTypes = miicsdg_list[['data_types']]
miicsdg_edgesMiicServer = miicsdg_list[['edges_miic_server']]
miic_sdg_miic = miicsdg_list[['miic']]

miicsdg::writeToFile(miicsdg_list , '~/testMIICSDG')
```

## Authors
- Sella Nadir
- Florent Guinot
- Nikita Lagrange
- Laurent-Philippe Albou
- Jonathan Desponds
- Hervé Isambert
