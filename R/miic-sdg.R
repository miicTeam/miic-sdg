#Sort in ascending order the degrees of the nodes
find_cycles = function(g) {
  Cycles = NULL
  for(v1 in igraph::V(g)) {
    if(igraph::degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = igraph::neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(igraph::all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}

degree_sort <- function(matrix){
  ln <- ncol(matrix)
  degree <- vector("list",length=ncol(matrix))
  names(degree) <- names(matrix)
  for(i in 1:ln)
  {
    count <- 0
    for(j in 1:ln)
    {
      vl <- matrix[i,j]
      if(vl==-2){
        count <- count+3
      }
      if(vl==1)
      {
        count <- count+2
      }
      if(vl==2)
      {
        count <- count+0.5
      }
    }
    degree[[i]] <- count
  }
  degree <- degree[order(unlist(degree),decreasing = FALSE)]
  return(degree)
}

#Convert miic matrix
convert_matrix <- function(matrix)
{
  nc <- ncol(matrix)
  new_matrix <- matrix(0,nc,nc)
  for(i in 1:nc)
  {
    for(j in 1:nc)
    {
      if(matrix[i,j]==2)
      {
        new_matrix[i,j] <- 1
      }
    }
  }
  nm <- names(matrix)
  row.names(new_matrix) <- nm
  colnames(new_matrix) <- nm
  return(new_matrix)
}

orient <- function(adj){
  cp_adj <- adj
  dg <- degree_sort(adj)
  nm_adj <- names(cp_adj)
  nm_degree <- names(dg)
  for(i in 1:length(dg))
  {
    if(dg[[i]]!=0)
    {
      ni <- nm_degree[i]
      nb <- which((cp_adj[ni,]==1)|(cp_adj[ni,]==6))
      if(length(nb)!=0)
      {
        for(j in nb)
        {
          nb_cycle <- c()
          nj <- nm_adj[j]
          adj1 <- cp_adj
          adj2 <- cp_adj
          adj1[ni,nj] <- -2
          adj1[nj,ni] <- 2
          adjc <- convert_matrix(adj1)
          G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
          c <- find_cycles(G)
          nb_cycle <- c(nb_cycle,length(c))
          adj2[ni,nj] <- 2
          adj2[nj,ni] <- -2
          adjc <- convert_matrix(adj2)
          G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
          c <- find_cycles(G)
          nb_cycle <- c(nb_cycle,length(c))
          argmin <- which.min(nb_cycle)
          if(argmin==1){
            cp_adj[ni,nj] <- -2
            cp_adj[nj,ni] <- 2
          }
          else{
            cp_adj[ni,nj] <- 2
            cp_adj[nj,ni] <- -2
          }
        }
      }
    }
  }
  return(cp_adj)
}

remove_cycle <- function(matrix)
{
  adjc <- convert_matrix(matrix)
  G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
  cycle <- find_cycles(G)
  nm <- colnames(matrix)
  while(length(cycle)!=0)
  {
    l <- length(cycle[[1]])-1
    nb_cycle <- c()
    idx <- c()
    for (i in 1:l)
    {
      x <- cycle[[1]][i]
      y <- cycle[[1]][i+1]
      matrix[x,y] <- -2
      matrix[y,x] <- 2
      adjc <- convert_matrix(matrix)
      G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
      cycl <- find_cycles(G)
      nb_cycle <- c(nb_cycle,length(cycl))
      idx <- c(idx,i)
      matrix[x,y] <- 2
      matrix[y,x] <- -2
    }
    argmin <- which.min(nb_cycle)
    nb_cycle <- sort(nb_cycle,index.return=TRUE)$ix
    ix <- idx[argmin]
    x <- cycle[[1]][ix]
    y <- cycle[[1]][ix+1]
    matrix[x,y] <- -2
    matrix[y,x] <- 2
    adjc <- convert_matrix(matrix)
    G <- igraph::graph_from_adjacency_matrix(adjc, mode = "directed")
    cycle <- find_cycles(G)
  }
  return(matrix)
}

convert_dag <- function(adj){
  nadj <- orient(adj)
  adjc <- convert_matrix(nadj)
  G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
  c <- find_cycles(G)

  while(length(c)!=0)
  {
    nadj <- remove_cycle(nadj)
    adjc <- convert_matrix(nadj)
    G <- igraph::graph_from_adjacency_matrix(adjc,mode = "directed")
    c <- find_cycles(G)
  }
  adjc <- convert_matrix(nadj)
  return(adjc)
}


createEdgesList = function(all.edges.summary, dag, state_order){

  if(is.na(state_order)){
    all.edges.summary$sign = NA
    all.edges.summary$partial_correlation = NA
  }


  all.edges.summary$ort_inferred[which(all.edges.summary$ort_inferred==6)]=1
  for(i in 1:nrow(all.edges.summary)){
    x = all.edges.summary$x[i]
    y = all.edges.summary$y[i]

    if(dag[which(rownames(dag)==x), which(colnames(dag)==y)] == 1){
      all.edges.summary$ort_inferred[i] = 2
      all.edges.summary$proba[i] = '0;1'
    }
    if(dag[which(rownames(dag)==y), which(colnames(dag)==x)] == 1){
      all.edges.summary$ort_inferred[i] = -2
      all.edges.summary$proba[i] = '1;0'
    }
  }

  return(all.edges.summary)
}


#' @title Write miicsdg files to a new folder. The function will create the folder you specified. These two files can be loaded at the page: https://miic.curie.fr/vis_NL.php
#' @param miicsdg_res the output of the miicddg algorithm
#' @param folder A new folder where to store the data
#'
#' @export
writeToFile <- function(miicsdg_res, folder){
  if(file.exists(folder))
    stop('Directory already exists!')

  dir.create(folder)
  write.table(miicsdg_res$edges_miic_server, file.path(folder, "summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(miicsdg_res$data_types, file.path(folder,"state_order.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  # write.table(miicsdg_res$orientations.prob,file.path(folder,"or_prob.tsv"),sep="\t",quote=FALSE, row.names=FALSE)

}


# create synthetic data trough the MIIC algorithm
#' @title Create synthetic data trough the MIIC algorithm
#' @param original_data A data frame containing the original data we want to mimic
#' @param n_synthetic_samples An integer representing the number of synthetic samples to be created. Default value: equal to the number of samples in the input data
#' @param method_DAG A string ('DFS', 'CPDAG_to_DAG', 'MIIC_to_DAG') corresponding to the method used to transform the MIIC network into a DAG
#' @param root_index An integer corresponding to the position of the feature in the original data frame that should be used as root in the DFS visit if DFS is chosen as method to build the DAG
#' @param ids_columns_discrete A vector representing the index of features that should be used as non-numerical
#' @param keep_continuous_variable_boundaries A boolean, if TRUE, original boundaries are kept for continuous variables
#' @param confidence_threshold  A floating point number corresponding to the threshold used to filter the less probable edges following the skeleton step. See the MIIC algorithm for more details
#' @param consistent A Boolean. If TRUE the method iterates over skeleton and orientation steps to ensure consistency of the network. See the MIIC algorithm for more details
#' @param verbose A boolean indicating if the execution should be documented with more informations
#'
#' @return A list containig 5 objects: \itemize{
#'                                     \item{\emph{synthetic_data}:}{ a data frame containing the synthetic data}
#'                                     \item{\emph{adjacency_matrix_DAG}:}{ an adjacency matrix representing the dag used for creating the data (to read in the sense row -> column)}
#'                                     \item{\emph{data_types}:}{ a data frame containing the information on the natture of variables (0 = discrete, 1 = continuous) useful for the visualization of the DAG on the miic web server}
#'                                     \item{\emph{edges_miic_server}:}{ the edges file that can be used for the visualization of the DAG on the miic web server. You can use the function miicsdg::writeToFile to save this file plus the data_types file}
#'                                     \item{\emph{miic}:}{ the object returned from running miic}
#'                                     }
#' @export
#'
#' @examples
#' library(miicsdg)
#' library(datasets)
#' data(iris)
#'
#' miicsdg_list = miicsdg::miicsdg(iris)
#'
#' miicsdg_synthetic_data = miicsdg_list[['synthetic_data']]
#' miicsdg_adjacencyMatrixDAG = miicsdg_list[['adjacency_matrix_DAG']]
#' miicsdg_dataTypes = miicsdg_list[['data_types']]
#' miicsdg_edgesMiicServer = miicsdg_list[['edges_miic_server']]
#' miic_sdg_miic = miicsdg_list[['miic']]
#'
#' miicsdg::writeToFile(miicsdg_list , '~/testMIICSDG')
#'
#'
miicsdg <- function(original_data, n_synthetic_samples = NA, method_DAG = 'MIIC_to_DAG', root_index = 1,
                                    ids_columnsDiscrete = NA, state_order = NA,
                                    keep_continuous_variable_boundaries = T,
                                    confidence_threshold = 1,
                                    consistent = FALSE, verbose = FALSE) {


  isNumberVector<-function(v, factorsAccepted = TRUE){
    values = stats::na.omit(unique(v))
    if(is.factor(v)){
      if(!factorsAccepted)
        return(FALSE)
      if(suppressWarnings(length(which(is.na(as.numeric(stats::na.omit(as.character(v)))))) > 0))
        return(FALSE)
    }
    if(length(which(!is.na(suppressWarnings(as.numeric(values))) == FALSE))>0)
      return(FALSE)
    return(TRUE)
  }

  originalSummary = NA
  n_dags_max = 1
  set.seed(root_index)
  originalDataAsGiven = original_data

  # set number of samples
  if(is.na(n_synthetic_samples))
    n_synthetic_samples = nrow(original_data)

  # if the list of discrete variables is given transform them to factors (in this way a numerical variable with many levels could become a factor)
  if(!is.na(ids_columnsDiscrete[1])){
    if(length(ids_columnsDiscrete) > 0){
      for(c in ids_columnsDiscrete){
        original_data[,c] = factor(original_data[,c], levels = unique(original_data[,c]))
      }
    }
  }

  # if the type of variables is not already given with the correct type
  originalData_copy = original_data

  if(is.na(state_order)){
    for(c in 1:ncol(original_data)){
      if(class(original_data[,c])=='character' | length(unique(original_data[,c])) <= 8){
        original_data[,c] = factor(original_data[,c], levels = unique(original_data[,c]))
      }
    }
  } else {
    for(c in 1:ncol(original_data)){
      pos_in_state_order = which(state_order$var_names == colnames(original_data)[c])
      if(state_order$var_type[pos_in_state_order] == 0){
        if(is.na(state_order$levels_increasing_order[pos_in_state_order])){
          original_data[,c] = factor(original_data[,c], levels = unique(original_data[,c]))
        } else {
          original_data[,c] = factor(original_data[,c], levels = strsplit(state_order$levels_increasing_order[pos_in_state_order], ',')[[1]])
        }
      }
    }
  }


  # create NA class for variables with NA as value in factors
  originalData_NA_factor = original_data

  for(c in 1:ncol(originalData_NA_factor)){
    if(class(originalData_NA_factor[,c]) == 'factor'){
      if(length(which(is.na(originalData_NA_factor[,c]))) > 0){
        levels(originalData_NA_factor[,c]) <- c(levels(originalData_NA_factor[,c]),"NA")
        originalData_NA_factor[which(is.na(originalData_NA_factor[,c])),c] = 'NA'
      }
    }
  }


  # create data types for MIIC execution
  dataType = c()

  for(c in 1:ncol(originalData_NA_factor)){
    if(length(unique(originalData_NA_factor[,c])) > 8 & isNumberVector(originalData_NA_factor[,c]) & class(originalData_NA_factor[,c]) != 'factor'){
      dataType = c(dataType, 'Continuous')
    } else {
      dataType = c(dataType, 'Discrete')
    }
  }

  # if the summary is not given then run miic
  if(class(originalSummary)=='logical'){
    set.seed(0)

    if(!verbose)
      sink("/dev/null")
    miic_res = NA

    if(consistent == TRUE)
      consistent_str = 'skeleton'
    else
      consistent_str = 'no'

    if(is.na(state_order)){
      dataType_MIIC = data.frame(var_names = colnames(original_data), var_type = dataType)
      dataType_MIIC$var_type = as.numeric(factor(dataType_MIIC$var_type, levels = c('Discrete','Continuous'), labels = c(0,1)))-1
      dataType_MIIC$levels_increasing_order = NA
    } else {
      dataType_MIIC = state_order
    }

    # if there is no confidence cut
    if(confidence_threshold == 1)
      miic_res = miic::miic(original_data, verbose = F, consistent = consistent_str, propagation = TRUE, state_order = dataType_MIIC)
    else
      miic_res = miic::miic(original_data, verbose = F, consistent = consistent_str, n_shuffles = 100,
                      conf_threshold = confidence_threshold, state_order = dataType_MIIC)

    #miic_res$all.edges.summary$proba = NA
    #miic_res$all.edges.summary$is_causal = NA

    if(!verbose)
      sink()
    originalSummaryR = miic_res$summary
    originalSummary = miic_res$summary
  }

  # if there are missing values create a complete version of the dataset using mice
  if(length(which(is.na(originalData_NA_factor))) > 0){
    sink("/dev/null")
    miceIMP = mice::mice(originalData_NA_factor, m = 1)
    originalData_NA_factor_complete = mice::complete(miceIMP )
    apply(originalData_NA_factor_complete,2,function(x){length(which(is.na(x)))})
    sink()
  } else{
    originalData_NA_factor_complete = original_data
  }


  # generate a dag from an unoriented network
  generateDAG = function(originalSummary, root){
    temporarySummary = originalSummary
    undirected_graph = igraph::graph_from_data_frame(temporarySummary[,c("x","y")], directed = FALSE, vertices = NULL)
    dfs_res = igraph::dfs(undirected_graph, root = root)
    ordering = names(as.list(dfs_res$order))
    for(r in 1:nrow(temporarySummary)){
      # if node y is seen before node x, swap them
      if(which(ordering==temporarySummary$y[r]) < which(ordering==temporarySummary$x[r])) {
        x = temporarySummary$x[r]
        temporarySummary$x[r] = temporarySummary$y[r]
        temporarySummary$y[r] = x
      }
    }
    return(temporarySummary)
  }

  get_ajdMat_from_data_frame = function(summary, original_data){
    dag = data.frame(matrix(0, nrow = length(colnames(original_data)), ncol = length(colnames(original_data))))
    colnames(dag) = colnames(original_data)
    rownames(dag) = colnames(original_data)
    for(i in 1:nrow(summary)){
      pos_x = which(colnames(dag) == summary[i, 1])
      pos_y = which(colnames(dag) == summary[i, 2])
      dag[pos_x, pos_y] = 1
    }

    return(dag)
  }

  # generate continuous distribution that follows the original data for orphan nodes
  generateOrphanDistributionContinuous <- function(synthetic_data, indexVariable){
    v = originalData_NA_factor[,indexVariable]
    v = v[which(!is.na(v))]
    pdf_of_data <- stats::density(v, from= 0, to=1, bw=0.1)
    y <- stats::rnorm(n_synthetic_samples, sample(v, size = length(v), replace = TRUE), pdf_of_data$bw)

    while (length(which(is.na(y))) > 0) {
      y[which(is.na(y))] = stats::rnorm(length(which(is.na(y))), sample(v, size = length(v), replace = TRUE), pdf_of_data$bw)
    }
    synthetic_data[,indexVariable] = y

    #plot(density(synthetic_data[,indexVariable]))
    return(synthetic_data)
  }

  # generate discrete distribution that follows the original data for orphan nodes
  generateOrphanDistributionDiscrete <- function(synthetic_data, indexVariable){
    probability = prop.table(table(originalData_NA_factor[,indexVariable]))
    synthetic_data[,indexVariable] = sample(x = names(probability), prob = probability, size = n_synthetic_samples, replace = TRUE)
    if(class(originalData_NA_factor[,indexVariable]) == 'integer')
      synthetic_data[,indexVariable] = as.integer(synthetic_data[,indexVariable])
    return(synthetic_data)
  }


  # run logistic or multinomial regression
  logistic_Regression_or_multinomial <- function(synthetic_data, node, node_parents, pos, mytype){
    howManyContinuous = 0
    for(node_parent in node_parents){
      pos2 = which(colnames(originalData_NA_factor)==node_parent)
      if(dataType[pos2]=='Continuous')
        howManyContinuous = howManyContinuous + 1
    }

    syntheticData_copy = synthetic_data
    originalData_NA_factor_copy = originalData_NA_factor

    if(howManyContinuous <= 2){
      # if there are between 0 and 2 continuous variables
      # discretize nodes
      for(node_parent in node_parents){
        pos2 = which(colnames(originalData_NA_factor)==node_parent)
        if(dataType[pos2]=='Continuous'){
          res_discretization = miic::discretizeMutual(x = originalData_NA_factor[,pos], y = originalData_NA_factor[,pos2])
          cutpoints = res_discretization$cutpoints2
          originalData_NA_factor_copy[,pos2] = cut(x = originalData_NA_factor[,pos2], breaks = cutpoints)


          cutpoints[length(cutpoints)] = max(synthetic_data[,pos2], na.rm = T)
          cutpoints[1] = min(synthetic_data[,pos2], na.rm = T) - 0.0001
          syntheticData_copy[,pos2] = cut(x = synthetic_data[,pos2], breaks = cutpoints)
          originalData_NA_factor_copy[,pos2] = factor(originalData_NA_factor_copy[,pos2], levels = levels( originalData_NA_factor_copy[,pos2]), labels = levels( syntheticData_copy[,pos2]))
        }
      }



      # make combinations of all parents values
      syntheticData_parents = as.data.frame(syntheticData_copy[, node_parents])
      colnames(syntheticData_parents) = node_parents
      syntheticData_parents_unique = as.data.frame(syntheticData_parents[!duplicated(syntheticData_parents), ])
      colnames(syntheticData_parents_unique) = node_parents

      # set row numbers
      syntheticData_parents_unique$N = 1:nrow(syntheticData_parents_unique)
      syntheticData_parents$row = 1:nrow(syntheticData_parents)
      rownames(syntheticData_parents_unique) = 1:nrow(syntheticData_parents_unique)

      originalData_NA_factor_parents = as.data.frame(originalData_NA_factor_copy[, node_parents])
      colnames(originalData_NA_factor_parents) = node_parents
      originalData_NA_factor_parents_unique = as.data.frame(originalData_NA_factor_parents[!duplicated(originalData_NA_factor_parents), ])
      colnames(originalData_NA_factor_parents_unique) = node_parents

      originalData_NA_factor_parents_unique$N = 1:nrow(originalData_NA_factor_parents_unique)


      originalData_NA_factor_parents$row = 1:nrow(originalData_NA_factor_parents)
      rownames(originalData_NA_factor_parents_unique) = 1:nrow(originalData_NA_factor_parents_unique)

      model = NA

      for(i in 1:nrow(syntheticData_parents_unique)){
        syntheticData_parents_i = merge(syntheticData_parents, syntheticData_parents_unique[i,])$row
        syntheticData_parents_i = sort(syntheticData_parents_i)

        originalData_NA_factor_parents_i = merge(originalData_NA_factor_parents, syntheticData_parents_unique[i,])$row
        originalData_NA_factor_parents_i = sort(originalData_NA_factor_parents_i)

        if(length(originalData_NA_factor_parents_i) > 5){
          probabylities = prop.table(table(originalData_NA_factor[originalData_NA_factor_parents_i, pos]))
          if (length(originalData_NA_factor_parents_i) == 0) {
            names_proba = names(probabylities)
            probabylities = rep(1/length(probabylities), length(probabylities))
            names(probabylities) = names_proba
          }

          synthetic_data[syntheticData_parents_i,pos] <- sample(names(probabylities), prob = probabylities, size = length(syntheticData_parents_i), replace = TRUE)

        } else {
          # there is no true data with this combination in the synthetic data, we need to go with prediction
          if(is.na(model[1])){
            if(mytype=='logistic'){
              data = tidyr::drop_na(originalData_NA_factor_copy[,c(node, node_parents)])
              data = data[, apply(data, 2, function(x){length(unique(x))})  > 1]
              data[,node] = droplevels(data[,node])

              model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents[which(node_parents %in% colnames(data))], collapse = ' + '))) ,
                             data = data, method = "rf", metric = "Kappa", na.action = stats::na.pass)
            } else {
              originalData_NA_factor[,pos] = factor(originalData_NA_factor[,pos], levels = unique(originalData_NA_factor[,pos]))
              data = tidyr::drop_na(originalData_NA_factor_copy[,c(node, node_parents)])
              data = data[, apply(data, 2, function(x){length(unique(x))})  > 1]
              data[,node] = droplevels(data[,node])
              model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents[which(node_parents %in% colnames(data))], collapse = ' + '))) ,
                             data = data,
                             method = "rf", metric = "Kappa", na.action = stats::na.pass)
            }
          }

          x <- stats::predict(model, newdata = syntheticData_copy[syntheticData_parents_i,c(node, node_parents)])
          synthetic_data[syntheticData_parents_i,pos] = as.character(x)
        }
      }
    } else {
      model = NA
      if(mytype=='logistic'){
        data = tidyr::drop_na(originalData_NA_factor[,c(node, node_parents)])
        data = dplyr::select_if(data, function(col) length(unique(col))>1)
        data[,node] = droplevels(data[,node])

        model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents[which(node_parents %in% colnames(data))], collapse = ' + '))) , data,
                       method = "rf", metric = "Kappa",na.action = stats::na.pass)

      } else {
        data = tidyr::drop_na(originalData_NA_factor[,c(node, node_parents)])
        data = dplyr::select_if(data, function(col) length(unique(col))>1)
        data[,node] = droplevels(data[,node])
        model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents[which(node_parents %in% colnames(data))], collapse = ' + '))) , data,
                       method = "rf", metric = "Kappa", na.action = stats::na.pass)
      }

      x <- stats::predict(model, newdata = synthetic_data)
      synthetic_data[,pos] = x
    }

    levels(synthetic_data[,node]) <- levels(originalData_NA_factor[,node])

    return(synthetic_data)
  }

  # run regression for continuous variable
  regression <- function(synthetic_data, node, node_parents, pos){
    # check if at least one predictor is continuous
    howManyContinuous = 0
    for(node_parent in node_parents){
      pos2 = which(colnames(originalData_NA_factor)==node_parent)
      if(dataType[pos2]=='Continuous')
        howManyContinuous = howManyContinuous + 1
    }

    syntheticData_copy = synthetic_data
    originalData_NA_factor_copy = originalData_NA_factor
    model = NA

    if(howManyContinuous < 2){
      if(howManyContinuous > 0){
        for(node_parent in node_parents){
          pos2 = which(colnames(originalData_NA_factor)==node_parent)
          if(dataType[pos2]=='Continuous'){
            res_discretization = miic::discretizeMutual(x = originalData_NA_factor[,pos], y = originalData_NA_factor[,pos2])
            originalData_NA_factor_copy[,pos2] = cut(x = originalData_NA_factor[,pos2], breaks = res_discretization$cutpoints2)

            cutpoints = res_discretization$cutpoints2
            cutpoints[length(cutpoints)] = max(synthetic_data[,pos2], na.rm = T)
            cutpoints[1] = min(synthetic_data[,pos2], na.rm = T) - 0.0001
            syntheticData_copy[,pos2] = cut(x = synthetic_data[,pos2], breaks = cutpoints)
            originalData_NA_factor_copy[,pos2] = factor(originalData_NA_factor_copy[,pos2], levels = levels( originalData_NA_factor_copy[,pos2]), labels = levels( syntheticData_copy[,pos2]))
          }
        }
      }
      # learn PDF and apply to subcategories
      # create combinations of parents
      syntheticData_parents = as.data.frame(syntheticData_copy[, node_parents])
      colnames(syntheticData_parents) = node_parents
      syntheticData_parents_unique = as.data.frame(syntheticData_parents[!duplicated(syntheticData_parents), ])
      colnames(syntheticData_parents_unique) = node_parents

      syntheticData_parents_unique$N = 1:nrow(syntheticData_parents_unique)


      rownames(syntheticData_parents_unique) = 1:nrow(syntheticData_parents_unique)
      syntheticData_parents$row = 1:nrow(syntheticData_parents)

      originalData_NA_factor_parents = as.data.frame(originalData_NA_factor_copy[, node_parents])
      colnames(originalData_NA_factor_parents) = node_parents
      originalData_NA_factor_parents_unique = as.data.frame(originalData_NA_factor_parents[!duplicated(originalData_NA_factor_parents), ])
      colnames(originalData_NA_factor_parents_unique) = node_parents

      originalData_NA_factor_parents_unique$N = 1:nrow(originalData_NA_factor_parents_unique)


      originalData_NA_factor_parents$row = 1:nrow(originalData_NA_factor_parents)
      rownames(originalData_NA_factor_parents_unique) = 1:nrow(originalData_NA_factor_parents_unique)
      #      originalData_NA_factor_parents$row = 1:nrow(originalData_NA_factor_parents)

      for(i in 1:nrow(syntheticData_parents_unique)){
        syntheticData_parents_i = merge(syntheticData_parents, syntheticData_parents_unique[i,])$row
        syntheticData_parents_i = sort(syntheticData_parents_i)

        originalData_NA_factor_parents_i = merge(originalData_NA_factor_parents, syntheticData_parents_unique[i,])$row
        originalData_NA_factor_parents_i = sort(originalData_NA_factor_parents_i)

        # if there is no enough information with this particular configuration then use standard regression techniques
        if(length(originalData_NA_factor_parents_i) - length(which(is.na(originalData_NA_factor[originalData_NA_factor_parents_i, pos]))) == 0){

          if(is.na(model[1])){
            model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents, collapse = '+'))) , data = tidyr::drop_na(originalData_NA_factor_copy[,c(node, node_parents)]), method = "rf")
          }
          x <- stats::predict(model, newdata = dplyr::select(syntheticData_copy[syntheticData_parents_i,], node_parents))

          synthetic_data[syntheticData_parents_i,pos] = x
        } else {
          # if the original data have a unique value
          if(length(unique(originalData_NA_factor_copy[originalData_NA_factor_parents_i,c(node)])) == 1){
            synthetic_data[syntheticData_parents_i,pos] = originalData_NA_factor_copy[originalData_NA_factor_parents_i[1],c(node)]
          } else {
            pdf_of_data <- stats::density(stats::na.omit(originalData_NA_factor[originalData_NA_factor_parents_i, pos]), from=0, to=1, bw=0.1)
            v =  stats::rnorm(length(syntheticData_parents_i), sample(originalData_NA_factor_complete[originalData_NA_factor_parents_i,pos],
                                                               size = length(syntheticData_parents_i), replace = TRUE), pdf_of_data$bw)

            while (length(which(is.na(v))) > 0) {
              v[which(is.na(v))] = stats::rnorm(length(which(is.na(v))), sample(originalData_NA_factor_complete[originalData_NA_factor_parents_i,pos],
                                                                                        size = length(syntheticData_parents_i), replace = TRUE), pdf_of_data$bw)
            }
            v[which(is.nan(v))] = NA
            synthetic_data[syntheticData_parents_i,pos] <- v
          }
        }
      }

    } else {
      ## classical RF
      X = tidyr::drop_na(originalData_NA_factor[,c(node, node_parents)])
      # y = X[,1]
      # X = dplyr::select(X, dplyr::one_of(node_parents))
      # rf3 = randomForest::randomForest(X,y,maxnodes=2400,ntree=50000)

      # rf3 = randomForest::randomForest(X,y)

      PR = dplyr::select(synthetic_data, dplyr::one_of(node_parents))
      # for(v in 1:ncol(X)){
      #   if(class(X[,v]) == 'factor'){
      #     PR[,v] = factor(PR[,v], levels = levels(X[,v]))
      #   }
      # }


      model <- caret::train(stats::as.formula(paste(node , ' ~ ', paste(node_parents, collapse = '+'))) , data = X, method = "rf")
      x <- stats::predict(model, newdata = PR)


      # x = stats::predict(rf3, PR)
      synthetic_data[,pos] = x
    }

    return(synthetic_data)
  }


  ## MAIN

  globalScore=c()
  globalStructure=list()

  summary= originalSummary

  # keep only existing edges
  if(length(which(summary[, 'ort_inferred'] == 0)) > 0)
    summary = summary[which(summary[, 'ort_inferred'] != 0), ]

  # if we want a DAG starting from the CPDAG from miic
  if(method_DAG == 'CPDAG_to_DAG'){
    temporarySummary = summary

    for(val in which(temporarySummary[,"ort_inferred"]==-2)){
      temporarySummary[val,"ort_inferred"] = 2
      x = temporarySummary[val,"x"]
      temporarySummary[val,"x"] = summary[val,"y"]
      temporarySummary[val,"y"] = x
    }

    # get indirected arcs that we need to orient for DAG creation
    idxUndirectedArcs = which(temporarySummary[,"ort_inferred"] == 1 | temporarySummary[,"ort_inferred"] == 6)

    if(length(idxUndirectedArcs) == 0){
      temporarySummaryPDAG = temporarySummary
    } else {

      # build PDAG
      temporarySummaryDirected = temporarySummary[-idxUndirectedArcs,]

      temporarySummaryPDAGRight = temporarySummary[idxUndirectedArcs,]
      temporarySummaryPDAGRight[, 'ort_inferred'] = 2

      temporarySummaryPDAGLeft = temporarySummary[idxUndirectedArcs,]
      temporarySummaryPDAGLeft[, 'ort_inferred'] = -2

      temporarySummaryPDAG = rbind(temporarySummaryDirected, temporarySummaryPDAGRight, temporarySummaryPDAGLeft)

      temporarySummaryPDAG_copy = temporarySummaryPDAG
      for(val in which(temporarySummaryPDAG[,"ort_inferred"]==-2)){
        temporarySummaryPDAG[val,"ort_inferred"] = 2
        x = temporarySummaryPDAG[val,"x"]
        temporarySummaryPDAG[val,"x"] = temporarySummaryPDAG_copy[val,"y"]
        temporarySummaryPDAG[val,"y"] = x
      }
    }

    gR = methods::new("graphNEL", nodes = unique(c(temporarySummaryPDAG$x,temporarySummaryPDAG$y)), edgemode = "directed")
    for(row in 1:nrow(temporarySummaryPDAG)){
      gR = graph::addEdge(temporarySummaryPDAG$x[row], temporarySummaryPDAG$y[row], gR, 5987)
    }


    # get DAG from pDAG
    set.seed(root_index)
    gDAG <- pcalg::pdag2dag(gR)

    g_bnlearn = bnlearn::as.bn(gDAG$graph)
    g_igraph = igraph::graph_from_data_frame(g_bnlearn$arcs[,c('from','to')], directed = TRUE, vertices = NULL)
    dag = get_ajdMat_from_data_frame(g_bnlearn$arcs[,c('from','to')], original_data)
  }

  if(method_DAG == 'DFS'){
    # create a DAG from a starting root node
    root = unique(c(summary$x,summary$y))[root_index]
    if(verbose)
      print(paste('root:', root))
    temporarySummary = generateDAG(originalSummary = summary, root = root)
    g_igraph = igraph::graph_from_data_frame(temporarySummary[,c("x","y")], directed = TRUE, vertices = NULL)
    dag = get_ajdMat_from_data_frame(temporarySummary[,c("x","y")], original_data)

  }

  if(method_DAG == 'MIIC_to_DAG'){
    adj = as.data.frame(miic_res$adj_matrix)
    dag <- convert_dag(adj)

    colnames(dag) = colnames(adj)
    rownames(dag) = rownames(adj)


    g_igraph <- igraph::graph_from_adjacency_matrix(dag, mode = "directed")
    igraph::is.dag(g_igraph)
  }

  distances = list()

  n_dags_found = 0
  results = list()


  # if we have not already tested n_dags_max DAGs
  if(n_dags_found < n_dags_max){

    # if the graph is a DAG
    if(igraph::is.dag(g_igraph)){
      if(verbose)
        print(paste('IS DAG:', igraph::is.dag(g_igraph)))


      # generate synthetic_data data frame
      synthetic_data = data.frame(matrix(NA, nrow = n_synthetic_samples, ncol = ncol(originalData_NA_factor)))
      colnames(synthetic_data) = colnames(originalData_NA_factor)

      plot(g_igraph, edge.arrow.size=0.2)

      g_bnlearn = bnlearn::as.bn(g_igraph)

      # get isolated nodes
      isolated = colnames(originalData_NA_factor)[which(!colnames(originalData_NA_factor) %in% g_bnlearn$arcs[,2] & !colnames(originalData_NA_factor) %in% g_bnlearn$arcs[,1])]

      # get orphan Nodes
      orphans = names(g_bnlearn$nodes)[which(!names(g_bnlearn$nodes) %in% g_bnlearn$arcs[,2])] # Nodes without parents.

      # generate data for orphan nodes
      for(node in c(isolated, orphans)){
        pos = which(colnames(originalData_NA_factor)==node)
        if(dataType[pos] == 'Continuous'){
          synthetic_data= generateOrphanDistributionContinuous(synthetic_data, pos)
        } else {
          if(dataType[pos] == 'Discrete'){
            synthetic_data = generateOrphanDistributionDiscrete(synthetic_data, pos)
          }
        }
        if(class(originalData_NA_factor[,pos]) == 'numeric')
          synthetic_data[,pos] = as.numeric(synthetic_data[,pos])
        if(class(originalData_NA_factor[,pos]) == 'factor')
          synthetic_data[,pos] = factor(synthetic_data[,pos], levels = levels(originalData_NA_factor[,pos]))
      }


      nodes = unique(c(names(g_bnlearn$nodes), isolated))
      simulation_finished = F

      # get nodes that have already evaluated parents
      list_of_parent_nodes = unique(c(isolated, orphans))

      while (!simulation_finished) {
        # loop on nodes
        for(node in nodes[!nodes %in% list_of_parent_nodes]){
          node_parents = g_bnlearn$arcs[,1][g_bnlearn$arcs[,2] == node]

          if(all(node_parents %in% list_of_parent_nodes)){
            if(verbose) print(paste0('Nodes already processed: ', round((length(list_of_parent_nodes) / length(nodes) * 100),2), '%'))
            if(verbose) print(paste0('Processing: ', node, '. It has ', length(node_parents), ' parents: ', paste(node_parents, collapse = ',')))

            # simulate data using parents
            pos = which(colnames(originalData_NA_factor)==node)

            # logistic regression
            if(length(unique(originalData_NA_factor[,pos])) == 2){
              synthetic_data = logistic_Regression_or_multinomial(synthetic_data, node, node_parents, pos, mytype = 'logistic')
            } else {

              # multinomial regression
              if(length(unique(originalData_NA_factor[,pos])) < 8 | !isNumberVector(originalData_NA_factor[,pos]) | dataType[pos] == 'Discrete'){
                synthetic_data = logistic_Regression_or_multinomial(synthetic_data = synthetic_data, node = node, node_parents = node_parents, pos = pos, mytype = 'multinomial')
              } else {
                if(dataType[pos] == 'Continuous'){
                  synthetic_data = regression(synthetic_data, node, node_parents, pos)
                }
              }
            }

            list_of_parent_nodes = c(list_of_parent_nodes, node)
          }

          if(class(originalData_NA_factor[,pos]) == 'numeric'){
            synthetic_data[,pos] = as.numeric(synthetic_data[,pos])

            # if original data is integer, keep integer
            if(all(na.omit(originalData_NA_factor[,pos]) == floor(na.omit(originalData_NA_factor[,pos]))) )
              synthetic_data[,pos] = round(synthetic_data[,pos])

            # set minimum and mazimum wrt original data
            if(keep_continuous_variable_boundaries){
              min = min(original_data[,pos],na.rm = T)
              max = max(original_data[,pos],na.rm = T)

              synthetic_data[which(synthetic_data[,pos] < min),pos] = min
              synthetic_data[which(synthetic_data[,pos] > max),pos] = max
            }
          }
          if(class(originalData_NA_factor[,pos]) == 'factor'){
            synthetic_data[,pos] = factor(synthetic_data[,pos], levels = levels(originalData_NA_factor[,pos]))
          }
        }
        simulation_finished = all(nodes %in% list_of_parent_nodes)
      }

      ## remove duplicated rows and keep the first nrow(originalData_NA_factor)
      synthetic_data = synthetic_data[!duplicated(synthetic_data), ]

      # put NAs as real NAs
      synthetic_data[synthetic_data=='NA'] = NA
      for(c in 1:ncol(synthetic_data)){
        if(class(synthetic_data[,c])=='factor'){
          levs = levels(synthetic_data[,c])
          if('NA' %in% levs){
            synthetic_data[,c] <- droplevels(synthetic_data[,c])
          }
        }
      }

      # this is a possible result, store it
      n_dags_found = n_dags_found + 1
      results[[n_dags_found]] = list()
      results[[n_dags_found]][['synthetic_data']] = synthetic_data
      results[[n_dags_found]][['dag']] = g_bnlearn
    } else {
      Cycles = NULL
      for(v1 in igraph::V(g_igraph)) {
        for(v2 in igraph::neighbors(g_igraph, v1, mode="out")) {
          Cycles = c(Cycles,
                     lapply(igraph::all_simple_paths(g_igraph, v2, v1, mode="out"), function(p) c(v1,p)))
        }
      }
    }
  }

  dataType_MIIC_server = dataType_MIIC
  dataType_MIIC$var_type = as.character(dataType_MIIC$var_type)
  dataType_MIIC$var_type[which(dataType_MIIC$var_type==0)] = 'Discrete'
  dataType_MIIC$var_type[which(dataType_MIIC$var_type==1)] = 'Continuous'

  res = list()
  res[['synthetic_data']] = results[[1]]$synthetic_data
  res[['adjacency_matrix_DAG']] = dag
  res[['data_types']] = dataType_MIIC_server
  res[['edges_miic_server']] = createEdgesList(all.edges.summary = miic_res$summary, dag = dag, state_order = state_order)
  res[['miic']] = miic_res


  return(res)
}
