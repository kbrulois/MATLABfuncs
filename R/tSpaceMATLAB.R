#' Run the matlab version of the tSpace algorithm
#' 
#' This is an interface to the matlab implementation of the tSpace algorithm described in 
#' Dermadi et al, modified for improved error handling. 
#' 
#' @param data A matrix-like object where cells correspond to rows and features correspond to columns
#' @param trajectories The number of trajectories to calculate.
#' @param k Number of nearest neighbors (passed to the k argument of wanderlust). If a vector of integers is provided (recommended), the tSpace calculation will be attempted on the smallest value in the vector and then reattempted on successively higher values of K should the algorithm fail. tSpace caclualtions that require an excessively high k value suggest the presence of outlier data points that should be removed from the input data for better performance.
#' @param l Number of nearest neighbors (passed to the l argumet of wanderlust). If set to 'auto' (default), l will be 2/3 of k (rounded to the nearest integer). If k is a vector, l must be a vector the same length and all corresponding values must be less than k.
#' @param metric A distance metric accepted by the dist function in matlab (wanderlust): euclidean, correlation, cityblock, chebyshev, minkowski, hamming, mahalanobis, jaccard, cosine or spearman.
#' @param graphs Number of graphs over which to average trajectories (wanderlust).
#' @param landmarks Number of landmarks (wanderlust). 
#' @param voting_scheme Arugment passed to wanderlust. 
#' @param label A character string to be incorporated into the tSpace column labels. Default is nothing.
#' @param matlab_version A character string specifying the active matlab version (e.g. "R2019a"). If set to "auto," the latest version of matlab will be used. Only works for mac users.
#' @param matlab_path The path to the MATLAB folder containing cyt and tSpace scripts.
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{tPCs}{A cell by tPC \code{matrix} of the PCA reduction of the trajectory data.}
#'   \item{traj}{A cell by trajectory \code{matrix} containing the trajectories.}
#'   \item{clusters}{A vector \code{vector} containing the cluster definitions used to initiate tSpace.}
#'   \item{pPCA}{A cell by pPCA \code{matrix} of the PCA reduction of the parameters.}
#' }
#'
#' @author Kevin Brulois
#' @export

tSpaceMATLAB = function(data, 
                     trajectories = 100,
                     k = c(20, 30, 50, 100),
                     l = "auto",
                     metric = "euclidean",
                     graphs = 5, 
                     landmarks = 20,
                     voting_scheme = "exponential",
                     label = "",
                     matlab_version = "auto",
                     matlab_path = "~/Documents/MATLAB") {
  
  on.exit({print(paste("removing temporary files"))
    try({file.remove(to.remove)
      Sys.unsetenv(c("tsp_distMetric",
                     "tsp_voting_scheme",
                     "path2tSpaceInput",
                     "path2tSpaceOutput"))})
    print(paste("done") )})
  
  tsp.dist.lut <- data.frame(ml = c("euclidean", "correlation", "cityblock", "chebyshev", "minkowski", "hamming", "mahalanobis", "jaccard", "cosine", "spearman"),
                             ml.abrv = c("eu", "cor", "manh", "cheb", "mink", "hamm", "maha", "jac", "cos", "sp"))
  
  getMatlabVersion <- function() {
    apps <- list.files("/Applications")
    matlab <- apps[grepl("MATLAB", apps)]
    matlab_ver <- sub("MATLAB_", "", matlab[gtools::mixedorder(matlab)][length(matlab)])
    sub(".app", "", matlab_ver)
  }
  
  if(matlab_version == "auto") {
    matlab_version <- getMatlabVersion()
  }
  matlab_command <- paste0("/Applications/MATLAB_", matlab_version, ".app/bin/matlab")
  
  path2tSpaceInput <- tempfile("tsp_in", fileext = ".csv")
  to.remove <- path2tSpaceInput
  
  path2tSpaceOutput <- tempfile("tsp_out", fileext = ".csv")
  to.remove <- c(to.remove, path2tSpaceOutput)
  
  Sys.setenv(tsp_distMetric = metric,
             tsp_voting_scheme = voting_scheme,
             path2tSpaceInput = path2tSpaceInput,
             path2tSpaceOutput = path2tSpaceOutput)
  
  data.table::fwrite(as.data.frame(data), file = path2tSpaceInput, row.names = FALSE)
  
  mat1 = "-nodisplay"
  mat2 = "-nosplash"
  mat3 = "-nodesktop"
  
  is.error <- TRUE
  iter <- 1
  k <- k[order(k)]
  
  pkg_ml_path <- system.file("exec", package="tSpaceMATLAB")
  
  while(is.error & iter <= length(k)) {
    is.error <- tryCatch({
      message("Trying tspace with k = ", k[iter])
      if(l == "auto") {l_param <- (2*k[iter]) %/% 3} else {l_param <- l[iter]}
      tSpaceParms = paste(k[iter], l_param, graphs, landmarks, 1, trajectories, 30, sep = ",")
      mat4 = paste0("-r", " ", "\'addpath(genpath(\"", matlab_path, "\"), genpath(\"", pkg_ml_path,"\"));", " ","try; tspace_ml(", tSpaceParms, ");", " ", "catch;" , " ", "end;", " ", "quit\'")
      args = c(mat1, mat2, mat3, mat4)
      output = system2(matlab_command, args= args, stdout=TRUE)
      #print(output)
      !file.exists(path2tSpaceOutput)
      
    }, error=function(e) TRUE)
    
    iter <- iter + 1
    
  }
  
  flug1 <- data.table::fread(path2tSpaceOutput, header = FALSE)
  
  flug1 <- as.matrix(flug1)
  rownames(flug1) <- rownames(data)
  
  tsp.p <- paste0(tsp.dist.lut[tsp.dist.lut$ml == metric, "ml.abrv"],
                  "T", trajectories,
                  "K", k[(iter-1)], 
                  "G", graphs,
                  "L", landmarks)
  
  L <- list(tPCs = flug1[,23:42], 
            traj = flug1[,43:ncol(flug1)], 
            clusters = flug1[,22], 
            pPCA = flug1[,2:21])
  
  if(!file.exists("~/.bcs.rds") | length(readRDS("~/.bcs.rds")) < 20) {
    bcs <- c(unique(replicate(10000, tolower(paste0(sample(LETTERS, 3, replace = FALSE), collapse = "")))),
             unique(replicate(10000, tolower(paste0(sample(LETTERS, 4, replace = FALSE), collapse = "")))))
    saveRDS(bcs, "~/.bcs.rds")
    
  }
  
  bcs <- readRDS("~/.bcs.rds")
  
  colnames(L$tPCs) <- c(paste0("tPC", label, "|", tsp.p, "|", bcs[1]), paste0(bcs[1], 2:20))
  colnames(L$traj) <- c(paste0("traj", label, "|", tsp.p, "|", bcs[3]), paste0(bcs[3], 2:trajectories))
  colnames(L$pPCA) <- c(paste0("pPC", label, "|", tsp.p, "|", bcs[4]), paste0(bcs[4], 2:20))
  
  saveRDS(bcs[-c(1:4)], "~/.bcs.rds")
  
  return(L)
  
}



