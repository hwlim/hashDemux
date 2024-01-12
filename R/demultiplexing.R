#' clustering-based sample demux
#'
#'
#' @param seurat_object Seurat object
#' @param assay Assay name in seurat_object with hash tag counts matrix (default: HTO)
#' @param expected_doublet_rate proportion of droplets expected to be doublets.
#' If value is NULL (default), it will be automatically estimated
#'
#' @return A Seurat object (obj) with three columns added to obj@meta.data. \code{sampleBC} column has
#'  tag assignment for each cell with doublets represented as "tag1_tag2". \code{classification} column has global
#'  classification information i.e. cells are classified as "Doublet", "Singlet"or "Negative".
#'  \code{confidence_score}a per-cell value between 0 and 1 with confidence score for the
#'  demultiplexing result
#'
#' @examples
#' # library(Seurat)
#' # library(dplyr)
#' # seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR",margin = 2)
#' # seurat_object = clustering_based_demux( seurat_object, assay = "HTO")
#'
#' @import Seurat
#' @import dplyr
#' @import doParallel
#' @importFrom foreach %dopar% %:% foreach
#' @export

clustering_based_demux = function(seurat_object, assay = "HTO", expected_doublet_rate = NULL,nCores=NULL){
  knns = seq(5,30,5)
  resolutions = c(1,2,3,4)
  if (is.null(nCores)) {
    doParallel::registerDoParallel(parallel::detectCores() -2)
  }else{
    doParallel::registerDoParallel(nCores)
  }
  results = foreach(k = knns) %:%
    foreach(resol = resolutions) %dopar% {
      res <- findMarkerTags(seurat_object, assay = assay,resol = resol,knn = k,logfc.threshold = 0.05)
      seurat_object$seurat_clusters = res[[2]]
      threshold = find_logfc(seurat_object,markers = res[[1]],logfc.thresholds = seq(0.05,0.5,0.05), expected_doublet_rate = expected_doublet_rate )

      filtered_markers = res[[1]] %>% dplyr::filter( avg_log2FC >= threshold)
      labels = label_clusters(seurat_object, filtered_markers )
      labels$sampleBC
    }

  names(results) = knns

  for (x in seq(length(results))) {
    names(results[[x]]) = resolutions
  }
  df = do.call(cbind, lapply(results, data.frame))

  predictions = apply(df, 1, function(x){
    tab = table(x)
    label = names(tab)[tab == max(tab)][1]
  })  %>% unlist()

  scores = apply(df, 1, function(x){
    tab = table(x)
    score = max(tab)/length(x)
  })  %>% unlist()

  seurat_object$sampleBC = predictions
  seurat_object$confidence_score = scores

  class= if_else(seurat_object$sampleBC == "Negative" , "Negative",
                 if_else(grepl("_",seurat_object$sampleBC) ,"Doublet" ,"Singlet"))
  class = factor(class)
  class = factor(class, levels = c( "Doublet","Negative", "Singlet" ))

  seurat_object$classification = class

  return(seurat_object)
}


# find marker tags per cluster
findMarkerTags2 <- function(seurat_object,assay = "HTO", resol = 1,
                                   logfc.threshold = 0.1 ,knn = 20 )
{
  DefaultAssay(seurat_object) = assay
  # build nearest-neighbor graph
  seurat_object <- Seurat::FindNeighbors(seurat_object, features = rownames(seurat_object), dims = NULL,
                                         assay = assay,k.param = knn )
  # cluster cells
  seurat_object <- Seurat::FindClusters(seurat_object,
                                        graph.name = paste0(assay,"_snn") ,
                                        resolution = resol )
  # find cluster markers
  markers <- Seurat::FindAllMarkers(seurat_object, only.pos = TRUE,pseudocount.use = 0.1, # compatible with Seurat5
                            assay = assay, verbose = F,
                            logfc.threshold = logfc.threshold  )
  return(list(markers, seurat_object$seurat_clusters))
}

findMarkerTags <- function(seurat_object,assay = "HTO", resol = 1,
                           logfc.threshold = 0.1 ,knn = 20 )
{
  DefaultAssay(seurat_object) = assay
  # build nearest-neighbor graph
  mtrx = seurat_object[[assay]]$data %>% t()
  #dist_mtrx = dist(mtrx %>% as.matrix() )
  nn = FindNeighbors(object = mtrx, k.param = knn)
  snn = nn$snn
  snn@assay.used = assay
  seurat_object@graphs[[paste0(assay,"_snn")]] = snn

  # cluster cells
  seurat_object <- Seurat::FindClusters(seurat_object,
                                        graph.name = paste0(assay,"_snn") ,
                                        resolution = resol )
  # find cluster markers
  markers <- Seurat::FindAllMarkers(seurat_object, only.pos = TRUE,pseudocount.use = 0.1, # compatible with Seurat5
                                    assay = assay, verbose = F,
                                    logfc.threshold = logfc.threshold  )
  return(list(markers, seurat_object$seurat_clusters))
}

## Assign a sample label to each cluster based on marker tags
##
## output: a data frame with two columns,
##        1. sampleBC: assign a tag to each Singlet, tag1_tag2 in case of Doublets and Negative otherwise
##        2. class: classify cells to either Singlet, Doublet or Negative
label_clusters <- function(seurat_object,markers)
{
  # marker expression plots
  if(nrow(markers) == 0){
    df = data.frame(cell_barcode = colnames(seurat_object), sampleBC = "Negative", class = "Negative")
    rownames(df) <- df$cell_barcode
    return(df)
  }

  top_markers <- markers %>% dplyr::filter( pct.1 == 1,  p_val_adj < 0.05)  %>%
    dplyr::group_by( cluster) %>%
    dplyr::top_n(2, wt =  avg_log2FC)

  ## assign cells to sample of origin
  cluster_barcode_dict <- top_markers %>% dplyr::select(c("cluster", "gene"))

  cell_cluster_dict <- seurat_object[["seurat_clusters"]]
  colnames(cell_cluster_dict) <- "cluster"
  cell_cluster_dict$cell_barcode <- rownames(cell_cluster_dict)

  cell_barcode_dict <- dplyr::left_join(cell_cluster_dict , cluster_barcode_dict , by = c("cluster") )

  cell_classification <- cell_barcode_dict %>% dplyr::group_by( cell_barcode) %>%
    dplyr::summarise(sampleBC = paste0(sort( gene), collapse = "_"))

  cell_classification$sampleBC[cell_classification$sampleBC == ""] = "Negative"

  cell_classification.global <- cell_classification %>%
    dplyr::mutate(class= if_else( sampleBC == "Negative" , "Negative",
                                 if_else(grepl("_", sampleBC) ,"Doublet" ,"Singlet")))

  cell_classification.global$class = factor(cell_classification.global$class,
                                            levels = c( "Doublet","Negative", "Singlet" ))
  rownames(cell_classification.global) <- cell_classification.global$cell_barcode

  cell_classification.global= cell_classification.global[colnames(seurat_object),]
  return( cell_classification.global)

}


# find logFC value that generate doublet rate close to the expected rate
find_logfc= function(seurat_object, markers,logfc.thresholds = seq(0.05,0.5,0.05), expected_doublet_rate = NULL){
  results = list()
  for(logfc.threshold in logfc.thresholds ) {
    filtered_markers = markers %>% dplyr::filter( avg_log2FC >= logfc.threshold)
    out <- label_clusters(seurat_object, filtered_markers)
    results[[as.character(logfc.threshold)]] = prop.table(table(out$class))
  }
  dfs <- lapply(results, data.frame)
  df = dplyr::bind_rows(dfs)
  colnames(df) = c("classification" , "proportion")
  df$logfc.threshold = rep(logfc.thresholds, each = 3)

  if(is.null(expected_doublet_rate)){
    expected_doublet_rate = 8e-6 * ncol(seurat_object) * (nrow(seurat_object)-1)  / nrow(seurat_object)
  }

  df2 = df[df$classification == "Doublet",]
  index = which.min(abs(df2$proportion - expected_doublet_rate))
  threshold = df2$logfc.threshold[index]

  return(threshold)
}

