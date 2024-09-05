


setClass("COV",
         slots = list(
           cov_data = "matrix",
           region = "list",
           sample = "list"
         )
)



#' get_cov_obj
#'
#' Get an object of the COV class.
#'
#' @param cov_data A numeric `matrix` or `data.frame` of the coverage data.
#' The row names are region names. The column names are sample names.
#' @param pheno_info A `data.frame` of phenotype data. The row names are sample name.
#' @param region_info A `data.frame` of region data. The row names are region name.
#'
#' The `cov_data` shouldn't contain region that is 0 in all sample, and the same for sample.
#'
#'
#' @export

get_cov_obj <- function(cov_data,
                        pheno_info = NULL,
                        region_info = NULL){

  cov_data <- check_pav_cov_data("cov_data", cov_data)

  if(!is.null(region_info)) region_info <- check_info(region_info, rownames(cov_data), "region", "region")
  if(!is.null(pheno_info)) pheno_info <- check_info(pheno_info, colnames(cov_data), "pheno", "sample")

  obj <- new("COV")

  obj@cov_data <- as.matrix(cov_data)
  obj@region <- as.list(region_info)
  obj@sample <- as.list(pheno_info)

  obj
}



