

setClass("PAV",
          slots = list(
           pav_data = "matrix",
           args = "list",
           region = "list",
           sample = "list"
         )
)


#' get_pav_obj
#'
#' Get an object of the PAV class.
#'
#' @param pav_data  A numeric `matrix` or `data.frame` of PAV data. `0` represent "absence" and `1` represent "presence".
#' The row names are region names. The column names are sample names.
#' @param pheno_info The phenotype data. The row names are sample names.
#' @param region_info Relevant information on target regions. The row names are region names.
#' @param add_softcore A logical value indicating whether to consider `soft-core` when determining the category of the region.
#' @param add_private A logical value indicating whether to consider `private` when determining category of the region.
#' @param softcore_loss_rate A numeric vector of loss rate.
#' @param use_binomial A logical value indicating whether to use binomial test.
#' @param softcore_p_value A numeric vector of p-value (binomial test). It only works when `use_binomial`=T
#'
#' `add_softcore` and `add_private` are logical values indicating whether to consider "softcore" and "private" when determining the category.
#'
#' If `add_softcore` is `TRUE`, the regions with loss rates not significantly larger than `softcore_loss_rate` will be considered as "softcore". Binomial tests (with a null hypothesis of loss rate < `softcore_loss_rate`) are carried out for each region if `use_binomial` is `T`. A p-vlaue below `softcore_p_value` means that this region is lost in a significant proportion and is a distributed region(loss rate > `softcore_loss_rate`).
#'
#' If `add_private` is `TRUE`, the regions present in only one sample will be considered as "private" region.
#'
#'
#' @importFrom methods new
#'
#' @export

get_pav_obj <- function(pav_data,
                        pheno_info = NULL,
                        region_info = NULL,
                        add_softcore = T,
                        add_private = T,
                        softcore_loss_rate = .1,
                        use_binomial = F,
                        softcore_p_value = .05){

  pav_data <- check_pav_cov_data("pav_data", pav_data)
  pav_data <- pav_data[order(rowSums(pav_data), decreasing = T), order(colSums(pav_data), decreasing = T)]

  if(!is.null(region_info)) region_info <- check_info(region_info, rownames(pav_data), "region", "region")
  if(!is.null(pheno_info)) pheno_info <- check_info(pheno_info, colnames(pav_data), "pheno", "sample")

  add_softcore <- match_logi("add_softcore", add_softcore)
  add_private <- match_logi("add_private", add_private)
  use_binomial <- match_logi("use_binomial", use_binomial)
  softcore_loss_rate <- match_num("softcore_loss_rate", softcore_loss_rate, 0, 1)
  softcore_p_value <- match_num("softcore_p_value", softcore_p_value, 0, 1)


  type_info_res <- get_type_info(ncol(pav_data), add_private, add_softcore, softcore_loss_rate, use_binomial, softcore_p_value)
  type_info <- type_info_res$type_info
  dist_sample_n <- type_info_res$dist_sample_n

  type <- get_region_type(pav_data, type_info)
  n <- as.integer(as.vector(rowSums(pav_data)))

  obj <- methods::new("PAV")

  obj@pav_data <- as.matrix(pav_data, dimnames = NULL)

  obj@args <- list(
    add_softcore = add_softcore,
    add_private = add_private,
    softcore_loss_rate = softcore_loss_rate,
    use_binomial = use_binomial,
    softcore_p_value = softcore_p_value,
    dist_n = dist_sample_n
  )

  obj@region <- list(
    name = rownames(pav_data),
    sample_n = n,
    type = type,
    info = as.list(region_info)
  )

  obj@sample <- list(
    name = colnames(pav_data),
    pheno = as.list(pheno_info)
  )

  return(obj)

}






