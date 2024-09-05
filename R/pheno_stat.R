

#' pheno_stat
#'
#' Perform phenotype association calculations
#'
#' @param pav_obj A PAV object.
#' @param pheno_name The name of phenotype used for calculation.
#' @param p_adjust_method The adjustment methods, pass to  \code{\link[stats]{p.adjust}}
#' @param pallel A logical value indicating whether to use parallel computing.
#' @param parallel_n The number of CPU cores used for parallel computing.
#'
#'
#' @export

pheno_stat <- function(pav_obj,
                      pheno_name,
                      p_adjust_method = "fdr",
                      pallel = F,
                      parallel_n = parallel::detectCores() - 1){

  check_class(pav_obj, "PAV")

  if(length(pav_obj@sample$pheno) == 0) {
    stop(pav_obj, " doesn't have `pheno_info`.")
  }
  pheno_data <- data.frame(pav_obj@sample$pheno)
  pheno_name <- match.arg(pheno_name, colnames(pheno_data), several.ok = T)
  pheno_data <- pheno_data[, pheno_name, drop = F]

  type = NULL

  region_data <- data.frame(region = pav_obj@region$name, type = pav_obj@region$type)
  pav_data <- pav_obj@pav_data[subset(region_data, type != "core")$region,]

  regions <- rownames(pav_data)
  phenos <- colnames(pheno_data)

  pheno_stat_p <- function(cur_pheno_data, cur_region_data){
    if(is.numeric(cur_pheno_data)){
      p_value <- tryCatch({
        presence_data <- cur_pheno_data[cur_region_data == 1]
        absence_data <- cur_pheno_data[cur_region_data == 0]
        wilcox_res <- stats::wilcox.test(presence_data, absence_data, exact = F)
        wilcox_res$p.value
      }, error = function(e) {
        NA
      })
    } else {
      p_value <- tryCatch({
        fisher_res <- stats::fisher.test(cur_pheno_data, cur_region_data, simulate.p.value = F)
        fisher_res$p.value
      }, error = function(e) {
        NA
      })
    }
    p_value
  }

  pheno_res <- do.call(
    rbind,
    lapply(phenos, function(pheno){
      if(pallel){
        snowfall::sfInit(parallel = TRUE, cpus = parallel_n )
        snowfall::sfExport("pheno_stat_p")
        snowfall::sfExport("regions")
        snowfall::sfExport("pheno_data")
        snowfall::sfExport("pheno")
        snowfall::sfExport("pav_data")
        res <- snowfall::sfLapply(regions, function(region){
          p_value <- pheno_stat_p(pheno_data[, pheno], pav_data[region, ])
          c(pheno, region, p_value)
        })
        snowfall::sfStop()
        res
      } else {
        res <- lapply(regions, function(region){
          p_value <- pheno_stat_p(pheno_data[, pheno], pav_data[region, ])
          c(pheno, region, p_value)
        })
      }
      data <- data.frame(t(as.data.frame(res)))
      rownames(data) <- NULL
      colnames(data) <- c("pheno", "region", "p_value")
      data$p_value <- as.numeric(data$p_value)
      data$p_adjusted <- stats::p.adjust(data$p_value, method = p_adjust_method)
      data
    })
  )
}






