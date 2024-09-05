

#' pav_pca
#'
#' `pav_pca()` will perform PCA analysis of PAV data using `prcomp()`. The `center`, `scale`, and `rank` will pass to `prcomp()`.
#'
#' @param pav_obj A PAV object.
#' @param center A logical value indicating whether the variables should be shifted to be
#' zero centered, pass to \code{\link[stats]{prcomp}}.
#' @param scale A logical value indicating whether the variables should be scaled to have
#' unit variance before the analysis takes place, pass to \code{\link[stats]{prcomp}}.
#' @param rank A number specifying the maximal rank, pass to \code{\link[stats]{prcomp}}.
#'
#' @param add_pheno_info A character string of `pheno_info` names.
#' @param pheno_info_color_list A list contains named vector of colors for `pheno_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#'
#' @param axis_text_size The size of tick labels on axis.
#' @param axis_title_size The size of axis title.
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_text_size The size of legend item labels.
#' @param legend_title_size The size of legend title.
#'
#' @importFrom ggplot2 ggplot geom_point labs aes xlab ylab theme_bw theme scale_color_manual
#'
#' @export

pav_pca <- function(pav_obj,
                    center = T,
                    scale = F,
                    rank = NULL,

                    add_pheno_info = NULL,
                    pheno_info_color_list = NULL,

                    axis_text_size = NULL,
                    axis_title_size = NULL,

                    legend_side = "right",
                    legend_text_size = NULL,
                    legend_title_size = NULL){

  check_class(pav_obj, "PAV")

  pav_data <- pav_obj@pav_data
  sample <- pav_obj@sample

  PC1 = PC2 = pheno = NULL

  pca_res <- stats::prcomp(t(pav_data), center = center, scale = scale, rank = rank)

  pca.var.per <- round(pca_res$sdev^2/sum(pca_res$sdev^2)*100, 2)

  p_data <- data.frame(sample = rownames(pca_res$x), PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])

  if(length(sample$pheno) == 0 || is.null(add_pheno_info)){
    p <- ggplot(p_data, aes(x = PC1, y = PC2)) + geom_point()
    pheno_col <- NULL
  } else {
    add_pheno_info <- match.arg(add_pheno_info, names(sample$pheno))
    pheno_col <- get_anno_palette(pheno_info_color_list, sample$pheno[add_pheno_info], num_re_func = F)
    pheno_data <- data.frame(sample = sample$name,
                            pheno = sample$pheno[[add_pheno_info]])
    p_data <- merge(p_data, pheno_data, by = "sample")

    p <- ggplot(p_data, aes(x = PC1, y = PC2, color = pheno)) +
      geom_point() + labs(color = add_pheno_info)
  }

  p <- p +
    xlab(paste("PC1(",pca.var.per[1],"%",")",sep=""))+
    ylab(paste("PC2(",pca.var.per[2],"%",")",sep=""))+
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size, color = "black", face = "bold"),
          axis.title = element_text(size = axis_title_size, color = "black", face = "bold"),
          legend.position = legend_side,
          legend.text = element_text(size = legend_text_size, color = "black", face = "bold"),
          legend.title = element_text(size = legend_title_size, color = "black", face = "bold"))

  if(!is.null(pheno_col)){
    if(is.numeric(sample$pheno[[add_pheno_info]])){
      p <- p +
        scale_color_gradientn(colours = pheno_col[[add_pheno_info]])
    }else{
      p <- p +
        scale_color_manual(values = pheno_col[[add_pheno_info]])
    }
  }

  print(p)
}








