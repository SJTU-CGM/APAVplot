

#' cov_heatmap
#'
#'
#'
#' @param cov_obj A COV object.
#' @param add_pheno_info A character vector of `pheno_info` names.
#' @param add_region_info A character vector of `region_info` names.
#'
#' @param cov_colors A vector of colors or a color mapping function, pass to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param region_info_color_list A list contains named vector of colors for `region_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
#' @param pheno_info_color_list A list contains named vector of colors for `pheno_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"), age = c("yellow", "red"))
#'
#' @param border A logical value or a string of color indicating whether draw border.
#'
#' @param cluster_rows A logical value indicating whether perform clustering on rows.
#' @param clustering_distance_rows Method of measuring distance when clustering on rows, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_rows Method to perform hierarchical clustering on rows, pass to \code{\link[stats]{hclust}}.
#' @param row_dend_side The position of the row dendrogram ("left", "right").
#' @param row_dend_width A \code{\link[grid]{unit}} object for the width of the row dendrogram.
#' @param row_sorted A vector of sorted row names. It only works when `cluster_rows = F`.
#'
#' @param show_row_names  A logical value indicating whether show row names.
#' @param row_names_side The position of row names ("left", "right").
#' @param row_names_size The size of row names.
#' @param row_names_rot The rotation of row names.
#'
#' @param cluster_columns A logical value indicating whether perform clustering on columns.
#' @param clustering_distance_columns Method of measuring distance when clustering on columns, pass to \code{\link[stats]{dist}}.
#' @param clustering_method_columns Method to perform hierarchical clustering on columns, pass to \code{\link[stats]{hclust}}.
#' @param column_dend_side The position of the column dendrogram ("top", "bottom").
#' @param column_dend_height A \code{\link[grid]{unit}} object for the height of the column dendrogram.
#' @param column_sorted A vector of sorted column names. It only works when `cluster_columns = F`.
#'
#' @param show_column_names A logical value indicating whether show column names.
#' @param column_names_side The position of column names ("top", "column").
#' @param column_names_size The size of column names.
#' @param column_names_rot The rotation of column names.
#'
#' @param anno_param_row_pheno A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#' @param anno_param_column_region A list contains parameters for the region annotation. These can be any of the following:
#' "show", "height", "border", "name_size", "name_rot" and "name_side".
#' @param anno_param_row_stat A list contains parameters for the stat annotation of rows. These can be any of the following:
#' "show", "width", "border", "title", "title_size", "title_rot", "title_side", "axis_side", "axis_at", "axis_labels" and "axis_labels_size".
#' @param anno_param_column_stat A list contains parameters for the stat annotation of columns. These can be any of the following:
#' "show", "height", "border", "title", "title_size", "title_rot", "title_side", "axis_side", "axis_at", "axis_labels" and "axis_labels_size".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#'
#' @param use_raster A logical value indicating whether render the heatmap body as a raster image, pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#'
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation rowAnnotation anno_density
#' @importFrom grid gpar
#'
#' @export

cov_heatmap <- function( cov_obj,
                         add_pheno_info = NULL,
                         add_region_info = NULL,

                         cov_colors = c("white", "#5680ae"),
                         region_info_color_list = NULL,
                         pheno_info_color_list = NULL,

                         border = T,

                         cluster_rows = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "complete",
                         row_dend_side = "left",
                         row_dend_width = grid::unit(5, "mm"),
                         row_sorted = c(),

                         show_row_names = F,
                         row_names_side = "left",
                         row_names_size = 10,
                         row_names_rot = 0,

                         cluster_columns = F,
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "complete",
                         column_dend_side = "top",
                         column_dend_height = grid::unit(5, "mm"),
                         column_sorted = c(),

                         show_column_names = F,
                         column_names_side = "bottom",
                         column_names_size = 10,
                         column_names_rot = 90,

                         anno_param_row_pheno = list(show = T,  width = 5, border = F,
                                                    name_size = NULL, name_rot = 90, name_side = "top"),
                         anno_param_column_region = list(show = T, height = 5, border = FALSE,
                                                       name_size = NULL, name_rot = 0, name_side = "right"),
                         anno_param_row_stat = list(show = F, width = 10, border = FALSE,
                                                    title = "Coverage", title_size = 10, title_side = "bottom", title_rot = 0,
                                                    axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),
                         anno_param_column_stat = list(show = F, height = 10, border = FALSE,
                                                       title = "Coverage", title_size = 10, title_side = "left", title_rot = 0,
                                                       axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),

                         legend_side = "right",
                         legend_title = "Coverage",
                         legend_title_size = NULL,
                         legend_text_size = NULL,
                         legend_grid_size =  grid::unit(4, "mm"),

                         use_raster = NULL

){

  check_class(cov_obj, "COV")

  cov_pheno <- cov_obj@sample
  if(length(cov_pheno) > 0 && !is.null(add_pheno_info)){
    add_pheno_info <- match.arg(add_pheno_info, names(cov_pheno), several.ok = T)
    cov_pheno <- cov_pheno[names(cov_pheno) %in% add_pheno_info]
  } else {
    cov_pheno <- NULL
  }

  cov_region <- cov_obj@region
  if(length(cov_region) > 0 && !is.null(add_region_info)){
    add_region_info <- match.arg(add_region_info, names(cov_region), several.ok = T)
    cov_region <- cov_region[names(cov_region) %in% add_region_info]
  } else {
    cov_region <- NULL
  }

  cov_data <- cov_obj@cov_data

  colors_pheno_info <- get_anno_palette(c(pheno_info_color_list, region_info_color_list), c(cov_pheno, cov_region))
  color_pheno <- colors_pheno_info[names(cov_pheno)]
  color_info <- colors_pheno_info[names(cov_region)]


  anno_param_row_pheno_def_args <- list(show = T,  width = 5, border = F,
                                       name_size = NULL, name_rot = 90, name_side = "top")
  anno_param_row_pheno <- merge_args(anno_param_row_pheno_def_args, anno_param_row_pheno)
  anno_param_column_region_def_args = list(show = T, height = 5, border = FALSE,
                                         name_size = NULL, name_rot = 0, name_side = "right")
  anno_param_column_region <- merge_args(anno_param_column_region_def_args, anno_param_column_region)
  anno_param_row_stat_def_args = list(show = T, width = 10, border = FALSE,
                                      title = "Coverage", title_size = 10, title_side = "bottom", title_rot = 0,
                                      axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_row_stat <- merge_args(anno_param_row_stat_def_args, anno_param_row_stat)
  anno_param_column_stat_def_args = list(show = T, height = 10, border = FALSE,
                                         title = "Coverage", title_size = 10, title_side = "left", title_rot = 0,
                                         axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_column_stat <- merge_args(anno_param_column_stat_def_args, anno_param_column_stat)


  if(length(cov_pheno) > 0){
    data_pheno <- as.data.frame(cov_pheno)

    if(anno_param_row_pheno$show){
      anno_right <- get_anno_row(data_pheno, color_pheno, anno_param_row_pheno)
    } else {
      anno_right <- NULL
    }
  } else {
    anno_right <- NULL
  }


  if(length(cov_region) > 0){
    data_info <-  as.data.frame(cov_region)

    if(anno_param_column_region$show){
      anno_bottom <- get_anno_column(data_info, color_info, anno_param_column_region)
    } else {
      anno_bottom <- NULL
    }
  } else {
    anno_bottom <- NULL
  }

  if(anno_param_row_stat$show){
    anno_left <- ComplexHeatmap::rowAnnotation(
      Coverage = ComplexHeatmap::anno_density(
        t(cov_data),
        border = anno_param_row_stat$border,
        width = grid::unit(anno_param_row_stat$width, 'mm'),
        axis_param = list(side = anno_param_row_stat$axis_side,
                          at = anno_param_row_stat$axis_at,
                          labels = anno_param_row_stat$axis_labels,
                          direction = "reverse",
                          gp = grid::gpar(fontsize = anno_param_row_stat$axis_labels_size, fontface = "bold"))
      ),
      annotation_label = anno_param_row_stat$title,
      annotation_name_side = anno_param_row_stat$title_side,
      annotation_name_rot = anno_param_row_stat$title_rot,
      annotation_name_gp = grid::gpar(fontsize = anno_param_row_stat$title_size, fontface = "bold")
    )
  } else {
    anno_left <- NULL
  }


  if(anno_param_column_stat$show){
    anno_top <- ComplexHeatmap::HeatmapAnnotation(
      Coverage = ComplexHeatmap::anno_density(
        t(cov_data),
        border = anno_param_column_stat$border,
        height = grid::unit(anno_param_column_stat$height, 'mm'),
        axis_param = list(side = anno_param_column_stat$axis_side,
                          at = anno_param_column_stat$axis_at,
                          labels = anno_param_column_stat$axis_labels,
                          gp = grid::gpar(fontsize = anno_param_column_stat$axis_labels_size, fontface = "bold"))
      ), annotation_label = anno_param_column_stat$title,
      annotation_name_side = anno_param_column_stat$title_side,
      annotation_name_rot = anno_param_column_stat$title_rot,
      annotation_name_gp = grid::gpar(fontsize = anno_param_column_stat$title_size, fontface = "bold")
    )
  } else {
    anno_top <- NULL
  }


  ht_main <- ComplexHeatmap::Heatmap(
    t(cov_data),
    name = "main",
    use_raster = use_raster,
    col = cov_colors,
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      title = legend_title,
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
      grid_height = legend_grid_size,
      grid_width = legend_grid_size
    ),
    border = border,
    cluster_rows = cluster_rows,
    clustering_distance_rows = clustering_distance_rows,
    clustering_method_rows = clustering_method_rows,
    cluster_columns = cluster_columns,
    clustering_distance_columns = clustering_distance_columns,
    clustering_method_columns = clustering_method_columns,
    row_dend_side = row_dend_side,
    row_dend_width = row_dend_width,
    column_dend_side = column_dend_side,
    column_dend_height = column_dend_height,
    show_row_names = ifelse(show_row_names & row_names_side == "left", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold"),
    show_column_names = show_column_names,
    column_names_rot = column_names_rot,
    column_names_gp = grid::gpar(fontsize = column_names_size, fontface = "bold"),
    column_names_side = column_names_side,
    left_annotation = anno_left,
    bottom_annotation = anno_bottom,
    top_annotation = anno_top
  )

  ht_right <- ComplexHeatmap::Heatmap(
    matrix(NA, ncol = 0, nrow = nrow(t(cov_data)),
           dimnames = list(rownames(t(cov_data)))),
    show_heatmap_legend = F,
    rect_gp =  grid::gpar(type = "none"),
    show_row_names = ifelse(show_row_names & row_names_side == "right", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold"),
    show_column_names = F,
    right_annotation = anno_right
  )

  lg_pheno <- get_legend(color_pheno, cov_pheno, legend_title_size, legend_text_size, legend_grid_size)
  lg_info <- get_legend(color_info, cov_region, legend_title_size, legend_text_size, legend_grid_size)

  ComplexHeatmap::draw(ht_main + ht_right,
                       main_heatmap = "main",
                       auto_adjust = FALSE,
                       heatmap_legend_list = c(lg_pheno, lg_info),
                       merge_legend = T,
                       heatmap_legend_side = legend_side)

}



