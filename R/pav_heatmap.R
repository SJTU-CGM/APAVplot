

#' pav_heatmap
#'
#' Plot a heatmap for a object of PAV class.
#'
#' @param pav_obj A PAV object.
#' @param region_type A vector of region types. These can be any of the following: "Core", "Softcore", "Distributed" and "Private".
#' @param add_pheno_info A character vector of `pheno_info` names.
#' @param add_region_info A character vector of `region_info` names.
#'
#' @param pav_colors  A named vector of colors for presence and absence. e.g. c(presence = "#5680ae", absence = "gray70")
#' @param type_colors A named vector of colors for types. e.g. c("distributed" = "red")
#' @param region_info_color_list A list contains named vector of colors for `region_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
#' @param pheno_info_color_list A list contains named vector of colors for `pheno_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"), age = c("yellow", "red"))
#'
#' @param border A logical value or a string of color indicating whether draw border.
#'
#' @param split_block A logical value indicating whether split columns based on region types.
#' @param block_name_size The size of block name.
#' @param block_name_rot The rotation of block name.
#'
#' @param cluster_rows A logical value indicating whether perform clustering on rows.
#' @param clustering_distance_rows Method of measuring distance when clustring on rows, pass to \code{\link[stats]{dist}}.
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
#' @param column_sorted A vector of sorted column names. It only works when `cluster_columns = F` and `split_block = F`.
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
#' @param use_raster Whether render the heatmap body as a raster image. pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @import ComplexHeatmap
#'
#' @export


pav_heatmap <- function(
    pav_obj,
    region_type,
    add_pheno_info = NULL,
    add_region_info = NULL,

    border = T,
    split_block = T,  # when sort is null.
    block_name_size = NULL,
    block_name_rot = 0,

    cluster_rows = F,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    row_dend_side = "left",
    row_dend_width = grid::unit(5, "mm"),
    row_sorted = c(),  #only work when `cluster_rows` = F

    show_row_names = F,
    row_names_side = "right",
    row_names_size = 10,
    row_names_rot = 0,

    cluster_columns = F,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    column_dend_side = "top",
    column_dend_height = grid::unit(5, "mm"),
    column_sorted = c(), #only work when `columncluster` = F

    show_column_names = F,
    column_names_side = "bottom",
    column_names_size = 10,
    column_names_rot = 90,

    anno_param_row_pheno = list(show = T,  width = 5, border = F,
                               name_size = NULL, name_rot = 90, name_side = "top"),
    anno_param_column_region = list(show = T, height = 5, border = FALSE,
                                  name_size = NULL, name_rot = 0, name_side = "right"),
    anno_param_row_stat = list(show = T, width = 10, bar_width = 1, border = FALSE,
                               title = "Presence\nNumber", title_size = 10, title_side = "bottom", title_rot = 0,
                               axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),
    anno_param_column_stat = list(show = T, height = 10, bar_width = 1, border = FALSE,
                                  title = "Presence\nNumber", title_size = 10, title_side = "left", title_rot = 0,
                                  axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8),

    pav_colors = c(presence = "#5680ae", absence = "gray70"),
    type_colors = NULL,
    region_info_color_list = NULL,
    pheno_info_color_list = NULL,

    legend_side = "right",
    legend_title = list(pav = "PAV", type = "Region"),
    legend_title_size = NULL,
    legend_text_size = NULL,
    legend_grid_size = grid::unit(4, "mm"),

    use_raster = NULL
){

  check_class(pav_obj, "PAV")

  pav_pheno <- pav_obj@sample$pheno
  if(length(pav_pheno) > 0 && !is.null(add_pheno_info)){
    add_pheno_info <- match.arg(add_pheno_info, names(pav_pheno), several.ok = T)
    pav_pheno <- pav_pheno[names(pav_pheno) %in% add_pheno_info]
  } else {
    pav_pheno <- NULL
  }

  pav_region_info <- pav_obj@region$info
  if(length(pav_region_info) > 0 && !is.null(add_region_info)){
    add_region_info <- match.arg(add_region_info, names(pav_region_info), several.ok = T)
    pav_region_info <- pav_region_info[names(pav_region_info) %in% add_region_info]
  } else {
    pav_region_info <- NULL
  }

  region_type <- match.arg(region_type, c("Core", "Softcore", "Distributed", "Private"), several.ok = TRUE)
  split_block <- match_logi("split_block", split_block)

  ## region_data
  regions_data <- data.frame(pav_obj@region[1:3], stringsAsFactors = F)
  if(length(pav_region_info) != 0){
    regions_data <- cbind(regions_data, data.frame(pav_region_info, stringsAsFactors = F))
  }
  rownames(regions_data) <- regions_data$name

  ## data_main
  data_main <- t(pav_obj@pav_data)

  if(!cluster_columns){
    if(length(column_sorted) == nrow(regions_data) && all(column_sorted %in% rownames(regions_data))){
      data_main <- data_main[, column_sorted]
      regions_data <- regions_data[column_sorted, ]
      split_block <- F
    } else {
      data_main <- data_main[, names(sort(colSums(data_main), decreasing = TRUE))]
      regions_data <- regions_data[colnames(data_main), ]
    }
  }

  type <- NULL
  regions_data <- subset(regions_data, type %in% region_type)
  data_main <- data_main[, rownames(regions_data)]

  ##sample_data
  samples_data <- data.frame(pav_obj@sample[1], stringsAsFactors = F)
  if(length(pav_pheno) != 0){
    samples_data <- cbind(samples_data, data.frame(pav_pheno, stringsAsFactors = F))
  }
  rownames(samples_data) <- samples_data$name

  if(!cluster_rows){
    if(length(row_sorted) == nrow(data_main) && all(row_sorted %in% rownames(data_main))){
      data_main <- data_main[row_sorted,]
      samples_data <- samples_data[row_sorted,]
    } else {
      data_main <- data_main[names(sort(rowSums(data_main), decreasing = TRUE)), , drop = F]
      samples_data <- samples_data[rownames(data_main), ]
    }
  }

  color_pav <- get_pav_palette(pav_colors)

  color_type <- get_type_palette(type_colors)[region_type]

  pheno_info_color_list_info <- get_anno_palette(c(pheno_info_color_list, region_info_color_list), c(pav_pheno, pav_region_info))
  color_pheno <- pheno_info_color_list_info[names(pav_pheno)]
  color_info <- pheno_info_color_list_info[names(pav_region_info)]

  ## anno param

  anno_param_row_pheno_def_args <- list(show = T,  width = 5, border = F,
                                       name_size = NULL, name_rot = 90, name_side = "top")
  anno_param_row_pheno <- merge_args(anno_param_row_pheno_def_args, anno_param_row_pheno)
  anno_param_column_region_def_args = list(show = T, height = 5, border = FALSE,
                                         name_size = NULL, name_rot = 0, name_side = "right")
  anno_param_column_region <- merge_args(anno_param_column_region_def_args, anno_param_column_region)
  anno_param_row_stat_def_args = list(show = T, width = 10, bar_width = 1, border = FALSE,
                                      title = "Presence\nNumber", title_size = 10, title_side = "bottom", title_rot = 0,
                                      axis_side = "bottom",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_row_stat <- merge_args(anno_param_row_stat_def_args, anno_param_row_stat)
  anno_param_column_stat_def_args = list(show = T, height = 10, bar_width = 1, border = FALSE,
                                         title = "Presence\nNumber", title_size = 10, title_side = "left", title_rot = 0,
                                         axis_side = "left",axis_at = NULL, axis_labels = NULL, axis_labels_size = 8)
  anno_param_column_stat <- merge_args(anno_param_column_stat_def_args, anno_param_column_stat)

  ## anno_left

  data_samplePN <- lapply(region_type, function(x){
    rowSums(data_main[, subset(regions_data, type == x)$name, drop = F])})

  if(anno_param_row_stat$show){
    anno_param_row_stat_color <- color_type
    anno_left <- ComplexHeatmap::rowAnnotation(
      PN = ComplexHeatmap::anno_barplot(
        data_samplePN, border = anno_param_row_stat$border,
        width = grid::unit(anno_param_row_stat$width, 'mm'),
        bar_width = anno_param_row_stat$bar_width,
        gp = grid::gpar(fill = anno_param_row_stat_color, col = NA),
        axis_param = list(side = anno_param_row_stat$axis_side,
                          at = anno_param_row_stat$axis_at,
                          labels = anno_param_row_stat$axis_labels,
                          direction = "reverse",
                          gp = grid::gpar(fontsize = anno_param_row_stat$axis_labels_size, fontface = "bold"))
      ), annotation_label = anno_param_row_stat$title,
      annotation_name_side = anno_param_row_stat$title_side,
      annotation_name_rot = anno_param_row_stat$title_rot,
      annotation_name_gp = grid::gpar(fontsize = anno_param_row_stat$title_size, fontface = "bold")
    )
  }else{
    anno_left <- NULL
  }

  ## anno_right

  if(length(pav_pheno) > 0){
    data_pheno <- samples_data[, names(pav_pheno), drop = F]

    if(anno_param_row_pheno$show){
      anno_right <- get_anno_row(data_pheno, color_pheno, anno_param_row_pheno)
    } else {
      anno_right <- NULL
    }
  } else {
    anno_right <- NULL
  }


  ## anno_top

  data_presenceN <- colSums(data_main)

  if(anno_param_column_stat$show){
    anno_param_column_stat_color <- color_type[match(regions_data$type, region_type)]
    anno_top <- ComplexHeatmap::HeatmapAnnotation(
      PN = ComplexHeatmap::anno_barplot(
        data_presenceN, bar_width = anno_param_column_stat$bar_width, border = anno_param_column_stat$border,
        height = grid::unit(anno_param_column_stat$height, 'mm'),
        gp = grid::gpar(fill = anno_param_column_stat_color, col = NA),
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

  #anno_bottom

  if(length(pav_region_info) > 0){
    data_info <- regions_data[, names(pav_region_info), drop = F]

    if(anno_param_column_region$show){
      anno_bottom <- get_anno_column(data_info, color_info, anno_param_column_region)
    } else {
      anno_bottom <- NULL
    }
  } else {
    anno_bottom <- NULL
  }


  lg <- list(
    ComplexHeatmap::Legend(
      labels = names(color_pav),
      legend_gp = grid::gpar(fill = color_pav),
      title = ifelse(is.null(legend_title$pav), "PAV", legend_title$pav),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
      grid_height = legend_grid_size,
      grid_width = legend_grid_size
    ), ComplexHeatmap::Legend(
      labels = region_type,
      legend_gp = grid::gpar(fill = color_type),
      title = ifelse(is.null(legend_title$type), "region", legend_title$type),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
      grid_height = legend_grid_size,
      grid_width = legend_grid_size
    ))

  lg_pheno <- get_legend(color_pheno, pav_pheno, legend_title_size, legend_text_size, legend_grid_size)
  lg_info <- get_legend(color_info, pav_region_info, legend_title_size, legend_text_size, legend_grid_size)

  if(split_block){
    column_split <- factor(regions_data$type, levels = region_type)
  }else {
    column_split <- NULL
  }

  ht_main <- ComplexHeatmap::Heatmap(
    data_main,
    name = "main",
    use_raster = use_raster,
    col = as.vector(color_pav[c("absence", "presence")]),
    show_heatmap_legend = F,
    column_split =column_split,
    column_gap = grid::unit(0, "mm"),
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
    column_names_rot = column_names_rot,
    column_names_gp = grid::gpar(fontsize = column_names_size, fontface = "bold"),
    column_names_side = column_names_side,
    column_title_rot = block_name_rot,
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = block_name_size, fontface = "bold"),
    show_column_names = show_column_names,
    cluster_column_slices = F,
    show_row_names = ifelse(show_row_names & row_names_side == "left", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold"),
    left_annotation = anno_left,
    bottom_annotation = anno_bottom,
    top_annotation = anno_top
  )

  ht_right <- ComplexHeatmap::Heatmap(
    matrix(NA,ncol = 0, nrow=nrow(data_main),
           dimnames = list(rownames(data_main))),
    show_heatmap_legend = F,
    rect_gp =  grid::gpar(type = "none"),
    show_row_names = ifelse(show_row_names & row_names_side == "right", T, F),
    row_names_rot = row_names_rot,
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold"),
    show_column_names = F,
    right_annotation = anno_right
  )

  ComplexHeatmap::draw(ht_main + ht_right,
                       main_heatmap = "main",
                       auto_adjust = FALSE,
                       heatmap_legend_list = c(lg,lg_pheno, lg_info),
                       merge_legend = T,
                       heatmap_legend_side = legend_side)

}













