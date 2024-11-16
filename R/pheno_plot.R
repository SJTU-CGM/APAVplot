

#' Draw a phenotype heat map
#'
#' Show the result of phenotype association in a heat map.
#'
#' @param pav_obj A PAV object.
#' @param pheno_stat_res The result from \code{\link[APAVplot]{pheno_stat}}.
#' @param add_region_info A character vector of `region_info` names.
#'
#' @param p_threshold The threshold of p_value/p_adjusted.
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param only_show_significant A logical value indicating whether only show p_value/p_adjusted that satisfies the condition.
#' @param flip A logical value indicating whether flip cartesian coordinates.
#'
#' @param p_colors A vector of colors or a color mapping function for p_value/p_adjusted, pass to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param na_col A string of color for NA values.
#' @param cell_border_color `NA` or a string of color for the border of cells.
#' @param region_info_color_list A list contains named vector of colors for `region_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
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
#' @param anno_param_region A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#' @import ComplexHeatmap
#'
#' @export

pheno_heatmap <- function(pav_obj,
                         pheno_stat_res,
                         add_region_info = NULL,

                         p_threshold = 0.01,
                         adjust_p = T,
                         only_show_significant = T,
                         flip = F,

                         p_colors = c("#B95758", "#f0d2d0"),
                         na_col = "gray",
                         cell_border_color = "white",
                         region_info_color_list = NULL,

                         cluster_rows = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "complete",
                         row_dend_side = "left",
                         row_dend_width = grid::unit(5, "mm"),
                         row_sorted = c(),

                         show_row_names = T,
                         row_names_side = "right",
                         row_names_size = 10,
                         row_names_rot = 0,

                         cluster_columns = F,
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "complete",
                         column_dend_side = "top",
                         column_dend_height = grid::unit(5, "mm"),
                         column_sorted = c(),

                         show_column_names = T,
                         column_names_side = "bottom",
                         column_names_size = 10,
                         column_names_rot = 90,

                         anno_param_region = list(show = T, width = 5, border = FALSE,
                                                name_size = NULL, name_rot = 90, name_side = "bottom"),

                         legend_side = "right",
                         legend_title = list(pav = "PAV", type = "region"),
                         legend_title_size = NULL,
                         legend_text_size = NULL,
                         legend_grid_size = grid::unit(4, "mm")){


  check_class(pav_obj, "PAV")

  check_pheno_stat_res(pheno_stat_res)

  if(adjust_p){
    pheno_stat_res$p <- pheno_stat_res$p_adjusted
  } else {
    pheno_stat_res$p <- pheno_stat_res$p_value
  }

  pheno_stat_res <-  pheno_stat_res[!is.na(pheno_stat_res$p), ,drop = F]

  p_data <- as.data.frame(data.table::dcast(data.table::data.table(pheno_stat_res),
                                            region ~ pheno, value.var = c("p")))
  rownames_p_data <- p_data[, 1]
  p_data <- p_data[, -1, drop = F]
  p_data <- apply(p_data, 2, as.numeric)
  rownames(p_data) <- rownames_p_data
  p_data <- p_data[apply(p_data, 1, function(x){min(x, na.rm = T)}) < p_threshold, ,drop = F]
  if(nrow(p_data) == 0){
    stop("No results after filtering.")
  }
  p_data <- p_data[order(apply(p_data, 1, function(x){min(x, na.rm = T)})), , drop = F]

  if(only_show_significant){
    p_data <- apply(p_data, 2, function(x){x[x > p_threshold] = NA; x})
  }

  p_data <- p_data[, colSums(p_data, na.rm = T) != 0, drop = F]

  ##

  if(length(pav_obj@region$info) > 0 && length(add_region_info) > 0){

    add_region_info <- match.arg(add_region_info, names(pav_obj@region$info), several.ok = T)

    region_info_data <- as.data.frame(pav_obj@region$info)
    rownames(region_info_data) <- pav_obj@region$name
    region_info_data <- region_info_data[rownames(p_data), add_region_info, drop = F]

    color_info <- get_anno_palette(region_info_color_list, as.list(region_info_data))

    anno_param_region_def_args = list(show = T, width = 5, border = FALSE,
                                    name_size = NULL, name_rot = 90, name_side = "bottom")
    anno_param_region <- merge_args(anno_param_region_def_args, anno_param_region)

    if(anno_param_region$show){
      if(flip){
        anno_param_region$name_side <- ifelse(anno_param_region$name_side == "top", "left", "right")
        anno_param_region$height <- anno_param_region$width
        anno_bottom <- get_anno_column(region_info_data, color_info, anno_param_region)
        anno_right <- NULL
      } else {
        anno_right <- get_anno_row(region_info_data, color_info, anno_param_region)
        anno_bottom <- NULL
      }
    } else {
      anno_right <- NULL
      anno_bottom <- NULL
      color_info <- NULL
    }
  } else {
    anno_right <- NULL
    anno_bottom <- NULL
    color_info <- NULL
  }

  if(flip) p_data <- t(p_data)
  ht <- ComplexHeatmap::Heatmap(p_data,
                                col = p_colors,
                                rect_gp = grid::gpar(col = cell_border_color),
                                name = ifelse(adjust_p, "p_adjusted", "p_value"),
                                right_annotation = anno_right,
                                bottom_annotation = anno_bottom,
                                na_col = na_col,
                                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
                                                            labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
                                                            grid_height = legend_grid_size, grid_width = legend_grid_size),
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
                                show_column_names = show_column_names,
                                column_names_rot = column_names_rot,
                                column_names_gp = grid::gpar(fontsize = column_names_size, fontface = "bold"),
                                column_names_side = column_names_side,
                                show_row_names = show_row_names,
                                row_names_rot = row_names_rot,
                                row_names_side = row_names_side,
                                row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold"))

  lg_info <- get_legend(color_info, region_info_data, legend_title_size, legend_text_size, legend_grid_size)

  ComplexHeatmap::draw(ht,
                       auto_adjust = FALSE,
                       heatmap_legend_list = lg_info,
                       merge_legend = T,
                       heatmap_legend_side = legend_side)
}


#' Draw a phenotype block plot
#'
#' Show the result of phenotype association of a certain phenotype in a block chart.
#'
#' @param pav_obj A PAV object.
#' @param pheno_stat_res The result from \code{\link[APAVplot]{pheno_stat}}.
#' @param pheno_name The name of phenotype used for grouping.
#' @param add_region_info A character vector of `region_info` names. The p_value/p_adjusted will also be added to `region_info`.
#'
#' @param p_threshold The threshold of p_value/p_adjusted.
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param only_show_significant A logical value indicating whether only show p_value/p_adjusted that satisfies the condition.
#' @param flip A logical value indicating whether flip cartesian coordinates.
#'
#' @param per_colors A vector of colors or a color mapping function for absence percentage, pass to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param na_col A string of color for NA values.
#' @param cell_border_color `NA` or a string of color for the border of cells.
#' @param region_info_color_list A list contains named vector of colors for `region_info` annotation.
#' e.g. list(source = c("reference" = "red", "novel" = "blue"), length = c("orange", "red"))
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
#' @param clustering_distance_columns Method of measuring distance when clustring on columns, pass to \code{\link[stats]{dist}}.
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
#' @param anno_param_region A list contains parameters for the phenotype annotation. These can be any of the following:
#' "show", "width", "border", "name_size", "name_rot" and "name_side".
#'
#' @param legend_side The position of legend ("top", "bottom", "right", "left").
#' @param legend_title The text for the legend title.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#' @param legend_grid_size A \code{\link[grid]{unit}} object for the size of legend grid.
#'
#' @export


pheno_block <- function(pav_obj,
                       pheno_stat_res,
                       pheno_name,
                       add_region_info = "p",

                       p_threshold = 0.01,
                       adjust_p = T,
                       only_show_significant = T,
                       flip = F,

                       per_colors = c("#d6deeb", "#376da1"),
                       na_col = "gray",
                       cell_border_color = "white",
                       region_info_color_list = NULL,

                       cluster_rows = T,
                       clustering_distance_rows = "euclidean",
                       clustering_method_rows = "complete",
                       row_dend_side = "left",
                       row_dend_width = grid::unit(5, "mm"),
                       row_sorted = c(),

                       show_row_names = T,
                       row_names_side = "right",
                       row_names_size = 10,
                       row_names_rot = 0,

                       cluster_columns = T,
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns = "complete",
                       column_dend_side = "top",
                       column_dend_height = grid::unit(5, "mm"),
                       column_sorted = c(),

                       show_column_names = T,
                       column_names_side = "bottom",
                       column_names_size = 10,
                       column_names_rot = 90,

                       anno_param_region = list(show = T, width = 5, border = FALSE,
                                              name_size = NULL, name_rot = 90, name_side = "bottom"),

                       legend_side = "right",
                       legend_title = list(pav = "PAV", type = "region"),
                       legend_title_size = NULL,
                       legend_text_size = NULL,
                       legend_grid_size = grid::unit(4, "mm")){

  check_class(pav_obj, "PAV")

  check_pheno_stat_res(pheno_stat_res)

  if(adjust_p){
    pheno_stat_res$p <- pheno_stat_res$p_adjusted
  } else {
    pheno_stat_res$p <- pheno_stat_res$p_value
  }

  pav_data <- pav_obj@pav_data
  check_obj_pheno(pav_obj)
  pheno_data <- as.data.frame(pav_obj@sample$pheno)
  rownames(pheno_data) <- pav_obj@sample$name
  pheno_name <- match.arg(pheno_name, unique(pheno_stat_res$pheno))

  pheno = p = NULL

  p_data <- subset(pheno_stat_res, pheno == pheno_name & p < p_threshold)
  p_data <- p_data[order(p_data$p, decreasing = F), ]
  if(nrow(p_data) == 0){
    stop("No results after filtering.")
    return(NULL)
  }
  p_pav_data <- pav_data[p_data$region, , drop = F]
  groups <- unique(pheno_data[[pheno_name]])

  groups_data <- as.data.frame(lapply(groups, function(x){
    curr_samples <- rownames(pheno_data)[pheno_data[[pheno_name]] == x]
    rowSums(p_pav_data[, curr_samples, drop = F]) / length(curr_samples)
  }))

  colnames(groups_data) <- paste0(names(table(pheno_data[[pheno_name]])),
                                  "(n=", table(pheno_data[[pheno_name]]), ")")

  ##

  region_data <- p_data[, c(ifelse(adjust_p, "p_adjusted", "p_value"), "region"), drop = F]
  if(length(pav_obj@region$info) > 0){
    region_data <- merge(region_data,
                       cbind(region = pav_obj@region$name, as.data.frame(pav_obj@region$info)),
                       by = "region")
  }
  rownames(region_data) <- region_data$region
  region_data <- region_data[, -1]

  if(length(add_region_info) > 0){

    add_region_info <- match.arg(add_region_info, colnames(region_data), several.ok = T)
    region_info_data <- region_data[rownames(groups_data), add_region_info, drop = F]

    color_info <- get_anno_palette(region_info_color_list, as.list(region_info_data))

    anno_param_region_def_args = list(show = T, width = 5, border = FALSE,
                                    name_size = NULL, name_rot = 90, name_side = "bottom")
    anno_param_region <- merge_args(anno_param_region_def_args, anno_param_region)

    if(anno_param_region$show){
      if(flip){
        anno_param_region$name_side <- ifelse(anno_param_region$name_side == "top", "left", "right")
        anno_param_region$height <- anno_param_region$width
        anno_bottom <- get_anno_column(region_info_data, color_info, anno_param_region)
        anno_right <- NULL
      } else {
        anno_right <- get_anno_row(region_info_data, color_info, anno_param_region)
        anno_bottom <- NULL
      }
    } else {
      anno_right <- NULL
      anno_bottom <- NULL
      color_info <- NULL
    }
  } else {
    anno_right <- NULL
    anno_bottom <- NULL
    color_info <- NULL
  }

  if(flip)  groups_data <- t(groups_data)

  ht <- ComplexHeatmap::Heatmap(as.matrix(groups_data),
                                name = "Presence(%)",
                                col = per_colors,
                                rect_gp = grid::gpar(col = cell_border_color),
                                right_annotation = anno_right,
                                bottom_annotation = anno_bottom,
                                na_col = na_col,
                                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
                                                            labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
                                                            grid_height = legend_grid_size, grid_width = legend_grid_size),
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
                                show_column_names = show_column_names,
                                column_names_rot = column_names_rot,
                                column_names_gp = grid::gpar(fontsize = column_names_size, fontface = "bold"),
                                column_names_side = column_names_side,
                                show_row_names = show_row_names,
                                row_names_rot = row_names_rot,
                                row_names_side = row_names_side,
                                row_names_gp = grid::gpar(fontsize = row_names_size, fontface = "bold")
  )

  lg_info <- get_legend(color_info, region_info_data, legend_title_size, legend_text_size, legend_grid_size)

  ComplexHeatmap::draw(ht,
                       auto_adjust = FALSE,
                       heatmap_legend_list = lg_info,
                       merge_legend = T,
                       heatmap_legend_side = legend_side)
}



#' Draw a Manhattan plot
#'
#' Show the result of phenotype association of a certain phenotype in a Manhattan plot.
#'
#' @param pav_obj A PAV object.
#' @param pheno_stat_res The result from \code{\link[APAVplot]{pheno_stat}}.
#' @param pheno_name The name of phenotype used for displaying.
#' @param chr The name in `region_info` denoting chromosomes.
#' @param bp The name in `region_info` denoting chromosomal position.
#'
#' @param adjust_p A logical value indicating whether adjust p_value.
#' @param highlight_top_n The top `n` points will be highlighted.
#' @param highlight_text_size The size of labels on highlight points.
#'
#' @param point_size The size of points.
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param x_text_angle The angle of tick labels.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @export

pheno_manhattan <- function(pav_obj,
                           pheno_stat_res,
                           pheno_name,
                           chr,
                           bp,
                           adjust_p = T,
                           # line_value = -log10(1e-05),

                           highlight_top_n = 5,
                           highlight_text_size = 4,
                           point_size = 1.5,
                           x_text_size = NULL,
                           x_text_angle = 0,
                           x_title_size = NULL,
                           y_text_size = NULL,
                           y_title_size = NULL
){

  check_class(pav_obj, "PAV")

  check_pheno_stat_res(pheno_stat_res)
  check_obj_pheno(pav_obj)
  pheno_name <- match.arg(pheno_name, names(pav_obj@sample$pheno))

  pheno = CHR = P = pos = CHR = REGION = NULL

  check_obj_region(pav_obj)
  region_list <- pav_obj@region
  bp <- match.arg(bp, names(region_list$info))
  if(!is.numeric(region_list$info[[bp]])){
    stop("the column of `bp` should be numerical")
  }
  d <- data.frame(REGION = region_list$name, CHR_name = region_list$info[[chr]],
		  BP = region_list$info[[bp]])
  z <- d$CHR_name
  zf <- as.numeric(as.factor(z)) + suppressWarnings(max(as.numeric(gsub("[chr|Chr]","", z)), na.rm = T))
  d$CHR <- unlist(lapply(1:length(z), function(x){
    n <- gsub("[chr|Chr]","", z[x])
    if(grepl("^[0-9]+$", n)){
      as.numeric(n)
    }else{
      zf[x]
    }
  }))
  d <- merge(d, subset(pheno_stat_res, pheno == pheno_name), by.x = "REGION", by.y = "region")
  if(adjust_p){
    d$P <- d$p_adjusted
  } else {
    d$P <- d$p_value
  }
  d$CHR <- as.numeric(d$CHR)
  d$BP <- as.numeric(d$BP)
  d <- d[order(d$CHR, d$BP),]
  d[is.na(d$CHR), "CHR"] <- "NA"
  d$pos <- 1:nrow(d)
  lengths <- unlist(lapply(unique(d$CHR), function(x){nrow(subset(d, CHR == x))}))
  ticks <- cumsum(lengths) - lengths/2

  n_CHR <- length(unique(d$CHR))
  CHR_colors <- rep(c("black","gray"),ceiling(length(unique(d$CHR))/2))[1:n_CHR]

  d2 <- subset(d, P <= sort(d$P)[highlight_top_n])

  ggplot(subset(d, !is.na(P)), aes(x = pos, y = -log10(P))) +
    geom_point(aes(color = factor(CHR, levels = unique(CHR))), size = point_size) +
    geom_point(data = d2, aes(x = pos, y = -log10(P)), color = "#B92427", size = point_size) +
    ggrepel::geom_text_repel(data = d2, aes(x = pos, y = -log10(P), label = REGION), size = highlight_text_size) +
    # geom_hline(yintercept = line_value, color = "red") +
    labs(x = "Chromosome", y = paste0("-log10(", ifelse(adjust_p, "p_adjusted", "p_value"), ")")) +
    theme_classic() +
    scale_color_manual(values = CHR_colors) +
    scale_x_continuous(expand = expansion(mult = c(.05,.05)), breaks = ticks, labels = unique(d$CHR_name)) +
    scale_y_continuous(expand = expansion(mult = c(.05,.05))) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(size = x_text_size, angle = x_text_angle, color = "black", face = "bold"),
          axis.title.x = element_text(size = x_title_size, color = "black", face = "bold"),
          axis.text.y = element_text(size = y_text_size, color = "black", face = "bold"),
          axis.title.y = element_text(size = y_title_size, color = "black", face = "bold")
    )
}



#' Draw a "PAV-phenotype" bar plot
#'
#' Show the presence/absence of a specified region in a specified phenotype group in a bar plot.
#'
#' @param pav_obj A PAV object.
#' @param pheno_name The name of phenotype.
#' @param region_name The name of region.
#' .
#' @param pav_colors A vector of colors for presence and absence.
#' @param bar_width A numeric vector giving the relative width of bars, ranging from 0 to 1.
#'
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#' @param legend_side The position of legend.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_fill_manual scale_y_continuous theme_classic theme
#'
#' @export

pheno_bar <- function(pav_obj,
                     pheno_name,
                     region_name,
                     pav_colors = c("gray70", "steelblue"),
                     bar_width = .8,
                     x_text_size = NULL,
                     x_title_size = NULL,
                     y_text_size = NULL,
                     y_title_size = NULL,
                     legend_side = "top",
                     legend_title_size = NULL,
                     legend_text_size = NULL
){

  check_class(pav_obj, "PAV")

  check_obj_pheno(pav_obj)
  pheno_name <- match.arg(pheno_name, names(pav_obj@sample$pheno))
  if(!is.character(pav_obj@sample$pheno[[pheno_name]])){
    stop("please select a character phenotype")
  }
  region_name <- match.arg(region_name, pav_obj@region$name)
  pav_data <- pav_obj@pav_data

  pheno = value = pav = NULL

  p_data <- data.table::melt(
    data.table::data.table(table(data.frame(pav = pav_data[region_name,],
                                            pheno = pav_obj@sample$pheno[[pheno_name]]))),
    id.vars = c("pav", "pheno"))

  ggplot(p_data, aes(x = pheno, y = value, fill = as.character(pav))) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    labs(title = region_name, fill = NULL, x = pheno_name, y = "Sample Number") +
    scale_fill_manual(values = pav_colors, breaks = c(0, 1), labels = c("Absence", "Presence")) +
    scale_y_continuous(expand = expansion(mult= c(0, .1))) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_text_size, color = "black", face = "bold"),
          axis.title.x = element_text(size = x_title_size, color = "black", face = "bold"),
          axis.text.y = element_text(size = y_text_size, color = "black", face = "bold"),
          axis.title.y = element_text(size = y_title_size, color = "black", face = "bold"),
          legend.position = legend_side,
          legend.title = element_text(size = legend_title_size, color = "black", face = "bold"),
          legend.text = element_text(size = legend_text_size, color = "black", face = "bold"))
}


#' Draw a "PAV-phenotype" violin plot
#'
#' Show specific phenotype data of the absence/presence group in a violin plot.
#'
#' @param pav_obj A PAV object.
#' @param pheno_name The name of phenotype.
#' @param region_name The name of region.
#' .
#' @param pav_colors A vector of colors for presence and absence.
#'
#' @param x_text_size The size of tick labels on x-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_text_size The size of tick labels on y-axis.
#' @param y_title_size The size of y-axis title.
#'
#' @param legend_side The position of legend.
#' @param legend_title_size The size of legend title.
#' @param legend_text_size The size of legend item labels.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter labs scale_color_manual theme_classic theme scale_x_discrete
#' @importFrom ggsignif geom_signif
#'
#' @export


pheno_violin <- function(pav_obj,
                        pheno_name,
                        region_name,
                        pav_colors = c("gray60", "steelblue"),
                        x_text_size = NULL,
                        x_title_size = NULL,
                        y_text_size = NULL,
                        y_title_size = NULL,
                        legend_side = "top",
                        legend_title_size = NULL,
                        legend_text_size = NULL){

  check_class(pav_obj, "PAV")

  check_obj_pheno(pav_obj)
  pheno_name <- match.arg(pheno_name, names(pav_obj@sample$pheno))
  if(!is.numeric(pav_obj@sample$pheno[[pheno_name]])){
    stop("please select a numerical phenotype.")
  }
  region_name <- match.arg(region_name, pav_obj@region$name)
  pav_data <- pav_obj@pav_data

  pav = pheno = NULL

  p_data <- data.frame(pav = as.character(pav_data[region_name,]),
                       pheno = pav_obj@sample$pheno[[pheno_name]])
  p_data <- p_data[!is.na(p_data$pheno),]

  ggplot(p_data, aes(x = pav, y = pheno, color = pav)) +
    geom_violin() + geom_jitter() +
    labs(title = region_name, color = NULL, y = pheno_name, x = "PAV") +
    ggsignif::geom_signif(comparisons = list(c("0","1")), color = "black") +
    scale_x_discrete(labels = c("0" = "Absence", "1" = "Presence")) +
    scale_y_continuous(expand = expansion(mult= c(0, .1))) +
    scale_color_manual(values = pav_colors, breaks = c("0", "1"), labels = c("Absence", "Presence")) +
    theme_classic() +
    theme(axis.text.x = element_text(size = x_text_size, color = "black", face = "bold"),
          axis.title.x = element_text(size = x_title_size, color = "black", face = "bold"),
          axis.text.y = element_text(size = y_text_size, color = "black", face = "bold"),
          axis.title.y = element_text(size = y_title_size, color = "black", face = "bold"),
          legend.position = legend_side,
          legend.title = element_text(size = legend_title_size, color = "black", face = "bold"),
          legend.text = element_text(size = legend_text_size, color = "black", face = "bold"))
}









