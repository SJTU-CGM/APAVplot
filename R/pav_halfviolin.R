



#' Show region number
#'
#' Plot a half-violin plot for a object of PAV class.
#'
#' @param pav_obj A PAV object.
#' @param violin_color A string of color for half-violin plot.
#' @param add_pheno_info A character string of `pheno_info` names.
#' @param pheno_info_color_list A list contains named vector of colors for `pheno_info` annotation.
#' e.g. list(gender = c("Male" = "green", "Female" = "red"))
#' @param y_title The text for the y-axis title.
#' @param x_text_size The size of tick labels on x-axis.
#' @param y_text_size The size of tick labels on y-axis.
#' @param x_title_size The size of x-axis title.
#' @param y_title_size The size of y-axis title.
#'
#' @importFrom ggplot2 geom_polygon geom_jitter scale_x_continuous theme_classic theme
#'
#' @export


pav_halfviolin <- function(pav_obj,
                     violin_color = "#7e9bc0",
                     add_pheno_info = NULL,
                     pheno_info_color_list = NULL,
                     y_title = "Target Region Number",
                     x_text_size  = NULL,
                     y_text_size = NULL,
                     x_title_size = NULL,
                     y_title_size = NULL){

  pav_data <- pav_obj@pav_data

  sample_data <-
    data.frame(
      sample = colnames(pav_data),
      number = colSums(pav_data)
    )

  sample <- pav_obj@sample
  if(!is.null(add_pheno_info) && length(pav_obj@sample$pheno) > 0){
    add_pheno_info <- match.arg(add_pheno_info, names(sample$pheno))
    if(is.numeric(sample$pheno[[add_pheno_info]])){
      stop("Only accepts discrete phenotype.")
    }
  } else {
    add_pheno_info <- NULL
  }

  loc = den = number = pheno = loc = x = NULL

  if(is.null(add_pheno_info)){
    p_data <- data.frame(loc = stats::density(sample_data$number)$x, den = stats::density(sample_data$number)$y)
    p_data$den <- p_data$den / max(p_data$den) /2
    p_data <- subset(p_data, loc >= min(sample_data$number) & loc <= max(sample_data$number))
    p_data <- rbind(p_data,
                    c(max(sample_data$number), 0),
                    c(min(sample_data$number), 0))
    p <- ggplot() +
      geom_polygon(data = p_data, aes(x = -den , y = loc), fill = violin_color) +
      geom_jitter(data = sample_data, aes(x = .25, y = number), width = .125) +
      scale_x_continuous(limits = c(-.5, .5)) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())

  } else {
    sample_data$pheno <- sample$pheno[[add_pheno_info]][match(sample$name, sample_data$sample)]
    pheno_list <- sort(unique(sample_data$pheno))
    p_data <- do.call(rbind, lapply(1:length(pheno_list), function(x){
      curr_pheno <- pheno_list[x]
      curr_data <- subset(sample_data, pheno == curr_pheno)
      if(nrow(curr_data) == 1) return(NULL)
      den_data <- data.frame(loc = stats::density(curr_data$number)$x, den = stats::density(curr_data$number)$y)
      den_data <- subset(den_data, loc >= min(curr_data$number) & loc <= max(curr_data$number))
      data.frame(loc = c(den_data$loc, rev(range(curr_data$number))),
                 den = c(den_data$den / max(den_data$den) / 2, 0, 0),
                 pheno = pheno_list[x],
                 x = x)
    }))
    sample_data$x <- match(sample_data$pheno, pheno_list)
    p <- ggplot() +
      geom_polygon(data = p_data, aes(x = -den + x, y = loc, fill = pheno)) +
      geom_jitter(data = sample_data, aes(x = x + .25, y = number, fill = pheno), width = .125) +
      scale_x_continuous(breaks = 1:length(pheno_list), labels = pheno_list) +
      labs(x = add_pheno_info) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = x_text_size, color = "black", face = "bold"),
            axis.title.x = element_text(size = x_title_size, color = "black", face = "bold"))

    color_pheno <- get_anno_palette(pheno_info_color_list, sample$pheno[add_pheno_info], "light")
    p <- p + scale_fill_manual(values = color_pheno[[add_pheno_info]])

  }

  p + labs(y = y_title) +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = y_text_size, color = "black", face = "bold"),
          axis.title.y = element_text(size = y_title_size, color = "black", face = "bold"))

}













