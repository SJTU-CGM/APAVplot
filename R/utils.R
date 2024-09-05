

match_color <- function(arg_name, def, input, len = NULL){
  if(!is.null(len)){
    input_n <- length(input)
    if(input_n < len){
      colors <- c(input, def[(input_n+1):len])
    }else{
      colors <- input[1:len]
    }
  }else{
    colors <- input
  }
  if(is.color(colors)){
    return(colors)
  }else{
    warning(paste0(arg_name, " does not seem like colors, use the default palette."));
    return(def)
  }

}


trunc_str <- function(strs, n){
  unlist(lapply(strs,
                function(x){
                  if(nchar(x) > n){
                    paste0(substring(x, 1, n), "...")
                  }else{
                    x
                  }
                }))
}


get_attr <- function(str, attr){
  res <- unlist(lapply(unlist(strsplit(str, ";")),
                       function(x){
                         l <- unlist(strsplit(x, "="))
                         if(l[1] == attr){
                           l[2]
                         }else{
                           NULL
                         }
                       }))
  if(is.null(res)){
    NA
  }else{
    res
  }
}

add_comma <- function(numbers){
  unlist(lapply(numbers, function(n){
    str <- rev(unlist(strsplit(as.character(n),"")))
    round <- floor(length(str)/3)
    chr <- c(unlist(lapply(1:round, function(i){c(str[((i-1)*3+1):(i*3)], ",")})),
             str[min((round*3+1),length(str)):length(str)])
    paste0(rev(chr), collapse = "")
  }))
}



check_info <- function(data, names, info_name, region_name){
  if(!is.data.frame(data)){
    stop("The `", info_name, "_info` should be a data frame.")
  }
  if(!all(names %in% rownames(data))){
    err_region <- setdiff(names, rownames(data))
    if(length(err_region) <= 5){
      err_region_message <- paste0(err_region, collapse = ",")
    } else {
      err_region_message <- paste0(paste0(err_region[1:5], collapse = ","), "...")
    }
    stop("Missing ", region_name, ": '", err_region_message, "' in ", info_name ,"_info.")
  } else {
    data[names, , drop = F]
  }
}


match_num <- function(arg_name, arg, min = -Inf, max = Inf){
  if(is.numeric(arg)){
    if(arg > min && arg < max){
      arg
    } else {
      stop(paste0(arg_name, "should less than ", max, " and greater than ", min, "."))
    }
  } else {
    stop(paste0(arg_name, " should be a numerical type."))
  }
}



get_type_info <- function(sample_n, private, softcore, prob, use_bin, p_value){
  if(softcore){
    if(use_bin){
      sample_dist_n <- dist_n(sample_n, prob, p_value)
    }else{
      sample_dist_n <- round(sample_n * prob)
    }
    type_info <- data.frame(
      sample_n = 1:sample_n,
      region_type = c(rep("Distributed",sample_n - 1 - sample_dist_n),
                    rep("Softcore",sample_dist_n),"Core")
    )
    dist_sample_n <- c(1, sample_n - 1 - sample_dist_n)
  } else {
    type_info <- data.frame(
      sample_n = 1:sample_n,
      region_type = c(rep("Distributed", sample_n - 1), "Core")
    )
    dist_sample_n <- c(1, sample_n - 1)
  }
  if(private){
    type_info[1, 2] <- "Private"
    dist_sample_n[1] <- 2
  }

  return(list(dist_sample_n = dist_sample_n, type_info = type_info))
}


#' @importFrom stats binom.test
dist_n <- function(n, prob, p_value){
  for(i in 1:n) {
    binom_res <- stats::binom.test(x = i, n = n, p = prob, alternative = 'greater')
    if(binom_res$p.value <= p_value){
      res <- i
      break
    }
  }
  return(res-1)
}



check_pav_cov_data <- function(data_name, data){
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("The `", data_name, "` should be a data frame or a matrix")
  } else if(!all(apply(data, 2, is.numeric))){
    stop("The `", data_name, "` should only contain numeric.")
  } else if(any(rowSums(data) == 0)){
    stop("Existing region in `", data_name, "` that is 0(absence) in all samples. Please remove.")
  } else {
    data
  }
}

get_region_type <- function(pav, type_info){
  type_data <- merge(
    data.frame(
      region = rownames(pav),
      col_sum = rowSums(pav)),
    type_info,
    by.x = "col_sum", all.x = TRUE,
    by.y = colnames(type_info)[1],
    sort = TRUE )
  type_data[match(rownames(pav), type_data$region), 3]
}






check_class <- function(obj, class_name){
  if(class(obj) != class_name){
    stop("Please input a `", class_name, "` object")
  }
}

match_logi <- function(arg_name, arg){
  if(is.logical(arg)){
    arg
  } else {
    stop(paste0(arg_name, " should be a logical type."))
  }
}

get_type_palette <- function(colors){
  groups <- c("Core", "Softcore", "Distributed", "Private")
  bc <- structure(c("#8FB4DC", "#debc58", "#B95758", "#A38CC2"),
                            names = groups)
  if(is.null(colors)){
    res <- bc
  }else{
    if(is.vector(colors) && is.color(colors)){
      if(!is.null(names(colors))){
        colors <- colors[intersect(groups, names(colors))]
      }else{
        if(length(colors) > length(groups)){
          colors <- colors[1:length(groups)]
        }
        names(colors) <- groups[1:length(colors)]
      }
      if(length(colors) == length(groups)){
        res <- colors
      }else{
        res <- merge_args(bc, colors)
      }
    }else{
      stop("Please input a vector of colors.")
    }
  }
  return(res)
}

get_gene_palette <- function(colors){
  groups <- c("transcript", "exon", "CDS", "five_prime_UTR", "three_prime_UTR")
  bc <- structure(c("gray70", "gray85", "#d77478", "#ecbf40", "#ecbf40"),
                  names = groups)
  if(is.null(colors)){
    res <- bc
  }else{
    if(is.vector(colors) && is.color(colors)){
      if(!is.null(names(colors))){
        colors <- colors[intersect(groups, names(colors))]
      }else{
        if(length(colors) > length(groups)){
          colors <- colors[1:length(groups)]
        }
        names(colors) <- groups[1:length(colors)]
      }
      if(length(colors) == length(groups)){
        res <- colors
      }else{
        res <- merge_args(bc, colors)
      }
    }else{
      stop("Please input a vector of colors.")
    }
  }
  return(res)
}

get_pav_palette <- function(colors){
  groups <- c("presence", "absence")
  bc <- structure(c("#5680ae", "gray85"),
                  names = groups)
  if(is.null(colors)){
    res <- bc
  }else{
    if(is.vector(colors) && is.color(colors)){
      if(!is.null(names(colors))){
        colors <- colors[intersect(groups, names(colors))]
      }else{
        if(length(colors) > length(groups)){
          colors <- colors[1:length(groups)]
        }
        names(colors) <- groups[1:length(colors)]
      }
      if(length(colors) == length(groups)){
        res <- colors
      }else{
        res <- merge_args(bc, colors)
      }
    }else{
      stop("Please input a vector of colors.")
    }
  }
  return(res)
}



get_palette <- function(data_list, mode, num_re_func = T){
  if(length(data_list) == 0) return(NULL)

  data_list_continue <- data_list[unlist(lapply(data_list, is.numeric))]
  data_list_discrete <- data_list[unlist(lapply(data_list, is.character))]

  all_n <- length(unique(unlist(data_list_discrete))) + length(data_list_continue)
  all_colors <- def_palette(all_n, mode)

  ## continue
  color_n <- length(data_list_continue)
  if(color_n > 0){
    colors <- all_colors[1:color_n]
    res_continue <- lapply(1:color_n, function(x){
      if(length(unique(data_list_continue[[x]])) == 1){
        structure(colors[x], names = unique(data_list_continue[[x]]))
      } else {
        if(num_re_func){
          circlize::colorRamp2(breaks = c(max(data_list_continue[[x]], na.rm = T),
                                          min(data_list_continue[[x]], na.rm = T)),
                               colors = c(colors[x], "gray90"))
        }else{
          c("gray90", colors[x])
        }
      }
    })
    names(res_continue) <- names(data_list_continue)
    res_continue
  } else {
    res_continue <- NULL
  }

  ## discrete
  if(length(data_list_discrete) > 0){
    named_colorN <- unlist(lapply(data_list_discrete, function(x){length(unique(x))}))
    colors <- all_colors[(color_n+1):all_n]
    Ns <- c(0, cumsum(as.vector(named_colorN)))
    res_discrete <- lapply(1:length(named_colorN), function(x){
      cur_col <- colors[(Ns[x]+1):Ns[x+1]]
      names(cur_col) <- unique(data_list_discrete[[x]])
      cur_col
    })
    names(res_discrete) <- names(named_colorN)
  } else {
    res_discrete <- NULL
  }

  return(c(res_continue, res_discrete))
}

get_anno_palette <- function(input_colors, data_list, mode = "def", num_re_func = T){
  def_colors <- get_palette(data_list, mode, num_re_func)
  res <- lapply(names(def_colors), function(x){
    x_def <- def_colors[[x]]
    x_input <- input_colors[[x]]
    x_data <- data_list[[x]]
    if(is.null(x_input)){
      x_def
    } else {
      if(is.numeric(x_data)){
        if(length(x_input) < 2){
          stop(x," should have at least two colors.")
        }
        if(num_re_func){
          circlize::colorRamp2(breaks = seq(min(x_data, na.rm = T),
                                            max(x_data, na.rm = T),
                                            length.out = length(x_input)),
                               colors = x_input)
        }else{
          x_input
        }
      } else {
        if(!is.null(names(x_input))){
          merge_args(x_def, x_input)
        } else {
          x_def
        }
      }
    }
  })
  names(res) <- names(def_colors)
  res
}



get_legend <- function(color_info, info_data, legend_title_size, legend_text_size, legend_grid_size){

  if(length(color_info) > 0){
    lg_info <- lapply(1:length(color_info), function(x){
      if(is.numeric(info_data[[names(color_info)[x]]]) && length(unique(info_data[[names(color_info)[x]]])) > 1){
        ComplexHeatmap::Legend(col_fun = color_info[[x]],
                               title = names(color_info)[x],
                               title_gp = grid::gpar(fontsize = legend_title_size,
                                                     fontface = "bold"),
                               labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
                               grid_height = legend_grid_size,
                               grid_width = legend_grid_size)
      } else {
        ComplexHeatmap::Legend(labels = names(color_info[[x]]),
                               legend_gp = grid::gpar(fill = color_info[[x]]),
                               title = names(color_info)[[x]],
                               title_gp = grid::gpar(fontsize = legend_title_size,
                                                     fontface = "bold"),
                               labels_gp = grid::gpar(fontsize = legend_text_size, fontface = "bold"),
                               grid_height = legend_grid_size,
                               grid_width = legend_grid_size)
      }

    })
  } else {
    lg_info <- NULL
  }
}

get_anno_row <- function(data_pheno, color_pheno, row_anno_pheno){
  anno <- NULL
  eval(parse(text = paste0(
    "anno <- ComplexHeatmap::rowAnnotation(",
    paste0(lapply(names(data_pheno), function(x){
      paste0(x, " = ComplexHeatmap::anno_simple(data_pheno[['", x, "']], border = row_anno_pheno$border,
             width = grid::unit(row_anno_pheno$width, 'mm'),
             col = structure(color_pheno[['", x, "']], names = names(color_pheno[['", x, "']])))")
    }), collapse = ","),
    ", annotation_name_side = row_anno_pheno$name_side, annotation_name_rot = row_anno_pheno$name_rot,
    annotation_name_gp = grid::gpar(fontsize = row_anno_pheno$name_size, fontface = 'bold')    )"
  )))
  return(anno)
}

get_anno_column <- function(data_info, color_info, column_anno_region){
  anno <- NULL
  eval(parse(text = paste0(
    "anno <- ComplexHeatmap::HeatmapAnnotation(",
    paste0(lapply(names(data_info), function(x){
      paste0(x, " = ComplexHeatmap::anno_simple(data_info[['", x, "']], border = column_anno_region$border,
               height = grid::unit(column_anno_region$height,'mm'),
               col = structure(color_info[['",x,"']], names = names(color_info[['",x,"']])))")
    }), collapse = ","),
    ", annotation_name_side = column_anno_region$name_side, annotation_name_rot = column_anno_region$name_rot,
      annotation_name_gp = grid::gpar(fontsize = column_anno_region$name_size, fontface = 'bold') )"
  )))
  return(anno)
}


check_pheno_stat_res <- function(data){
  if(!all(colnames(data) == c("pheno", "region", "p_value", "p_adjusted"))){
    stop("please input result from function `pheno_stat`.")
  }
}


check_obj_pheno <- function(pav_obj){
  if(length(pav_obj@sample$pheno) == 0)
    stop("can't find `pheno_info` in pav_obj.")
}

check_obj_region <- function(pav_obj){
  if(length(pav_obj@region$info) == 0)
    stop("can't find `region_info` in pav_obj.")
}




merge_args <- function(def_args, args){
  return(c(args[names(args) %in% names(def_args)], def_args[!names(def_args) %in% names(args)]))
}

is.color <- function(colors){
  all(unlist(lapply(colors, function(x){
    x %in% colors() || setequal(grep("#[0-9a-fA-F]{6}", x), 1)
  })))
}

def_palette <- function(n, mode = "def"){
  if(n <= 10){
    if(mode == "dark"){
      colors <- c("#a73133", "#376da1", "#e38e28", "#298022", "#7e4fb0",
                  "#e9ba2c", "#b0437a", "#994c21", "#0d7c6f", "#7e306e")
    }else if(mode == "light"){
      colors <- c("#cb7871", "#7e9bc0", "#f2b372", "#76aa6a", "#aa88cb",
                  "#f5d07a", "#ce83a4", "#c08566", "#6da79c", "#aa749c")
    }else{
      colors <- c("#b7514d", "#5880ae", "#ea9e4a", "#4d9242", "#9a71c1",
                  "#efc351", "#bd608c", "#aa643e", "#408e82", "#914e82")
    }
    return(colors[1:n])
  }else{
    if(mode == "dark"){
      colors <- unlist(lapply(seq(0, 1, length.out=n+1),
                              function(x){grDevices::hsv(x, 0.7, 0.65)}))
    }else{
      colors <- unlist(lapply(seq(0, 1, length.out=n+1),
                              function(x){grDevices::hsv(x, 0.45, 0.85)}))
    }
    return(colors[1:n])
  }
}

get_color <- function(colors, groups){
  bc <- structure(def_palette(length(groups)), names = groups)
  if(is.null(colors)){
    res <- bc
  }else{
    if(is.vector(colors) && is.color(colors)){
      if(!is.null(names(colors))){
        colors <- colors[intersect(groups, names(colors))]
      }else{
        if(length(colors) > length(groups)){
          colors <- colors[1:length(groups)]
        }
        names(colors) <- groups[1:length(colors)]
      }
      if(length(colors) == length(groups)){
        res <- colors
      }else{
        res <- merge_args(bc, colors)
      }
    }else{
      stop("Please input a vector of colors.")
    }
  }
  return(res)
}

