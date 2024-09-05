

#' Coverage Data
#'
#' @format ## `cov_data`
#' A data frame with 906 rows and 111 columns:
"cov_data"


#' PAV Data
#'
#' "1" represents presence and "0" represents absence.
#'
#' @format ## `pav_data`
#' A data frame with 906 rows and 111 columns:
"pav_data"


#' Phenotype Data
#'
#'
#' @format ## `pheno_info_data`
#' A data frame with 111 rows and 5 columns:
#' \describe{
#'   \item{Genetic_sex}{Genetic sex}
#'   \item{DNA_source}{DNA source}
#'   \item{Region}{Region}
#'   \item{Country}{Country}
#'   \item{Coverage_mean}{Average of coverage}
#' }
#' @source <https://www.nature.com/articles/nature18964>
"pheno_info_data"


#' Gene information
#'
#'
#' @format ## `gene_info_data`
#' A data frame with 906 rows and 5 columns:
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{length}{Length}
#'   \item{chr_n}{The Number of Chromosome}
#' }
"gene_info_data"


#' Simulation result
#'
#'
#' @format ## `sim_res`
#' A data frame with 1110 rows and 5 columns:
#' \describe{
#'   \item{Round}{Simulation round}
#'   \item{SampleN}{The number of selected samples}
#'   \item{Core}{The number of core regions}
#'   \item{Pan}{The number of all regions}
#'   \item{Delta}{The number of novel regions}
#' }
"sim_res"


#' Simulation result with grouped samples
#'
#'
#' @format ## `sim_group_res`
#' A data frame with 1110 rows and 6 columns:
#' \describe{
#'   \item{Round}{Simulation round}
#'   \item{SampleN}{The number of selected samples}
#'   \item{Core}{The number of core regions}
#'   \item{Pan}{The number of all regions}
#'   \item{Delta}{The number of novel regions}
#'   \item{Group}{The group of samples}
#' }
"sim_group_res"


#' Coverage profile of elements
#'
#' @format ## `ele_cov`
#' A data frame with 11 rows and 15 columns:
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start}{Start position}
#'   \item{End}{End position}
#'   \item{Length}{Length}
#'   \item{Annotation}{Annotation of elements}
#'   \item{sampleN}{Coverage}
#' }
"ele_cov"


#' PAV profile of elements
#'
#' @format ## `ele_pav`
#' A data frame with 11 rows and 15 columns:
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start}{Start position}
#'   \item{End}{End position}
#'   \item{Length}{Length}
#'   \item{Annotation}{Annotation of elements}
#'   \item{sampleN}{PAV}
#' }
"ele_pav"


#' Phenotype data
#'
#' @format ## `ele_pheno`
#' A data frame with 10 rows and 3 columns
"ele_pheno"


#' GFF data
#'
#' @format ## `ele_gff`
#' A data frame with 17 rows and 9 columns
"ele_gff"


#' Depth data
#'
#' @format ## `ele_depth`
#' A data frame with 3034 rows and 12 columns
"ele_depth"


