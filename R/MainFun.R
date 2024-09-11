################################################################################
# check_qc ---------------------------------------------------------------------

#' @title check_qc_data
#' @author Zhiwei Zhou
#' @param path working path. It should contain two folders: "poolQC", "ms2". "poolQC" contains
#' @param column 'c18', 'hilic'
#' @param polarity "positive", 'negative'
#' @importFrom magrittr %>%
#' @importFrom data.table %between%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importClassesFrom SpectraTools 'SpectraData'
#' @importClassesFrom massdataset 'ms2_data'
#' @importClassesFrom massdataset 'mass_dataset'
#' @import DoddLabDatabase
#' @export
#' @example
#' check_qc_data(path = '~/Project/00_IBD_project/Data/20230808_data_quality_check/data_quality_check/',
#'               column = 'c18',
#'               polarity = 'positive')
#
# check_qc_data(path = '~/Project/00_IBD_project/Data/20230808_data_quality_check/data_quality_check/',
#               column = 'c18',
#               polarity = 'positive')

check_qc_data <- function(path = '.',
                          column = c('c18', 'hilic'),
                          polarity = c('positive', 'negative')) {

  # browser()
  column <- match.arg(column)
  polarity <- match.arg(polarity)

  if ('object.RData' %in% list.files(path = path)) {
    load(file.path(path, 'object.RData'), envir = environment())
  } else {
    # process the qc data
    object <- process_qc_data(path = path,
                              column = column,
                              polarity = polarity)
  }


  # get summary of qc data
  object <- object %>% massdataset::mutate_mean_intensity()
  object <- object %>% massdataset::mutate_rsd()
  summary_table_qc <- get_summary_object(object) %>%
    sjmisc::rotate_df() %>%
    tibble::rownames_to_column(var = 'parameter') %>%
    dplyr::rename(value = V1)

  # get ISTDs in qc sample
  summary_table_istd <- get_summary_istd(object = object,
                                         column = column,
                                         polarity = polarity)

  # get biological metabolite in qc sample
  summary_table_biomet <- get_summary_biomet(object = object,
                                             column = column,
                                             polarity = polarity)

  summary_table <- list('poolQC' = summary_table_qc, 'ISTD' = summary_table_istd, 'BiologicalMet' = summary_table_biomet)
  writexl::write_xlsx(summary_table, path = file.path(path, 'Quality_Check_Summary.xlsx'), format_headers = FALSE)

  # generate plots
  cat(crayon::blue('Export TIC & EIC plots...\n'))
  data_files <- list.files(path = path, pattern = 'mzML|mzxml', recursive = TRUE)
  idx <- which(dirname(data_files) != 'ms2')
  msdata <- RaMS::grabMSdata(file.path(path, data_files[idx]))

  # TIC plot
  temp_tic_plot <- plot_tic(msdata = msdata,
                            mode = 'both')
  # ISTD EIC plot
  temp_istd_plot <- lapply(seq_along(summary_table_istd$name_istd), function(i){
    plot_eic(msdata = msdata,
             target_mz = summary_table_istd$mz_istd[i],
             target_metabolite = summary_table_istd$name_istd[i],
             mz_tol = 10,
             mode = 'both')
  })
  temp_istd_plot <- temp_istd_plot %>% do.call(c, .)


  # biomet EIC plot
  temp_biomet_plot <- lapply(seq_along(summary_table_biomet$name_biomet), function(i){
    plot_eic(msdata = msdata,
             target_mz = summary_table_biomet$mz_biomet[i],
             target_metabolite = summary_table_biomet$name_biomet[i],
             mz_tol = 10,
             mode = 'both')
  })
  temp_biomet_plot <- temp_biomet_plot %>% do.call(c, .)

  # RSD plot
  temp_rsd_plot1 <- massqc::massqc_rsd_plot(object = object)
  temp_rsd_plot2 <- massqc::massqc_cumulative_rsd_plot(object = object, rsd_cutoff = 30, color = "black")
  # temp_rsd_plot <- patchwork::wrap_plots(temp_rsd_plot1, temp_rsd_plot2, nrow = 1)
  temp_rsd_plot <- list(temp_rsd_plot1, temp_rsd_plot2)


  # ms2 plot
  load(file.path(path, '/01_metabolite_annotation/00_intermediate_data/ms2_result'), envir = environment())
  load(file.path(path, '/01_metabolite_annotation/00_intermediate_data/annot_table'), envir = environment())

  temp_ms2_plot <- lapply(seq_along(summary_table_biomet$name_biomet), function(i){
    cat(i, ' ')

    if (is.na(summary_table_biomet$variable_id[i])) {
      return(NULL)
    }

    plot_ms2_biomet(variable_id = summary_table_biomet$variable_id[i],
                    metabolite_name = summary_table_biomet$name_biomet[i],
                    annot_table = annot_table,
                    ms2_result = ms2_result)
  })

  # export plots
  list_plot <- c(temp_tic_plot, temp_rsd_plot, temp_istd_plot, temp_biomet_plot, temp_ms2_plot)

  ggplot2::ggsave(
    filename = file.path(path, 'Quality_Check_Plot.pdf'),
    plot = gridExtra::marrangeGrob(list_plot, nrow = 1, ncol = 2),
    width = 15, height = 9
  )


  cat(crayon::blue('Quality Check Job: Done...\n\n'))

  cat(crayon::silver('--------------------------------------------------\n'))
  cat(crayon::green('Brief Summary: \n\n'))
  cat(crayon::silver('--------------------------------------------------\n'))
  cat(crayon::underline('PoolQC sample result:\n'))
  print(summary_table_qc %>% knitr::kable())
  cat('**********************\n\n')
  cat(crayon::underline('ISTD result:\n'))
  print(summary_table_istd %>% knitr::kable())
  cat('**********************\n\n')
  cat(crayon::underline('Biological metabolite result:\n'))
  print(summary_table_biomet %>% knitr::kable())
}



################################################################################
# process_qc_data --------------------------------------------------------------

#' @title process_qc_data
#' @author Zhiwei Zhou
#' @param path '.'
#' @param column 'c18', 'hilic'. Default: 'c18'
#' @param polarity 'positive', 'negative'. Default: "positive"
#' @importFrom magrittr %>%
#' @importFrom data.table %between%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importClassesFrom SpectraTools 'SpectraData'
#' @importClassesFrom massdataset 'ms2_data'
#' @importClassesFrom massdataset 'mass_dataset'
#' @import DoddLabDatabase
#' @export
#' @example
#' process_qc_data(path = '~/Project/00_IBD_project/Data/20230808_data_quality_check/data_quality_check/',
#'                column = 'c18',
#'                polarity = 'positive')

# process_qc_data(path = '~/Project/00_IBD_project/Data/20230808_data_quality_check/data_quality_check/',
#                column = 'c18',
#                polarity = 'positive')


process_qc_data <- function(path = '.',
                            column = c('c18', 'hilic'),
                            polarity = c('positive', 'negative')){
  column <- match.arg(column)
  polarity <- match.arg(polarity)

  # raw data procssing
  cat(crayon::blue('Peak picking for QC samples ...\n'))
  parameter_set <- DoddLabRawMS::initialize_raw_parameter_class(column = column)
  parameter_set@para_peak_grouping$bw <- 2

  DoddLabRawMS::process_raw_data(parameter_set = parameter_set,
                                 path = path)

  cat(crayon::blue('Extract and Align MS2 for QC samples ...\n'))
  require(DoddLabDatabase)
  parameter_set_annotation <- DoddLabMetID::initialize_annotation_parameter_class(path = path,
                                                                                  lib = 'dodd',
                                                                                  ce = '20',
                                                                                  column = column,
                                                                                  polarity = polarity,
                                                                                  is_rt_score = TRUE,
                                                                                  is_ms2_score = TRUE)

  DoddLabMetID::annotate_metabolite(parameter_set_annotation = parameter_set_annotation)

  cat(crayon::blue('Creat tidymass object for following analysis\n'))

  peak_table <- readr::read_csv(file.path(path, 'data.csv'), show_col_types = FALSE)

  # sample name
  sample_id <- colnames(peak_table)[-c(1:3)]
  sample_info <- tibble::tibble(sample_id = sample_id,
                                injection.order = seq(length(sample_id)),
                                class = 'Subject',
                                group = 'PoolQC')

  # expression data
  expression_data <- peak_table %>%
    dplyr::select(-c(name:rt)) %>%
    as.data.frame()

  # variable info
  variable_info <- peak_table %>%
    dplyr::rename(variable_id = name) %>%
    dplyr::select(variable_id:rt) %>%
    as.data.frame()

  rownames(expression_data) <- variable_info$variable_id

  # creat mass_data object
  object <- massdataset::create_mass_dataset(expression_data = expression_data,
                                             sample_info = sample_info,
                                             variable_info = variable_info)


  # load ms2 data
  load(file.path(path, '01_metabolite_annotation/00_intermediate_data/ms2_data_combined'))

  temp_ms2_data <-
    new(
      Class = "ms2_data",
      column = column,
      polarity = parameter_set_annotation@para_general$polarity,
      variable_id = names(ms2_data_combined$spec),
      ms2_spectrum_id = names(ms2_data_combined$spec),
      ms2_mz = ms2_data_combined$info$mz,
      ms2_rt = ms2_data_combined$info$rt,
      ms2_file = ms2_data_combined$info$filename,
      ms2_spectra = ms2_data_combined$spec,
      mz_tol = parameter_set_annotation@para_ms2_match$mz_tol_combine_ms1_ms2,
      rt_tol = parameter_set_annotation@para_ms2_match$rt_tol_combine_ms1_ms2
    )

  temp_ms2_data = list(name = temp_ms2_data)
  temp_name <- paste(sort(unique(basename(ms2_data_combined$info$filename))), collapse = ";")
  names(temp_ms2_data) <- temp_name
  object@ms2_data <- temp_ms2_data
  save(object, file = file.path(path, 'object.RData'))

  cat(crayon::blue('Job Done!\n'))

  return(object)
}



################################################################################
# get_summary_object -----------------------------------------------------------

#' @title get_summary_object
#' @author Zhiwei Zhou
#' @param object mass
#' @export

# object <- object %>% massdataset::mutate_mean_intensity()
# object <- object %>% massdataset::mutate_rsd()
# get_summary_object(object)

get_summary_object <- function(object){
  num_sample <- massdataset::get_sample_number(object)
  num_feature <- massdataset::get_variable_number(object)
  num_ms2 <- object@ms2_data[[1]]@ms2_spectra %>% length()
  ms2_percentage <- round(num_ms2/num_feature*100, 2)

  if (!('rsd' %in% colnames(object@variable_info))) {
    object <- object %>%
      massdataset::mutate_rsd()
  }

  rsd_percentage <- round(sum(object@variable_info$rsd <= 30)/num_feature * 100, 2)

  result_table <- tibble::tibble(num_sample = num_sample,
                                 num_feature = num_feature,
                                 num_ms2 = num_ms2,
                                 ms2_percentage = ms2_percentage,
                                 rsd_percentage = rsd_percentage)

  return(result_table)
}



################################################################################
# get_summary_istd -------------------------------------------------------------

# get_summary_istd(object = object, column = 'c18', polarity = 'positive')

get_summary_istd <- function(object,
                             column,
                             polarity,
                             mz_ppm = 10,
                             rt_second = 25) {

  # browser()

  switch (column,
    'c18' = {
      mode <- ifelse(polarity == 'positive', 'c18_pos', 'c18_neg')
    },
    'hilic' = {
      mode <- ifelse(polarity == 'positive', 'hilic_pos', 'hilic_neg')
    }
  )

  cat(crayon::blue('Extract internal standard result...\n'))

  require(DoddLabTool)
  istd_result <- DoddLabTool::match_istd(object = object,
                                         mode = mode,
                                         mz_ppm = mz_ppm,
                                         rt_second = rt_second)


  temp_istd_object <- object %>%
    massdataset::activate_mass_dataset(what = 'variable_info') %>%
    massdataset::filter(variable_id %in% istd_result$variable_id) %>%
    massdataset::arrange(match(variable_id, istd_result$variable_id))

  temp_istd_int <- apply(temp_istd_object@expression_data, 1, mean)
  temp_rsd <- apply(temp_istd_object@expression_data, 1, function(x){
    round(sd(x)/mean(x)*100, 2)
  })
  idx <- match(istd_result$variable_id, names(temp_istd_int))
  istd_result <- istd_result %>%
    mutate(int = temp_istd_int[idx],
           rsd = temp_rsd[idx])

  return(istd_result)

  cat(crayon::blue('Done...\n'))
}



# get_summary_biomet -------------------------------------------------------------

#' @title get_summary_biomet
#' @author Zhiwei Zhou
#' @param object tidymass object
#' @param column 'hilic',  'c18'
#' @param polarity "positive", 'negative'
#' @param mz_ppm 10 ppm in defalut
#' @param rt_second 25 s in defalut
#' @export
#' @examples
#' load('~/Project/00_IBD_project/Data/20230808_data_quality_check/B0002/hilic_pos/object.RData')
#' get_summary_biomet(object = object, polarity = 'positive', column = 'hilic')

# get_summary_biomet(object = object, column = 'c18', polarity = 'positive')
# load('~/Project/00_IBD_project/Data/20230808_data_quality_check/B0002/hilic_pos/object.RData')
# get_summary_biomet(object = object, polarity = 'positive', column = 'hilic')

get_summary_biomet <- function(object,
                               column,
                               polarity,
                               mz_ppm = 10,
                               rt_second = 25) {

  # browser()

  switch (column,
          'c18' = {
            mode <- ifelse(polarity == 'positive', 'c18_pos', 'c18_neg')
          },
          'hilic' = {
            mode <- ifelse(polarity == 'positive', 'hilic_pos', 'hilic_neg')
          }
  )

  cat(crayon::blue('Extract real metabolite result...\n'))

  require(DoddLabTool)
  biomet_result <- DoddLabTool::match_biomet(object = object,
                                             mode = mode,
                                             mz_ppm = mz_ppm,
                                             rt_second = rt_second)


  temp_biomet_object <- object %>%
    massdataset::activate_mass_dataset(what = 'variable_info') %>%
    massdataset::filter(variable_id %in% biomet_result$variable_id) %>%
    massdataset::arrange(match(variable_id, biomet_result$variable_id))

  temp_biomet_int <- apply(temp_biomet_object@expression_data, 1, mean)
  temp_rsd <- apply(temp_biomet_object@expression_data, 1, function(x){
    round(sd(x)/mean(x)*100, 2)
  })
  idx <- match(biomet_result$variable_id, names(temp_biomet_int))
  biomet_result <- biomet_result %>%
    mutate(int = temp_biomet_int[idx],
           rsd = temp_rsd[idx])

  return(biomet_result)

  cat(crayon::blue('Done...\n'))
}




################################################################################
# plot_tic ---------------------------------------------------------------------

plot_tic <- function(msdata = NULL,
                     data_files,
                     mode = c('both', 'overlapped', 'separated')){
  if (is.null(msdata)) {
    msdata <- RaMS::grabMSdata(data_files, grab_what = 'TIC')
  }

  # # keep the segment 1
  # temp_msdata <- lapply(unique(msdata$TIC$filename), function(x){
  #   temp_data <- msdata$TIC %>%
  #     dplyr::filter(filename == x)
  #
  #   temp_data %>% dplyr::slice(seq(1, nrow(temp_data), 3))
  # }) %>% dplyr::bind_rows()
  if (mode == 'overlapped') {
    temp_plot <- ggplot2::ggplot(msdata$TIC) +
      ggplot2::geom_line(ggplot2::aes(x=rt, y=int, color=filename), lwd = 1) +
      ggplot2::xlab('Retention time (minute)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle('TIC') +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(axis.text.y = element_text(hjust = 0.5, angle = 0),
                     legend.position = 'none')

    return(temp_plot)
  }

  if (mode == 'separated') {
    temp_plot <- ggplot2::ggplot(msdata$TIC) +
      ggplot2::geom_line(aes(x=rt, y=int, color=filename)) +
      ggplot2::facet_wrap(facets = 'filename', ncol = 1) +
      ggplot2::scale_y_continuous(n.breaks = 3) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle('TIC') +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(axis.text.y = element_text(hjust = 0.5, angle = 0),
                     legend.position = 'none')

    return(temp_plot)
  }

  if (mode == 'both') {
    temp_plot1 <- ggplot2::ggplot(msdata$TIC) +
      ggplot2::geom_line(ggplot2::aes(x=rt, y=int, color=filename), lwd = 1) +
      ggplot2::xlab('Retention time (minute)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle('TIC') +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(axis.text.y = element_text(hjust = 0.5, angle = 0),
                     legend.position = 'none')

    temp_plot2 <- ggplot2::ggplot(msdata$TIC) +
      ggplot2::geom_line(aes(x=rt, y=int, color=filename)) +
      ggplot2::facet_wrap(facets = 'filename', ncol = 1) +
      ggplot2::scale_y_continuous(n.breaks = 3) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      # ggplot2::ggtitle('TIC') +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(axis.text.y = element_text(hjust = 0.5, angle = 0),
                     legend.position = 'none')

    # temp_plot <- patchwork::wrap_plots(temp_plot1, temp_plot2, nrow = 1)
    temp_plot <- list(temp_plot1, temp_plot2)

    return(temp_plot)
  }

}



################################################################################
# plot_eic ---------------------------------------------------------------------

# target_mz <- 171.1177

plot_eic <- function(msdata = NULL,
                    data_files,
                    target_mz,
                    target_metabolite = '', # title
                    mz_tol = 10, # ppm
                    mode = c('both', 'overlapped', 'separated')
){
  # browser()
  if (is.null(msdata)) {
    msdata <- RaMS::grabMSdata(data_files, grab_what = 'MS1')
  }
  require(data.table)
  require(RaMS)

  mode <- match.arg(mode)

  # target_data <- msdata$MS1[mz %between% RaMS::pmppm(target_mz, ppm = mz_tol)]

  temp_mz_tol <- RaMS::pmppm(target_mz, ppm = mz_tol)
  temp_idx <- msdata$MS1$mz %>% RaMS::between(lower = temp_mz_tol[1], upper = temp_mz_tol[2]) %>% which()
  target_data <- msdata$MS1[temp_idx,]

  if (mode == 'overlapped') {
    temp_plot <- ggplot2::ggplot(target_data) +
      ggplot2::geom_line(ggplot2::aes(x=rt, y=int, color=filename)) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle(target_metabolite) +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(legend.position = c(0.8, 0.8))

    return(temp_plot)
  }

  if (mode == 'separated') {
    temp_plot <- ggplot2::ggplot(target_data, ggplot2::aes(x=rt, y=int, color=filename)) +
      ggplot2::geom_line() +
      facet_wrap(facets = 'filename', ncol = 1) +
      scale_y_continuous(n.breaks = 3) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle(target_metabolite) +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(legend.position = 'none')

    return(temp_plot)
  }

  if (mode == 'both') {
    temp_plot1 <- ggplot2::ggplot(target_data) +
      ggplot2::geom_line(ggplot2::aes(x=rt, y=int, color=filename)) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle(target_metabolite) +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(legend.position = c(0.8, 0.8))

    temp_plot2 <- ggplot2::ggplot(target_data, ggplot2::aes(x=rt, y=int, color=filename)) +
      ggplot2::geom_line() +
      facet_wrap(facets = 'filename', ncol = 1) +
      scale_y_continuous(n.breaks = 3) +
      ggplot2::xlab('RT (min)') +
      ggplot2::ylab('Intensity') +
      ggplot2::ggtitle(target_metabolite) +
      ZZWtool::ZZWTheme() +
      ggplot2::theme(legend.position = 'none')

    # temp_plot <- patchwork::wrap_plots(temp_plot1, temp_plot2, nrow = 1)
    temp_plot <- list(temp_plot1, temp_plot2)

    return(temp_plot)
  }
}

#
# temp_plot <- lapply(istd_peak$variable_id, function(x){
#   object_c18_pos %>%
#     intensity_plot(variable_id = x,
#                    color_by = 'class',
#                    order_by = 'injection.order',
#                    interactive = FALSE) +
#     scale_x_discrete(labels = NULL) +
#     ZZWtool::ZZWTheme() +
#     theme(axis.text.y.left = element_text(hjust = 0.5)) +
#     ggtitle(paste(x, 'Before intra-batch normalization')) +
#     geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed')
# })
#
# ggsave(
#   filename = "./02_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization.pdf",
#   plot = gridExtra::marrangeGrob(temp_plot, nrow=1, ncol=1),
#   width = 15, height = 9
# )

################################################################################
# plot_ms2_biomet --------------------------------------------------------------

#' @title plot_ms2_biomet
#' @author Zhiwei Zhou
#' @param variable_id
#' @param metabolite_name
#' @param annot_table intermidate data from DoddLabMetID
#' @param ms2_result intermidate data from DoddLabMetID
#' @export
#' @examples
#' variable_id <-  'M166T279'
#' metabolite_name <- 'L-Phenylalanine'
#' plot_ms2_biomet(variable_id = 'M166T279',
#'                 metabolite_name = 'L-Phenylalanine',
#'                 annot_table = annot_table,
#'                 ms2_result = ms2_result)


# variable_id <-  'M166T279'
# metabolite_name <- 'L-Phenylalanine'
# plot_ms2_biomet(variable_id = 'M166T279',
#                 metabolite_name = 'L-Phenylalanine',
#                 annot_table = annot_table,
#                 ms2_result = ms2_result)

plot_ms2_biomet <- function(variable_id,
                            metabolite_name,
                            annot_table,
                            ms2_result){
  temp_result <- annot_table %>%
    dplyr::filter(feature_name == variable_id) %>%
    dplyr::filter(name == metabolite_name) %>%
    dplyr::slice(1)

  if (nrow(temp_result) < 1) {
    return(NULL)
  }

  temp_feature_name <- temp_result$feature_name
  temp_feature_mz <- temp_result$mz
  temp_feature_rt <- temp_result$rt
  temp_cpd_id <- temp_result$id
  temp_cpd_name <- temp_result$name
  temp_adduct <- temp_result$adduct
  temp_mz_error <- temp_result$mz_error
  temp_rt_error <- temp_result$rt_error
  temp_ms2_score_forward <- temp_result$msms_score_forward
  temp_ms2_score_reverse <- temp_result$msms_score_reverse
  temp_ms2_match_fragment <- temp_result$msms_matched_frag
  temp_mz_lib <- temp_result$mz_lib
  temp_rt_lib <- temp_result$rt_lib
  temp_smiles <- temp_result$smiles
  temp_inchikey <- temp_result$inchikey

  temp_ms2_obj <- which(names(ms2_result) == temp_feature_name) %>%
    ms2_result[[.]]
  temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
    temp_ms2_obj@matchedFragments[[.]]

  text <- paste(c(
    paste0('Feature: ', temp_feature_name),
    paste0('Feature m/z: ', round(temp_feature_mz, 4)),
    paste0('Feature RT: ', round(temp_feature_rt, 1)),
    paste0('Compound ID: ', temp_cpd_id),
    paste0('Compound name: ', temp_cpd_name),
    paste0('Adduct: ', temp_adduct),
    paste0('SMILES: ', temp_smiles),
    paste0('InChIKey: ', temp_inchikey),
    paste0('m/z lib: ', round(temp_mz_lib, 4)),
    paste0('m/z error: ', temp_mz_error),
    paste0('RT lib: ', temp_rt_lib),
    paste0('RT error: ', temp_rt_error),
    paste0('MS2 score forward: ', temp_ms2_score_forward),
    paste0('MS2 score reverse: ', temp_ms2_score_reverse),
    paste0('MS2 matched frag.: ', temp_ms2_match_fragment)
  ),
  collapse = '\n')

  suppressMessages(
    temp_plot <- DoddLabMetID::plot_id_ms2(obj_spec = temp_ms2_obj) +
      ggplot2::scale_colour_manual(
        name = 'Attribute',
        labels= c(paste0('Experiment'),
                  'Unmatched fragments',
                  paste0('Library ')),
        values = c(
          'experiment' = 'black',
          'library' = 'red',
          'frag_unmatch' = 'gray'
        )
      ) +
      ggplot2::scale_shape_manual(
        name = 'Label',
        labels= c('matched' = "Matched",
                  'unmatched' = "Unmatched"),
        values = c(
          'matched' = 16,
          'unmatched' = 4
        )
      ) +
      ggplot2::ggtitle(label = paste0(temp_feature_name,
                                      ': ', temp_cpd_name,
                                      ' (ID: ', temp_cpd_id,
                                      ')')) +
      ZZWtool::ZZW_annotate_text2(label = text, x = 0, y = 1) +
      ggplot2::theme(legend.position = c(0.85, 0.85),
                     title = ggplot2::element_text(vjust = 0.5))
  )

  return(temp_plot)


}
