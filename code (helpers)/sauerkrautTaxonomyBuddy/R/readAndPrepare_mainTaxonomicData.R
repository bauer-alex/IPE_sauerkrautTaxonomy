
#' Read and prepare all data of the taxonomic stool measurements
#' 
#' This function wraps the functions for taxonomy, blood marker and metabolome
#' data preparation and the main experiment data preparation and creates one
#' \code{tse} experiment data object for the main measurements in the experiment.
#' Additionally, it creates a taxonomic tree, calculates clr transformed relative
#' abundance values and calculates different statistics (richness, evenness, etc.)
#' on the stool samples and stores them in the \code{tse} object.
#' 
#' @param path_absTaxonomyData,path_relTaxonomyData Path to the absolute and
#' relative abundance taxonomy files \code{merged_abundance_table_counts.tsv}
#' and \code{merged_abundance_table_relab.tsv}
#' @inheritParams readAndPrepare_bloodMarkerData
#' @inheritParams readAndPrepare_experimentData
#' @inheritParams readAndPrepare_taxonomy
#' 
#' @import checkmate dplyr mia TreeSummarizedExperiment
#' @importFrom vegan vegdist
#' @export
#' 
#' @return \code{tse} data object
#' 
readAndPrepare_mainTaxonomicData <- function(path_absTaxonomyData,
                                             path_relTaxonomyData,
                                             path_groupInfoRdata,
                                             path_participantRdata,
                                             path_samplesIDxlsx,
                                             path_dietInfoCsv,
                                             path_nutrientInfoCsv,
                                             path_bloodMetabolomeXlsx,
                                             path_stoolMetabolomeXlsx,
                                             path_T4questXlsx,
                                             path_stoolInfoXlsx,
                                             path_bodyMeasuresXlsx,
                                             path_sampleLookupXlsx,
                                             path_krautLookupXlsx,
                                             path_bloodMarkerBostonRdata,
                                             path_bloodMarkerZentrallaborRdata,
                                             path_bloodMarkerBostonXlsx,
                                             path_bloodMarkerZentrallaborXlsx,
                                             aggregation_level             = "species",
                                             remove_unobservedBacteria     = TRUE,
                                             exclude_rareStoolBacteria     = TRUE,
                                             exclude_medicatedParticipants = TRUE,
                                             exclude_sickParticipants) {
  
  checkmate::assert_file(path_absTaxonomyData,              extension = "tsv")
  checkmate::assert_file(path_relTaxonomyData,              extension = "tsv")
  checkmate::assert_file(path_groupInfoRdata,               extension = "Rdata")
  checkmate::assert_file(path_participantRdata,             extension = "Rdata")
  checkmate::assert_file(path_samplesIDxlsx,                extension = "xlsx")
  checkmate::assert_file(path_dietInfoCsv,                  extension = "csv")
  checkmate::assert_file(path_nutrientInfoCsv,              extension = "csv")
  checkmate::assert_file(path_bloodMetabolomeXlsx,          extension = "XLSX")
  checkmate::assert_file(path_stoolMetabolomeXlsx,          extension = "XLSX")
  checkmate::assert_file(path_T4questXlsx,                  extension = "xlsx")
  checkmate::assert_file(path_stoolInfoXlsx,                extension = "xlsx")
  checkmate::assert_file(path_bodyMeasuresXlsx,             extension = "xlsx")
  checkmate::assert_file(path_sampleLookupXlsx,             extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,              extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerBostonRdata,       extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerZentrallaborRdata, extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerBostonXlsx,        extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerZentrallaborXlsx,  extension = "xlsx")
  checkmate::assert_choice(aggregation_level, choices = c("species","strain"))
  checkmate::assert_logical(remove_unobservedBacteria,     len = 1)
  checkmate::assert_logical(exclude_rareStoolBacteria,     len = 1)
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  
  # read all data -----------------------------------------------------------
  # basic experiment data
  dat_exp <- readAndPrepare_experimentData(path_groupInfoRdata           = path_groupInfoRdata,
                                           path_participantRdata         = path_participantRdata,
                                           path_samplesIDxlsx            = path_samplesIDxlsx,
                                           path_dietInfoCsv              = path_dietInfoCsv,
                                           path_nutrientInfoCsv          = path_nutrientInfoCsv,
                                           path_bloodMetabolomeXlsx      = path_bloodMetabolomeXlsx,
                                           path_stoolMetabolomeXlsx      = path_stoolMetabolomeXlsx,
                                           path_T4questXlsx              = path_T4questXlsx,
                                           path_stoolInfoXlsx            = path_stoolInfoXlsx,
                                           path_bodyMeasuresXlsx         = path_bodyMeasuresXlsx,
                                           path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                           path_krautLookupXlsx          = path_krautLookupXlsx,
                                           exclude_medicatedParticipants = exclude_medicatedParticipants,
                                           exclude_sickParticipants      = exclude_sickParticipants)
  
  # absolute and relative abundance taxonomy data
  # Note: The 'exclude_rareStoolBacteria' argument is only implemented for relative abundance data.
  #       Accordingly, the absolute abundance data is filtered after applying the 'readAndPrepare_taxonomy()' function.
  datList_absTax <- readAndPrepare_taxonomy(path_taxonomyData             = path_absTaxonomyData,
                                            path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                            path_krautLookupXlsx          = path_krautLookupXlsx,
                                            aggregation_level             = aggregation_level,
                                            remove_unobservedBacteria     = remove_unobservedBacteria,
                                            exclude_rareStoolBacteria     = FALSE, # functionality only available for relative abundance data
                                            exclude_medicatedParticipants = exclude_medicatedParticipants,
                                            exclude_sickParticipants      = exclude_sickParticipants)
  datList_relTax <- readAndPrepare_taxonomy(path_taxonomyData             = path_relTaxonomyData,
                                            path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                            path_krautLookupXlsx          = path_krautLookupXlsx,
                                            aggregation_level             = aggregation_level,
                                            remove_unobservedBacteria     = remove_unobservedBacteria,
                                            exclude_rareStoolBacteria     = exclude_rareStoolBacteria,
                                            exclude_medicatedParticipants = exclude_medicatedParticipants,
                                            exclude_sickParticipants      = exclude_sickParticipants)
  dat_absTax    <- datList_absTax$dat_stoolSamples
  dat_relTax    <- datList_relTax$dat_stoolSamples
  dat_absTax    <- dat_absTax[which(dat_absTax[[aggregation_level]] %in% dat_relTax[[aggregation_level]]),]
  dat_lookupTax <- datList_absTax$dat_lookupStool
  
  # blood marker data
  datList_blood <- readAndPrepare_bloodMarkerData(path_bloodMarkerBostonRdata       = path_bloodMarkerBostonRdata,
                                                  path_bloodMarkerZentrallaborRdata = path_bloodMarkerZentrallaborRdata,
                                                  path_bloodMarkerBostonXlsx        = path_bloodMarkerBostonXlsx,
                                                  path_bloodMarkerZentrallaborXlsx  = path_bloodMarkerZentrallaborXlsx,
                                                  path_sampleLookupXlsx             = path_sampleLookupXlsx,
                                                  path_krautLookupXlsx              = path_krautLookupXlsx,
                                                  exclude_medicatedParticipants     = exclude_medicatedParticipants,
                                                  exclude_sickParticipants          = exclude_sickParticipants)
  dat_blood <- datList_blood$main_data
  
  
  # properly format the data objects ----------------------------------------
  # prepare the assays, i.e. matrices for absolute and relative abundances
  matrix_absAbd           <- dat_absTax %>% select(starts_with("sample")) %>% as.matrix()
  rownames(matrix_absAbd) <- dat_absTax[[aggregation_level]]
  matrix_relAbd           <- dat_relTax %>% select(starts_with("sample")) %>% as.matrix()
  rownames(matrix_relAbd) <- dat_relTax[[aggregation_level]]
  
  # prepare colData, i.e. information on the individual samples
  dat_blood <- dat_blood %>% select(-participant_id, -timepoint)
  dat_col <- dat_exp %>% 
    arrange(sample_id) %>% 
    mutate(sample_id = paste0("sample_", sample_id)) %>% 
    as.data.frame() %>% 
    dplyr::left_join(dat_blood, by = "part_time")
  row.names(dat_col) <- dat_col$sample_id
  dat_col            <- dat_col %>% select(-sample_id)
  
  # prepare rowData, i.e. (taxonomic) information on the species / bacteria strains
  dat_row <- dat_absTax %>% 
    select(starts_with(aggregation_level)) %>% 
    as.data.frame()
  rownames(dat_row) <- dat_row[[aggregation_level]]
  colnames(dat_row) <- gsub(pattern     = paste0(aggregation_level, "_"),
                            replacement = "",
                            colnames(dat_row))
  
  
  
  # create a tse data object ------------------------------------------------
  tse <- TreeSummarizedExperiment(assays  = list(absAbd = matrix_absAbd,
                                                 relAbd = matrix_relAbd),
                                  colData = dat_col,
                                  rowData = dat_row)
  
  
  
  # add a taxonomic tree ----------------------------------------------------
  tse <- tse %>% mia::addTaxonomyTree()
  
  
  
  # add clr transformation for abundance values -----------------------------
  relAbd_mat  <- assay(tse, "relAbd") 
  pseudocount <- relAbd_mat[relAbd_mat > 0] %>% min()
  
  tse <- tse %>% transformAssay(assay.type  = "relAbd",
                                method      = "clr", 
                                pseudocount = pseudocount,
                                name        = "clr")
  
  
  
  # estimate some taxonomy variation statistics -----------------------------
  tse <- tse %>% mia::estimateRichness(  assay.type = "relAbd", index = c("hill", "observed"),              name = c("richness_hill", "richness_observed"))
  tse <- tse %>% mia::estimateEvenness(  assay.type = "relAbd", index = c("pielou", "simpson_evenness"),    name = c("evenness_pielou", "evenness_simpson"))
  tse <- tse %>% mia::estimateDiversity( assay.type = "relAbd", index = c("shannon", "inverse_simpson"),    name = c("diversity_shannon", "diversity_invSimpson"))
  tse <- tse %>% mia::estimateDominance( assay.type = "relAbd", index = c("dbp", "core_abundance"),         name = c("dominance_dbp", "dominance_coreAbundance"))
  tse <- tse %>% mia::estimateDiversity( assay.type = "relAbd", index = "log_modulo_skewness",              name = "rarity_logModuloSkewness")
  tse <- tse %>% mia::estimateDivergence(assay.type = "relAbd", reference = "median", FUN = vegan::vegdist, name = "divergence_toMedian")
  
  
  return(tse)
}
