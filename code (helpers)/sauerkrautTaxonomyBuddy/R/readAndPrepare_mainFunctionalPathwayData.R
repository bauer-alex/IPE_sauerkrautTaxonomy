
#' Read and prepare all data of the functional pathway stool measurements
#' 
#' This function wraps the functions for pathway, blood marker and metabolome
#' data preparation and the main experiment data preparation and creates one
#' \code{tse} experiment data object for the main measurements in the experiment.
#' 
#' @param data_type One of \code{c("stool","kraut")}, indicating which samples
#' should be prepared and returned. Defaults to \code{"stool"}.
#' @inheritParams readAndPrepare_bloodMarkerData
#' @inheritParams readAndPrepare_experimentData
#' @inheritParams readAndPrepare_functionalData
#' 
#' @import checkmate dplyr mia TreeSummarizedExperiment
#' @importFrom vegan vegdist
#' @export
#' 
#' @return \code{tse} data object
#' 
readAndPrepare_mainFunctionalPathwayData <- function(path_KEGG_KO,
                                                     path_KEGG_pathways,
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
                                                     data_type                     = "stool",
                                                     exclude_medicatedParticipants = TRUE,
                                                     exclude_sickParticipants) {
  
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
  checkmate::assert_file(path_KEGG_KO,                      extension = "txt")
  checkmate::assert_file(path_KEGG_pathways,                extension = "txt")
  checkmate::assert_file(path_sampleLookupXlsx,             extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,              extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerBostonRdata,       extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerZentrallaborRdata, extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerBostonXlsx,        extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerZentrallaborXlsx,  extension = "xlsx")
  checkmate::assert_choice(data_type, choices = c("stool","kraut"))
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  
  # read all data -----------------------------------------------------------
  # basic experiment data
  if (data_type == "stool") {
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
  } else { # data_type == "kraut"
    
    dat_lookupList <- readAndPrepare_lookupData(path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                                path_krautLookupXlsx          = path_krautLookupXlsx,
                                                exclude_medicatedParticipants = exclude_medicatedParticipants,
                                                exclude_sickParticipants      = exclude_sickParticipants)
    dat_lookupKraut <- dat_lookupList$dat_lookupKraut
  }
  
  # absolute and relative abundance functional analysis data
  datList_tpmFun <- readAndPrepare_functionalData(path_KEGG_KO             = path_KEGG_KO,
                                                  path_KEGG_pathways       = path_KEGG_pathways,
                                                  path_sampleLookupXlsx    = path_sampleLookupXlsx,
                                                  path_krautLookupXlsx     = path_krautLookupXlsx,
                                                  abundance_unit           = "TPM",
                                                  exclude_sickParticipants = exclude_sickParticipants)
  datList_relFun <- readAndPrepare_functionalData(path_KEGG_KO             = path_KEGG_KO,
                                                  path_KEGG_pathways       = path_KEGG_pathways,
                                                  path_sampleLookupXlsx    = path_sampleLookupXlsx,
                                                  path_krautLookupXlsx     = path_krautLookupXlsx,
                                                  abundance_unit           = "relAbd",
                                                  exclude_sickParticipants = exclude_sickParticipants)
  dat_tpmFun <- datList_tpmFun[[paste0("datPathways_", data_type, "Samples")]]
  dat_relFun <- datList_relFun[[paste0("datPathways_", data_type, "Samples")]]
  
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
  matrix_tpmAbd <- dat_tpmFun %>% select(starts_with("sample")) %>% as.matrix()
  matrix_relAbd <- dat_relFun %>% select(starts_with("sample")) %>% as.matrix()
  
  # prepare colData, i.e. information on the individual samples
  if (data_type == "stool") {
    dat_blood <- dat_blood %>% select(-participant_id, -timepoint)
    dat_col   <- dat_exp %>% 
      arrange(sample_id) %>% 
      mutate(sample_id = paste0("sample_", sample_id)) %>% 
      as.data.frame() %>% 
      dplyr::left_join(dat_blood, by = "part_time")
  
  } else { # data_type == "kraut"
    dat_col <- dat_lookupKraut %>% mutate(sample_id = paste0("sample_", sample_id)) %>% as.data.frame()
  }
  row.names(dat_col) <- dat_col$sample_id
  dat_col            <- dat_col %>% select(-sample_id)
  
  # prepare rowData, i.e. information on the individual functions
  dat_row           <- dat_tpmFun %>% 
    select(otu_id, starts_with("taxonomy")) %>% 
    as.data.frame()
  rownames(dat_row) <- dat_row$otu_id
  dat_row           <- dat_row %>% select(-otu_id)
  
  
  
  # create a tse data object ------------------------------------------------
  tse <- TreeSummarizedExperiment(assays  = list(tpmAbd = matrix_tpmAbd,
                                                 relAbd = matrix_relAbd),
                                  colData = dat_col,
                                  rowData = dat_row)
  
  return(tse)
}
