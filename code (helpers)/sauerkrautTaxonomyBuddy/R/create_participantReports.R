
#' Create one measurement summary pdf report for each participant
#' 
#' This function calls the R package internal Quarto report once for each
#' participant and stores the rendered pdf reports in the specified directory.
#' 
#' @param path_sauerkrautBuddyDir Path to the \code{sauerkrautBuddy} directory
#' as downloaded from GitHub.
#' @param path_targetDir Path to the target directory where the individual
#' participant pdf reports should be stored.
#' @inheritParams readAndPrepare_mainTaxonomicData
#' 
#' @import checkmate dplyr
#' @importFrom quarto quarto_render
#' @export
#' 
create_participantReports <- function(path_sauerkrautBuddyDir,
                                      path_targetDir,
                                      path_groupInfoRdata,
                                      path_participantRdata,
                                      path_samplesIDxlsx,
                                      path_absTaxonomyData,
                                      path_relTaxonomyData,
                                      path_sampleLookupXlsx,
                                      path_krautLookupXlsx,
                                      path_bloodMarkerBostonRdata,
                                      path_bloodMarkerZentrallaborRdata,
                                      path_bloodMetabolomeXlsx,
                                      path_stoolMetabolomeXlsx,
                                      path_bodyMeasuresXlsx) {
  
  checkmate::assert_directory(path_sauerkrautBuddyDir)
  checkmate::assert_directory(path_targetDir)
  checkmate::assert_file(path_groupInfoRdata,               extension = "Rdata")
  checkmate::assert_file(path_participantRdata,             extension = "Rdata")
  checkmate::assert_file(path_samplesIDxlsx,                extension = "xlsx")
  checkmate::assert_file(path_absTaxonomyData,              extension = "tsv")
  checkmate::assert_file(path_relTaxonomyData,              extension = "tsv")
  checkmate::assert_file(path_sampleLookupXlsx,             extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,              extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerBostonRdata,       extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerZentrallaborRdata, extension = "Rdata")
  checkmate::assert_file(path_bloodMetabolomeXlsx,          extension = "XLSX")
  checkmate::assert_file(path_stoolMetabolomeXlsx,          extension = "XLSX")
  checkmate::assert_file(path_bodyMeasuresXlsx,             extension = "xlsx")
  
  
  # read participant information
  dat <- readAndPrepare_experimentData(path_groupInfoRdata      = path_groupInfoRdata,
                                       path_participantRdata    = path_participantRdata,
                                       path_samplesIDxlsx       = path_samplesIDxlsx,
                                       path_bloodMetabolomeXlsx = path_bloodMetabolomeXlsx,
                                       path_stoolMetabolomeXlsx = path_stoolMetabolomeXlsx,
                                       path_bodyMeasuresXlsx    = path_bodyMeasuresXlsx)
  
  
  
  # render the Quarto document once for each participant --------------------
  for (id in levels(dat$participant_id)) {
    
    # create a list with all parameters to pass to the Quarto document
    params_list <- list(participant_id                    = id,
                        path_groupInfoRdata               = path_groupInfoRdata,
                        path_participantRdata             = path_participantRdata,
                        path_samplesIDxlsx                = path_samplesIDxlsx,
                        path_absTaxonomyData              = path_absTaxonomyData,
                        path_relTaxonomyData              = path_relTaxonomyData,
                        path_sampleLookupXlsx             = path_sampleLookupXlsx,
                        path_krautLookupXlsx              = path_krautLookupXlsx,
                        path_bloodMarkerBostonRdata       = path_bloodMarkerBostonRdata,
                        path_bloodMarkerZentrallaborRdata = path_bloodMarkerZentrallaborRdata,
                        path_bloodMetabolomeXlsx          = path_bloodMetabolomeXlsx,
                        path_stoolMetabolomeXlsx          = path_stoolMetabolomeXlsx,
                        path_targetXlsx                   = file.path(path_targetDir, paste0("Sauerkrautstudie - Daten Proband ",id,".xlsx")))
    
    # render the Quarto document
    target_file <- paste0("Sauerkrautstudie - Daten Proband ", id, ".pdf")
    quarto::quarto_render(input          = file.path(path_sauerkrautBuddyDir, "quarto/participant_report.qmd"),
                          output_format  = "pdf",
                          output_file    = target_file,
                          execute_params = params_list)
    
    # move the Quarto document to the target directory
    file.rename(from = target_file,
                to   = file.path(path_targetDir, target_file))
  }
}