
#' Helper function to define objects in the global environment stating the paths to all data files
#' 
#' The paths are all set to the data repository with the relative path
#' '../data/', which is the path necessary for running the main code files
#' from the code appendix for the paper.
#' 
#' @export
#' 
define_dataPaths <- function() {
  
  # main experiment data on participant level
  assign("path_groupInfoRdata",   "../data/1a - group info.Rdata", envir = .GlobalEnv)
  assign("path_participantRdata", "../data/1b - participant data.Rdata", envir = .GlobalEnv)
  assign("path_samplesIDxlsx",    "../data/1c - sample IDs.xlsx", envir = .GlobalEnv)
  assign("path_dietInfoCsv",      "../data/1d - food frequency questionnaire.csv", envir = .GlobalEnv)
  assign("path_nutrientInfoCsv",  "../data/1e - food frequency questionnaire (nutrients).csv", envir = .GlobalEnv)
  assign("path_T4questXlsx",      "../data/1f - opinion questionnaire.xlsx", envir = .GlobalEnv)
  assign("path_stoolInfoXlsx",    "../data/1g - blood markers.xlsx", envir = .GlobalEnv)
  assign("path_bodyMeasuresXlsx", "../data/1h - body measurements.xlsx", envir = .GlobalEnv)
  
  # symptom diary data
  assign("path_symptomDiaryData", "../data/1i - symptom diaries.xlsx", envir = .GlobalEnv)
  
  # stool sample taxonomy data
  assign("path_absTaxonomyData",  "../data/2a - stool taxonomy (absolute abundances).tsv", envir = .GlobalEnv)
  assign("path_relTaxonomyData",  "../data/2b - stool taxonomy (relative abundances).tsv", envir = .GlobalEnv)
  assign("path_sampleLookupXlsx", "../data/2c - meta information blood and stool samples.xlsx", envir = .GlobalEnv)
  assign("path_krautLookupXlsx",  "../data/2d - meta information sauerkraut samples.xlsx", envir = .GlobalEnv)
  
  # stool sample functional data
  assign("path_KEGGKO_funData",       "../data/3a - KEGG genes.txt", envir = .GlobalEnv)
  assign("path_KEGGPathways_funData", "../data/3b - KEGG pathways.txt", envir = .GlobalEnv)
  
  # blood marker data
  assign("path_bloodMarkerBostonRdata",       "../data/4a - blood data part 1.Rdata", envir = .GlobalEnv)
  assign("path_bloodMarkerZentrallaborRdata", "../data/4b - blood data part 2.Rdata", envir = .GlobalEnv)
  assign("path_bloodMarkerBostonXlsx",        "../data/4c - blood data part 1 including duplicates.xlsx", envir = .GlobalEnv)
  assign("path_bloodMarkerZentrallaborXlsx",  "../data/4d - blood data part 2 including duplicates.xlsx", envir = .GlobalEnv)
  
  # metabolome data
  assign("path_bloodMetabolomeXlsx", "../data/5a - blood SCFAs.xlsx", envir = .GlobalEnv)
  assign("path_stoolMetabolomeXlsx", "../data/5b - stool SCFAs.xlsx", envir = .GlobalEnv)
  
}
