
#' Read and prepare the symptom diary data
#' 
#' Each subject in the study filled out a symptom diary to track potential
#' irregularities during the study period. Read the corresponding Excel file
#' containing all subjects and prepare the dataset.
#' 
#' @param path Absolute or relative path to the Excel file
#' 
#' @import checkmate dplyr
#' @importFrom readxl read_excel
#' @export
#' 
readAndPrepare_symptomDiaries <- function(path) {
  
  checkmate::assert_file(path, extension = "xlsx")
  
  
  dat <- readxl::read_excel(path, sheet = 1, skip = 1, na = c("0","999"),
                            .name_repair = "unique_quiet")
  
  
  
  # data preparation --------------------------------------------------------
  # rename variables
  dat <- dat %>% 
    rename(subject = "Probanden-ID")
  
  # TODO
  # date values in several columns (e.g. column "SA1aa") are not read correctly,
  # but as numbers. In case we want to work with this information as date values,
  # manual transformation is necessary. Excel numbers behind dates encode the
  # number of days since 1899-12-31.
  # To see how such a conversion can look like, see openxlsx::convertToDate
  
  
  return(dat)
}