
#' Read and prepare the blood marker data
#' 
#' Since the blood marker data were already previously prepared by Virginie,
#' this function mainly reads Virginie's final \code{.Rdata} data files (for the
#' two laboratories in Boston and the 'Zentrallabor' in Germany), merges them
#' and performs some cosmetic preparations of the data.
#' Additionally, the function uses the raw data Excel files to prepare duplicate
#' datasets, and then returns these and the main measurements as a named list.
#' 
#' @param path_bloodMarkerBostonRdata Path to the file
#' \code{BloodData_Boston.Rdata}
#' @param path_bloodMarkerZentrallaborRdata Path to the file
#' \code{BloodData_Zentrallabor.Rdata}
#' @param path_bloodMarkerBostonXlsx Path to the file
#' \code{23-019 Michels 10 sent 111423(final serum results).xlsx}, for preparing
#' the data on duplicate measurements.
#' @param path_bloodMarkerZentrallaborXlsx Path to the file
#' \code{Probensammlung_mit Zentrallab_Duplikate.xlsx}, for preparing the data
#' on duplicate measurements.
#' @inheritParams readAndPrepare_experimentData
#' @inheritParams readAndPrepare_lookupData
#' 
#' @import checkmate dplyr
#' @importFrom readxl read_excel
#' @export
#' 
readAndPrepare_bloodMarkerData <- function(path_bloodMarkerBostonRdata,
                                           path_bloodMarkerZentrallaborRdata,
                                           path_bloodMarkerBostonXlsx,
                                           path_bloodMarkerZentrallaborXlsx,
                                           path_sampleLookupXlsx,
                                           path_krautLookupXlsx,
                                           exclude_medicatedParticipants = TRUE,
                                           exclude_sickParticipants) {
  
  checkmate::assert_file(path_bloodMarkerBostonRdata,       extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerZentrallaborRdata, extension = "Rdata")
  checkmate::assert_file(path_bloodMarkerBostonXlsx,        extension = "xlsx")
  checkmate::assert_file(path_bloodMarkerZentrallaborXlsx,  extension = "xlsx")
  checkmate::assert_file(path_sampleLookupXlsx,             extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,              extension = "xlsx")
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  # load the main measurements data
  load(path_bloodMarkerBostonRdata)
  dat_boston <- data_merged %>% select(-Sample_ID)
  load(path_bloodMarkerZentrallaborRdata)
  dat_zl     <- data_merged %>% select(-Gruppe)
  rm(data_merged)
  
  # read the duplicate measurements data
  datD_boston  <- readxl::read_excel(path_bloodMarkerBostonXlsx, sheet = "Results")
  sheet_names  <- c("T1", "T2", "T3", "T4", "T5")
  datD_zl_list <- lapply(sheet_names, function(sheet) {
    readxl::read_excel(path_bloodMarkerZentrallaborXlsx, sheet = sheet, skip = 1) %>% 
      mutate(timepoint = sheet) %>% 
      suppressMessages()
  })
  names(datD_zl_list) <- sheet_names
  
  # read the lookup data for duplicate Boston measurement information
  dat_lookupList <- readAndPrepare_lookupData(path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                              path_krautLookupXlsx          = path_krautLookupXlsx,
                                              exclude_medicatedParticipants = exclude_medicatedParticipants,
                                              exclude_sickParticipants      = exclude_sickParticipants)
  
  
  
  # merge main data ---------------------------------------------------------
  dat <- dat_zl %>% 
    dplyr::full_join(dat_boston, by = c("Part_Time","Participant_ID","Timepoint"))
  
  
  
  # final main data preparation ---------------------------------------------
  # make all column names lowercase
  colnames(dat) <- colnames(dat) %>% tolower()
  
  # sort dataset and order columns
  dat <- dat %>%
    arrange(participant_id, timepoint) %>% 
    select(participant_id, timepoint, part_time, tnfr2, il6, scrp, rage, lbp,
           fabp2, zonulin, glucose, fructosamine, insulin)
  
  # format variables
  dat <- dat %>% 
    mutate(participant_id = factor(participant_id),
           timepoint      = factor(timepoint))
  
  
  
  # prepare duplicate data --------------------------------------------------
  lookup_dupl <- dat_lookupList$dat_lookupBlood_BostonDuplicates
  
  # Boston data
  datD_boston <- datD_boston %>% 
    dplyr::rename(dupl_sampleID = "SK-serum/1.000µl",
                  tnfr2         = "TNFR2 (pg/mL)",
                  rage          = "sRAGE (pg/mL)",
                  il6           = "IL-6 (pg/mL)",
                  fabp2         = "FABP2 (pg/mL)",
                  zonulin       = "Zonulin (ng/mL)",
                  lbp           = "LBP (ng/mL)") %>% 
    mutate(dupl_sampleID = gsub(dupl_sampleID, pattern = "SK-", replacement = "")) %>% 
    filter(dupl_sampleID %in% lookup_dupl$dupl_sampleID) %>% 
    dplyr::left_join(lookup_dupl, by = "dupl_sampleID") %>% 
    select(orig_subjectTime, dupl_sampleID, everything()) %>% 
    mutate(tnfr2   = as.numeric(tnfr2),
           rage    = as.numeric(rage),
           il6     = as.numeric(il6),
           fabp2   = case_when(fabp2 == ">5000.000" ~ "6000",
                               TRUE                 ~ fabp2),
           fabp2   = as.numeric(fabp2),
           zonulin = as.numeric(zonulin),
           lbp     = as.numeric(lbp))
  
  # Zentrallabor data
  datD_zl_list <- datD_zl_list %>% lapply(function(x) {
    x <- x %>% select(1, (ncol(x) - 4):ncol(x))
    colnames(x) <- c("participant_id", "glucose", "scrp", "fructosamine", "insulin", "timepoint")
    x %>% 
      mutate(glucose      = as.character(glucose),
             scrp         = as.character(scrp),
             fructosamine = as.character(fructosamine),
             insulin      = as.character(insulin))
  })
  datD_zl <- datD_zl_list %>% 
    dplyr::bind_rows() %>% 
    filter(participant_id != "Participant_ID",
           !is.na(glucose) | !is.na(scrp) | !is.na(fructosamine) | !is.na(insulin)) %>% 
    mutate(orig_subjectTime = paste0(participant_id, "_", timepoint)) %>% 
    select(orig_subjectTime, glucose, scrp, fructosamine, insulin) %>% 
    mutate(glucose      = as.numeric(glucose),
           scrp         = case_when(is.na(scrp)    ~ NA_character_,
                                    scrp == "<0,1" ~ "0.05",
                                    TRUE           ~ scrp),
           scrp         = as.numeric(scrp),
           fructosamine = as.numeric(fructosamine),
           insulin      = as.numeric(insulin))
  
  
  
  # remove medicated participants -------------------------------------------
  if (exclude_medicatedParticipants) {
    dat <- dat %>% 
      filter(!(participant_id %in% c("SK005", "SK093"))) %>% 
      mutate(participant_id = droplevels(participant_id))
    
    datD_boston <- datD_boston %>% 
      filter(!grepl("SK005", orig_subjectTime),
             !grepl("SK093", orig_subjectTime))
    datD_zl <- datD_zl %>% 
      filter(!grepl("SK005", orig_subjectTime),
             !grepl("SK093", orig_subjectTime))
  }
  
  
  
  # remove sick participants ------------------------------------------------
  # taxonomy paper
  biased_stoolMeasurements <- c("SK026_T2", "SK064_T3", "SK070_T1", "SK073_T5",
                                "SK085_T2", "SK086_T1", "SK086_T2", "SK086_T3",
                                "SK086_T4", "SK086_T5", "SK090_T2", "SK092_T3")
  
  # blood paper
  biased_bloodMeasurements <- c("SK007_T1", "SK009_T3", "SK035_T2", "SK050_T1",
                                "SK053_T5", "SK064_T3", "SK070_T1", "SK073_T3",
                                "SK073_T5", "SK080_T2", "SK085_T2", "SK086_T1",
                                "SK086_T2", "SK086_T3", "SK086_T4", "SK086_T5",
                                "SK094_T3", "SK094_T5")
  
    
  if (exclude_sickParticipants == "taxonomy paper") {
    dat <- dat %>% 
      filter(!(part_time %in% biased_stoolMeasurements)) %>% 
      mutate(participant_id = droplevels(participant_id))
    
    datD_boston <- datD_boston %>% filter(!(orig_subjectTime %in% biased_stoolMeasurements))
    datD_zl     <- datD_zl     %>% filter(!(orig_subjectTime %in% biased_stoolMeasurements))
    
  } else if (exclude_sickParticipants == "blood paper") {
    dat <- dat %>% 
      filter(!(part_time %in% biased_bloodMeasurements)) %>% 
      mutate(participant_id = droplevels(participant_id))
    
    datD_boston <- datD_boston %>% filter(!(orig_subjectTime %in% biased_bloodMeasurements))
    datD_zl     <- datD_zl     %>% filter(!(orig_subjectTime %in% biased_bloodMeasurements))
  }
  
  
  list("main_data"                   = dat,
       "duplicates_BostonData"       = datD_boston,
       "duplicates_ZentrallaborData" = datD_zl) %>% 
    return()
}
