
#' Read and prepare the sauerkraut study experiment data
#' 
#' This function basically wraps Virginies code to create one table comprising
#' all the study participants' metadata. On top, variable naming is done more
#' consistently with this function.
#' 
#' @param path_groupInfoRdata Path to the \code{Gruppe_Info.Rdata} file
#' @param path_participantRdata Path to the \code{DATA_participant.Rdata} file
#' @param path_samplesIDxlsx Path to the \code{Samples_ID.xlsx} file
#' @param path_dietInfoCsv Path to the \code{FFQ2_nondietary_data.csv} file
#' @param path_nutrientInfoCsv Path to the \code{FFQ2_nutrintake.csv} file
#' @param path_bloodMetabolomeXlsx Path to the \code{FREI-01-24TASA CDT.XLSX} file
#' @param path_stoolMetabolomeXlsx Path to the \code{FREI-0302-21TASA CDT.XLSX} file
#' @param path_T4questXlsx Path to the \code{Dateneingabe Meinungsfragebögen Sauerkrautstudie.xlsx} file
#' @param path_stoolInfoXlsx Path to the \code{Probensammlung_mit Zentrallab_Duplikate.xlsx} file
#' @param path_bodyMeasuresXlsx Path to the \code{Koerpermesswerte.xlsx} file
#' @inheritParams readAndPrepare_lookupData
#' 
#' @import checkmate dplyr
#' @importFrom readxl read_excel
#' @importFrom tidyr pivot_wider
#' @export
#' 
#' @return data.frame
#' 
readAndPrepare_experimentData <- function(path_groupInfoRdata,
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
                                          exclude_medicatedParticipants = TRUE,
                                          exclude_sickParticipants) {
  
  checkmate::assert_file(path_groupInfoRdata,      extension = "Rdata")
  checkmate::assert_file(path_participantRdata,    extension = "Rdata")
  checkmate::assert_file(path_samplesIDxlsx,       extension = "xlsx")
  checkmate::assert_file(path_dietInfoCsv,         extension = "csv")
  checkmate::assert_file(path_nutrientInfoCsv,     extension = "csv")
  checkmate::assert_file(path_bloodMetabolomeXlsx, extension = "XLSX")
  checkmate::assert_file(path_stoolMetabolomeXlsx, extension = "XLSX")
  checkmate::assert_file(path_stoolMetabolomeXlsx, extension = "XLSX")
  checkmate::assert_file(path_T4questXlsx,         extension = "xlsx")
  checkmate::assert_file(path_stoolInfoXlsx,       extension = "xlsx")
  checkmate::assert_file(path_bodyMeasuresXlsx,    extension = "xlsx")
  checkmate::assert_file(path_sampleLookupXlsx,    extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,     extension = "xlsx")
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  load(path_groupInfoRdata)
  load(path_participantRdata)
  dat_sample        <- readxl::read_excel(path_samplesIDxlsx, sheet = 1,
                                          .name_repair = "unique_quiet")
  dat_diet          <- read.csv(path_dietInfoCsv)
  dat_nutrients     <- read.csv(path_nutrientInfoCsv)
  dat_bloodMeta     <- readxl::read_excel(path_bloodMetabolomeXlsx, na = "< 0")
  dat_stoolMeta     <- readxl::read_excel(path_stoolMetabolomeXlsx)
  dat_T4quest       <- readxl::read_excel(path_T4questXlsx, sheet = 1, skip = 1, na = "999")
  dat_stoolInfoList <- lapply(c("T1","T2","T3","T4","T5"), function(tp) {
    readxl::read_excel(path_stoolInfoXlsx, sheet = tp, skip = 2) %>% 
      suppressMessages() %>% 
      mutate(timepoint = tp)
  })
  dat_bodyMeasures  <- readxl::read_excel(path_bodyMeasuresXlsx, sheet = 1, na = "N/A") %>% suppressMessages()
  dat_lookupList    <- readAndPrepare_lookupData(path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                                 path_krautLookupXlsx          = path_krautLookupXlsx,
                                                 exclude_medicatedParticipants = exclude_medicatedParticipants,
                                                 exclude_sickParticipants      = exclude_sickParticipants)
  dat_lookupStool   <- dat_lookupList$dat_lookupStool
  dat_lookupKraut   <- dat_lookupList$dat_lookupKraut
  
  
  
  # individual data preparation ---------------------------------------------
  # sample data preparation
  dat_sample <- as.data.frame(dat_sample)
  dat_sample <- dat_sample[, 1:4]
  
  
  # diet data
  dat_diet <- dat_diet %>% 
    dplyr::rename(Participant_ID = subj_id,
                  diet_veggie    = fleischlos) %>% 
    mutate(diet_veggie = factor(diet_veggie, levels = c("1", "3", "2"),
                                labels = c("omnivore", "vegetarian", "vegan")),
           diet_veggie = droplevels(diet_veggie)) %>% 
    select(Participant_ID, diet_veggie)
  
  
  # nutrient intake data
  dat_nutrients <- dat_nutrients %>% 
    dplyr::rename(Participant_ID         = subj_id,
                  dailyConsumption_fiber = ZB,
                  dailyConsumption_fat   = ZF,
                  dailyConsumption_water = ZW) %>% 
    select(Participant_ID, starts_with("daily"))
  
  
  # blood metabolome data
  dat_bloodMeta <- dat_bloodMeta %>% 
    select("Client Sample ID", Analyte, "Results (ng/mL)") %>% 
    dplyr::rename(sample_id = "Client Sample ID",
                  results   = "Results (ng/mL)") %>% 
    mutate(results = case_when(is.na(results) ~ 0.01,
                               TRUE           ~ results),
           Analyte = tolower(Analyte),
           Analyte = gsub("2-", "", Analyte),
           Analyte = gsub(" ", "", Analyte),
           Analyte = gsub("acid", "Acid", Analyte),
           Analyte = paste0("blood_", Analyte)) %>% 
    tidyr::pivot_wider(values_from = results, id_cols = sample_id, names_from = Analyte)
  
  
  # stool metabolome data
  dat_stoolMeta <- dat_stoolMeta %>% 
    select("Client Sample ID", Analyte, Results) %>% 
    dplyr::rename(sample_id = "Client Sample ID") %>% 
    mutate(sample_id = as.numeric(sample_id),
           Analyte   = tolower(Analyte),
           Analyte   = gsub("2-", "", Analyte),
           Analyte   = gsub(" ", "", Analyte),
           Analyte   = gsub("acid", "Acid", Analyte),
           Analyte   = paste0("stool_", Analyte)) %>% 
    tidyr::pivot_wider(values_from = Results, id_cols = sample_id, names_from = Analyte)
  
  
  # T4 questionnaire data
  dat_T4quest <- dat_T4quest %>% 
    dplyr::rename(Participant_ID          = "ID-Nr.",
                  goodTaste_fresh         = "SA4a",
                  goodTaste_past          = "SA4b",
                  byeffects_fresh         = "SA6aa",
                  byeffectsOthers_fresh   = "SA6ab",
                  byeffects_past          = "SA6ba",
                  byeffectsOthers_past    = "SA6bb",
                  posEffects_fresh        = "SA7caa",
                  posEffects_past         = "SA7cba",
                  easyDiet_fresh          = "SA8aa",
                  easyDiet_past           = "SA8ba",
                  easyNotEatingOtherFoods = "SA9a",
                  favoriteSauerkrautType  = "SA10",
                  futureDailySK_fresh     = "SA11a",
                  futureDailySK_past      = "SA11d",
                  participant_inMenopause = "SA13") %>% 
    mutate(flatulence_afterFresh      = case_when(grepl("2", byeffects_fresh) | grepl("Blähung", byeffectsOthers_fresh) |
                                                    grepl("Darmwind", byeffectsOthers_fresh) ~ "yes",
                                                  TRUE                                       ~ "no"),
           flatulence_afterFresh      = factor(flatulence_afterFresh),
           flatulence_afterPast       = case_when(grepl("2", byeffects_past) | grepl("Blähung", byeffectsOthers_past) |
                                                    grepl("Darmwind", byeffectsOthers_past) ~ "yes",
                                                  TRUE                                      ~ "no"),
           flatulence_afterPast       = factor(flatulence_afterPast),
           betterDigestion_afterFresh = case_when(grepl("1", posEffects_fresh) ~ "yes",
                                                  TRUE                         ~ "no"),
           betterDigestion_afterFresh = factor(betterDigestion_afterFresh),
           betterDigestion_afterPast  = case_when(grepl("1", posEffects_past) ~ "yes",
                                                  TRUE                        ~ "no"),
           betterDigestion_afterPast  = factor(betterDigestion_afterPast),
           byeffects_fresh            = case_when(is.na(goodTaste_fresh)     ~ NA_character_,
                                                  is.na(byeffects_fresh)     ~ "0",
                                                  TRUE                       ~ byeffects_fresh),
           byeffects_fresh            = factor(byeffects_fresh,
                                               levels = c("0","1","1+2","2","2+3","1+3","3","1+2+3"),
                                               labels = c("no adverse reaction",
                                                          "stomachache",
                                                          "stomachache & flatulence",
                                                          "flatulence",
                                                          "flatulence & diarrhea",
                                                          "stomachache & diarrhea",
                                                          "diarrhea",
                                                          "all of the above")),
           byeffects_past            = case_when(is.na(goodTaste_past)     ~ NA_character_,
                                                 is.na(byeffects_past)     ~ "0",
                                                 TRUE                       ~ byeffects_past),
           byeffects_past             = factor(byeffects_past,
                                               levels = c("0","1","1+2","2","2+3","1+3","3","1+2+3"),
                                               labels = c("no adverse reaction",
                                                          "stomachache",
                                                          "stomachache & flatulence",
                                                          "flatulence",
                                                          "flatulence & diarrhea",
                                                          "stomachache & diarrhea",
                                                          "diarrhea",
                                                          "all of the above")),
           posEffects_fresh           = case_when(is.na(goodTaste_fresh)     ~ NA_character_,
                                                  is.na(posEffects_fresh)    ~ "0",
                                                  TRUE                       ~ posEffects_fresh),
           posEffects_fresh           = factor(posEffects_fresh,
                                               levels = c("0","1","1+2","2","2+3","1+3","3","1+2+3"),
                                               labels = c("no positive effects",
                                                          "better digestion",
                                                          "better digestion & overall condition",
                                                          "better overall condition",
                                                          "better overall condition & less hunger",
                                                          "better digestion & less hunger",
                                                          "less hunger",
                                                          "all of the above")),
           posEffects_past            = case_when(is.na(goodTaste_fresh)     ~ NA_character_,
                                                  is.na(posEffects_past)     ~ "0",
                                                  TRUE                       ~ posEffects_past),
           posEffects_past            = factor(posEffects_past,
                                               levels = c("0","1","1+2","2","2+3","1+3","3","1+2+3"),
                                               labels = c("no positive effects",
                                                          "better digestion",
                                                          "better digestion & overall condition",
                                                          "better overall condition",
                                                          "better overall condition & less hunger",
                                                          "better digestion & less hunger",
                                                          "less hunger",
                                                          "all of the above")),
           favoriteSauerkrautType     = factor(favoriteSauerkrautType, levels = c(1, 2, 3),
                                               labels = c("fresh sauerkraut", "pasteurized sauerkraut", "no difference")),
           futureDailySK_fresh        = factor(futureDailySK_fresh, levels = c(1, 2),
                                               labels = c("yes, I could see myself eating fresh sauerkraut daily", "no")),
           futureDailySK_past         = factor(futureDailySK_past, levels = c(1, 2),
                                               labels = c("yes, I could see myself eating pasteurized sauerkraut daily", "no")),
           participant_inMenopause    = grepl("3", participant_inMenopause) %>% factor(levels = c(FALSE, TRUE), labels = c("no", "yes"))) %>% 
    rowwise() %>% 
    mutate(id_num = strsplit(Participant_ID, split = " ")[[1]][2],
           id_num = case_when(nchar(id_num) == 1 ~ paste0("00", id_num),
                              nchar(id_num) == 2 ~ paste0("0",  id_num),
                              TRUE               ~ as.character(id_num)),
           Participant_ID = paste0("SK", id_num)) %>% 
    select(-starts_with("SA"), -id_num)
  
  
  # stool info data
  dat_stoolInfoList <- lapply(dat_stoolInfoList, function(dat) {
    dat %>% 
      mutate(stoolType_healthy = case_when(is.na(StoolType)                ~ NA_character_,
                                           StoolType >= 3 & StoolType <= 4 ~ "yes",
                                           TRUE                            ~ "no"),
             stoolType_healthy = factor(stoolType_healthy)) %>% 
      select(Participant_ID, timepoint, Collection_date, bowel_mvnt_freq, StoolType, stoolType_healthy, PH_value) %>% 
      dplyr::rename(sample_collectionDate = Collection_date,
                    bowelMovement_freq    = bowel_mvnt_freq,
                    stool_type            = StoolType,
                    ph_value              = PH_value) %>% 
      mutate(bowelMovement_freq = factor(bowelMovement_freq, levels = c("> 3 x taeglich", "3 x taeglich", "2 x taeglich", "1 x taeglich",
                                                                        "1 x alle 2 Tage", "1 x alle 3 Tage", "< 1 x alle 3 Tage")),
             ph_value = case_when(is.na(ph_value)          ~ NA_character_,
                                  ph_value == "zu trocken" ~ NA_character_,
                                  TRUE                     ~ ph_value),
             ph_value = ph_value %>% as.numeric())
  })
  dat_stoolInfo <- dat_stoolInfoList %>%
    dplyr::bind_rows() %>% 
    arrange(Participant_ID, timepoint)
  
  
  # body measures data
  # Note: We currently only use blood pressure data from this data file.
  #       Blood pressure was measured thrice at each timepoint. We discard the
  #       first one (due to potential biases, e.g. due to excitement) and take
  #       the average of the remaining two measurements.
  dat_bodyMeasures <- dat_bodyMeasures %>% 
    select("ID-Nummer", contains("Blutdruck")) %>% 
    select(-contains("_1")) %>% 
    dplyr::rename(Participant_ID = "ID-Nummer") %>% 
    filter(!is.na(Participant_ID),
           substr(Participant_ID, 1, 2) == "SK") %>% 
    mutate(across(starts_with("Blutdruck"), function(x) { ifelse(grepl("/", x), x, NA) } ))
  colnames(dat_bodyMeasures) <- colnames(dat_bodyMeasures) %>%
    gsub(pattern = "Blutdruck ", replacement = "bloodPressure_") %>% 
    gsub(pattern = "Sart",       replacement = "T0")
  dat_bodyMeasures <- dat_bodyMeasures %>% 
    rowwise() %>% 
    mutate(systolic_t0_2  = strsplit(bloodPressure_T0_2, split = "/")[[1]][1],
           systolic_t0_3  = strsplit(bloodPressure_T0_3, split = "/")[[1]][1],
           systolic_t1_2  = strsplit(bloodPressure_T1_2, split = "/")[[1]][1],
           systolic_t1_3  = strsplit(bloodPressure_T1_3, split = "/")[[1]][1],
           systolic_t2_2  = strsplit(bloodPressure_T2_2, split = "/")[[1]][1],
           systolic_t2_3  = strsplit(bloodPressure_T2_3, split = "/")[[1]][1],
           systolic_t3_2  = strsplit(bloodPressure_T3_2, split = "/")[[1]][1],
           systolic_t3_3  = strsplit(bloodPressure_T3_3, split = "/")[[1]][1],
           systolic_t4_2  = strsplit(bloodPressure_T4_2, split = "/")[[1]][1],
           systolic_t4_3  = strsplit(bloodPressure_T4_3, split = "/")[[1]][1],
           systolic_t5_2  = strsplit(bloodPressure_T5_2, split = "/")[[1]][1],
           systolic_t5_3  = strsplit(bloodPressure_T5_3, split = "/")[[1]][1],
           diastolic_t0_2 = strsplit(bloodPressure_T0_2, split = "/")[[1]][2],
           diastolic_t0_3 = strsplit(bloodPressure_T0_3, split = "/")[[1]][2],
           diastolic_t1_2 = strsplit(bloodPressure_T1_2, split = "/")[[1]][2],
           diastolic_t1_3 = strsplit(bloodPressure_T1_3, split = "/")[[1]][2],
           diastolic_t2_2 = strsplit(bloodPressure_T2_2, split = "/")[[1]][2],
           diastolic_t2_3 = strsplit(bloodPressure_T2_3, split = "/")[[1]][2],
           diastolic_t3_2 = strsplit(bloodPressure_T3_2, split = "/")[[1]][2],
           diastolic_t3_3 = strsplit(bloodPressure_T3_3, split = "/")[[1]][2],
           diastolic_t4_2 = strsplit(bloodPressure_T4_2, split = "/")[[1]][2],
           diastolic_t4_3 = strsplit(bloodPressure_T4_3, split = "/")[[1]][2],
           diastolic_t5_2 = strsplit(bloodPressure_T5_2, split = "/")[[1]][2],
           diastolic_t5_3 = strsplit(bloodPressure_T5_3, split = "/")[[1]][2],
           across(contains("stolic"), as.numeric),
           bloodPressure_avgSys_t0 = mean(c(systolic_t0_2, systolic_t0_3), na.rm = TRUE),
           bloodPressure_avgSys_t1 = mean(c(systolic_t1_2, systolic_t1_3), na.rm = TRUE),
           bloodPressure_avgSys_t2 = mean(c(systolic_t2_2, systolic_t2_3), na.rm = TRUE),
           bloodPressure_avgSys_t3 = mean(c(systolic_t3_2, systolic_t3_3), na.rm = TRUE),
           bloodPressure_avgSys_t4 = mean(c(systolic_t4_2, systolic_t4_3), na.rm = TRUE),
           bloodPressure_avgSys_t5 = mean(c(systolic_t5_2, systolic_t5_3), na.rm = TRUE),
           bloodPressure_avgDia_t0 = mean(c(diastolic_t0_2, diastolic_t0_3), na.rm = TRUE),
           bloodPressure_avgDia_t1 = mean(c(diastolic_t1_2, diastolic_t1_3), na.rm = TRUE),
           bloodPressure_avgDia_t2 = mean(c(diastolic_t2_2, diastolic_t2_3), na.rm = TRUE),
           bloodPressure_avgDia_t3 = mean(c(diastolic_t3_2, diastolic_t3_3), na.rm = TRUE),
           bloodPressure_avgDia_t4 = mean(c(diastolic_t4_2, diastolic_t4_3), na.rm = TRUE),
           bloodPressure_avgDia_t5 = mean(c(diastolic_t5_2, diastolic_t5_3), na.rm = TRUE),
           across(contains("bloodPressure"), function(x) { ifelse(is.nan(x), NA_real_, x) })) %>% 
    ungroup() %>% 
    select(Participant_ID, starts_with("bloodPressure_avg"))
  dat_bodyMeasures_long <- data.frame(Participant_ID       = rep(dat_bodyMeasures$Participant_ID, times = length(0:5)),
                                      timepoint            = rep(c("T0","T1","T2","T3","T4","T5"), each = nrow(dat_bodyMeasures)),
                                      bloodPressure_avgSys = c(dat_bodyMeasures$bloodPressure_avgSys_t0,
                                                               dat_bodyMeasures$bloodPressure_avgSys_t1,
                                                               dat_bodyMeasures$bloodPressure_avgSys_t2,
                                                               dat_bodyMeasures$bloodPressure_avgSys_t3,
                                                               dat_bodyMeasures$bloodPressure_avgSys_t4,
                                                               dat_bodyMeasures$bloodPressure_avgSys_t5),
                                      bloodPressure_avgDia = c(dat_bodyMeasures$bloodPressure_avgDia_t0,
                                                               dat_bodyMeasures$bloodPressure_avgDia_t1,
                                                               dat_bodyMeasures$bloodPressure_avgDia_t2,
                                                               dat_bodyMeasures$bloodPressure_avgDia_t3,
                                                               dat_bodyMeasures$bloodPressure_avgDia_t4,
                                                               dat_bodyMeasures$bloodPressure_avgDia_t5))
  
  
  
  # merge all datasets ------------------------------------------------------
  dat <- dat_sample %>% 
    dplyr::left_join(Gruppe_Info,           by = "Participant_ID") %>% 
    merge(PartiInfo_T0[,c(1:3,6,9)],        by = "Participant_ID", all = TRUE) %>% 
    dplyr::left_join(dat_diet,              by = "Participant_ID") %>% 
    dplyr::left_join(dat_nutrients,         by = "Participant_ID") %>% 
    dplyr::left_join(dat_bloodMeta,         by = c("Sample_ID" = "sample_id")) %>% 
    dplyr::left_join(dat_stoolMeta,         by = c("Sample_ID" = "sample_id")) %>% 
    dplyr::left_join(dat_T4quest,           by = "Participant_ID") %>% 
    dplyr::left_join(dat_stoolInfo,         by = c("Participant_ID" = "Participant_ID", "Timepoint" = "timepoint")) %>% 
    dplyr::left_join(dat_bodyMeasures_long, by = c("Participant_ID" = "Participant_ID", "Timepoint" = "timepoint"))
  dat <- dat[order(dat$Part_Time),]
  dat <- dataPreparation_enrichExperimentData(dat)
  
  rownames(dat) <- dat[, "Sample_ID"]
  
  
  
  # final data preparation --------------------------------------------------
  # rename variables
  dat <- dat %>% 
    dplyr::rename(participant_id   = Participant_ID,
                  sample_id        = Sample_ID,
                  part_time        = Part_Time,
                  timepoint        = Timepoint,
                  group            = Gruppe,
                  gender           = Geschlecht,
                  age              = Alter,
                  bmi_t0           = BMI_T0,
                  hipWaistRatio_t0 = HipWaist_T0,
                  treatment        = Treatment,
                  intervention     = Intervention,
                  period           = Period,
                  inter_treat      = Inter_treat)
  
  # order columns and consistently format the sample id's
  dat <- dat %>%
    select(sample_id, participant_id, timepoint, part_time, group, treatment,
           intervention, period, inter_treat, sample_collectionDate, everything()) %>% 
    mutate(sample_id = case_when(nchar(sample_id) == 1 ~ paste0("00", sample_id),
                                 nchar(sample_id) == 2 ~ paste0("0",  sample_id),
                                 TRUE                  ~ as.character(sample_id)))
  
  
  
  # remove medicated and / or sick participants -----------------------------
  dat <- dat %>% 
    filter(part_time %in% dat_lookupStool$subject_time) %>% 
    mutate(participant_id = droplevels(participant_id))
  
    
  return(dat %>% as_tibble())
}



#' Helper function to further enrich the experiment data
#' 
#' This function basically wraps Virginies function \code{add_meta()} to
#' perform some data preparation steps for the participants' metadata forming
#' the experiment data.
#' Currently only called from within \code{\link{readAndPrepare_experimentData}}.
#' 
#' @param dat data.frame comprising the experiment data
#' 
#' @import checkmate
#' 
#' @return data.frame
#' 
dataPreparation_enrichExperimentData <- function(dat){
  
  checkmate::assert_data_frame(dat)
  
  
  dat <- dat %>% 
    mutate(Treatment    = case_when(Timepoint %in% c("T1","T3","T5") ~ "Baseline",
                                    Timepoint %in% c("T2","T4")      ~ "After",
                                    TRUE                             ~ NA_character_),
           Intervention = case_when(Timepoint %in% c("T1","T2") & Gruppe == "A" ~ "Fresh",
                                    Timepoint %in% c("T3","T4") & Gruppe == "A" ~ "Pasteurized",
                                    Timepoint %in% c("T1","T2") & Gruppe == "B" ~ "Pasteurized",
                                    Timepoint %in% c("T3","T4") & Gruppe == "B" ~ "Fresh",
                                    Timepoint == "T5"                           ~ "FollowUp",
                                    TRUE                                        ~ NA_character_),
           Intervention = factor(Intervention, levels = c("Fresh", "Pasteurized", "FollowUp")),
           Period       = case_when(Timepoint %in% c("T1", "T2") ~ "P1",
                                    Timepoint %in% c("T3", "T4") ~ "P2",
                                    Timepoint == "T5"            ~ "P3",
                                    TRUE                         ~ NA_character_))
  
  dat <- dat %>% 
    mutate(Inter_treat    = paste(Treatment, Intervention, sep = "."),
           Inter_treat    = factor(Inter_treat, levels = c("Baseline.Fresh", "After.Fresh", "Baseline.Pasteurized", "After.Pasteurized", "Baseline.FollowUp"), ordered = TRUE),
           Treatment      = factor(Treatment,   levels = c("Baseline", "After"),          ordered = TRUE),
           Timepoint      = factor(Timepoint,   levels = c("T1", "T2", "T3", "T4", "T5"), ordered = TRUE),
           Gruppe         = factor(Gruppe,      levels = c("A", "B"), labels = c("Fresh first", "Pasteurized first")),
           Participant_ID = factor(Participant_ID),
           Period         = factor(Period,      levels = c("P1", "P2", "P3"), ordered = TRUE))
  
  
  return(dat)
}



#' Prepare the stool and sauerkraut sample lookup datasets
#' 
#' This function prepares the lookup datasets both for the stool samples,
#' blood samples and the samples taken from the sauerkraut jars.
#' 
#' @param path_sampleLookupXlsx Path to the \code{Metabolon_Blutproben_Stuhlproben_intern.xlsx}
#' file which contains one sheet \code{Sample Manifest} summarizing all relevant
#' sample information, including the information which stool and blood samples were
#' analyzed twice in the lab for us to measure the measurement variability.
#' @param path_krautLookupXlsx Path to the \code{Sauerkrautproben für Sequenzierung.xlsx}
#' file which contains lookup information on the sauerkraut samples.
#' @param exclude_medicatedParticipants Indicator if the two participants SK005
#' and SK093 that were on medication during the whole study should be removed
#' from the data. Defaults to TRUE.
#' @param exclude_sickParticipants One of \code{c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper")}.
#' If one of the latter is specified, all measurements biased by a recent sickness
#' are excluded. Since the taxonomy paper (mainly stool measurements) and the blood paper
#' (blood measurements) analyze different biosamples, the exclusion list differs
#' between the two papers. For the stool type analysis, a few more measurements
#' are excluded compared to the taxonomy paper.
#' 
#' Note that medicated or sick participants' measurements' information is not
#' excluded from the returned blood Boston duplicate data since this is only
#' a secondary lookup table used for matching original with duplicate measurements
#' measured by the Boston laboratory, and not used for filtering other data.
#' 
#' @import checkmate dplyr
#' @importFrom readxl read_excel
#' @export
#' 
#' @return List of three data.frames
#' 
readAndPrepare_lookupData <- function(path_sampleLookupXlsx,
                                      path_krautLookupXlsx,
                                      exclude_medicatedParticipants = TRUE,
                                      exclude_sickParticipants) {
  
  checkmate::assert_file(path_sampleLookupXlsx,            extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,             extension = "xlsx")
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  # read raw data
  dat_stool       <- readxl::read_excel(path         = path_sampleLookupXlsx,
                                        sheet        = "Sample Manifest",
                                        skip         = 13,
                                        .name_repair = "unique_quiet")
  dat_blood_part1 <- readxl::read_excel(path         = path_sampleLookupXlsx,
                                        sheet        = "Duplikat Serum Boston 500µl",
                                        .name_repair = "unique_quiet")
  dat_blood_part2 <- readxl::read_excel(path         = path_sampleLookupXlsx,
                                        sheet        = "Duplikat Serum Boston 1000µl",
                                        .name_repair = "unique_quiet")
  dat_kraut       <- readxl::read_excel(path         = path_krautLookupXlsx,
                                        sheet        = "Tabelle1",
                                        skip         = 24,
                                        col_names    = c("sample_id", "sample_description"))
  
  
  
  # format the stool lookup table -------------------------------------------
  dat_stool <- dat_stool %>% 
    rename(sample_id          = "Sample Number",
           subject_time       = "Unique Tube Label ID*",
           id_duplicateSample = "Duplikat Stuhl alk. Braunschweig") %>% 
    mutate(is_duplicateSample = grepl("_D", subject_time) %>% factor(levels = c(FALSE, TRUE), labels = c("no", "yes")),
           subject_time       = gsub("_D",  "",   subject_time),
           subject_time       = gsub("SK-", "SK", subject_time),
           subject_time       = gsub("-T",  "_T", subject_time)) %>% 
    filter(!is.na(subject_time)) %>% 
    mutate(sample_id = case_when(nchar(sample_id) == 1 ~ paste0("00", sample_id),
                                 nchar(sample_id) == 2 ~ paste0("0",  sample_id),
                                 TRUE                  ~ as.character(sample_id)),
           id_duplicateSample = case_when(id_duplicateSample == "war schon 1x  aufgetaut!!!" ~ NA_character_,
                                          nchar(id_duplicateSample) == 1                     ~ paste0("00", id_duplicateSample),
                                          nchar(id_duplicateSample) == 2                     ~ paste0("0",  id_duplicateSample),
                                          TRUE                                               ~ as.character(id_duplicateSample))) %>% 
    select(sample_id, subject_time, is_duplicateSample, id_duplicateSample)
  
  
  
  # format the blood lookup table -------------------------------------------
  dat_blood_part1 <- dat_blood_part1 %>% 
    dplyr::rename(orig_subjectTime  = "Unique Tube Label ID*",
                  dupl_sampleID     = "Duplikat Serum/Boston 500µl") %>% 
    select(orig_subjectTime, dupl_sampleID) %>% 
    mutate(orig_subjectTime  = gsub(orig_subjectTime, pattern = "SK-", replacement = "SK"),
           orig_subjectTime  = gsub(orig_subjectTime, pattern = "_D",  replacement = ""),
           orig_subjectTime  = gsub(orig_subjectTime, pattern = "-",   replacement = "_"),
           dupl_sampleID = as.character(dupl_sampleID))
  dat_blood_part2 <- dat_blood_part2 %>% 
    dplyr::rename(orig_subjectTime  = "Unique Tube Label ID*",
                  dupl_sampleID     = "Duplikat Serum Boston 1000µl") %>% 
    select(orig_subjectTime, dupl_sampleID) %>% 
    mutate(orig_subjectTime = gsub(orig_subjectTime, pattern = "SK-", replacement = "SK"),
           orig_subjectTime = gsub(orig_subjectTime, pattern = "-",   replacement = "_"))
  
  dat_bloodBostonDupl <- dat_blood_part1 %>% 
    dplyr::bind_rows(dat_blood_part2)
  
  
  
  # format the sauerkraut lookup table --------------------------------------
  dat_kraut <- dat_kraut %>% 
    mutate(sample_glass = case_when(grepl("frisch", sample_description)             ~ "fresh glass 1",
                                    grepl("past",   sample_description)             ~ "past. glass 1",
                                    sample_description %in% c("2.1a","2.2a","2.3a") ~ "fresh glass 2",
                                    sample_description %in% c("2.1b","2.2b","2.3b") ~ "past. glass 2",
                                    sample_description %in% c("3.1a","3.2a","3.3a") ~ "fresh glass 3",
                                    sample_description %in% c("3.1b","3.2b","3.3b") ~ "past. glass 3"),
           sample_glass = factor(sample_glass)) %>% 
    select(sample_id, sample_glass, sample_description)
  
  
  
  # mark and exclude medicated participants ---------------------------------
  dat_stool <- dat_stool %>% 
    mutate(subject_isMedicated = case_when(grepl("SK005", subject_time) ~ "yes",
                                           grepl("SK086", subject_time) ~ "yes",
                                           grepl("SK093", subject_time) ~ "yes",
                                           TRUE                         ~ "no"),
           subject_isMedicated = factor(subject_isMedicated))
  
  # exclude measurements from medicated participants
  if (exclude_medicatedParticipants) {
    dat_stool <- dat_stool %>% 
      filter(subject_isMedicated == "no")
  }
  
  
  
  # mark and exclude sick participants --------------------------------------
  # taxonomy paper
  biased_stoolMeasurements <- c("SK073_T5")
  
  # stool type analysis
  biased_stoolTypeMeasurements <- c("SK026_T2", "SK064_T3", "SK070_T1", "SK073_T5",
                                    "SK085_T2", "SK090_T2", "SK092_T3")
  # blood paper
  biased_bloodMeasurements <- c("SK007_T1", "SK009_T3", "SK035_T2", "SK050_T1",
                                "SK053_T5", "SK064_T3", "SK070_T1", "SK073_T3",
                                "SK073_T5", "SK080_T2", "SK085_T2", "SK094_T3", "SK094_T5")
  
  # add information to the dataset
  dat_stool <- dat_stool %>% 
    mutate(subject_stoolBiasedBySickness     = case_when(subject_time %in% biased_stoolMeasurements ~ "yes",
                                                         TRUE                                       ~ "no"),
           subject_stoolBiasedBySickness     = factor(subject_stoolBiasedBySickness),
           subject_stoolTypeBiasedBySickness = case_when(subject_time %in% biased_stoolTypeMeasurements ~ "yes",
                                                         TRUE                                           ~ "no"),
           subject_stoolTypeBiasedBySickness = factor(subject_stoolTypeBiasedBySickness),
           subject_bloodBiasedBySickness     = case_when(subject_time %in% biased_bloodMeasurements ~ "yes",
                                                         TRUE                                       ~ "no"),
           subject_bloodBiasedBySickness     = factor(subject_bloodBiasedBySickness))
  
      
  # exclude measurements from sick participants
  if (exclude_sickParticipants == "taxonomy paper") {
    dat_stool <- dat_stool %>% filter(subject_stoolBiasedBySickness == "no")
    
  } else if (exclude_sickParticipants == "blood paper") {
    dat_stool <- dat_stool %>% filter(subject_bloodBiasedBySickness == "no")
    
  } else if (exclude_sickParticipants == "stool type analysis") {
    dat_stool <- dat_stool %>% filter(subject_stoolTypeBiasedBySickness == "no")
    
  }
  
  
  
  # return results ----------------------------------------------------------
  list("dat_lookupStool"                  = dat_stool,
       "dat_lookupBlood_BostonDuplicates" = dat_bloodBostonDupl,
       "dat_lookupKraut"                  = dat_kraut) %>% 
    return()
}
