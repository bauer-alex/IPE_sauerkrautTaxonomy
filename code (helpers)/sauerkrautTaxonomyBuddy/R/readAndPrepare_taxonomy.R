
#' Read and prepare the absolute or relative abundance taxonomy data
#' 
#' Taxonomy data were collected for two types of data: the stool samples of the
#' individual participants, and some sauerkraut samples. This function prepares
#' tables for the absolute or relative abundances of both. The function returns
#' a list of five tables, two for stool samples (original samples and some
#' duplicate samples to check measurement variability), one for the sauerkraut
#' samples and two lookup tables (one each for the stool and sauerkraut samples).
#' 
#' @param path_taxonomyData Path to either the \code{merged_abundance_table_counts.tsv}
#' or the \code{merged_abundance_table_relab.tsv} file.
#' @param aggregation_level One of \code{c("species","strain")}, indicating if
#' the datasets should be aggregated on species level or on bacteria strain level.
#' Defaults to \code{"species"}.
#' @param remove_unobservedBacteria Indicator if bacteria species / strains without
#' any appearance should be removed from the stool and sauerkraut datasets.
#' Defaults to TRUE.
#' @param exclude_rareStoolBacteria Indicator if bacteria species / strains in
#' the stool dataset should be excluded if they don't have a minimum relative
#' abundance of 0.01% in at least 10% of all samples. Defaults to TRUE.
#' @inheritParams readAndPrepare_lookupData
#' 
#' @import checkmate dplyr
#' @importFrom readxl read_excel
#' @export
#' 
#' @return List of five data.frames
#' 
readAndPrepare_taxonomy <- function(path_taxonomyData,
                                    path_sampleLookupXlsx,
                                    path_krautLookupXlsx,
                                    aggregation_level             = "species",
                                    remove_unobservedBacteria     = TRUE,
                                    exclude_rareStoolBacteria     = TRUE,
                                    exclude_medicatedParticipants = TRUE,
                                    exclude_sickParticipants) {
  
  checkmate::assert_file(path_taxonomyData,     extension = "tsv")
  checkmate::assert_file(path_sampleLookupXlsx, extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,  extension = "xlsx")
  checkmate::assert_choice(aggregation_level,   choices = c("species","strain"))
  checkmate::assert_logical(remove_unobservedBacteria,     len = 1)
  checkmate::assert_logical(exclude_rareStoolBacteria,     len = 1)
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  dat             <- read.table(path_taxonomyData, header = TRUE, na.strings = "-")
  dat_lookupList  <- readAndPrepare_lookupData(path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                               path_krautLookupXlsx          = path_krautLookupXlsx,
                                               exclude_medicatedParticipants = exclude_medicatedParticipants,
                                               exclude_sickParticipants      = exclude_sickParticipants)
  dat_lookupStool <- dat_lookupList$dat_lookupStool
  dat_lookupKraut <- dat_lookupList$dat_lookupKraut
  
  
  
  # organize the classification information in individual columns -----------
  # limit the data to species information or keep it on bacteria strain level
  class_list <- dat$clade_name %>% strsplit("\\|")
  rows_obs   <- which(sapply(class_list, length) == ifelse(aggregation_level == "species", 7, 8))
  dat        <- dat %>% slice(rows_obs)
  class_list <- class_list[rows_obs]

  if (aggregation_level == "species") {  
    
    dat <- dat %>% 
      select(-clade_name, -clade_taxid) %>% 
      mutate(species_kingdom = sapply(class_list, function(x) { x[1] }),
             species_phylum  = sapply(class_list, function(x) { x[2] }),
             species_class   = sapply(class_list, function(x) { x[3] }),
             species_order   = sapply(class_list, function(x) { x[4] }),
             species_family  = sapply(class_list, function(x) { x[5] }),
             species_genus   = sapply(class_list, function(x) { x[6] }),
             species         = sapply(class_list, function(x) { x[7] })) %>% 
      mutate(across(starts_with("species"), function(x) { gsub(".__", "", x) })) %>% 
      mutate(across(starts_with("species"), factor)) %>% 
      select(starts_with("species"), everything())
    
  } else { # aggregation_level == "strain"
    dat <- dat %>% 
      select(-clade_name, -clade_taxid) %>% 
      mutate(strain_kingdom = sapply(class_list, function(x) { x[1] }),
             strain_phylum  = sapply(class_list, function(x) { x[2] }),
             strain_class   = sapply(class_list, function(x) { x[3] }),
             strain_order   = sapply(class_list, function(x) { x[4] }),
             strain_family  = sapply(class_list, function(x) { x[5] }),
             strain_genus   = sapply(class_list, function(x) { x[6] }),
             strain_species = sapply(class_list, function(x) { x[7] }),
             strain         = sapply(class_list, function(x) { x[8]})) %>% 
      mutate(across(starts_with("strain"), function(x) { gsub(".__", "", x) })) %>% 
      mutate(across(starts_with("strain"), factor)) %>% 
      select(starts_with("strain"), everything())
  }
  
  
  
  # separate the stool from the sauerkraut samples --------------------------
  # delete the chicken sample, it only being a control sample of the laboratory
  col_chicken <- which(grepl("chicken", colnames(dat)))
  dat         <- dat %>% select(-all_of(col_chicken))
  
  # simplify the column names
  cols_samples <- which(grepl("Sauerkraut", colnames(dat)))
  sample_ids   <- sapply(cols_samples, function(x) { strsplit(colnames(dat)[x], "\\_")[[1]][6] })
  colnames(dat)[cols_samples] <- paste0("sample_", sample_ids)
  
  # ensure proper sorting of the sample columns
  cols_sample <- colnames(dat)[grepl("sample_", colnames(dat))]
  dat <- dat %>% 
    select(starts_with("species"), starts_with("strain"), sort(cols_sample))
  
  # separate the stool samples from the sauerkraut samples
  col_firstKrautSample <- which(grepl("sample_508", colnames(dat)))
  dat_stool <- dat %>% select(all_of(1:(col_firstKrautSample - 1)))
  dat_kraut <- dat %>% select(starts_with(aggregation_level),
                              all_of(col_firstKrautSample:ncol(dat)))
  
  
  
  # remove rare stool bacteria species / strains from the data --------------
  if (exclude_rareStoolBacteria) {
    
    if (max(dat_stool$sample_001) > 100) { # absolute abundance data
      warning("Not excluding rare stool bacteria, since this functionality is currently only applicable to relative abundance data.")
      
    } else { # relative abundance data
      
      message("Only retaining stool bacteria with a minimum abundance of 0.01% in at least 10% of samples...")
      
      n_samples <- colnames(dat_stool)[grepl("sample", colnames(dat_stool))] %>% length()
      
      dat_stool <- dat_stool %>% 
        mutate(n_samplesWithRelevantAbundance = rowSums(across(starts_with("sample"), function(x) { ifelse(x > 0.01, 1, 0) } ))) %>% 
        filter(n_samplesWithRelevantAbundance >= 0.1 * n_samples)
    }
  }
  
  
  
  # move the duplicated stool measurements to another dataset ---------------
  sampleIds_dupl <- dat_lookupStool %>% 
    filter(!is.na(id_duplicateSample)) %>% 
    pull(id_duplicateSample) %>% 
    sort()
  
  # split 'dat_stool' in two, placing the duplicate samples in an extra dataset
  dat_stoolDuplicates <- dat_stool %>%
    select(starts_with(aggregation_level), all_of(paste0("sample_", sampleIds_dupl)))
  dat_stool           <- dat_stool %>%
    select(-all_of(paste0("sample_", sampleIds_dupl)))
  
  
  
  # remove medicated and / or sick participants -----------------------------
  # main stool dataset
  sample_cols <- colnames(dat_stool)[grepl("sample_", colnames(dat_stool))]
  cols_toDrop <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  dat_stool   <- dat_stool %>% select(-cols_toDrop)

  # duplicates stool dataset
  sample_cols         <- colnames(dat_stoolDuplicates)[grepl("sample_", colnames(dat_stoolDuplicates))]
  cols_toDrop         <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  dat_stoolDuplicates <- dat_stoolDuplicates %>% select(-cols_toDrop)
  
  
  
  # remove unobserved species / strains from the data ------------------------
  if (remove_unobservedBacteria) {
    # remove the unobserved bacteria from the stool, the stool duplicate and the sauerkraut samples
    dat_list <- lapply(list(dat_stool, dat_stoolDuplicates, dat_kraut), function(dat) {
      
      dat %>%
        mutate(row_sum    = rowSums(across(starts_with("sample"))),
               row_isZero = (row_sum == 0)) %>%
        filter(!row_isZero) %>% 
        select(-row_sum, -row_isZero)
      
    })
    
    # save the cleaned datasets
    dat_stool           <- dat_list[[1]]
    dat_stoolDuplicates <- dat_list[[2]]
    dat_kraut           <- dat_list[[3]]
  }
  
  
  
  # final preparations ------------------------------------------------------
  dat_list <- list("dat_stoolSamples"          = dat_stool           %>% as_tibble(),
                   "dat_stoolDuplicateSamples" = dat_stoolDuplicates %>% as_tibble(),
                   "dat_krautSamples"          = dat_kraut           %>% as_tibble(),
                   "dat_lookupStool"           = dat_lookupStool,
                   "dat_lookupKraut"           = dat_lookupKraut)
  
  return(dat_list)
}
