
#' Prepare the functional datasets
#' 
#' This function reads the KEGG-KO and KEGG-pathway datasets and prepares them.
#' 
#' @param path_KEGG_KO Path to the file \code{KO-sum-wide.txt}
#' @param path_KEGG_pathways Path to the file \code{kegg_otu_table.txt}
#' @param abundance_unit One of \code{c("TPM", "relAbd")}, specifying if the
#' returned abundances should be in unit TPM (transcripts per million reads
#' mapped) or returned as relative abundances (i.e., TPM divided by one million).
#' @inheritParams readAndPrepare_lookupData
#' 
#' @import checkmate dplyr stringr
#' @importFrom readxl read_excel
#' @export
#' 
#' @return List of data.frames
#' 
readAndPrepare_functionalData <- function(path_KEGG_KO,
                                          path_KEGG_pathways,
                                          path_sampleLookupXlsx,
                                          path_krautLookupXlsx,
                                          abundance_unit = "TPM",
                                          exclude_medicatedParticipants = TRUE,
                                          exclude_sickParticipants) {
  
  checkmate::assert_file(path_KEGG_KO,       extension = "txt")
  checkmate::assert_file(path_KEGG_pathways, extension = "txt")
  checkmate::assert_file(path_sampleLookupXlsx, extension = "xlsx")
  checkmate::assert_file(path_krautLookupXlsx,  extension = "xlsx")
  checkmate::assert_choice(abundance_unit, choices = c("TPM", "relAbd"))
  checkmate::assert_logical(exclude_medicatedParticipants, len = 1)
  checkmate::assert_choice(exclude_sickParticipants, choices = c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper"), null.ok = TRUE)
  
  
  dat_KO           <- read.table(path_KEGG_KO,       header = TRUE)
  dat_PW           <- read.table(path_KEGG_pathways, header = FALSE) # automatic header parsing doesn't work
  header_PW        <- readLines(path_KEGG_pathways, n = 1) %>% strsplit(split = "\t")
  colnames(dat_PW) <- header_PW[[1]]
  dat_lookupList   <- readAndPrepare_lookupData(path_sampleLookupXlsx         = path_sampleLookupXlsx,
                                                path_krautLookupXlsx          = path_krautLookupXlsx,
                                                exclude_medicatedParticipants = exclude_medicatedParticipants,
                                                exclude_sickParticipants      = exclude_sickParticipants)
  dat_lookupStool  <- dat_lookupList$dat_lookupStool
  dat_lookupKraut  <- dat_lookupList$dat_lookupKraut
  
  
  
  # base preparations -------------------------------------------------------
  ### KO
  # delete the chicken sample, it only being a control sample of the laboratory
  col_chicken <- dat_KO %>% colnames() %>% str_which("chicken")
  dat_KO      <- dat_KO %>% select(-all_of(col_chicken))
  
  # simplify the column names
  colnames(dat_KO)[2:ncol(dat_KO)] <- colnames(dat_KO)[2:ncol(dat_KO)] %>% 
    sapply(function(x) { strsplit(x, "\\_")[[1]][6] }) %>% 
    paste0("sample_", .)
  
  
  ### pathways
  # delete the chicken sample, it only being a control sample of the laboratory
  col_chicken <- dat_PW %>% colnames() %>% str_which("chicken")
  dat_PW      <- dat_PW %>% select(-all_of(col_chicken))
  
  # split the taxonomic information
  dat_PW <- dat_PW %>% 
    mutate(taxonomy_level1 = sapply(taxonomy,        function(x) { strsplit(x, split = ";")[[1]][2] %>% gsub("_", " ", .) }, USE.NAMES = FALSE),
           taxonomy_level2 = sapply(taxonomy,        function(x) { strsplit(x, split = ";")[[1]][3] %>% gsub("_", " ", .) }, USE.NAMES = FALSE),
           taxonomy_level3 = sapply(taxonomy,        function(x) { strsplit(x, split = ";")[[1]][4] %>% gsub("_", " ", .) }, USE.NAMES = FALSE),
           taxonomy_level3 = sapply(taxonomy_level3, function(x) { strsplit(x, split = " \\(")[[1]][1] },                    USE.NAMES = FALSE),
           taxonomy_level3 = case_when(taxonomy_level3 == "" ~ "unnamed",
                                       TRUE                  ~ taxonomy_level3),
           across(starts_with("taxonomy"), factor)) %>% 
    select(-taxonomy)
  
  # simplify the column names and reorder columns
  dat_PW <- dat_PW %>%
    dplyr::rename(otu_id = "#OTU ID") %>% 
    mutate(otu_id = toupper(otu_id)) %>% 
    select(otu_id, starts_with("taxonomy"), everything())
  colnames(dat_PW)[5:ncol(dat_PW)] <- colnames(dat_PW)[5:ncol(dat_PW)] %>% 
    sapply(function(x) { strsplit(x, "\\_")[[1]][6] }) %>% 
    paste0("sample_", .)
  
  
  
  # normalize the TPM abundances --------------------------------------------
  if (abundance_unit == "relAbd") {
    ### KO
    dat_KO <- dat_KO %>% 
      mutate(across(starts_with("sample"), function(x) { x / 1000000 }))
    
    ### pathways
    dat_PW <- dat_PW %>% 
      mutate(across(starts_with("sample"), function(x) { x / 1000000 }))
  }
  
  
  
  # separate the stool from the sauerkraut samples ---------------------------
  ### KO
  # separate the stool samples from the sauerkraut samples
  col_firstKrautSample <- which(grepl("sample_508", colnames(dat_KO)))
  dat_KOstool          <- dat_KO %>% select(all_of(1:(col_firstKrautSample - 1)))
  dat_KOkraut          <- dat_KO %>% select(KO, all_of(col_firstKrautSample:ncol(dat_KO)))
  
  
  ### pathways
  # separate the stool samples from the sauerkraut samples
  col_firstKrautSample <- which(grepl("sample_508", colnames(dat_PW)))
  dat_PWstool          <- dat_PW %>% select(all_of(1:(col_firstKrautSample - 1)))
  dat_PWkraut          <- dat_PW %>% select(otu_id, starts_with("taxonomy"), all_of(col_firstKrautSample:ncol(dat_PW)))
  
  
  
  # move the duplicated stool measurements to another dataset ----------------
  sampleIds_dupl <- dat_lookupStool %>% 
    filter(!is.na(id_duplicateSample)) %>% 
    pull(id_duplicateSample) %>% 
    sort()
  
  ### KO
  # split 'dat_KOstool' in two, placing the duplicate samples in an extra dataset
  dat_KOstoolDuplicates <- dat_KOstool %>%
    select(KO, all_of(paste0("sample_", sampleIds_dupl)))
  dat_KOstool           <- dat_KOstool %>%
    select(-all_of(paste0("sample_", sampleIds_dupl)))
  
  
  ### pathways
  # split 'dat_PWstool' in two, placing the duplicate samples in an extra dataset
  dat_PWstoolDuplicates <- dat_PWstool %>%
    select(otu_id, starts_with("taxonomy"), all_of(paste0("sample_", sampleIds_dupl)))
  dat_PWstool           <- dat_PWstool %>%
    select(-all_of(paste0("sample_", sampleIds_dupl)))
  
  
  
  # exclude medicated and / or sick participants ------------------------------
  ### KO
  # main stool dataset
  sample_cols <- colnames(dat_KOstool)[grepl("sample_", colnames(dat_KOstool))]
  cols_toDrop <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  if (length(cols_toDrop) > 0) { 
    dat_KOstool <- dat_KOstool %>% select(-all_of(cols_toDrop))
  }
  
  # duplicates stool dataset
  sample_cols         <- colnames(dat_KOstoolDuplicates)[grepl("sample_", colnames(dat_KOstoolDuplicates))]
  cols_toDrop         <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  if (length(cols_toDrop) > 0) {
    dat_KOstoolDuplicates <- dat_KOstoolDuplicates %>% select(-cols_toDrop)
  }
  
  
  ### pathways
  # main stool dataset
  sample_cols <- colnames(dat_PWstool)[grepl("sample_", colnames(dat_PWstool))]
  cols_toDrop <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  if (length(cols_toDrop) > 0) { 
    dat_PWstool <- dat_PWstool %>% select(-all_of(cols_toDrop))
  }
  
  # duplicates stool dataset
  sample_cols         <- colnames(dat_PWstoolDuplicates)[grepl("sample_", colnames(dat_PWstoolDuplicates))]
  cols_toDrop         <- sample_cols[!(sample_cols %in% paste0("sample_", dat_lookupStool$sample_id))]
  if (length(cols_toDrop) > 0) {
    dat_PWstoolDuplicates <- dat_PWstoolDuplicates %>% select(-cols_toDrop)
  }
  
  
  
  # final preparations -------------------------------------------------------
  dat_list <- list("datKO_stoolSamples"                = dat_KOstool           %>% as_tibble(),
                   "datKO_stoolDuplicateSamples"       = dat_KOstoolDuplicates %>% as_tibble(),
                   "datKO_krautSamples"                = dat_KOkraut           %>% as_tibble(),
                   "datPathways_stoolSamples"          = dat_PWstool           %>% as_tibble(),
                   "datPathways_stoolDuplicateSamples" = dat_PWstoolDuplicates %>% as_tibble(),
                   "datPathways_krautSamples"          = dat_PWkraut           %>% as_tibble(),
                   "dat_lookupStool"                   = dat_lookupStool,
                   "dat_lookupKraut"                   = dat_lookupKraut)
  
  return(dat_list)
  
}