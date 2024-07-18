
#' Calculate Bray-Curtis dissimilarities and (un)weighted UniFrac distances
#' 
#' Based on a \code{TreeSummarizedExperiment} object and the relative
#' abundances contained in it this function calculates the pairwise Bray-Curtis
#' dissimilarities (based on function \code{vegan::vegdist}) and the pairwise
#' (un)weighted UniFrac distances (based on function \code{mia::calculateUnifrac})
#' between all samples, and returns them in form of a list of two data frames,
#' comprising the intra-personal distances and the inter-personal distances.
#' 
#' @param tse Object of class \code{TreeSummarizedExperiment}
#' @param assay_relAbd Character name of the assay in the \code{tse} object that
#' contains relative abundances, as a base for calculating the Bray-Curtis
#' dissimiliraties.
#' 
#' @import checkmate dplyr SingleCellExperiment
#' @importFrom vegan vegdist
#' @importFrom mia calculateUnifrac
#' @export
#' 
#' @return List of two data.frames
#' 
calculate_distances <- function(tse, assay_relAbd) {
  
  checkmate::assert_class(tse, "TreeSummarizedExperiment")

  
  # retrieve the sample dataset from the tse object
  dat <- tse %>% colData() %>% as.data.frame()
  
  
  
  # calculate dissimilarities and distances ---------------------------------
  # Bray-Curtis dissimilarities
  distMatrix_BC  <- tse %>% assay(assay_relAbd) %>% t() %>% vegan::vegdist(method = "bray")
  
  # weighted UniFrac distances
  distMatrix_wUF <- tse %>% mia::calculateUnifrac(assay.type = "relAbd", weighted = TRUE, normalized = FALSE)
  
  # unweighted UniFrac distances
  distMatrix_uUF <- tse %>% mia::calculateUnifrac(assay.type = "relAbd", weighted = FALSE)
  
  
  # create long dataframes from the distance matrices
  dat_BC <- distMatrix_BC %>% as.matrix() %>% as.data.frame() %>% 
    mutate(sample_t1 = row.names(.)) %>%
    pivot_longer(-sample_t1, names_to = "sample_t2", values_to = "brayCurtis")
  dat_wUF <- distMatrix_wUF %>% as.matrix() %>% as.data.frame() %>% 
    mutate(sample_t1 = row.names(.)) %>%
    pivot_longer(-sample_t1, names_to = "sample_t2", values_to = "uniFrac_weighted")
  dat_uUF <- distMatrix_uUF %>% as.matrix() %>% as.data.frame() %>% 
    mutate(sample_t1 = row.names(.)) %>%
    pivot_longer(-sample_t1, names_to = "sample_t2", values_to = "uniFrac_unweighted")
  
  
  
  # create the base datasets ------------------------------------------------
  partID_vec <- dat$participant_id %>% unique() %>% as.character() %>% sort()
  t_vec      <- dat$timepoint      %>% unique() %>% as.character() %>% sort()
  
  # intra-personal base comparison dataset
  tDat_intraPersonal <- expand.grid(index_1 = 1:length(t_vec),
                                          index_2 = 1:length(t_vec)) %>% 
    filter(index_2 > index_1) %>% 
    arrange(index_1) %>% 
    mutate(timepoint_1 = t_vec[index_1],
           timepoint_2 = t_vec[index_2])
  results_intraPersonal <- data.frame(participant_id = partID_vec %>% rep(each = nrow(tDat_intraPersonal)),
                                      timepoint_1    = tDat_intraPersonal$timepoint_1 %>% rep(times = length(partID_vec)),
                                      timepoint_2    = tDat_intraPersonal$timepoint_2 %>% rep(times = length(partID_vec))) %>% 
    mutate(comparison = paste(timepoint_1, "vs.", timepoint_2))
  
  # inter-personal base comparison dataset
  idDat_interPersonal <- expand.grid(index_1 = 1:length(partID_vec),
                                     index_2 = 1:length(partID_vec)) %>% 
    filter(index_2 > index_1) %>% 
    arrange(index_1) %>% 
    mutate(participantID_1 = partID_vec[index_1],
           participantID_2 = partID_vec[index_2])
  results_interPersonal <- data.frame(timepoint       = t_vec %>% rep(each = nrow(idDat_interPersonal)),
                                      participantID_1 = idDat_interPersonal$participantID_1 %>% rep(times = length(t_vec)),
                                      participantID_2 = idDat_interPersonal$participantID_2 %>% rep(times = length(t_vec))) %>% 
    mutate(comparison = paste(participantID_1, "vs.", participantID_2))
  
  
  # add the sample names to the datasets
  lookup_dat <- dat %>% 
    mutate(partID_t = paste(participant_id, timepoint),
           sample   = row.names(.)) %>% 
    select(partID_t, sample)
  row.names(lookup_dat) <- NULL
  
  results_intraPersonal <- results_intraPersonal %>% 
    mutate(partID_t1 = paste(participant_id, timepoint_1),
           partID_t2 = paste(participant_id, timepoint_2)) %>% 
    dplyr::left_join(lookup_dat, by = c("partID_t1" = "partID_t")) %>% 
    dplyr::rename(sample_t1 = sample) %>% 
    dplyr::left_join(lookup_dat, by = c("partID_t2" = "partID_t")) %>% 
    dplyr::rename(sample_t2 = sample) %>% 
    select(-starts_with("partID"))
  
  results_interPersonal <- results_interPersonal %>% 
    mutate(partID_t1 = paste(participantID_1, timepoint),
           partID_t2 = paste(participantID_2, timepoint)) %>% 
    dplyr::left_join(lookup_dat, by = c("partID_t1" = "partID_t")) %>% 
    dplyr::rename(sample_t1 = sample) %>% 
    dplyr::left_join(lookup_dat, by = c("partID_t2" = "partID_t")) %>% 
    dplyr::rename(sample_t2 = sample) %>% 
    select(-starts_with("partID"))
  
  
  # remove rows with timepoints that we have no data for
  results_intraPersonal <- results_intraPersonal %>% 
    filter(!is.na(sample_t1) & !is.na(sample_t2))
  results_interPersonal <- results_interPersonal %>% 
    filter(!is.na(sample_t1) & !is.na(sample_t2))
  
  
  
  # add the dissimilarities and distances to the dataset --------------------
  results_intraPersonal <- results_intraPersonal %>% 
    dplyr::left_join(dat_BC,  by = c("sample_t1", "sample_t2")) %>% 
    dplyr::left_join(dat_wUF, by = c("sample_t1", "sample_t2")) %>% 
    dplyr::left_join(dat_uUF, by = c("sample_t1", "sample_t2"))
  results_interPersonal <- results_interPersonal %>% 
    dplyr::left_join(dat_BC,  by = c("sample_t1", "sample_t2")) %>% 
    dplyr::left_join(dat_wUF, by = c("sample_t1", "sample_t2")) %>% 
    dplyr::left_join(dat_uUF, by = c("sample_t1", "sample_t2"))
  
  
  
  # final adjustments -------------------------------------------------------
  # intra-personal: also categorize the timepoint comparisons (also by adding intervention information)
  results_intraPersonal <- results_intraPersonal %>% 
    mutate(partID_t2 = paste0(participant_id, "_", timepoint_2)) %>% 
    dplyr::left_join(dat %>% select(part_time, intervention, group, age, gender, bmi_t0), by = c("partID_t2" = "part_time")) %>% 
    mutate(comparison_group = case_when(comparison %in% c("T1 vs. T3", "T1 vs. T5", "T3 vs. T5")                    ~ "baselines",
                                        comparison %in% c("T1 vs. T2", "T3 vs. T4") & intervention == "Fresh"       ~ "fresh intervention",
                                        comparison %in% c("T1 vs. T2", "T3 vs. T4") & intervention == "Pasteurized" ~ "pasteurized intervention",
                                        TRUE                                                                        ~ "other")) %>% 
    select(participant_id, group, timepoint_1, timepoint_2, comparison,
           comparison_group, brayCurtis, uniFrac_weighted, uniFrac_unweighted,
           age, gender, bmi_t0)
  
  # inter-personal
  results_interPersonal <- results_interPersonal %>% 
    mutate(timepoint_group = case_when(timepoint %in% c("T1", "T3", "T5") ~ "baselines",
                                       TRUE                               ~ "interventions")) %>% 
    select(timepoint, timepoint_group, participantID_1, participantID_2, comparison,
           brayCurtis, uniFrac_weighted, uniFrac_unweighted)
  
  
  results_list <- list("datDist_intraPersonal" = results_intraPersonal,
                       "datDist_interPersonal" = results_interPersonal)
  
  return(results_list)
}
