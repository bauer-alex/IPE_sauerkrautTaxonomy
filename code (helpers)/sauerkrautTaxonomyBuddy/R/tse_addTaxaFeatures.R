
#' Add different taxa features to a tse object
#' 
#' This function calculates the following features based on an object \code{tse}
#' with class \code{TreeSummarizedExperiment} and adds them as additional
#' columns in \code{colData(tse)}:
#' 1. relative abundance and rel. abundance change (after the intervention) of
#' L. paracasei
#' 2. dominant (= most abundant) species, genus, etc.
#' 3. relative abundance of the dominant species, genus, etc.
#' 4. species, genus, etc. with the biggest growth and loss
#' 
#' @param tse Object of class \code{TreeSummarizedExperiment}
#' @param assay_relAbd Character name of the assay in the \code{tse} object that
#' contains relative abundances.
#' 
#' @import checkmate dplyr tidyr SummarizedExperiment
#' @importFrom mia addPerSampleDominantFeatures
#' @export
#' 
#' @return \code{tse} object with additional columns in \code{colData(tse)}
#' 
tse_addTaxaFeatures <- function(tse, assay_relAbd) {
  
  checkmate::assert_class(tse, classes = "TreeSummarizedExperiment")
  checkmate::assert_choice(assay_relAbd, choices = SummarizedExperiment::assayNames(tse))
  
  
  # preparations ------------------------------------------------------------
  # perform all analyses for all taxonomic hierarchy levels
  hierarchies <- c("phylum", "class", "order", "family", "genus", "species")
  
  # save some information in separate objects, for easier access
  dat_relAbd       <- tse %>% assay(assay_relAbd) %>% as.data.frame()
  dat_row          <- tse %>% rowData() %>% as.data.frame()
  dat_col          <- tse %>% colData() %>% as.data.frame()
  samples_vec      <- colnames(dat_relAbd)
  participants_vec <- dat_col$participant_id %>% unique() %>% as.character() %>% sort()
  
  # create a long version of dat_relAbd
  dat_relAbd_long <- dat_relAbd %>% 
    mutate(species = row.names(.) %>% gsub(pattern = "species:", replacement = "")) %>% 
    tidyr::pivot_longer(-species, names_to = "sample_id", values_to = "relAbd") %>% 
    arrange(sample_id, species)
  
  
  
  # calculate the L. paracasei features -------------------------------------
  # add the main fresh sauerkraut species' abundance values to the data
  mat_relAbd                     <- tse %>% assay("relAbd")
  mat_row                        <- mat_relAbd[grep("Lacticaseibacillus_paracasei", row.names(mat_relAbd)),]
  colData(tse)$relAbd_Lparacasei <- unname(mat_row)[match(row.names(colData(tse)), names(mat_row))]
  
  # calculate the change in L. paracasei after the intervention
  dat_LpChange <- colData(tse) %>% 
    as.data.frame() %>% 
    select(participant_id, treatment, intervention, relAbd_Lparacasei) %>% 
    filter(intervention != "FollowUp") %>% 
    arrange(participant_id, intervention, treatment) %>% 
    group_by(participant_id, intervention) %>% 
    summarize(treatment                = "After",
              relAbd_Lparacasei_change = diff(relAbd_Lparacasei)) %>% 
    ungroup()
  
  matchingID_colDat <- paste(colData(tse)$participant_id, colData(tse)$treatment, colData(tse)$intervention)
  matchingID_newDat <- paste(dat_LpChange$participant_id, dat_LpChange$treatment, dat_LpChange$intervention)
  colData(tse)$relAbd_Lparacasei_change <- dat_LpChange$relAbd_Lparacasei_change[match(matchingID_colDat, matchingID_newDat)]
  
  
  
  # calculate other features ------------------------------------------------
  for (hierarchy in hierarchies) {
    
    # 1) dominant species, genus, etc. ----------------------------------------
    tse <- mia::addPerSampleDominantFeatures(x          = tse,
                                             rank       = hierarchy,
                                             name       = paste0("dominant_", hierarchy),
                                             assay_name = assay_relAbd)
  
  
    
    # 2) relative abundance of the dominant species, genus, etc. --------------
    
    # save colData(tse) again as separate object, for easier handling
    dat_col          <- tse %>% colData() %>% as.data.frame()
    
    # create a dataset that contains the list of species in the most dominant
    # feature of every sample
    dat_sampleSpecies_list <- lapply(samples_vec, function(x) {
      
      # extract the dominant feature for the current sample in the current hierarchy
      dominant_feature <- dat_col %>% 
        filter(row.names(.) == x) %>% 
        pull(paste0("dominant_", hierarchy))
      
      # extract all species that form this feature
      species_vec <- dat_row %>% 
        filter(.[[hierarchy]] == dominant_feature) %>% 
        pull(species) %>% 
        as.character()
      
      # return a data.frame with the sample and its species
      data.frame(sample_id = x,
                 species   = species_vec)
    })
    dat_sampleSpecies <- dat_sampleSpecies_list %>% dplyr::bind_rows()
    
    
    # add the relative abundances to the dataset
    dat_sampleSpecies_relAbd <- dat_sampleSpecies %>% 
      dplyr::left_join(dat_relAbd_long, by = c("sample_id", "species"))
    
    # aggregate these relative abundances to get the relative abundance for the
    # whole dominant feature
    dat_dom_relAbd <- dat_sampleSpecies_relAbd %>% 
      group_by(sample_id) %>% 
      summarize(relAbd = sum(relAbd))
    
    # add the data to colData(tse)
    colData(tse)$relAbd_dom <- dat_dom_relAbd$relAbd[match(row.names(colData(tse)), dat_dom_relAbd$sample_id)]
    colnames(colData(tse))[ncol(colData(tse))] <- paste0("dominant_", hierarchy, "_relAbd")
    
    
    
    # 3) species, genus, etc. with the highest win / loss during the intervention --------
    # create one dataset containing all relative abundances and experiment information
    dat_col_toJoin <- dat_col %>% 
      mutate(sample_id = row.names(.)) %>% 
      select(participant_id, sample_id, treatment, intervention)
    
    dat_base <- dat_relAbd_long %>% 
      dplyr::left_join(dat_col_toJoin, by = "sample_id") %>% 
      filter(intervention != "FollowUp")
    
    # ensure proper sorting of the data
    dat_base <- dat_base %>% arrange(participant_id, intervention, treatment)
    
    # aggregate species' abundances
    if (hierarchy != "species") {
      # add the hierarchy information to the data
      dat_base <- dat_base %>% dplyr::left_join(dat_row, by = "species")
      
      # rename the hierarchy column for easier handling
      colnames(dat_base)[colnames(dat_base) == hierarchy] <- "hierarchy"
      
      # summarize over the hierarchy
      dat_base <- dat_base %>% 
        group_by(participant_id, intervention, treatment, hierarchy) %>% 
        summarize(relAbd = sum(relAbd)) %>% 
        ungroup()
    }
    
    # consistently name the relevant species, genus, etc. column
    colnames(dat_base)[colnames(dat_base) %in% c("species", "hierarchy")] <- "taxonomy"
    
    # per person, for each species, genus, etc., calculate the change in each intervention
    dat_diff <- dat_base %>% 
      group_by(participant_id, intervention, taxonomy) %>% 
      summarize(treatment = "After",
                change    = diff(relAbd)) %>% 
      ungroup()
    
    # only keep the biggest changes
    dat_changes <- dat_diff %>% 
      group_by(participant_id, intervention) %>% 
      filter(change %in% c(min(change), max(change))) %>% 
      ungroup() %>% 
      mutate(change_type = case_when(change > 0 ~ "maxGrowth",
                                     change < 0 ~ "maxLoss",
                                     TRUE       ~ NA_character_)) %>% 
      tidyr::pivot_wider(id_cols = c("participant_id", "intervention", "treatment"),
                         names_from = "change_type", values_from = c("change", "taxonomy"))
    
    # reorder and properly name the columns
    dat_changes <- dat_changes %>% 
      select(participant_id, intervention, treatment, taxonomy_maxGrowth, change_maxGrowth, taxonomy_maxLoss, change_maxLoss)
    
    # add the data to colData(tse)
    matchingID_colDat <- paste(colData(tse)$participant_id, colData(tse)$treatment, colData(tse)$intervention)
    matchingID_newDat <- paste(dat_changes$participant_id, dat_changes$treatment, dat_changes$intervention)
    colData(tse)$taxonomy_maxGrowth <- dat_changes$taxonomy_maxGrowth[match(matchingID_colDat, matchingID_newDat)]
    colData(tse)$change_maxGrowth   <- dat_changes$change_maxGrowth[match(  matchingID_colDat, matchingID_newDat)]
    colData(tse)$taxonomy_maxLoss   <- dat_changes$taxonomy_maxLoss[match(  matchingID_colDat, matchingID_newDat)]
    colData(tse)$change_maxLoss     <- dat_changes$change_maxLoss[match(    matchingID_colDat, matchingID_newDat)]
    colnames(colData(tse))[ncol(colData(tse)) - 3] <- paste0(hierarchy, "_maxGrowth")
    colnames(colData(tse))[ncol(colData(tse)) - 2] <- paste0(hierarchy, "_maxGrowth_growth")
    colnames(colData(tse))[ncol(colData(tse)) - 1] <- paste0(hierarchy, "_maxLoss")
    colnames(colData(tse))[ncol(colData(tse))]     <- paste0(hierarchy, "_maxLoss_loss")
  }
  
  
  
  # return tse object -------------------------------------------------------
  return(tse)
}
