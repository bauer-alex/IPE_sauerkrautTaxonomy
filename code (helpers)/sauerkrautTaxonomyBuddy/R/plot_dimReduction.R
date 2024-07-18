
#' Plot the results of an estimated dimension reduction method
#' 
#' This function was initially based on the function \code{plot_ordination}
#' written by Virginie (which has nothing to do with the equally named function
#' from the \code{phyloseq} package), but was substantially thinned out by
#' dropping some of its optional functionalities.
#' 
#' @param tse Object of class \code{TreeSummarizedExperiment}
#' @param dimReduction_name Character name of a dimension reduction technique
#' result stored in the \code{tse} object. The results was usually stored in the
#' \code{tse} object by calling a function like \link[scater]{runMDS}.
#' @param assay_relAbd Character name of the assay in the \code{tse} object that
#' contains relative abundances.
#' @param var_participantID Character name of the participant id variable in
#' \code{colData(tse)}. Defaults to \code{"participant_id"}.
#' @param var_timepoint Character name of the timepoint variable in
#' \code{colData(tse)}, containing the indicator at which discrete timepoint
#' in the (crossover) study design the measuremtn was taken. Defaults to
#' \code{"timepoint"}.
#' @param var_group Character name of the group variable in \code{colData(tse)}.
#' @param var_interTreat Character name of the combined intervention-treament
#' variable in \code{colData(tse)}, containing values like
#' \code{"Baseline.Condition1"}, \code{"FollowUp.Condition1"}, etc.
#' @param taxRank_dominance Character name of the taxonomic rank for which the
#' dominant bacteria type should be plotted. Defaults to \code{"family"}.
#' @param taxRank_topNfocus Optional number specifying the number of taxa
#' categories that should be highlighted. The 'taxRank_topNfocus' most frequent
#' categories will be highlighted and all other categories will be grayed out.
#' Currently limited to maximum 13.
#' @param covars_categorical,covars_metric Optional vectors of variable names
#' of \code{colData(tse)} for each of which one extra plot should be created
#' and returned as a separate element in the output list.
#' 
#' @import checkmate dplyr mia SingleCellExperiment
#' @importFrom mia taxonomyRanks
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment assayNames
#' @export
#' 
#' @return Named list of ggplot objects
#' 
plot_dimReduction <- function(tse, dimReduction_name, assay_relAbd,
                              var_participantID  = "participant_id",
                              var_timepoint      = "timepoint",
                              var_group          = "group",
                              var_interTreat     = "inter_treat",
                              taxRank_dominance  = "family",
                              taxRank_topNfocus  = NULL,
                              covars_categorical = NULL,
                              covars_metric      = NULL) {
  
  checkmate::assert_class(tse, "TreeSummarizedExperiment")
  checkmate::assert_choice(dimReduction_name,  choices = SingleCellExperiment::reducedDimNames(tse))
  checkmate::assert_choice(assay_relAbd,       choices = assayNames(tse))
  checkmate::assert_choice(var_participantID,  choices = colData(tse) %>% colnames())
  checkmate::assert_choice(var_timepoint,      choices = colData(tse) %>% colnames())
  checkmate::assert_choice(var_group,          choices = colData(tse) %>% colnames())
  checkmate::assert_choice(var_interTreat,     choices = colData(tse) %>% colnames())
  checkmate::assert_choice(taxRank_dominance,  mia::taxonomyRanks(tse))
  checkmate::assert_number(taxRank_topNfocus,  lower = 1, upper = 13, null.ok = TRUE)
  checkmate::assert_subset(covars_categorical, choices = colData(tse) %>% colnames())
  checkmate::assert_subset(covars_metric,      choices = colData(tse) %>% colnames())
  
  
  # prepare the tse object --------------------------------------------------
  # extract and add different bacteria statistics to the data
  tse <- tse %>% tse_addTaxaFeatures(assay_relAbd = assay_relAbd)
  
  # calculate the total abundance  
  colData(tse)$Total_abundance  <- colSums(assays(tse)[[assay_relAbd]])
  
  
  
  # create a dataset for plotting -------------------------------------------
  dat <- as.data.frame(cbind(reducedDim(tse, dimReduction_name), colData(tse)))

  # rename some variables of the hierarchy of interest for easier handling
  colnames(dat)[colnames(dat) == paste0("dominant_", taxRank_dominance)]  <- "dominant_feature"
  colnames(dat)[colnames(dat) == paste0(taxRank_dominance, "_maxGrowth")] <- "taxonomy_highestGrowth"
  colnames(dat)[colnames(dat) == paste0(taxRank_dominance, "_maxLoss")]   <- "taxonomy_highestLoss"
  
  # make the species / genus / etc. labels a bit nicer
  dat <- dat %>% mutate(dominant_feature = gsub(pattern = "_", replacement = " ", x = dominant_feature))
  
  # sort the dominant / highest growth / highest loss taxa levels by frequency,
  # and add the relative frequency to the category labels
  taxDom_tab       <- dat$dominant_feature       %>% table() %>% sort(decreasing = TRUE) %>% prop.table() %>% as.data.frame() %>% dplyr::rename(Level = ".") %>% 
    mutate(Freq_label = paste0("(", 100 * round(Freq, 3), "%)"),
           Label      = paste(Level, Freq_label, sep = " "))
  taxMaxGrowth_tab <- dat$taxonomy_highestGrowth %>% table() %>% sort(decreasing = TRUE) %>% prop.table() %>% as.data.frame() %>% dplyr::rename(Level = ".") %>% 
    mutate(Freq_label = paste0("(", 100 * round(Freq, 3), "%)"),
           Label      = paste(Level, Freq_label, sep = " "))
  taxMaxLoss_tab   <- dat$taxonomy_highestLoss   %>% table() %>% sort(decreasing = TRUE) %>% prop.table() %>% as.data.frame() %>% dplyr::rename(Level = ".") %>% 
    mutate(Freq_label = paste0("(", 100 * round(Freq, 3), "%)"),
           Label      = paste(Level, Freq_label, sep = " "))
  dat <- dat %>% 
    mutate(dominant_feature       = factor(dominant_feature,       levels = taxDom_tab$Level,       labels = taxDom_tab$Label),
           taxonomy_highestGrowth = factor(taxonomy_highestGrowth, levels = taxMaxGrowth_tab$Level, labels = taxMaxGrowth_tab$Label),
           taxonomy_highestLoss   = factor(taxonomy_highestLoss,   levels = taxMaxLoss_tab$Level,   labels = taxMaxLoss_tab$Label))
  
  # ensure consistent column names
  if (!all(colnames(dat)[1:2] == c("V1", "V2"))) {
    if (grepl("PCA",  dimReduction_name)) dat <- dplyr::rename(dat, V1 = PC1,   V2 = PC2)
    if (grepl("tSNE", dimReduction_name)) dat <- dplyr::rename(dat, V1 = TSNE1, V2 = TSNE2)
    if (grepl("UMAP", dimReduction_name)) dat <- dplyr::rename(dat, V1 = UMAP1, V2 = UMAP2)
  }
  
  # name the central variables, for easier handling in dplyr and ggplot
  colnames(dat)[colnames(dat) == var_participantID] <- "participant_id"
  colnames(dat)[colnames(dat) == var_timepoint]     <- "timepoint"
  colnames(dat)[colnames(dat) == var_group]         <- "group"
  colnames(dat)[colnames(dat) == var_interTreat]    <- "inter_treat"
  
  # make the 'inter_treat' labels a bit nicer
  dat <- dat %>% mutate(inter_treat = case_when(inter_treat == "Baseline.Fresh" ~ "Fresh (baseline)",
                                                inter_treat == "After.Fresh"    ~ "Fresh (intervention)",
                                                inter_treat == "Baseline.Pasteurized" ~ "Pasteurized (baseline)",
                                                inter_treat == "After.Pasteurized"    ~ "Pasteurized (intervention)",
                                                inter_treat == "Baseline.FollowUp"    ~ "Follow-up"),
                        inter_treat = factor(inter_treat, levels = c("Fresh (baseline)", "Fresh (intervention)",
                                                                     "Pasteurized (baseline)", "Pasteurized (intervention)",
                                                                     "Follow-up")))
  
  
  dat_participants    <- dat %>% filter(timepoint == "T2")
  dat_participants$V2 <- dat_participants$V2 - (IQR(dat_participants$V2) / 15)
  
  
  
  # prepare plot settings ---------------------------------------------------
  # x and y axis labels, potentially with explained variations
  dimReduction_attributes <- tse %>% SingleCellExperiment::reducedDim(dimReduction_name) %>% attributes() %>% names()
  if ("eig" %in% dimReduction_attributes) {
    
    eigenvalues <- tse %>% SingleCellExperiment::reducedDim(dimReduction_name) %>% attr("eig")
    evalues_rel <- eigenvalues / sum(eigenvalues[eigenvalues > 0])          
    x_lab       <- paste("Dimension 1 (", round(100 * evalues_rel[[1]],1), "%", ")", sep = "")
    y_lab       <- paste("Dimension 2 (", round(100 * evalues_rel[[2]],1), "%", ")", sep = "")
    
  } else {
    x_lab <- "Dimension 1"
    y_lab <- "Dimension 2"
  }
  
  # shapes of timepoint points
  shape_vector <- c(1, 16, 2, 17, 3, 4:15)[1:length(unique(dat$timepoint))]
  
  # color vector for most dominant taxa categories
  if (!is.null(taxRank_topNfocus)) {
    
    n_dominantTaxa  <- dat$dominant_feature       %>% unique() %>% length()
    n_maxGrowthTaxa <- dat$taxonomy_highestGrowth %>% unique() %>% length()
    n_maxLossTaxa   <- dat$taxonomy_highestLoss   %>% unique() %>% length()
    
    # define available color vectors
    base_colors    <- RColorBrewer::brewer.pal(n = min(8, taxRank_topNfocus), name = "Dark2")
    further_colors <- c("firebrick4", "dodgerblue", "dodgerblue4", "black", "pink")
    
    # define the final color vector for the focused categories
    if (taxRank_topNfocus > length(base_colors)) {
      base_colors <- c(base_colors, further_colors[1:(taxRank_topNfocus - 8)])
    }
    
    # construct the final color vector
    taxDominant_colors  <- c(base_colors, rep(gray(0.8), max(0, n_dominantTaxa  - taxRank_topNfocus)))
    taxMaxGrowth_colors <- c(base_colors, rep(gray(0.8), max(0, n_maxGrowthTaxa - taxRank_topNfocus)))
    taxMaxLoss_colors   <- c(base_colors, rep(gray(0.8), max(0, n_maxLossTaxa   - taxRank_topNfocus)))
  }
  
  
  
  # create plots ------------------------------------------------------------
  # create a base plot
  gg_basePlot <- ggplot() + 
    geom_point(data = dat, aes(x = V1, y = V2), size = 2, color = "gray20", alpha = 0.5) +  
    ggtitle(dimReduction_name) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # plot where the measurements are differentiated by the abundance of L. paracasei
  gg_Lparacasei_byIntervention <- dat %>% 
    filter(intervention != "FollowUp") %>% 
    ggplot(aes(x = V1, y = V2, color = relAbd_Lparacasei, shape = treatment)) +
    geom_point(size = 3, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    scale_shape_manual(values = shape_vector[c(2,4)]) +
    ggtitle("Measurements by the relative abundance of L. paracasei", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab) +
    scale_color_gradient2("L. paracasei rel. abundance",
                          mid = "gray95", high = "dodgerblue4", midpoint = 0)
  
  
  # ... same plot with categorized values
  gg_Lparacasei_byIntervention_cat <- dat %>% 
    mutate(relAbd_Lparacasei_cat = case_when(is.na(relAbd_Lparacasei)  ~ NA_character_,
                                             relAbd_Lparacasei >= 0.2  ~ ">= 20%",
                                             relAbd_Lparacasei >= 0.1  ~ "[10%, 20%)",
                                             relAbd_Lparacasei >= 0.03 ~ "[3%, 10%)",
                                             relAbd_Lparacasei >= 0.01 ~ "[1%, 3%)",
                                             TRUE                      ~ "< 1%"),
           relAbd_Lparacasei_cat = factor(relAbd_Lparacasei_cat, levels = c("< 1%", "[1%, 3%)", "[3%, 10%)", "[10%, 20%)", ">= 20%"))) %>% 
    filter(intervention != "FollowUp") %>% 
    ggplot(aes(x = V1, y = V2, color = relAbd_Lparacasei_cat, shape = treatment)) +
    geom_point(size = 3, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    scale_shape_manual(values = shape_vector[c(2,4)]) +
    ggtitle("Measurements by the relative abundance of L. paracasei - categorized", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab) +
    scale_color_manual("L. paracasei rel. abundance", values = c("gray95", "dodgerblue1", "dodgerblue4", "firebrick1", "firebrick4"))
  
  
  # plot where the measurements are differentiated by the abundance change of L. paracasei
  gg_LparacaseiChange_byIntervention <- dat %>% 
    filter(intervention != "FollowUp",
           treatment    == "After") %>% 
    ggplot(aes(x = V1, y = V2, color = relAbd_Lparacasei_change)) +
    geom_point(size = 3, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle("Measurements by the relative abundance change of L. paracasei", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab) +
    scale_color_gradient2("L. paracasei rel. abundance change",
                          label = function(x) { paste(round(100 * x, 0), "PP") },
                          low = "firebrick4", mid = "gray95", high = "dodgerblue4")
  
  # ... same plot with categorized values
  gg_LparacaseiChange_byIntervention_cat <- dat %>% 
    mutate(relAbd_Lparacasei_change_cat = case_when(is.na(relAbd_Lparacasei_change)  ~ NA_character_,
                                                    relAbd_Lparacasei_change >= 0.2  ~ ">= 20 PP",
                                                    relAbd_Lparacasei_change >= 0.1  ~ "[10 PP, 20 PP)",
                                                    relAbd_Lparacasei_change >= 0.05 ~ "[5 PP, 10 PP)",
                                                    relAbd_Lparacasei_change >= 0.01 ~ "[1 PP, 5 PP)",
                                                    TRUE                      ~ "< 1 PP"),
           relAbd_Lparacasei_change_cat = factor(relAbd_Lparacasei_change_cat, levels = c("< 1 PP", "[1 PP, 5 PP)", "[5 PP, 10 PP)", "[10 PP, 20 PP)", ">= 20 PP"))) %>% 
    filter(intervention != "FollowUp",
           treatment    == "After") %>% 
    ggplot(aes(x = V1, y = V2, color = relAbd_Lparacasei_change_cat)) +
    geom_point(size = 3, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle("Measurements by the relative abundance change of L. paracasei - categorized", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab) +
    scale_color_manual("L. paracasei rel. abundance change", values = c("gray95", "dodgerblue1", "dodgerblue4", "firebrick1", "firebrick4"))
  
  
  # plot where the measurements are differentiated by the most abundant taxa
  gg_dominantTaxa <- ggplot(dat) + 
    geom_point(aes(x = V1, y = V2, color = dominant_feature), size = 5, alpha = 0.7, stroke = 0) +  
    ggtitle(label    = paste("Measurements by most dominant", taxRank_dominance),
            subtitle = paste(dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$dominant_feature)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_dominantTaxa <- if (is.null(taxRank_topNfocus)) {
    gg_dominantTaxa + scale_color_manual(paste("dominant", taxRank_dominance, "(share of dominated microbiome samples)"))
  } else {
    gg_dominantTaxa + scale_color_manual(paste("dominant", taxRank_dominance, "(share of dominated microbiome samples)"),
                                         values = taxDominant_colors)
  }

  
  # ... same plot differentiated by age
  gg_dominantTaxa_byAge <- dat %>% 
    mutate(age_cat = case_when(age >= 50 ~ "age >= 50",
                               age >= 30 ~ "age [30, 50)",
                               TRUE      ~ "age < 30"),
           age_cat = factor(age_cat, levels = c("age < 30", "age [30, 50)", "age >= 50"))) %>% 
    ggplot(aes(x = V1, y = V2, color = dominant_feature)) + 
    geom_point(size = 2, alpha = 0.7, stroke = 1.5) +  
    facet_wrap(~ age_cat) +
    ggtitle(label    = paste("Measurements by most dominant", taxRank_dominance),
            subtitle = paste(dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$dominant_feature)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_dominantTaxa_byAge <- if (is.null(taxRank_topNfocus)) {
    gg_dominantTaxa_byAge + scale_color_manual(taxRank_dominance)
  } else {
    gg_dominantTaxa_byAge + scale_color_manual(taxRank_dominance, values = taxDominant_colors)
  }
  
  
  # plot where the measurements are differentiated by the taxa with highest growth, for fresh sauerkraut
  gg_maxGrowthTaxa_fresh <- dat %>% 
    filter(treatment    == "After",
           intervention == "Fresh") %>% 
    ggplot(aes(x = V1, y = V2, color = taxonomy_highestGrowth)) + 
    geom_point(size = 2, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle(label    = paste("Measurements by", taxRank_dominance, "with highest growth"),
            subtitle = paste("Fresh sauerkraut -", dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$taxonomy_highestGrowth)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_maxGrowthTaxa_fresh <- if (is.null(taxRank_topNfocus)) {
    gg_maxGrowthTaxa_fresh + scale_color_manual(taxRank_dominance)
  } else {
    gg_maxGrowthTaxa_fresh + scale_color_manual(taxRank_dominance, values = taxMaxGrowth_colors)
  }
  
  
  # plot where the measurements are differentiated by the taxa with highest growth, for pasteurized sauerkraut
  gg_maxGrowthTaxa_past <- dat %>% 
    filter(treatment    == "After",
           intervention == "Pasteurized") %>% 
    ggplot(aes(x = V1, y = V2, color = taxonomy_highestGrowth)) + 
    geom_point(size = 2, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle(label    = paste("Measurements by", taxRank_dominance, "with highest growth"),
            subtitle = paste("Pasteurized sauerkraut -", dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$taxonomy_highestGrowth)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_maxGrowthTaxa_past <- if (is.null(taxRank_topNfocus)) {
    gg_maxGrowthTaxa_past + scale_color_manual(taxRank_dominance)
  } else {
    gg_maxGrowthTaxa_past + scale_color_manual(taxRank_dominance, values = taxMaxGrowth_colors)
  }
  
  
  # plot where the measurements are differentiated by the taxa with highest loss, for fresh sauerkraut
  gg_maxLossTaxa_fresh <- dat %>% 
    filter(treatment    == "After",
           intervention == "Fresh") %>% 
    ggplot(aes(x = V1, y = V2, color = taxonomy_highestLoss)) + 
    geom_point(size = 2, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle(label    = paste("Measurements by", taxRank_dominance, "with highest loss"),
            subtitle = paste("Fresh sauerkraut -", dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$taxonomy_highestLoss)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_maxLossTaxa_fresh <- if (is.null(taxRank_topNfocus)) {
    gg_maxLossTaxa_fresh + scale_color_manual(taxRank_dominance)
  } else {
    gg_maxLossTaxa_fresh + scale_color_manual(taxRank_dominance, values = taxMaxLoss_colors)
  }
  
  
  # plot where the measurements are differentiated by the taxa with highest loss, for pasteurized sauerkraut
  gg_maxLossTaxa_past <- dat %>% 
    filter(treatment    == "After",
           intervention == "Pasteurized") %>% 
    ggplot(aes(x = V1, y = V2, color = taxonomy_highestLoss)) + 
    geom_point(size = 2, alpha = 0.7, stroke = 1.5) +
    facet_wrap(~ intervention) +
    ggtitle(label    = paste("Measurements by", taxRank_dominance, "with highest loss"),
            subtitle = paste("Pasteurized sauerkraut -", dimReduction_name, "\n",
                             ifelse(is.null(taxRank_topNfocus) || taxRank_topNfocus >= length(unique(dat$taxonomy_highestLoss)),
                                    paste(taxRank_dominance, "categories sorted by frequency"),
                                    paste("focus on the top", taxRank_topNfocus, taxRank_dominance, "categories")))) +
    xlab(x_lab) + ylab(y_lab)
  # add color scale
  gg_maxLossTaxa_past <- if (is.null(taxRank_topNfocus)) {
    gg_maxLossTaxa_past + scale_color_manual(taxRank_dominance)
  } else {
    gg_maxLossTaxa_past + scale_color_manual(taxRank_dominance, values = taxMaxLoss_colors)
  }
  
  
  # plot where the measurements are differentiated by experiment groups
  gg_groups <- ggplot(dat, aes(x = V1, y = V2, color = group)) + 
    geom_point(aes(shape = timepoint), size = 2, alpha = 0.7, stroke = 1.5) + 
    stat_ellipse() +
    scale_shape_manual(values = shape_vector) +
    scale_color_brewer("experimental group", palette = "Dark2") +
    ggtitle("Measurements by experimental groups", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # plot where the measurements are differentiated by their underlying intervention times
  gg_interventionTimes <- ggplot(dat, aes(x = V1, y = V2, color = inter_treat)) + 
    geom_point(aes(shape = timepoint), size = 2, stroke = 1.5) + 
    stat_ellipse() +
    scale_shape_manual(values = shape_vector) +
    scale_color_brewer("intervention time", palette = "Paired") +
    ggtitle("Measurements by intervention times", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab)
  
  # same plot, without the timepoint shapes
  gg_interventionTimes_v2 <- dat %>% 
    filter(intervention != "FollowUp") %>%
    ggplot(aes(x = V1, y = V2, color = inter_treat)) + 
    geom_point(size = 3, stroke = 0) + 
    stat_ellipse(size = 0.9) +
    scale_color_manual(NULL, values = c("#A6CEE3","#1F78B4","#FDBF6F","#FF7F00")) +
    ggtitle("Measurements by intervention times", subtitle = dimReduction_name) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # plot where each participant's measurements are connected by a polygon to visualize variation
  gg_participantVariation <- ggplot(dat, aes(x = V1, y = V2)) +
    geom_point(aes(shape = timepoint, color = inter_treat), size = 2, alpha = 0.7, stroke = 1.5) + 
    geom_polygon(aes(fill = participant_id), alpha = 0.3, show.legend = FALSE) +
    scale_shape_manual(values = shape_vector) +
    geom_text(data = dat_participants, aes(label = participant_id, color = participant_id), check_overlap = FALSE, fontface = "bold") +
    guides(color = "none") + 
    ggtitle(label    = "Measurements by participants",
            subtitle = paste(dimReduction_name, "\nEach participant's measurements are connected by a polygon")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just in minimal form, for publication
  gg_participantVariation_minimal <- ggplot(dat, aes(x = V1, y = V2)) +
    geom_point(aes(color = inter_treat), size = 2, alpha = 0.7) + 
    geom_polygon(aes(fill = participant_id), alpha = 0.3, show.legend = FALSE) +
    guides(color = "none") + 
    ggtitle(label    = "Measurements by participants",
            subtitle = paste(dimReduction_name, "\nEach participant's measurements are connected by a polygon")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just without T5
  gg_participantVariation_withoutT5 <- dat %>% 
    filter(timepoint != "T5") %>% 
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = timepoint, color = inter_treat), size = 2, alpha = 0.7, stroke = 1.5) + 
    geom_polygon(aes(fill = participant_id), alpha = 0.3, show.legend = FALSE) +
    scale_shape_manual(values = shape_vector) +
    geom_text(data = dat_participants, aes(label = participant_id, color = participant_id), check_overlap = FALSE, fontface = "bold") +
    guides(color = "none") + 
    ggtitle(label    = "Measurements by participants, excluding T5",
            subtitle = paste(dimReduction_name, "\nEach participant's measurements are connected by a polygon")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just splitted by interventions (and with lines instead of polygons)
  gg_participantVariation_byIntervention <- dat %>% 
    filter(intervention %in% c("Fresh", "Pasteurized")) %>% 
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = treatment), size = 2, alpha = 0.7, stroke = 1.5, color = "gray50") + 
    geom_line(aes(color = participant_id), alpha = 0.5, show.legend = FALSE) +
    geom_text(data = dat %>% filter(intervention != "FollowUp", timepoint %in% c("T2","T4")),
              aes(label = participant_id, color = participant_id), check_overlap = FALSE,
              hjust = -0.2, size = 3, fontface = "bold") +
    facet_wrap(~ intervention) +
    scale_shape_manual(values = shape_vector) +
    guides(color = "none") + 
    ggtitle(label    = "Measurements by participants",
            subtitle = paste(dimReduction_name, "\nEach participant's measurements are connected by a line")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just colored by the participants' initial Shannon diversity
  gg_participantVariation_byInt_div <- dat %>% 
    filter(intervention %in% c("Fresh", "Pasteurized")) %>% 
    arrange(participant_id, timepoint) %>% 
    group_by(participant_id, intervention) %>% 
    mutate(divShannon_initial = first(diversity_shannon)) %>% 
    ungroup() %>% 
    mutate(divShannon_initial_cat = case_when(divShannon_initial >= quantile(divShannon_initial, probs = .66) ~ "high diversity",
                                              divShannon_initial >= quantile(divShannon_initial, probs = .33) ~ "mid diversity",
                                              TRUE ~ "low diversity"),
           divShannon_initial_cat = factor(divShannon_initial_cat, levels = c("low diversity", "mid diversity", "high diversity"))) %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(color = divShannon_initial_cat, shape = treatment), size = 2, alpha = 0.7, stroke = 1.5) + 
    geom_line(aes(color = divShannon_initial_cat, group = participant_id), alpha = 0.5, show.legend = FALSE) +
    facet_wrap(~ intervention) +
    scale_shape_manual(values = shape_vector) +
    scale_color_manual(values = c("firebrick3", "lightblue", "dodgerblue3")) +
    ggtitle(label    = "Measurements by participants",
            subtitle = paste(dimReduction_name, "\nEach participant's measurements are connected by a line, colored by the initial (T1 or T3) measurement")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just splitted by size of change in observed richness
  gg_participantVariation_byRichness <- dat %>% 
    filter(intervention %in% c("Fresh", "Pasteurized")) %>% 
    group_by(intervention, participant_id) %>% 
    filter(n() == 2) %>% 
    mutate(richness_change = abs(diff(richness_observed))) %>% 
    ungroup() %>% 
    group_by(participant_id) %>% 
    mutate(richness_change = mean(richness_change)) %>% 
    ungroup() %>% 
    mutate(high_richnessChange = case_when(richness_change > median(richness_change) ~ "top 50%",
                                           TRUE                                      ~ "bottom 50%")) %>% 
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = treatment, color = high_richnessChange), size = 2, alpha = 0.7, stroke = 1.5) + 
    geom_line(aes(group = participant_id, color = high_richnessChange), alpha = 0.5, show.legend = FALSE) +
    facet_wrap(~ high_richnessChange) +
    scale_shape_manual(values = shape_vector) +
    scale_color_brewer(palette = "Dark2") +
    guides(color = "none") +
    ggtitle(label    = "Measurements by change in observed richness",
            subtitle = paste(dimReduction_name, "\nParticipants are categorized based on how much their observed richness changed over both interventions\nEach participant's measurements are connected by a line")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  # same plot, just splitted by size of change in Shannon diversity
  gg_participantVariation_byDiversity <- dat %>% 
    filter(intervention %in% c("Fresh", "Pasteurized")) %>% 
    group_by(intervention, participant_id) %>% 
    filter(n() == 2) %>% 
    mutate(diversity_change = abs(diff(diversity_shannon))) %>% 
    ungroup() %>% 
    group_by(participant_id) %>% 
    mutate(diversity_change = mean(diversity_shannon)) %>% 
    ungroup() %>% 
    mutate(high_diversityChange = case_when(diversity_change > median(diversity_change) ~ "top 50%",
                                            TRUE                                        ~ "bottom 50%")) %>% 
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = treatment, color = high_diversityChange), size = 2, alpha = 0.7, stroke = 1.5) + 
    geom_line(aes(group = participant_id, color = high_diversityChange), alpha = 0.5, show.legend = FALSE) +
    facet_wrap(~ high_diversityChange) +
    scale_shape_manual(values = shape_vector) +
    scale_color_brewer(palette = "Dark2") +
    guides(color = "none") +
    ggtitle(label    = "Measurements by change in Shannon diversity",
            subtitle = paste(dimReduction_name, "\nParticipants are categorized based on how much their Shannon diversity changed over both interventions\nEach participant's measurements are connected by a line")) +
    xlab(x_lab) + ylab(y_lab)
  
  
  
  # create the first plot list ----------------------------------------------
  gg_list <- list("gg_basePlot"                            = gg_basePlot,
                  "gg_Lparacasei_byIntervention"           = gg_Lparacasei_byIntervention,
                  "gg_Lparacasei_byIntervention_cat"       = gg_Lparacasei_byIntervention_cat,
                  "gg_LparacaseiChange_byIntervention"     = gg_LparacaseiChange_byIntervention,
                  "gg_participantVariation_byInt_div"      = gg_participantVariation_byInt_div,
                  "gg_dominantTaxa"                        = gg_dominantTaxa,
                  "gg_dominantTaxa_byAge"                  = gg_dominantTaxa_byAge,
                  "gg_maxGrowthTaxa_fresh"                 = gg_maxGrowthTaxa_fresh,
                  "gg_maxGrowthTaxa_pasteurized"           = gg_maxGrowthTaxa_past,
                  "gg_maxLossTaxa_fresh"                   = gg_maxLossTaxa_fresh,
                  "gg_maxLossTaxa_pasteurized"             = gg_maxLossTaxa_past,
                  "gg_groups"                              = gg_groups,
                  "gg_interventionTimes"                   = gg_interventionTimes,
                  "gg_interventionTimes_v2"                = gg_interventionTimes_v2,
                  "gg_participantVariation"                = gg_participantVariation,
                  "gg_participantVariation_minimal"        = gg_participantVariation_minimal,
                  "gg_participantVariation_withoutT5"      = gg_participantVariation_withoutT5,
                  "gg_participantVariation_byIntervention" = gg_participantVariation_byIntervention,
                  "gg_participantVariation_byDiversity"    = gg_participantVariation_byDiversity,
                  "gg_participantVariation_byRichness"     = gg_participantVariation_byRichness,
                  "gg_participantVariation_byDiversity"    = gg_participantVariation_byDiversity)
  
  
  
  # create plots for optional covariates ------------------------------------
  if (!is.null(covars_categorical) | !is.null(covars_metric)) {
    
    covar_vec <- c(covars_categorical, covars_metric)
    
    for (x in covar_vec) {
      
      # rename the variable for easier handling
      dat_covar <- dat
      colnames(dat_covar)[colnames(dat_covar) == x] <- "x"
      
      # hack part 1: specific handling of the gender variable, to give it correct labels
      # (somewhat messy hack, I know...)
      if (x == "gender") {
        dat_covar <- dat_covar %>% 
          mutate(x = factor(x, levels = c("M", "W"), labels = c("male", "female")))
      }
      
      # plot
      gg <- ggplot(dat_covar, aes(x = V1, y = V2, color = x)) +
        geom_point(size = 3, alpha = 0.9, stroke = 0) +
        ggtitle(paste("Measurements by", x), subtitle = dimReduction_name) +
        xlab(x_lab) + ylab(y_lab)
      
      # variable type specifics
      if (x %in% covars_categorical) {
        gg <- gg + stat_ellipse()
        
        # hack part 2: change the gender variable label
        x <- ifelse(x == "gender", "sex", x)
        
        if (length(unique(dat$x)) <= 8) { gg <- gg + scale_color_brewer(x, palette = "Dark2") }
        
      } else { # metric variable
        gg <- gg + scale_color_continuous(x, type = "viridis")
      }
      
      
      # add the plot to the results list
      gg_list[[length(gg_list) + 1]]  <- gg
      names(gg_list)[length(gg_list)] <- paste0("covar_", x)
    }
  }
  
  
  
  # return plot list --------------------------------------------------------
  return(gg_list)
  
}
