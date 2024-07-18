
#' Boxplots with participant development lines for crossover studies
#' 
#' This function uses boxplots and one line between the boxplots per participant
#' to visualize the development of some feature in a crossover study.
#' Additionally, it returns a summary table of the visualized data.
#' 
#' @param dat data.frame containing all the variables stated as argument
#' @param y_var Character name of the target variable
#' @param intervention_var Character name of the intervention group variable.
#' The variable in the dataset should have exactly two unique values.
#' @param t_var Character name of the measurement time variable. The variable
#' in the dataset should have exactly two unique values.
#' @param subject_var Character name of a subject id variable
#' @param cutOff_relevantChange Numeric, non-negative cut-off parameter used for
#' differentiating the three categories \code{c("pos. change", "no change", "neg. change")}.
#' Defaults to 0.01.
#' @param y_isRelAbundance Indicator if the \code{y_var} variable contains
#' relative abundances, in which case the plot labels are slightly altered.
#' Defaults to FALSE.
#' @param title Optional character title for the plot
#' @param print_summaryTable Indicator if a summary table of the individual
#' changes of the participants in \code{y_var} should be returned in a joint
#' list with the ggplot object. Defaults to \code{TRUE}.
#' 
#' @import checkmate dplyr ggplot2
#' @export
#' 
#' @return ggplot2 object
#' 
plot_crossoverBoxplots <- function(dat, y_var, intervention_var, t_var,
                                   subject_var,
                                   cutOff_relevantChange = 0.01,
                                   y_isRelAbundance      = FALSE,
                                   title                 = NULL,
                                   return_summaryTable   = TRUE) {
  
  checkmate::assert_data_frame(dat)
  checkmate::assert_choice(y_var,            choices = colnames(dat))
  checkmate::assert_choice(t_var,            choices = colnames(dat))
  checkmate::assert_choice(intervention_var, choices = colnames(dat))
  checkmate::assert_choice(subject_var,      choices = colnames(dat))
  checkmate::assert_vector(dat[[t_var]]            %>% unique(), len = 2)
  checkmate::assert_vector(dat[[intervention_var]] %>% unique(), len = 2)
  checkmate::assert_number(cutOff_relevantChange, lower = 0)
  checkmate::assert_logical(y_isRelAbundance,     len = 1)
  checkmate::assert_character(title,              len = 1)
  checkmate::assert_logical(return_summaryTable,  len = 1)
  
  
  # data preparation --------------------------------------------------------
  # rename all used variables, for ease of use
  dat <- dat %>% 
    rename(y            = y_var,
           t            = t_var,
           intervention = intervention_var,
           subject      = subject_var)
  
  # create a variable that indicates if a person's 'y_var' value fell or sunk
  change_labels <- c("neg." = paste0("neg. change <= ", cutOff_relevantChange, 
                                     ifelse(y_isRelAbundance, " PP", "")),
                     "no"   = paste0("change < ", cutOff_relevantChange, 
                                     ifelse(y_isRelAbundance, " PP", "")),
                     "pos." = paste0("pos. change >= ", cutOff_relevantChange, 
                                     ifelse(y_isRelAbundance, " PP", "")))
  dat <- dat %>% 
    arrange(subject, intervention, t) %>% 
    group_by(subject, intervention) %>% 
    filter(n() == 2) %>% 
    mutate(y_diff         = diff(y),
           y_changeFactor = y[2] / y[1]) %>% 
    ungroup() %>% 
    mutate(y_change = case_when(y_diff <= -1*cutOff_relevantChange ~ change_labels["neg."],
                                y_diff <  cutOff_relevantChange    ~ change_labels["no"],
                                TRUE                               ~ change_labels["pos."]),
           y_change = factor(y_change, levels = unname(change_labels)))
  
  
  
  # summary table -----------------------------------------------------------
  if (return_summaryTable) {
    results_tab <- dat %>% 
      group_by(intervention, subject) %>% 
      slice(1) %>% 
      ungroup() %>% 
      group_by(intervention) %>% 
      summarize(share_posChange  = sum(y_change == change_labels["pos."]) / n(),
                share_noChange   = sum(y_change == change_labels["no"])   / n(),
                share_negChange  = sum(y_change == change_labels["neg."]) / n(),
                median_posChange = median(y_diff[y_change == change_labels["pos."]]),
                median_noChange  = median(y_diff[y_change == change_labels["no"]]),
                median_negChange = median(y_diff[y_change == change_labels["neg."]]),
                median_posChangeFactor = median(y_changeFactor[y_change == change_labels["pos."]]),
                median_noChangeFactor  = median(y_changeFactor[y_change == change_labels["no"]]),
                median_negChangeFactor = median(y_changeFactor[y_change == change_labels["neg."]])) %>% 
      mutate(across(starts_with("share"),  function(x) { paste0(round(100 * x, 1), "%") }),
             across(contains("Factor"),    function(x) { ifelse(is.na(x), "-", paste0(ifelse(x >= 1, "+", ""), round(100 * (x - 1)), "%")) }),
             median_posChange = case_when(is.na(median_posChange) ~ "-",
                                          TRUE                    ~ paste0(round(median_posChange, 2), ifelse(y_isRelAbundance, " PP", ""))),
             median_noChange  = case_when(is.na(median_noChange)  ~ "-",
                                          TRUE                    ~ paste0(round(median_noChange, 2),  ifelse(y_isRelAbundance, " PP", ""))),
             median_negChange = case_when(is.na(median_negChange) ~ "-",
                                          TRUE                    ~ paste0(round(median_negChange, 2), ifelse(y_isRelAbundance, " PP", ""))))
  }
  
  
  
  # plot --------------------------------------------------------------------
  # prepare the color vector for the plot depending on which categories appear
  # in the data
  col_vector <- c(if (change_labels["neg."] %in% unique(dat$y_change)) { "firebrick3" } else { NULL },
                  if (change_labels["no"]   %in% unique(dat$y_change)) { "gray80"     } else { NULL },
                  if (change_labels["pos."] %in% unique(dat$y_change)) { "turquoise3" } else { NULL })
  
  gg <- dat %>% 
    ggplot(aes(x = t, y = y)) +
    geom_boxplot() +
    geom_line(aes(group = subject, col = y_change), alpha = .3) +
    facet_grid(cols = vars(intervention)) +
    scale_y_continuous(labels = function(x) { paste0(x, ifelse(y_isRelAbundance, "%", "")) } ) +
    ggtitle(title) + ylab(y_var) +
    scale_color_manual("change", values = col_vector) +
    theme(axis.title.x = element_blank(),
          plot.title   = element_text(hjust = 0.5))
  
  
  
  # return results ----------------------------------------------------------
  if (return_summaryTable) {
    return(list("ggplot_object" = gg,
                "results_table" = results_tab))
  } else {
    return(gg)
  }
}
