#' Encode DNA Sequences into Position-Specific Features
#'
#' Converts DNA sequences into a data frame of categorical
#' nucleotide features suitable for machine learning models.
#'
#' @param dna_strings A character vector containing DNA sequences
#'   of equal length, using the nucleotides A, T, C and G.
#'
#' @return A data frame with one column per nucleotide position.
#'   Each column is a factor with levels A, T, C and G.
#'
#' @examples
#' dna_encoding(c("GGACA", "AGACT"))
#'
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A Modification for Multiple RNA Sites
#'
#' Uses a trained random forest classifier to predict the
#' probability and modification status of multiple RNA sites
#' based on sequence and biological features.
#'
#' @param ml_fit A trained random forest classification model.
#' @param feature_df A data frame containing the required
#'   features: gc_content, RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold A numeric probability threshold
#'   between 0 and 1 used to classify m6A modification.
#'   Defaults to 0.5. Probabilities greater than this threshold
#'   are classified as Positive.
#'
#' @return A data frame containing the original input features
#'   and two additional columns:
#'   \describe{
#'     \item{predicted_m6A_prob}{Predicted probability of m6A modification.}
#'     \item{predicted_m6A_status}{Predicted status: Positive or Negative.}
#'   }
#'
#' @examples
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#'
#' example_data <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' ))
#'
#' prediction_multiple(
#'   ml_fit = rf_fit,
#'   feature_df = example_data
#' )
#'
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- stats::predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}
#' Predict m6A Modification for a Single RNA Site
#'
#' Predicts the probability and modification status of a
#' single RNA site using a trained random forest classifier.
#' The function uses sequence and biological features and
#' internally calls prediction_multiple().
#'
#' @param ml_fit A trained random forest classification model.
#' @param gc_content Numeric GC content of the input sequence.
#' @param RNA_type Character string specifying the RNA type:
#'   mRNA, lincRNA, lncRNA, or pseudogene.
#' @param RNA_region Character string specifying the RNA region:
#'   CDS, intron, 3'UTR, or 5'UTR.
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to the splice junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character string containing a five-nucleotide
#'   DNA sequence composed of A, T, C, and G.
#' @param positive_threshold Numeric classification threshold
#'   between 0 and 1. Defaults to 0.5.
#'
#' @return A named character vector containing:
#'   \describe{
#'     \item{predicted_m6A_prob}{Predicted m6A probability.}
#'     \item{predicted_m6A_status}{Predicted classification:
#'       Positive or Negative.}
#'   }
#'
#' @examples
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#'
#' prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = 0.5,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 150,
#'   distance_to_junction = 50,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA"
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
