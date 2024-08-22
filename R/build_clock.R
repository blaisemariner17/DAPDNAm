#'generate region metadata from regions of interest given provided annotation files
#'
#' @param test_samples test samples for elastic net
#' @param alph alpha to use for glmnet
#' @param metaData sample metaData
#' @param region_metaData region metaData
#' @param perc_meth imputed percent methylation matrix
#' @return Function returns a list of predicted age as a column to the input metaData, region weights as a column to the input region_metaData, and clock results as a list
#' @export build_clock

build_clock <- function(test_samples,
                        alph = 0.5,
                        metaData,
                        region_metaData,
                        perc_meth
                        ) {
  test_samples <- order(test_samples)

  # Read in meta info with known chronological ages/sex and technical variables.
  all_info <- metaData

  lids <- colnames(perc_meth)
  meta <- all_info[,c("lid_pid", "Age_at_sample","dog_id", "Cohort", "Breed_size", "Sex")]

  n_regions <- nrow(perc_meth)

  # Convert the 'Age_at_sample' column to numeric
  meta$Age_at_sample <- as.numeric(meta$Age_at_sample)

  # Read in data - This should be the imputed df you are working with
  epi <- perc_meth

  #filter epi for dogs you have meta data for
  epi <- epi[, colnames(epi) %in% meta$lid_pid]

  # Reorder the columns in 'epi' to match the order of 'meta$lid'
  meta <- meta[order(meta$lid_pid),]
  epi <- epi[,order(colnames(epi))]

  # Remove test subject(s)
  train <- epi[, !colnames(epi) %in% test_samples]
  test <- epi[, colnames(epi) %in% test_samples]

  # Create a vector of training and test ages for elastic net model construction
  trainage <- meta$Age_at_sample[!rownames(meta) %in% test_samples]
  testage <- meta$Age_at_sample[rownames(meta) %in% test_samples]
  test_id <- meta$dog_id[rownames(meta) %in% test_samples]
  test_lid <- meta$lid_pid[rownames(meta) %in% test_samples]

  #### Elastic-net model building ####

  # Using N-fold internal CV, train the elastic net model using the training data
  # Note with larger sample sizes, N-fold internal CV becomes intractable
  # model <- cv.glmnet(train, trainage, nfolds = X, alpha = alph, standardize = F)

  for (row_ in 1:nrow(train)){
    if (sum(is.na(train[row_,])) > 0){
      stop(paste("Row", row_, "has NAs. Please impute before running this function."))
    }
  }

  ## Transpose the train matrix to be samples x features
  model <- glmnet::cv.glmnet(t(train), trainage, nfolds = 10, alpha = alph, standardize = FALSE)

  # Predict age using the test sample from parameters that minimized MSE during internal CV
  predicted <- predict(model, newx = t(test), s = "lambda.min")

  rownames(predicted) <- test_samples

  # Calculate mean squared error (MSE)
  mse <- mean((predicted - testage) ^ 2)

  median_ae <- MLmetrics::MedianAE(predicted, testage)
  mae <- MLmetrics::MAE(predicted, testage)
  r2 <- MLmetrics::R2_Score(predicted, testage)

  # Extract weights for this model
  weights <- unlist(coef(model, lambda = "lambda.min"))[, 1]

  meta <- data.frame(meta[meta$lid_pid %in% rownames(predicted),])

  for (lid_pid in meta$lid_pid){
    meta$predicted[meta$lid_pid == lid_pid] <- predicted[rownames(predicted) == lid_pid,1]
  }

  region_metaData <- region_metaData[rownames(region_metaData) %in% names(weights),]
  weights <- weights[names(weights) %in% rownames(region_metaData)]

  region_metaData <- region_metaData[order(rownames(region_metaData)),]
  weights <- weights[order(names(weights))]

  if (all(names(weights) == rownames(region_metaData))){
    reg_meta <- cbind(weights, region_metaData)
  } else {
    stop("regions used in clock do not match rownames(region_metaData)")
  }

  return_ <- list()
  return_[["metaData"]] <- meta
  return_[["region_metaData"]] <- reg_meta
  return_[["clock_results"]] <- data.frame("alpha" = alph,
                                           "MSE" = mse,
                                           "MAE" = mae,
                                           "Median_ae" = median_ae,
                                           "R2" = r2,
                                           "n_regions" = n_regions,
                                           "training_samples" = length(colnames(perc_meth)[!colnames(perc_meth) %in% test_samples]),
                                           "testing_samples" = length(test_samples))

  return(return_)
}
