#'generate region metadata from regions of interest given provided annotation files
#'
#' @param dog_id the dog_id for testing
#' @param alph alpha to use for glmnet
#' @param metaData sample metaData
#' @param perc_meth imputed percent methylation matrix
#' @return Function returns a list of predicted age as a column to the input metaData, region weights as a column to the input region_metaData, and clock results as a list
#' @export build_clock_leave_one_out

build_clock_leave_one_out <- function(dog_id,
                        alph,
                        metaData,
                        perc_meth
                        ) {
  # Read in meta info with known chronological ages/sex and technical variables.
  all_info <- metaData

  lid_pids <- colnames(perc_meth)
  meta <- all_info[,c("lid_pid", "Age_at_sample","dog_id", "Cohort", "Breed_size", "Sex")]

  n_regions <- nrow(perc_meth)

  # Convert the 'Age_at_sample' column to numeric
  meta$Age_at_sample <- as.numeric(meta$Age_at_sample)

  # Read in data - This should be the imputed df you are working with
  epi <- perc_meth

  #filter epi for dogs you have meta data for and the metadata for the perc meth you have data for
  epi <- epi[, colnames(epi) %in% meta$lid_pid]
  meta <- meta[meta$lid_pid %in% colnames(epi),]

  # Reorder the columns in 'epi' to match the order of 'meta$lid'
  meta <- meta[order(meta$lid_pid),]
  epi <- epi[,order(colnames(epi))]

  # Remove test subject(s)
  # SAMP indexes from 1 to N samples
  test_lid_pid <- metaData$lid_pid[metaData$dog_id == dog_id]
  train <- epi[, colnames(epi)!=test_lid_pid]
  test <- epi[, colnames(epi)==test_lid_pid]

  # Create a vector of training and test ages for elastic net model construction
  trainage <- meta$Age_at_sample[colnames(epi)!=test_lid_pid]
  testage <- meta$Age_at_sample[colnames(epi)==test_lid_pid]
  test_id <- meta$dog_id[colnames(epi)==test_lid_pid]
  test_lid <- meta$lid_pid[colnames(epi)==test_lid_pid]

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

  return(predicted[1,1])
}
