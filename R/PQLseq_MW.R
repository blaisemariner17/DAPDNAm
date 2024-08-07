#'modified PQLseq function by Marina Watowich to output covariate results
#' #### Updated MW Feb 11, 2019
#'
#' @param RawCountDataSet methylation mat
#' @param Covariates covariates mat
#' @param numCore numer of cores
#' @param RelatednessMatrix Relatedness Matrix
#' @return BMM mixed-model result
#' @export pqlseq_MW
#'

########################################################################################################################

pqlseq_MW <- function(RawCountDataSet, Phenotypes=NULL, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL,
                   fit.model="PMM", fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, numCore=1,
                   filtering=TRUE, verbose=FALSE, ...) {
  # specify the number of cores we want to use
  if(numCore > 1){
    if(numCore>detectCores()){warning("PQLseq:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()-1}
  }

  registerDoParallel(numCore)

  # cl <- makeCluster(numCore)
  # registerDoParallel(cl,cores=numCore)
  # on.exit(stopCluster(cl))

  # filtering genes/sites
  if (filtering & fit.model == "PMM"){
    unfilterIdx <- apply(RawCountDataSet, 1, function(x){length(x[x>5])>=2} )
    CountData   <- RawCountDataSet[unfilterIdx,]
  }else{
    CountData   <- RawCountDataSet
  }
  rm(RawCountDataSet)

  numVar <- dim(CountData)[1]
  numIDV <- dim(CountData)[2]

  # remove the intercept
  if(length(unique(Covariates[,1])) == 1){
    Covariates<- Covariates[,-1]
  }

  if(is.null(Covariates)){
    numCov <- 0
  }else{
    numCov     <- dim(Covariates)[2]
    Covariates <- as.matrix(Covariates)
  }

  cat(paste("## number of total individuals: ", numIDV,"\n"))
  cat(paste("## number of total genes/sites: ", numVar,"\n"))
  cat(paste("## number of adjusted covariates: ", numCov,"\n"))

  if(is.null(Phenotypes)){
    numPhe 		<- 0
  }else{
    Phenotypes 	<- as.matrix(Phenotypes)
    numPhe 		<- 1
  }# end fi

  CountData  <- as.matrix(CountData)
  # Phenotypes <- as.matrix(Phenotypes)


  if(is.null(RelatednessMatrix)){
    stop("PQLseq::please input relatedness matrix!")
  }else{
    RelatednessMatrix <- as.matrix(RelatednessMatrix)
    scalerM           <- diag(numIDV)-(rep(1,numIDV)%*%t(rep(1,numIDV)))/numIDV
    eig               <- eigen(RelatednessMatrix)
    eigval            <- eig$value
    eigvector         <- eig$vectors
    if(any(eigval<1e-10)){
      warning("PQLseq::the relatedness matrix is singular, it has been modified!")
      RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T,maxit=500)$mat)
    }
    rm(scalerM)
    rm(eig)
    rm(eigval)
    rm(eigvector)
  }

  RelatednessMatrix <- list(RelatednessMatrix, diag(numIDV))

  #***********************************#
  #       Poisson Mixed Model         #
  #***********************************#
  if(fit.model == "PMM"){
    cat("# fitting Poisson mixed model ... \n")
    if(is.null(LibSize)){
      LibSize <- apply(CountData, 2, sum)
      LibSize <- as.matrix(LibSize)
    }else{
      LibSize <- as.matrix(t(LibSize))
    }


    # do parallel using foreach function
    iVar   <- NULL
    resPMM <-foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA

      if(numPhe==0 && numCov==0){
        model0 <- try(glm(formula = as.numeric(CountData[iVar,]) ~ 1 + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = as.numeric(CountData[iVar,]) ~ 1 + offset(log(LibSize)), na.action = na.omit)),rownames(model.frame(formula = as.numeric(CountData[iVar,]) ~ 1 + offset(log(LibSize)), na.action = na.pass)))
      }else if(numPhe !=0 && numCov==0){
        model0 <- try(glm(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.pass)))
      }else{
        model0 <- try(glm(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.pass)))
      }

      if(verbose) {cat(paste("NO. Gene = ",iVar,"\n"))}

      tmpRelatednessMatrix <- RelatednessMatrix
      if(class(tmpRelatednessMatrix) == "matrix") {
        tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
      }else {
        for(ik in seq_len(length(tmpRelatednessMatrix)) ) {tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]}
      }

      names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

      if(class(model0)[1]!="try-error"){
        # t1 <- system.time(model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix)))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))
      }else{
        model1 <- NULL
      }

      if(!is.null(model1)&(class(model1)!="try-error")){
        if(verbose){cat(paste("PQLseq::PMM::tau = ", model1$theta,"\n"))}
        numAnalysis <- length(idx)

        # beta        <- model1$coefficients[length(model1$coefficients)]
        # se_beta     <- sqrt(diag(model1$cov)[length(model1$coefficients)])
        # pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)

        all_beta        <- model1$coefficients[-1]
        all_se_beta     <- sqrt(diag(model1$cov))[-1]
        all_pvalue 		<- pchisq((all_beta/all_se_beta)^2, 1, lower.tail = F)


        sigma2      <- model1$theta[2]+model1$theta[3]
        h2          <- model1$theta[2]/(sigma2)
        LL 			<- model1$LL
        tau1        <- model1$theta[2]
        tau2        <- model1$theta[3]
        converged   <- model1$converged
      }else{converged <- FALSE}

      res <- data.frame(numAnalysis, matrix(all_beta,nrow=1), matrix(all_se_beta,nrow=1),
                        matrix(all_pvalue,nrow=1),h2, sigma2, LL, converged)

      # res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta,
      # 				  pvalue = pvalue, h2 = h2, sigma2 = sigma2, LL=LL,
      # 				  converged = converged)


    }# end for iVar, parallel
    rm(iVar)
    closeAllConnections()

    colnames(resPMM) <- c("numIDV",
                          c(paste0("beta_covariates",1:numCov),"beta_predictor"),
                          c(paste0("se_covariates",1:numCov),"se_predictor"),
                          c(paste0("pval_covariates",1:numCov),"pval_predictor"),
                          "h2","sigma2","LL","converged")
    rownames(resPMM) <- rownames(CountData)
    return(resPMM)
  }# end PMM
  #***********************************#
  #       Binomial Mixed Model        #
  #***********************************#
  if(fit.model == "BMM"){
    cat("# fitting binomial mixed model ... \n")
    if(is.null(LibSize)){
      stop("PQLseq::BMM::ERROR: please input the LibSize (total counts) file!!")
    }else{
      LibSize <- as.matrix(LibSize)
    }

    ratio               <- CountData/LibSize
    ratio[is.na(ratio)] <- 0
    flag                <- ratio>1.0
    sumflag             <- apply(flag,1, sum)
    idx                 <- which(sumflag>0)

    if (length(idx)>0){
      CountData <- CountData[-idx,]
      LibSize   <- LibSize[-idx,]
    }else{
      CountData <- CountData
      LibSize   <- LibSize
    }

    numVar <- dim(CountData)[1]
    numIDV <- dim(CountData)[2]
    iVar   <- NULL

    # do parallel #this is the area to change if want to make a df, which allows different classes of columns (matrix forces all columns to be same class)
    resBMM <- foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
      all_beta 	<- all_se_beta <- all_pvalue <- rep(NA,numCov+1)
      if(verbose){cat(paste("NO. Gene/Site = ",iVar,"\n"))}
      if(sum(dim(LibSize)==dim(CountData)) != 2){
        stop("PQLseq::BMM::ERROR: the dimensions of read counts and total read counts do not match!")
      }

      LibSize <- as.matrix(LibSize)

      if(numCov == 0){
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,])
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.pass)))
      }else{
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,] )
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.pass)))
      }

      model0$numTotal <- LibSize[iVar,idx]
      model0$numSucc  <- CountData[iVar,idx]

      redflag <- FALSE
      for( ierr in c(2:dim(model.matrix(model0))[2])){
        if(length(unique(model.matrix(model0)[,ierr])) == 1){
          warning(paste("PQLseq::BMM::the ",ierr-1,"-th column of covariates are the same for gene/site ",rownames(CountData)[iVar],"!",sep = "") )
          redflag <- TRUE
        }
      }
      if(!redflag){

        tmpRelatednessMatrix <- RelatednessMatrix
        if(class(tmpRelatednessMatrix) == "matrix") {
          tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
        }else {
          for(ik in seq_len(length(tmpRelatednessMatrix)) ) {
            tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]
          }
        }
        names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

        # t1 <- system.time(model1 <- try( PQLseq.fit(model0, tmpRelatednessMatrix) ))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))

        if(class(model1) != "try-error"&!is.null(model1)){
          if(verbose){cat(paste("PQLseq::BMM::tau = ", model1$theta,"\n"))}
          numAnalysis <- length(idx)


          # beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
          # se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )

          all_beta        <- model1$coefficients[-1]
          all_se_beta     <- sqrt(diag(model1$cov))[-1]
          all_pvalue 		<- pchisq((all_beta/all_se_beta)^2, 1, lower.tail = F)

          pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
          sigma2      <- model1$theta[2]+model1$theta[3]
          h2          <- model1$theta[2]/(sigma2)
          LL 	      <- model1$LL
          AIC         <- -2*(LL) + 2*(length(model1$coefficients[-1]))
          tau1        <- model1$theta[2]
          tau2        <- model1$theta[3]
          converged   <- model1$converged
        }else{converged <- FALSE}

        # res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta,
        # 			pvalue = pvalue, h2 = h2, sigma2 = sigma2,
        # 			converged = converged)

        res <- data.frame(numAnalysis, matrix(all_beta,nrow=1), matrix(all_se_beta,nrow=1),
                          matrix(all_pvalue,nrow=1),h2, sigma2, LL, AIC, converged)  #MW add AIC

        # res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta,
        # 				  pvalue = pvalue, h2 = h2, sigma2 = sigma2, LL=LL,
        # 				  converged = converged)


      }# end for iVar, parallel

    }
    rm(iVar)


    closeAllConnections()
    rownames(resBMM) <- rownames(CountData)
    colnames(resBMM) <- c("numIDV",
                          c(paste0("beta_",colnames(Covariates)),"beta_predictor"), #MW changed these names to be more intuitive
                          c(paste0("se_",colnames(Covariates)),"se_predictor"),
                          c(paste0("pval_",colnames(Covariates)),"pval_predictor"),
                          "h2","sigma2","LL", "AIC", "converged") #MW added AIC
    return(resBMM)
  }# end BMM

}# end function PQLseq


##########################################################
#           	   PQLseq FIT FUNCTION					 #
##########################################################

PQLseq.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {

  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  # if((method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
  if(method.optim == "AI") {
    fixtau.old 	<- rep(0, length(RelatednessMatrix)+1)
    # to use average information method to fit alternative model
    model1 		<- PQLseq.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new 	<- 1*(model1$theta < 1.01 * tol)

    while(any(fixtau.new != fixtau.old)) {
      fixtau.old <- fixtau.new
      model1 	<- PQLseq.AI(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(model1$theta < 1.01 * tol)
    }
  }else{
    model1 <- NULL
  }
  return(model1)
}

##########################################################
#       PQLseq FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################

PQLseq.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {

  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
  numIDV <- length(y)
  offset <- model0$offset
  if(is.null(offset)) {offset <- rep(0, numIDV)}

  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  D <- mu.eta/sqrt(model0$family$variance(mu))

  if(family$family %in% c("binomial")){
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
    mu <- model0$numTotal*mu
  }

  Y <- eta - offset + (y - mu)/mu.eta
  X <- model.matrix(model0)
  alpha <- model0$coef

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  numK <- length(RelatednessMatrix)
  idxtau <- which(fixtau == 0)
  numK2 <- sum(fixtau == 0)
  if(numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)

    H <- tau[1]*diag(1/D^2)
    for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

    Hinv 	<- chol2inv(chol(H))
    HinvX 	<- crossprod(Hinv, X)
    XHinvX 	<- crossprod(X, HinvX)

    P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

    if(class(P) == "try-error"){
      stop("Error in P matrix calculation!")
    }

    PY <- crossprod(P, Y)
    tau0 <- tau
    for(ik in 1:numK2) {
      if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
      else {
        PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
        tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
      }
    }
  }

  for (iter in seq_len(maxiter)) {
    alpha0 	<- alpha
    tau0 	<- tau
    model1 	<- PQLseq.AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

    tau <- as.numeric(model1$tau)
    cov <- as.matrix(model1$cov)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta) + offset


    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta

    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
    if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

      iter <- maxiter
      break
    }
  }

  converged <- ifelse(iter < maxiter, TRUE, FALSE)
  res <- y - mu


  # D is updated , so the P should also be updated
  # P <- model1$P


  H <- tau[1]*diag(1/D^2)
  for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

  Hinv 	<- chol2inv(chol(H))
  HinvX 	<- crossprod(Hinv, X)
  XHinvX 	<- crossprod(X, HinvX)

  P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))


  ## calculate the likelihood
  PY 	<- crossprod(P, Y)
  YPY <- crossprod(Y, PY)

  # ## REML
  # LL 	<- try(1/2*length(coefficients)*log(2*pi)-1/2*(logDET(diag(D^2))+ logDET(H)+ logDET(XHinvX) + YPY))

  ## ML
  LL 	<- try(-1/2*(logDET(diag(D^2))+ logDET(H)+ logDET(XHinvX) + YPY))

  if(class(LL)=="try-error"){
    cat("summary of D:\n")
    LL <- NA
  }

  return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, LL=LL, P = P, residuals = res, cov = cov, converged = converged))

}# end function


#########################################
#             CODE END                  #
#########################################


## Log Determinant

logDET <- function(mat){
  eig_val <- eigen(mat)$values
  ld 		<- sum(log(eig_val))
  return(ld)
}






