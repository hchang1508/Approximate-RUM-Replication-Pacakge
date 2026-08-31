# Functions to conduct empirical analysis on Fishing data.
# The following functions require the auxiliary functions in Fish_functions.R.
#
# Sample splitting, mlogit preparation, aggregation, and the six aggregate-share
# model likelihoods that tbl:in-sample_d1 and tbl:out-sample estimate.

library(dplyr)
library(fastDummies)
library(mlogit)
library(Rsolnp)


split.sample <- function(dat, K = 2,index,data_path){
  # Input
  #     dat: micro data
  #          It has to include the 'chid' variable to identify individuals.
  #     K: how many split samples do we use?
  # Output
  #     samples.list: the list of training and test samples.
  #                   The k-th test sample is the k-th subsample.
  
  ids <- sort(unique(dat[, "chid"]))
  set.seed(index); training.idx <- sample(seq(1:K), length(ids), replace = T)
  
  # Store samples in a list
  samples.list <- list()
  for (k in 1:K) {
    idx <- which(training.idx == k)
    dat.test          <- dat[ which(dat[, "chid"] %in% idx), ]
    dat.train         <- dat[ which(dat[, "chid"] %in% setdiff(ids, idx)), ] 
    samples.list[[k]] <- list(dat.train, dat.test)
    names(samples.list[[k]]) <- c('train', 'test')
    
    name_temp_train = paste0(data_path,'data.train.',index,'.',k,'.csv')
    write.csv(dat.train,name_temp_train) #save kth split training data
 
    name_temp_test = paste0(data_path,'data.test.',index,'.',k,'.csv')
    write.csv(dat.test,name_temp_test) #save kth split test data   
    
  }
  return(samples.list)
}


# Function to prepare data in appropriate (mlogit) format and with alternative and boat-type fixed effects
prepare.mlogit.dat <- function(dat, pol.label, id.var = 'chid', avg.char = T){
  # Input
  #     dat:      micro (individual-level) or macro (alternative-level) data. Typically it will be micro data.
  #     pol.label: names of characteristics included in the logit model
  #     id.var:   individual ID ('chid' for micro data and 'type' for macro data)
  #     avg.char: use average characteristics for each individual? 
  #               (you can set avg.char = T even with micro data)
  #     NOTE: input data is assumed to be of d=1.
  # Output
  #     dat:      micro or macro data in mlogit format and with alternative fixed effects
  
  CS        <- unique(dat$alt) %>% as.character()
  # Create dummies 
  dat <- dummy_cols(dat, select_columns = 'alt', remove_first_dummy = T) # Prevent the matrix from becoming singular
  # Convert the sample to an mlogit.data object.
  dat <- mlogit.data(data = dat, shape = 'long', choice = 'mode', alt.levels = CS, id.var = id.var) %>% as.data.frame()
  dat[, grep('^idx', colnames(dat))] <- NULL
  if (avg.char){
    for (char in pol.label) {
      for (alt in CS) {
        dat[which(dat$alt == alt), char] <- mean(dat[which(dat$alt == alt), char])
      }
    }
  }
  dat$alt <- as.character(dat$alt)
  return(dat)
}



# Function to generate macro data for representability analysis
gen.macro.data <- function(dat, pol.label){
  # Input
  #     dat: micro data
  #          It has to have variables alt (alternatives), X1.0 and X0.1 (two alternative characteristics).
  # Output
  #     dat.avg: macro (average) data
  #              It has an artificial variable named 'type' for use of some functions
  
  dat.macro <- data.frame(alt = unique(dat$alt), chid = 1)
  for (v in pol.label){
    dat$v   <- dat[, v]
    dat.aux <- group_by(dat, alt) %>% 
      summarize(., v.avg = mean(v)) %>% as.data.frame()  
    dat.macro[, v] <- dat.aux$v.avg
  }
  dat.macro$alt <- as.character(dat.macro$alt)
  return(dat.macro)
}



# Function to conduct regressions using alternative models with macro data
run.alt.regressions.macro <- function(dat.train, dat.test, pol.label, emp.CP.train, emp.CP.test, spec){
  # Input
  #     dat.train:    micro data (with average characteristics) for estimation
  #     dat.test:     micro data (with average characteristics) for testing
  #     pol.label:       names of characteristics included in the logit model
  #     emp.CP.train: in-sample empirical CP
  #     emp.CP.test:  out-of-sample empirical CP
  #     spec:         specification (multinomial logit, nested logit with chater and non-charter nests,
  #                   and nested logit with boat and non-boat nests)
  # Output
  #     output:       the list of param (estimated parameters), pred.in.macro (in-sample prediction), 
  #                   loss.in.macro (in-sample prediction error), pred.out.macro (out-of-sample prediction),
  #                   and loss.out.macro (out-of-sample prediction error)
  
  CS        <- unique(dat.train$alt)
  n.alt     <- length(CS)
  reg.formula.base <- Reduce(function(x, y) paste(x, y, sep = ' + '), c(pol.label, 0))
  reg.formula.base <- paste('mode', reg.formula.base, sep = ' ~ ') %>% as.formula()
  
  # Regression and prediction
  if (spec == 'mnl_aggre'){
    
    dat.train.macro <- gen.macro.data(dat.train, pol.label)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label)
    
    ll <- function(param) {
      logit.cp <-compute.logitCP(dat = dat.train.macro, pol.label = pol.label, b=param,cs=CS)
      ll.val   <- -sum(emp.CP.train * log(logit.cp)) # minus log-likelihood
      return(ll.val)
    } 
    
    sol   <- solnp(pars = c(rep(0,2)), fun = ll, 
                   LB = c(rep(-1000,2)), 
                   UB = c(rep(1000,2)),control=c(inner.iter=1600,tol=1e-8 ) )
    param <- sol[["pars"]]
    
    names(param) <- pol.label
    
    # In-sample prediction
    pred.in.avg <- compute.logitCP(dat = dat.train.macro, pol.label = pol.label, b=param,cs=CS)
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- compute.logitCP(dat = dat.test.macro, pol.label = pol.label, b=param,cs=CS)
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
    
  } else if (spec == 'nest_charter_aggre'){
    
    dat.train <- prepare.mlogit.dat(dat.train, pol.label)    # Replace individual characteristics by averages
    dat.test  <- prepare.mlogit.dat(dat.test, pol.label)
    
    # Use macro data because the mlogit package does not work
    dat.train.macro <- gen.macro.data(dat.train, pol.label)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label)
    nests  <- list('charter' = c('charter'), 'other' = c('beach', 'boat', 'pier'))
    n.char <- length(pol.label)
    # Define the objective function
    ll    <- function(param){
      cp       <- pred.nest.macro(dat = dat.train.macro, pol.label = pol.label, b = param[1:n.char], 
                                  lambda = param[n.char+1], nests = nests)
      ll.val   <- -sum(emp.CP.train * log(cp)) # MINUS log-likelihood
      return(ll.val)
    } 
    # MINIMIZE the negative of the log-likelihood
    ## Note: "Solution not reliable....Problem Inverting Hessian" when we start with lambda = 1.
    sol   <- solnp(pars = c(rep(0, n.char), 0.5), fun = ll, 
                   LB = c(-100, rep(-1000, n.char-1), 0), UB = c(1000, rep(1000, n.char-1), 1),control=c(inner.iter=1600,tol=1e-6))
    param <- sol[["pars"]]
    names(param) <- c(pol.label, 'nest')
    
    print('charter')
    print(param)
    # In-sample prediction
    pred.in.avg <- pred.nest.macro(dat = dat.train.macro, pol.label = pol.label,
                                   b = param[1:n.char], lambda = param[n.char+1], nests = nests)
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- pred.nest.macro(dat = dat.test.macro, pol.label = pol.label,
                                    b = param[1:n.char], lambda = param[n.char+1], nests = nests)
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
    
  } else if (spec == 'nest_boat_aggre'){
    
    dat.train <- prepare.mlogit.dat(dat.train, pol.label)    # Replace individual characteristics by averages
    dat.test  <- prepare.mlogit.dat(dat.test, pol.label)
    
    # Use macro data because the mlogit package does not work
    dat.train.macro <- gen.macro.data(dat.train, pol.label)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label)
    nests <- list('boat' = c('boat'), 'other' = c('beach', 'charter', 'pier'))
    n.char <- length(pol.label)
    
    # Define the objective function
    ll    <- function(param){
      cp       <- pred.nest.macro(dat = dat.train.macro, pol.label = pol.label, b = param[1:n.char], 
                                  lambda = param[n.char+1], nests = nests)
      ll.val   <- -sum(emp.CP.train * log(cp)) # MINUS log-likelihood
      return(ll.val)
    }
    
    # MINIMIZE the negative of the log-likelihood
    ## Note: "Solution not reliable....Problem Inverting Hessian" when we start with lambda = 1.
    sol   <- solnp(pars = c(rep(0, n.char), 0.5), fun = ll, 
                   LB = c(-100, rep(-1000, n.char-1), 0), UB = c(1000, rep(1000, n.char-1), 1),control=c(inner.iter=1600,tol=1e-6))
    
    param <- sol[["pars"]]
    names(param) <- c(pol.label, 'nest')
    print('boat')
    print(param)
    # In-sample prediction
    pred.in.avg <- pred.nest.macro(dat = dat.train.macro, pol.label = pol.label,
                                   b = param[1:n.char], lambda = param[n.char+1], nests = nests)
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- pred.nest.macro(dat = dat.test.macro, pol.label = pol.label,
                                    b = param[1:n.char], lambda = param[n.char+1], nests = nests)
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
    
  } else if (spec == 'random_coef_aggre'){
    
    dat.train <- prepare.mlogit.dat(dat.train, pol.label)    # Replace individual characteristics by averages
    dat.test  <- prepare.mlogit.dat(dat.test, pol.label)
    
    # Use macro data because the mlogit package does not work
    dat.train.macro <- gen.macro.data(dat.train, pol.label)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label)
    n.char <- length(pol.label)
    # Define the objective function
    if ('X_boattype' %in% pol.label){
      rc = c(X1.0 = 'ln', X0.1 = 'ln', X_boattype = 'n')  
    } else {
      rc = c(X1.0 = 'ln', X0.1 = 'ln')  
    }
    ll <- function(param) {
      logit.cp <- pred.rc.macro(dat = dat.train.macro, pol.label = pol.label, rc = rc,
                                mu = param[1:n.char], sigma.vec = param[(n.char+1):(2*n.char)])
      ll.val   <- -sum(emp.CP.train * log(logit.cp)) # minus log-likelihood
      return(ll.val)
    }
    # MINIMIZE the negative of the log-likelihood
    ## Initial value is taken to make the parameter value close to the original estimate
    ## "Solution not reliable....Problem Inverting Hessian" if we set the initial sd to zero.
    sol   <- solnp(pars = c(-4, rep(0, n.char-1), rep(0.1, n.char)), fun = ll, 
                   LB = c(rep(-1000, n.char), rep(0, n.char)), 
                   UB = c(rep(1000, n.char), rep(1000, n.char)) ,control=c(inner.iter=1600,tol=1e-6))
    param <- sol[["pars"]]
    names(param) <- c(paste0('mu.', names(rc)), paste0('sd.', names(rc)))
    # In-sample prediction
    mu.label    <- grep('^mu\\.', names(param), value = T)
    sigma.label <- grep('^sd\\.', names(param), value = T)
    pred.in.avg <- pred.rc.macro(dat = dat.train.macro, pol.label = pol.label, rc = rc,
                                 mu = param[mu.label], sigma.vec = param[sigma.label])
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- pred.rc.macro(dat = dat.test.macro, pol.label = pol.label, rc = rc,
                                  mu = param[mu.label], sigma.vec = param[sigma.label])
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
    
  } else if (spec == 'mnl_FE_aggre'){
    
    fe.label        <- grep('^alt_', colnames(dat.train), value = T)
    pol.label       <- setdiff(pol.label, c('X_boattype')) 
    n.char          <- length(pol.label)
    pol.label.fe    <- c(pol.label, fe.label) # augmented characteristics 
    dat.train.macro <- gen.macro.data(dat.train, pol.label.fe)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label.fe)
    
    ll <- function(param) {
      logit.cp <-compute.logitCP(dat = dat.train.macro, pol.label = pol.label.fe, b=param,cs=CS)
      ll.val   <- -sum(emp.CP.train * log(logit.cp)) # minus log-likelihood
      return(ll.val)
    } 
    
    sol   <- solnp(pars = c(rep(0,5)), fun = ll, 
                   LB = c(rep(-1000,5)), 
                   UB = c(rep(1000,5)),control=c(inner.iter=1600,tol=1e-6))
    param <- sol[["pars"]]
    
    names(param) <- pol.label
    # In-sample prediction
    pred.in.avg <- compute.logitCP(dat = dat.train.macro, pol.label = pol.label.fe, b=param,cs=CS)
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- compute.logitCP(dat = dat.test.macro, pol.label = pol.label.fe, b=param,cs=CS)
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
    
  } else if (spec == 'random_coef_FE_aggre'){
    
    fe.label        <- grep('^alt_', colnames(dat.train), value = T)
    pol.label       <- setdiff(pol.label, c('X_boattype')) 
    n.char          <- length(pol.label)
    pol.label.fe    <- c(pol.label, fe.label) # augmented characteristics 
    dat.train.macro <- gen.macro.data(dat.train, pol.label.fe)
    dat.test.macro  <- gen.macro.data(dat.test, pol.label.fe)
    
    if ('X_boattype' %in% pol.label){
      rc = c(X1.0 = 'ln', X0.1 = 'ln', X_boattype = 'n')  
    } else {
      rc = c(X1.0 = 'ln', X0.1 = 'ln')  
    }
    ll <- function(param) {
      logit.cp <- pred.rc.macro(dat = dat.train.macro, pol.label = pol.label.fe, rc = rc,
                                mu = param[1:n.char], sigma.vec = param[(n.char+1):(2*n.char)],fc.val=param[(2*n.char+1):(2*n.char+3)])
      ll.val   <- -sum(emp.CP.train * log(logit.cp)) # minus log-likelihood
      return(ll.val)
    } 
    
    sol   <- solnp(pars = c(-4, rep(0, n.char-1), rep(0.1, n.char),rep(0,3)), fun = ll, 
                   LB = c(rep(-1000, n.char), rep(0, n.char),rep(-1000, 3)), 
                   UB = c(rep(1000, n.char), rep(1000, n.char),rep(1000, 3)) ,control=c(inner.iter=1600,tol=1e-6))
    param <- sol[["pars"]]

    names(param) <- c(paste0('mu.', names(rc)), paste0('sd.', names(rc)))
    # In-sample prediction
    mu.label    <- grep('^mu\\.', names(param), value = T)
    sigma.label <- grep('^sd\\.', names(param), value = T)
    pred.in.avg <- pred.rc.macro(dat = dat.train.macro, pol.label = pol.label.fe, rc = rc,
                                 mu = param[mu.label], sigma.vec = param[sigma.label],fc.val=param[(2*n.char+1):(2*n.char+3)])
    loss.in.avg <- dist.prob(emp.CP.train, pred.in.avg)
    # Out-of-sample prediction
    pred.out.avg <- pred.rc.macro(dat = dat.test.macro, pol.label = pol.label.fe, rc = rc,
                                 mu = param[mu.label], sigma.vec = param[sigma.label],fc.val=param[(2*n.char+1):(2*n.char+3)])
    loss.out.avg <- dist.prob(emp.CP.test, pred.out.avg)
  }
  # Output
  output <- list(param, pred.in.avg, loss.in.avg, pred.out.avg, loss.out.avg)
  names(output) <- c('param', 'pred.in', 'loss.in', 'pred.out', 'loss.out')
  return(output)
}



# Macro nested-logit prediction
pred.nest.macro <- function(dat, pol.label, b, lambda, nests){
  # Input
  #     dat:        macro data
  #     pol.label:  names of characteristics included in the logit model
  #     b:          coefficients
  #     lambda:     nest parameter
  #     nests:      the list of nests
  # Output
  #     cp.vec:     the vector of choice probabilities
  
  CS <- unique(dat$alt)
  n.alt     <- length(CS)
  
  # Denominator of nest probabilities
  nest.prob.num.list <- lapply( nests, function(nest) compute.expsum(dat, pol.label, b, lambda, nest) ) 
  nest.prob.denom    <- Reduce('+', nest.prob.num.list)    
  
  # Within-nest choice probabilities and nest probabilities
  cp.vec        <- rep(NA, n.alt)
  names(cp.vec) <- CS
  for (idx.nest in 1:length(nests)) {
    # Within-nest choice probabilitiy
    nest     <- nests[[idx.nest]]
    logit.cp <- compute.logitCP(dat.input = dat, pol.label, b/lambda, 
                                cs = nest, b.scale = 1)    # logit CP within the nest
    # Nest probability
    nest.p   <- nest.prob.num.list[[idx.nest]]/nest.prob.denom
    # Choice probability
    cp <- logit.cp * nest.p
    cp.vec[nest] <- cp
  }
  
  return(cp.vec)
}



# Auxiliary function used to compute the probabilities of choosing from a nest
compute.expsum <- function(dat, pol.label, b, lambda, nest){
  # Input
  #     dat:       macro data
  #     pol.label: names of characteristics included in the logit model
  #     b:         coefficients
  #     lambda:    nest parameter
  #     nest:      nest (a set of alternatives)
  # Output
  #     expsum:    the sum of exponential terms used to compute nest probabilities
  X         <- as.matrix(dat[which(dat$alt %in% nest), pol.label]) 
  if (length(nest) == 1){
    X <- as.numeric(X)
  }
  lin.idx   <- X %*% b/lambda 
  expsum    <- sum(exp(lin.idx))^lambda
  return(expsum)
}



# Macro random-coefficient prediction
pred.rc.macro <- function(dat, pol.label, rc = NULL, mu = NULL, sigma.vec = NULL, fc.val = NULL, fe = F, n.sim = 100) {
  # Input
  #     dat:       macro data
  #     pol.label: names of characteristics included in the logit model
  #     rc:        named vector to indicate the distribution of random coefficients
  #                (e.g., rc = c(X1.0 = 'ln', X0.1 = 'n') where 'ln' denotes log-normal distribution
  #                 and 'n' denotes normal distribution)
  #     mu:        a vector of the mean of coefficients (assumed to be in the same order as rc)
  #     sigma.vec: a vector of the standard deviation of coefficients (assumed to be in the same order as rc)
  #     fc.val:    the values of fixed coefficients in the same order of characteristics vectors
  #                (assumed to be in the same order of pol.label)
  #     fe:        do we use alternative fixed effects?
  #                TRUE: use only fixed effects. FALSE: use only characteristics 
  #                FEs are assumed to be named as "alt_<the name of the alternative>"
  #     n.sim:     the number of simulations
  # Output
  #     logit.cp:  the vector of choice probabilities
  
  set.seed(1)
  CS     <- unique(dat$alt)
  n.alt  <- length(CS)
  n.char <- length(pol.label)
  
  # Simulate (and fix) the set of coefficients
  rc.label  <- names(rc)
  b.sim     <- matrix(nrow = n.sim, ncol = length(pol.label))
  colnames(b.sim) <- pol.label
  for (idx.char in 1:length(rc)) {
    char <- rc.label[idx.char]
    b.sim[, char] <- rnorm(n.sim, mu[idx.char], sigma.vec[idx.char]) 
    if (rc[char] == 'ln'){ # If log-normal
      b.sim[, char] <- exp(b.sim[, char])
    }
  }
  # Fill in the fixed coefficients
  fc.label <- setdiff(pol.label, rc.label)    # Names of variables with fixed coefficients
  if (length(fc.label) > 0){
    for (idx.char in 1:length(fc.label)) {
      # The column index (idx represents the idx-th regressor with fixed coefficients, which is the idx.col-th regressor)
      char <- fc.label[idx.char]    
      b.sim[, char] <- fc.val[idx.char]
    } 
  }
  # Predict CP
  logit.sim <- matrix(nrow = n.sim, ncol = n.alt)
  for (i in 1:n.sim) {
    b  <- b.sim[i,]
    logit.sim[i,] <- compute.logitCP(dat.input = dat, pol.label, b, cs = CS, b.scale = 1) # logit CP
  }
  logit.cp <- apply(logit.sim, 2, FUN = function(x) mean(x, na.rm = T))
  names(logit.cp) <- CS
  
  return(logit.cp)
}



# Function to conduct the whole empirical analysis
do.empirical.analysis <- function(dat, K, pol.label, specify.scale = NULL,r,data_path){
  # Input
  #     dat:           micro data
  #     K:             we do K-fold cross-validation
  #     pol.label:     names of characteristics
  #     specify.scale: scale of separating vectors. 
  #                    specify.scale = NULL means we search for and use the smallest scale for perfect in-sample fit.
  # Output
  #     output:        results (results of empirical exercises) and samples (split samples used)
  
  cat('\n')
  cat('########################\n')
  cat(sprintf('%s-fold cross-validation', K))
  cat('\n########################') 
  cat('\n')
  
  # Split the sample
  samples <- split.sample(dat, K = K,r,data_path)
  
  
  # to store results
  cv.results.list <- list()
  for (k in 1:K) {
    # Indicate the sample number
    cat('\n')
    cat('#########\n')
    cat(sprintf('Sample %s', k))
    cat('\n#########') 
    cat('\n')
    
    # Take the k-th sample and generate average characteristics with d = 1
    dat.train <- samples[[k]][['train']]    
    dat.test  <- samples[[k]][['test']]

    
    dat.train <- prepare.mlogit.dat(dat.train, pol.label)    # Replace individual characteristics by averages
    dat.test  <- prepare.mlogit.dat(dat.test, pol.label)    
    
    
    dat.train.avg=gen.macro.data(dat.train,pol.label )
    dat.test.avg=gen.macro.data(dat.test,pol.label )
 
    # In-sample fit
    ## Compute the empirical choice probability
    emp.CP.train <- compute.emp.CP(dat.train)
    
    # Out-of sample fit
    ## Out-of-sample empirical choice probabilities
    emp.CP.test <- compute.emp.CP(dat.test)
    ## Out-of-sample predictions using estimates above and test data
    
    
    # Store the results in lists
    pred.in <- list()    # list of in-sample predictions
    loss.in <- list()    # list of in-sample prediction errors
    pred.out <- list()   # list of out-of-sample predictions
    loss.out <- list()
    
    
    #######################
    # Alternative methods #
    #######################
    dat.train[,c(4,5)] <- samples[[k]][['train']][,c(4,5)] #replace average characeristic with individual characteristics
    dat.test[,c(4,5)] <- samples[[k]][['test']][,c(4,5)] #replace average characeristic with individual characteristics
    
    # Estimate each model and evaluate its in-sample and out-of-sample prediction performance
    alt.specs <- c('mnl_aggre', 'nest_charter_aggre', 'nest_boat_aggre', 'random_coef_aggre',
                   'mnl_FE_aggre','random_coef_FE_aggre')
    regs      <- list()
    for (spec in alt.specs) {
      # Regression and prediction
      regs[[spec]]     <- run.alt.regressions.macro(dat.train, dat.test, 
                                                    pol.label, emp.CP.train, emp.CP.test, spec)
      pred.in[[spec]]  <- regs[[spec]][['pred.in']]
      loss.in[[spec]]  <- regs[[spec]][['loss.in']]
      pred.out[[spec]] <- regs[[spec]][['pred.out']]
      loss.out[[spec]] <- regs[[spec]][['loss.out']]
    }
    
    ########################################
    # Store the results from the k-th sample
    ########################################
    cv.results.list[[k]] <- list(pred.in, loss.in, pred.out, loss.out)
    names(cv.results.list[[k]]) <- c('in-sample_pred', 'in-sample_loss', 
                                     'out-of-sample_pred', 'out-of-sample_loss')
  }
  output <- list(cv.results.list, samples)
  names(output) <- c('results', 'samples')
  return(output)
}
