##########################################################################################

#File name: grfit_function

#Purpose: The function of proposed Rule ensemble method 

#Author: Ke Wan

#R Version: R-4.2.1

#Input : data and the parameters of proposed method

#Output: object of "grfit"

#Required R packages: rpart_4.1.19, grpreg_3.4.0

#########################################################################################
library(rpart)
library(grpreg)
#----------------------------------------------------------------------------------------
#The function of proposed methods
#----------------------------------------------------------------------------------------
grfit <- function(
  data,				#input data frame (y, g, x)
  meandepth = 2, 		#the mean depth of each tree-base function 
  learnrate = 0.01,   #shrinkage rate for each boosting steps    
  ntrees = 500, 		#number of trees
  weight = rep(1,nrow(data)),  # inverse probability weight for each individual                         
  sampfrac = NULL,    #the fraction size of training sample for rule generation
  wins = 0.025,		#the quantile of winsorized
  nfolds = 10         #nfold cross-validation for estimating the parameter of group lasso
){
  
  data <- data.frame(data)
  y_names <- colnames(data)[1]			#outcome name
  g_names <- colnames(data)[2]			#treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		#covariates names
  
  #Calculate the transformed outcome and create the transformed outcome data
  y <- data[,y_names]
  g <- data[,g_names]
  x <- data[,x_names]
  
  data.learn <- data.frame(y = y,g = g,x)
  
  #create the sample fraction for training the model
  n <- nrow(data)
  if(is.null(sampfrac)){
    size <- min(n/2,100 + 6*sqrt(n))
  }else{
    size <- sampfrac*n
  }			   					
  subsample <- list()
  subsample <- mapply(function(i)
    #bootstrap
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    }
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
    ,i=1:ntrees,SIMPLIFY=FALSE)
  
  #------------------------------------------------------------------------------------------
  #Step 1: Rule Generation using GBT
  #------------------------------------------------------------------------------------------
  rules <- c() #Initialize the rule 
  data.learn <- data.learn #Initialize the training data
  
  #initialize the memory function eta
  eta0 <- mean(data.learn$y)
  eta <- rep(eta0, length(data$y))
  
  #pseudo residual
  data.learn$y <- y - eta0
  
  #the number of treminal nodes for each tree
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1/(2^meandepth - 2))), base = 2))
  
  #GBM(gradient boosting tree) 
  
  for(i in 1:ntrees){
    
    tree.control <- rpart::rpart.control(maxdepth = maxdepth[i])
    
    #built tree base learner
    tree <- rpart::rpart(y~., control = tree.control,data = data.learn[subsample[[i]], ],weight = weight[subsample[[i]]])
    
    #extract the rules of built tree
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
    paths <- paths[-1]
    rules <- c(rules,paths)
    
    #upgrade the memory function eta
    eta <- eta + learnrate * predict(tree, newdata = data.learn)
    
    #update the pseudo residual 
    data.learn$y <- (y - eta)
  }
  
  #Sub_function 1: remove the duplicate and complement rules
  fit.rule <- remove.duplicate.complement(data,rules)
  
  #create the set of rules for ensemble
  rules <- fit.rule$rules										#	all rules
  rules_pred_all <- rules[grep("g",rules)] 					#	predict rules
  rules_pred <- gsub(" & g< 0.5","",rules_pred_all)
  rules_pred <- gsub("g< 0.5 & ","",rules_pred)  
  rules_pred <- gsub(" & g>=0.5","",rules_pred)
  rules_pred <- gsub("g>=0.5 & ","",rules_pred)
  rules_pred <- unique(rules_pred)
  rules_prog <- rules[!is.element(rules,rules_pred_all)]		#	prognostic rules
  
  if (length(rules_pred)!= 0){
    
    rulevars <- rule.predict(data,rules)
    rulevars_pred <- rule.predict(data,rules_pred)
    rulevars_prog <- rule.predict(data,rules_prog) 
    
    #Sub_function 2 : transform rules into binary values
    colnames(rulevars) <- rules
    colnames(rulevars_pred) <- rules_pred
    colnames(rulevars_prog) <- rules_prog
    rules_pred_nums <- ncol(rulevars_pred)
    rules_prog_nums <- ncol(rulevars_prog)
    
  }else {	
    
    rulevars <- rule.predict(data,rules)
    rulevars_prog <- rule.predict(data,rules_prog) 
    
    colnames(rulevars) <- rules
    colnames(rulevars_prog) <- rules_prog
    rules_prog_nums <- ncol(rulevars_prog)
  }	
  
  #------------------------------------------------------------------------------------------
  #Step 2: Rule ensemble
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------
  #"winsorize" and normalize the linear terms (Freideman and Popescu, 2008)
  #------------------------------------------------------------------------
  
  #create the "winsorize" version of linear terms 
  data.win <- winsorize(data,wins = wins) #Sub_function 3 : the winsorized function
  cut_point <- data.win[[2]]
  data.win <- data.win[[1]]
  
  dat.linear <- data.win[,-c(1,2)]
  linear.nums <- ncol(dat.linear)
  
  # Rulefit group lasso
  linear_id <- 1:linear.nums
  prog_id <- (length(linear_id) + 1):(length(linear_id) + length(rules_prog))
  pred_id0 <- (length(linear_id) + length(rules_prog) + 1):(length(linear_id) + length(rules_prog) + length(rules_pred))
  pred_id1 <- (length(linear_id) + length(rules_prog) + 1):(length(linear_id) + length(rules_prog) + length(rules_pred))
  
  
  if (length(rules_pred)!= 0){
    
    group_id <- c(linear_id,prog_id,pred_id0,pred_id1)
    
    linear.names <- colnames(x)
    prog.names <- rules_prog
    
    pred.names0 <- paste0(rules_pred,"_c")
    pred.names1 <- paste0(rules_pred,"_t")
    g.names <- c(linear.names,prog.names,pred.names0,pred.names1)
    
    rulevars_pred0 <- rulevars_pred*(data$g-1)/2
    rulevars_pred1 <- rulevars_pred*(data$g)/2
    
    g.var <- cbind(dat.linear,rulevars_prog,rulevars_pred0,rulevars_pred1)
    colnames(g.var) <- g.names
    
    fit.cv <- grpreg::cv.grpreg(g.var, data$y, group = group_id, penalty="grLasso",nfolds = nfolds)
    lambda <- fit.cv$lambda.min
    gr.fit <- grpreg::grpreg(g.var, data$y, group = group_id, penalty="grLasso",lambda=lambda)
    
    # prognostic effect
    coef <- coef(gr.fit)[-1]
    coef.prg <- coef[c(linear_id,prog_id)]
    names(coef.prg) <- rules_prog
    coef.prg <- coef.prg[coef.prg!=0]
    
    # predictive effect/treatment effect
    coef.trt0 <- coef[grep("_c",names(coef))]
    coef.trt1 <- coef[grep("_t",names(coef))]
    names(coef.trt0) <- rules_pred
    names(coef.trt1) <- rules_pred
    coef.trt <- (coef.trt1 + coef.trt0)/2
    coef.trt <- coef.trt[coef.trt!=0]
    
  }else{
    
    group_id <- c(linear_id,prog_id)
    
    linear.names <- colnames(x)
    prog.names <- rules_prog
    
    g.var <- cbind(dat.linear,rulevars_prog)
    g.names <- c(linear.names,prog.names)
    colnames(g.var) <- g.names
    
    fit.cv <- grpreg::cv.grpreg(g.var, data$y, group = group_id, penalty="grLasso",nfolds = nfolds)
    lambda <- fit.cv$lambda.min
    gr.fit <- grpreg::grpreg(g.var, data$y, group = group_id, penalty="grLasso",lambda=lambda)
    
    # prognostic effect
    coef <- coef(gr.fit)[-1]
    coef.prg <- coef[c(linear_id,prog_id)]
    names(coef.prg) <- rules_prog
    coef.prg <- coef.prg[coef.prg!=0]
    #add 20221208
    coef.trt <- 0
  }	
  
  rule_all <- list(rules_prog = rules_prog,rules_pred = rules_pred)
  coef_all <- coef(gr.fit)
  
  res <- list(model = gr.fit,		#	the model of proposed method
              mu = coef.prg,				#	estimsted mean effect
              tau = coef.trt, 			#	estimated HTE			
              rules = rule_all 			#	all the created rules in rules generation 
  )
  
  return(res)
}			


#--------------------------------------------------------------------------------------------------
#Sub_function 1: remove the duplicate and complement rules (refer the code of R package "pre" )
#--------------------------------------------------------------------------------------------------

remove.duplicate.complement <- function(data,rules,duplicate = TRUE, complement = TRUE){
  
  #remove the duplicate rules
  if(duplicate == TRUE){
    
    #calculate the binary value for each rules (F2)
    rulevars <- rule.predict(data,rules)
    
    #remove the rules which have same binary value
    duplicates <- duplicated(rulevars,MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[,!duplicates,drop = FALSE]
    rules <- rules[!duplicates]
    
  } else {
    
    duplicates.removed <- NULL
    
  }
  
  #remove the complement rules
  if(complement == TRUE){
    
    # find columns with equal variance to reduce the number of comparisons
    vars <- apply(rulevars, 2, var_bin) # The function to calculate the variance
    vars_distinct <- lapply(
      unique(vars), function(x) { 
        idx <- which(is_almost_eq(x, vars)) # The function to check near equality
        list(var = x, n = length(idx), idx = idx)
      }) # 
    
    #store the complements
    complements <- logical(ncol(rulevars))
    
    #identify the complements
    for (va in vars_distinct) {
      
      #if no variance same to this columns then next columns
      if(va$n < 2L) next
      
      idx <- va$idx
      idx <- setdiff(idx, which(complements))
      if(length(idx) < 2)next
      
      n_idx <- length(idx)
      for(j in 1:(n_idx - 1)){
        if (complements[idx[j]])next
        this_val <- rulevars[, idx[j]]
        is_compl <- which(apply(rulevars[, idx[(j + 1):n_idx], drop = FALSE], 2, function(x) all(x != this_val))) + j
        if (length(is_compl) > 0)
          complements[idx[is_compl]] <- TRUE
      }
    }
    
    complements <- which(complements)
    complements.removed <- rules[complements]
    
    if (length(complements) > 0)
      rules <- rules[-complements]
    rulevars <- rulevars[,-complements,drop = FALSE]
    
  } else {
    
    complements.removed <- NULL
    
  }
  
  if (length(duplicates.removed)!= 0|length(complements.removed) != 0) {
    return(list(rules = rules, rulevars = rulevars,
                duplicates.removed = duplicates.removed,
                complements.removed = complements.removed))
  }
  
  list(rules = rules, rulevars = rulevars)
}

#standard error of rule terms
var_bin <- function(x) {
  p <- mean(x)
  p*(1L-p)
}

#check near equality
is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(x), length(x) == 1L)
  x_abs <- abs(x)
  xy <- if (x_abs > tolerance) {abs(x - y) / x_abs} else {abs(x - y)}
  xy <= tolerance
}

#---------------------------------------------------------------------------------------------------
#Sub_function 2: transform the rules into binary value (0: against the rules; 1: obey the rules)
#---------------------------------------------------------------------------------------------------	
rule.predict <- function(data,rules){
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr,data)
  colnames(x) <- names(rules)
  x <- apply(x,2,as.numeric)
  return(x)
}

#-----------------------------------------------------------------
#Sub_function 3: the winsorized function
#-----------------------------------------------------------------	
winsorize <- function(data,wins = 0.025){
  
  data.all <- data
  data <- data[,-c(1:2)]
  cut.point <- matrix(0,2,ncol(data))
  data <- sapply(1:ncol(data),function(j,data,wins){
    if(length(unique(data[,j]))> 1){		
      lb <- quantile(data[,j],prob = wins)
      ub <- quantile(data[,j],prob = 1 - wins)
      data[which(data[,j] < lb),j] <- lb
      data[which(data[,j] > ub),j] <- ub
      data.w <- data[,j]
    }else{
      data.w <- data[,j]
    }		
    return(data.w)
  },data,wins)
  
  lb.point <- apply(data,2,min)
  ub.point <- apply(data,2,max)
  win_cut <- rbind(lb = lb.point,ub = ub.point)
  
  colnames(data) <- colnames(data.all)[-c(1:2)]	
  data <- data.frame(y=data.all$y,g=data.all$g,data)	
  return(list(data = data,win_cut = win_cut))
}


#-------------------------------------------------------------------
#Sub_function 4: Predict function
#-------------------------------------------------------------------
predict.grfit <- function(object,data){
  
  #transform the rules into binary values based on the input data frame 
  rules_prog <- object$rules[[1]]
  rules_pred <- object$rules[[2]]
  
  if(length(rules_pred)!= 0){
    
    rulevars_pred <- rule.predict(data,rules_pred)
    rulevars_prog <- rule.predict(data,rules_prog) 
    
    #Sub_function 2 : transform rules into binary values
    colnames(rulevars_pred) <- rules_pred
    colnames(rulevars_prog) <- rules_prog
    rules_pred_nums <- ncol(rulevars_pred)
    rules_prog_nums <- ncol(rulevars_prog)
    
    #normalized the covariates of input data based on the standard deviation for that of training data
    dat.linear <- data[,-c(1:2)]
    
    linear.names <- colnames(data[,-c(1:2)])
    prog.names <- rules_prog
    
    pred.names0 <- paste0(rules_pred,"_c")
    pred.names1 <- paste0(rules_pred,"_t")
    g.names <- c(linear.names,prog.names,pred.names0,pred.names1)
    
    rulevars_pred0 <- rulevars_pred*(data$g - 1)/2
    rulevars_pred1 <- rulevars_pred*(data$g)/2
    
    data.can0 <- data.can1 <- data
    data.can0$g <- rep(0,nrow(data))
    data.can1$g <- rep(1,nrow(data))
    rulevars_pred00 <- rulevars_pred*(data.can0$g - 1)/2
    rulevars_pred01 <- rulevars_pred*(data.can0$g)/2
    rulevars_pred11 <- rulevars_pred*(data.can1$g)/2
    rulevars_pred10 <- rulevars_pred*(data.can1$g - 1)/2
    
    g.var <- cbind(dat.linear,rulevars_prog,rulevars_pred0,rulevars_pred1)
    g.var0 <- cbind(dat.linear,rulevars_prog,rulevars_pred00,rulevars_pred01)
    g.var1 <- cbind(dat.linear,rulevars_prog,rulevars_pred10,rulevars_pred11)
    
    colnames(g.var) <- g.names
    
    y <- predict(object$model,as.matrix(g.var),type = "link")
    trt0 <- predict(object$model,as.matrix(g.var0),type = "link")
    trt1 <- predict(object$model,as.matrix(g.var1),type = "link")
    
  }else{
    
    rulevars_prog <- rule.predict(data,rules_prog)
    colnames(rulevars_prog) <- rules_prog
    rules_prog_nums <- ncol(rulevars_prog)
    
    #normalized the covariates of input data based on the standard deviation for that of training data
    dat.linear <- data[,-c(1:2)]
    
    linear.names <- colnames(data[,-c(1:2)])
    prog.names <- rules_prog
    
    g.names <- c(linear.names,prog.names)
    
    g.var <- cbind(dat.linear,rulevars_prog)
    
    colnames(g.var) <- g.names
    
    y <- predict(object$model,as.matrix(g.var),type = "link")
    trt0 <- rep(0,nrow(g.var))
    trt1 <- rep(0,nrow(g.var))
  }
  
  return(list(y = y, trt0 = trt0,trt1 = trt1))
}

#-------------------------------------------------------------------------------------------------
#Sub_function : transform the rules into binary value (0: against the rules; 1: obey the rules)
#--------------------------------------------------------------------------------------------------	
rule.predict <- function(data,rules){
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr,data)
  colnames(x) <- names(rules)
  x <- apply(x,2,as.numeric)
  return(x)
}
