.run_grplasso = function(file, dnaseKnots, dnaseThres, chrlist, index, num_require){

  data_all <- vector('list', length(chrlist))
  for(i in 1:length(chrlist)){
    chr.i <- chrlist[i]
    data_all[[i]] <- .fast.read.table(paste(chr.i, file, sep='_'), 3)
  }
  names(data_all) <- chrlist
  .grplasso_difflambda(dnaseKnots, dnaseThres, data_all, index, num_require)
}


.grplasso_difflambda <- function(dnaseKnots, dnaseThres, data_all, index, num_require){
  dnaseThres[1] <- 1
  # b: piecewice linear B-spline bases (matrix)
  b <- bs(dnaseThres, degree=1, knots=dnaseKnots) 
  chrlist <- names(data_all)
  N <- length(chrlist)
  temp <- NULL
  #combine all the chromosome together-length(data_all) = number of chromosome
  for(i in 1:length(data_all)){
    temp <- rbind(temp, t(data_all[[i]]))
  }
  #
  data_all_combined <- .aggregate_data(temp)
  data_id <- t(apply(as.matrix(data_all_combined[, 1]), 1, function(x){strsplit(as.character(x), '_')[[1]]}))
  b_data <- b[match(data_id[, 1], 0:(length(dnaseThres)-1)), ] #order the dnase spline value in the order of dnaseThres
  # Construct Histone part of the model matrix
  aa <- data_id[, -1]
  aa <- as.data.frame(aa)
  d <- model.matrix(~., data=aa) #Design matrix for histone
  # y: average ChIP read count
  y <- as.numeric(data_all_combined[, 3])/as.numeric(data_all_combined[, 2])
  # data_x the design matrix for DNase and Histone
  data_x <- cbind(b_data, d)
  # Maximal value of penalty parameter lambda
  lambda_max <- lambdamax(x=data_x, y=log(y+1e-4), index=index, model=LinReg(), weights = as.numeric(data_all_combined[, 2]), lambda=lambda, center=F, standardize=F)
  lambda <- lambda_max*unique(c(seq(from=1, to=0.1, length.out=100), seq(0.1, 0.01, length.out=10)))
  # Select most informative Histone
  m <- grplasso(x=data_x, y=log(y+1e-4), index=index, model=LinReg(), weights = as.numeric(data_all_combined[, 2]), lambda=lambda[1], center=F, standardize=F, control=grpl.control(trace=0))
  choose <- which(m$coefficients!=0)
  i <- 2
  while(length(which(m$coefficients!=0))<=(num_require) & i<length(lambda)){
    m <- grplasso(x=data_x, y=log(y+1e-4), index=index, model=LinReg(), weights = as.numeric(data_all_combined[,2]), lambda=lambda[i], center=F, standardize=F, control=grpl.control(trace=0))
    i <- i+1
    choose_temp <- which(m$coefficients != 0)
    if(length(setdiff(choose_temp,choose)) != 0){
      choose=c(choose, setdiff(choose_temp, choose))
    }    
  }
  # choose: chosen parameters from group lasso
  #choose=choose[order(choose)]
  # remove spline part and intercept to have only the selected histones

  choose <- choose[!choose %in% 1:5]
  seq <- seq(from=1, to=length(choose), by=2)
  choose <- choose[seq]
  choose <- (choose-6)/2+1
  return(list(choose=choose, data_all_combined=data_all_combined))
}


.aggregate_data = function(temp){
  temp <- as.data.frame(temp)
  colnames(temp) <- c('V1','V2','V3')
  data_aggr_all1 <- aggregate(as.numeric(as.character(V2))~., data=temp[, -3], FUN=sum)
  data_aggr_all2 <- aggregate(as.numeric(as.character(V3))~., data=temp[, -2], FUN=sum)
  data_aggr <- cbind(as.character(data_aggr_all1[, 1]), data_aggr_all1[, 2], data_aggr_all2[, 2])
  colnames( data_aggr) <- c('id', 'total positions', 'total chip')
  return(data_aggr)
}


