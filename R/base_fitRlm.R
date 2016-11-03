.fitRlm = function(data_all_combined, filterout, object){
  dnaseKnots <- object[['dnaseKnots']]
  dnaseThres <- object[['dnaseThres']]
  dnaseThres[1] <- 1
  # b: piecewice linear B-spline bases (matrix)
  b <- bs(dnaseThres, degree=1, knots=dnaseKnots)
  # data_id: group indicator of histone summary statistics (35_0_0 -> 35 0 0)
  if(object[['dataNum']] > 1){ 
    data_id <- t(apply(as.matrix(data_all_combined[, 1]), 1, function(x){strsplit(as.character(x), '_')[[1]]}))
    }else{
    data_id <- as.matrix(data_all_combined[, 1], ncol=1)
  }
  if(length(filterout) == 0){     
    # b_data: Spline part of the model of the Dnase read counts
    b_data <- b[match(data_id[, 1], 0:(length(dnaseThres)-1)),]
    # y: average ChIP read count
    y <- as.numeric(data_all_combined[, 3])/as.numeric(data_all_combined[, 2])
    data_x <- cbind(as.data.frame(b_data), as.data.frame(data_id[, -1]))
    y <- log(y+1e-4)
    data_fit <- cbind(as.data.frame(y), as.data.frame(b_data), as.data.frame(data_id[, -1]))
    colnames(data_fit) <- c('y', paste('v', 1:(dim(data_fit)[2]-1), sep=''))
    b_data_n <- dim(b_data)[2]
    h_data_n <- object[['dataNum']]-1
    # xvars: "model matrix"
    xvars <- paste('v', 1:(b_data_n+h_data_n), sep='')
    t <- paste(xvars, collapse=" + ")
    model <- paste('y~', t, sep='')
    data_fit <- cbind(as.data.frame(data_fit), weights=as.numeric(data_all_combined[, 2]))
    # Fit model
    m <- lm(model, data = data_fit, weights=weights)
    y <- cbind(data_all_combined[, 1], exp(m$fitted.values))
  }else{
    ##filter out not important variables
    select_histone <- (1:(object[['dataNum']]-1))[-filterout]
    tt <- data_id[, 1]
    for(i in select_histone){
      tt <- paste(tt, data_id[, 1+i], sep='_')
    }
    data_all_new <- data_all_combined
    data_all_new[, 1] <- tt
    data_all_combined <- .aggregate_data(data_all_new)
    data_id <- t(apply(as.matrix(data_all_combined[, 1]), 1, function(x){strsplit(as.character(x),'_')[[1]]}))
    b_data <- b[match(data_id[, 1], 0:(length(dnaseThres)-1)), ]
    y <- as.numeric(data_all_combined[, 3])/as.numeric(data_all_combined[, 2])
    y <- log(y+1e-4)
    data_fit <- cbind(as.data.frame(y), as.data.frame(b_data), as.data.frame(data_id[, -1]))
    colnames(data_fit) <- c('y', paste('v', 1:(dim(data_fit)[2] - 1), sep=''))
    b_data_n <- dim(b_data)[2]
    h_data_n <- object[['dataNum']] - 1 - length(filterout)
    xvars <- paste('v', 1:(b_data_n + h_data_n), sep = '')
    t <- paste(xvars, collapse = " + ")
    model <- paste('y~', t, sep = '')
    data_fit <- cbind(as.data.frame(data_fit), weights = as.numeric(data_all_combined[, 2]))
    m <- lm(model, data = data_fit, weights = weights)
    y <- cbind(data_all_combined[, 1], exp(m$fitted.values))
  }
  return(y)
}
