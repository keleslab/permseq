.quantileFunction = function(out, p1){
  total <- out[dim(out)[1], 2]
  back <- out[dim(out)[1], 1]
  out <- out[-dim(out)[1], ]
  out <- out[order(out[, 1], decreasing = F), ]

  out_cum <- c(0,cumsum(out[, 2])) + back
  status1 <- status2 <- 0
  j <- floor(total * p1)
  g <- total * p1 - j

  x1 <- out[min(which(out_cum >= j)) - 1, 1]
  x2 <- out[min(which(out_cum >= j + 1)) - 1, 1]

  if(g == 0 & j%%2 == 0){
    r <- 0
  }else{
    r <- 1
  }
  q1 <- (1 - r) * x1 + r * x2
  if(length(q1) == 0){
    q1 <- 2
  }
  return(q1)
}
