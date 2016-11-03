#fast read table function
.fast.read.table = function(file, nrow){
  # read contents
  countC <- scan( file = file, what = character(), sep = "\t", skip = 0 )
  countC1 <- matrix( countC, nrow, byrow = TRUE )
  return(countC1)
} 
