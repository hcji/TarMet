peakshape <- function(x, p, w1, w2=NULL, asym=FALSE, normalize=TRUE){
  if (is.null(w2)) {
    w2 <- w1
  }
  
  y <- rep(0, length(x))
  
  # Gaussian
  if (!asym) {
    y <- exp(-((x-p)/(0.6005612*w1))^2)
  }
  
  # Half-Gaussian & Half Lorentzian
  if (asym) {
    y[x<p] <- exp(-((x[x<p]-p)/(0.6005612*w1))^2)
    y[x>=p] <- 1/(1+((x[x>=p]-p)/(0.5*w2))^2)
  }
  
  if (normalize){
    y <- y/sum(y)
  }

  return(y)
}

conv <- function(y, W1 = seq(1,32,1), W2=NULL, asym=FALSE) {
  if (is.null(W2)) {
    W2 <- W1
  }
  
  res <- lapply(seq_along(W1), function(i) {
    w1 <- W1[i]
    w2 <- W2[i]
    point <- min(length(y), 10*w1, 10*w2)
    d <- -point : point
    m <- peakshape(d, p=0, w1=w1, w2=w2, asym=asym, normalize=TRUE)
    l <- (length(d)-1)/2
    r <- convolve(y, m, type='open')[(l+1):(length(y)+l)]
    r
  })
  
  res <- do.call(cbind, res)
  colnames(res) <- W1
  return(res)
}

matcor <- function(A,B){
  a <- A-mean(A)
  b <- B-mean(B)
  return(sum(a*b)/sqrt(sum(a*a)*sum(b*b)))
}

