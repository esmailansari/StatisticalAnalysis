# selecting prime numbers below n
prime = function(n)
{
  n = as.integer(n)
  if(n > 1e8) stop("n too large")
  primes = rep(TRUE, n)
  primes[1] = FALSE
  last.prime = 2L
  fsqr = floor(sqrt(n))
  while (last.prime <= fsqr)
  {
    primes[seq.int(2L*last.prime, n, last.prime)] = FALSE
    sel = which(primes[(last.prime+1):(fsqr+1)])
    if(any(sel)){
      last.prime = last.prime + min(sel)
    }else last.prime = fsqr+1
  }
  which(primes)
}

# ------------------------------------------------
# Creating van de Corput method
vdc <- function(b,N) {
 #b is not vector here.
 out = rep(0,N)
 numd = ceiling(log(N)/log(b))
 bb = 1/b^(1:numd)
 a = NULL;
 out[1] = 0;
 
 for (i in 2:N){
   a = nbe(a,b)
   fa = rev(a)
   out[i] = sum(fa*bb[1:length(a)])
 }
 
 return(out)
}

# ----------------------------------------------------
#next b-ary number estimation
nbe <- function(a,b) {
  
  numd = length(a);
  na = a;
  carry = TRUE; 
  if (!is.null(a)) {
  for (i in numd:1) {
 
    if (carry) {
      
      if (a[i] == (b-1)) {na[i]=0
      } else {
        na[i] = a[i] + 1
        carry = FALSE
      }    
    } 
  }
  }
  if (carry) na = c(1,na)
  return(na)
}

# ---------------------------------------------
# creating Halton sequence 
halton <- function(b,N) {
  
  dim = length(b)
  out = matrix(rep(0,N*dim),c(N,dim))
  for (i in 1:dim) 
    out[,i] = vdc(b[i],N)
  return(out)
  
}


# --------------------------------------------------

hammersley <- function(b,N){
  
  dim = length(b) 
  out = matrix(rep(0,N*dim),c(N,dim))
  out[2:N,2:dim] = halton(b[1:dim-1],N-1)
  out[,1] = (0:(N-1))/N
  return(out)
  
}