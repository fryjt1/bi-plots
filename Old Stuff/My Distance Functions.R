library(GenSA)
library(flexclust)

Cos.Dist = function(X.mat)
{
  XXt = X.mat%*%t(X.mat)
  lengths = sqrt(rowSums(X.mat^2))
  temp = sweep(XXt,2,lengths,FUN="/")
  sim = sweep(temp,1,lengths,FUN="/")
  return(1-sim)
}

Euclidean.Dist = function(X.mat)
{
  xtx = rowSums(X.mat^2)
  XXt = X.mat%*%t(X.mat)
  temp = sweep(-2*XXt,2,xtx,FUN="+")
  temp2 = sweep(temp,1,xtx,FUN="+")
  temp2[which(abs(temp2)<1e-7)] = 0
  dist = sqrt(temp2)
  
  return(dist)
}

Pythag.Dist = function(X.mat)
{
  xtx = rowSums(X.mat^2)
  XXt = X.mat%*%t(X.mat)
  temp = sweep(-2*XXt,2,xtx,FUN="+")
  dist = sweep(temp,1,xtx,FUN="+")
  
  return(dist)
}

PCA.Dist = function(X.mat)
{
  dist = X.mat%*%t(X.mat)
  return(dist)
}

Cos.Dist2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  Xy = X.mat%*%y.vec
  y.length = sqrt(sum(y.vec^2))
  X.lengths = matrix(sqrt(rowSums(X.mat^2)),ncol=1)
  sim = (Xy/X.lengths)/y.length
  
  return(1-sim)
}

Euclidean.Dist2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  xtx = matrix(rowSums(X.mat^2),ncol=1)
  Xty = X.mat%*%y.vec
  yty = sum(y.vec^2)
  dist = sqrt(xtx - 2*Xty + yty)
  
  return(dist)
}

Pythag.Dist2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  xtx = matrix(rowSums(X.mat^2),ncol=1)
  Xty = X.mat%*%y.vec
  yty = sum(y.vec^2)
  dist = xtx - 2*Xty + yty
  
  return(dist)
}

PCA.Dist2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  dist = X.mat%*%y.vec
  return(dist)
}

Min.Euclidean2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  ata = sum(y.vec^2)
  ctc = rowSums(X.mat^2)
  cta = X.mat%*%y.vec
  dist = sqrt(ctc - (cta)^2/ata)
  
  return(dist)
}

Min.Pythag2 = function(X.mat,y.vec)
{
  y.vec = matrix(y.vec,ncol=1)
  ata = sum(y.vec^2)
  ctc = rowSums(X.mat^2)
  cta = X.mat%*%y.vec
  dist = ctc - (cta)^2/ata
  
  return(dist)
}

Calc.Min.Dist = function(X.matrix,axis.vector) # This uses a sequence...not a closed form answer
{
  my.seq = seq(-5,5,length=100)
  my.axis = matrix(axis.vector,length(axis.vector),1)
  test.seq = apply(my.axis,1,function(c) c*my.seq)
  
  all.dist = apply(X.matrix,1,function(c) dist2(matrix(c,1,ncol(X.matrix)),test.seq))
  min.dist = apply(all.dist,2,min)
  
  return(matrix(min.dist,ncol=1))
}

Euclidean.MDS = function(X.matrix)
{
  n = nrow(X.matrix)
  x.start = c(rnorm(2*n))
  HD.dist = Euclidean.Dist(X.matrix)
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(-4,2*n), upper = rep(4,2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Cos.MDS = function(X.matrix)
{
  n = nrow(X.matrix)
  x.start = c(rnorm(2*n))
  HD.dist = Cos.Dist(X.matrix)
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(-4,2*n), upper = rep(4,2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Euclidean.Stress = function(x,X,Z,a)
{
  a = matrix(a,ncol=1)
  x = matrix(x,ncol=1)
  HD.dist = Euclidean.Dist2(X,a)
  LD.dist = Euclidean.Dist2(Z,x)
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}

PCA.Stress = function(x,X,Z,a)
{
  a = matrix(a,ncol=1)
  x = matrix(x,ncol=1)
  HD.dist = X%*%a
  LD.dist = Z%*%x
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}