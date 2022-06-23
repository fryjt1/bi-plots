library(GenSA)
library(flexclust)
library(lhs)

W.Euclidean = function(w,X.mat)
{
  X.curl = X.mat%*%diag(w)
  D = as.matrix(dist(X.curl))
  return(D)
}

W.Euclidean2 = function(w,X.mat,y.vec)
{
  X.mat = X.mat%*%diag(w)
  y.vec = matrix(y.vec*w,ncol=1)
  
  xtx = matrix(rowSums(X.mat^2),ncol=1)
  Xty = X.mat%*%y.vec
  yty = sum(y.vec^2)
  dist = sqrt(xtx - 2*Xty + yty)
  
  return(dist)
}
  
W.Euclidean.stress = function(x,X.mat,HD.dist)
{
  D.curl = W.Euclidean(x,X.mat)
  stress = sum((HD.dist-D.curl)^2)/sum(HD.dist^2)
  return(stress)
}


Cos.Dist = function(X.mat)
{
  XXt = X.mat%*%t(X.mat)
  lengths = sqrt(rowSums(X.mat^2))
  temp = sweep(XXt,2,lengths,FUN="/")
  sim = sweep(temp,1,lengths,FUN="/")
  return(1-sim)
}

Acos.dist = function(X.mat)
{
  n = nrow(X.mat)
  D = matrix(0,n,n)
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      D[i,j] <- D[j,i] <- acos(sum(sqrt(X[i,]*X[j,])))
    }
  }
  return(D)
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

Euclidean.MDS.SA = function(X.matrix)
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

General.MDS.SA = function(D,bounds=c(-4,4))
{
  n = nrow(D)
  x.start = c(rnorm(2*n))
  HD.dist = D
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(bounds[1],2*n), upper = rep(bounds[2],2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Sqrt.Manhattan.MDS.SA = function(X.matrix)
{
  n = nrow(X.matrix)
  x.start = c(rnorm(2*n))
  HD.dist = sqrt(as.matrix(dist(X.matrix,method="manhattan")))
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = as.matrix(dist(Z.mat,method="euclidean"))
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(-10,2*n), upper = rep(10,2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Minkowski.MDS.SA = function(X.matrix,my.power)
{
  n = nrow(X.matrix)
  x.start = c(rnorm(2*n))
  HD.dist = as.matrix(dist(X.matrix,method="minkowski",p=my.power))
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = as.matrix(dist(Z.mat,method="euclidean"))
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(-8,2*n), upper = rep(8,2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Manhattan.MDS.SA = function(X.matrix)
{
  n = nrow(X.matrix)
  x.start = c(rnorm(2*n))
  HD.dist = as.matrix(dist(X.matrix,method="manhattan"))
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = as.matrix(dist(Z.mat))
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  my.bound = 12
  temp = GenSA(lower = rep(-my.bound,2*n), upper = rep(my.bound,2*n), fn = my.fun,
               control=list(verbose=FALSE,smooth=TRUE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Euclidean.MDS = function(X.matrix,nstart=50)
{
  n = nrow(X.matrix)
  HD.dist = Euclidean.Dist(X.matrix)
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  if(nstart>1)
  {
   starts = rbind(8*randomLHS(n=nstart-1,k=2*n)-4,c(cmdscale(dist(X.matrix),k=2)))
  }else
  {
    starts = matrix(c(cmdscale(dist(X.matrix),k=2)),nrow=1)
  }
  
  sols = matrix(NA,nrow=nrow(starts),ncol=2*n)
  mins = rep(NA,nrow(starts))
  
  for(k in 1:nrow(starts))
  {
    opt = optim(par=starts[k,],fn=my.fun,method="L-BFGS-B",
                lower=rep(-4,2),upper=rep(4,2))
    sols[k,] = opt$par
    mins[k] = opt$value
  }
  best.sol = sols[which.min(mins),]
  
  Z = matrix(best.sol,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Manhattan.MDS = function(X.matrix,nstart=50)
{
  n = nrow(X.matrix)
  HD.dist = as.matrix(dist(X.matrix,method="manhattan"))
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  if(nstart>1)
  {
    starts = rbind(16*randomLHS(n=nstart-1,k=2*n)-8,c(cmdscale(dist(X.matrix),k=2)),c(Z.paper))
  }else
  {
    starts = matrix(c(cmdscale(dist(X.matrix),k=2)),nrow=1)
  }
  
  sols = matrix(NA,nrow=nrow(starts),ncol=2*n)
  mins = rep(NA,nrow(starts))
  
  for(k in 1:nrow(starts))
  {
    opt = optim(par=starts[k,],fn=my.fun,method="L-BFGS-B",
                lower=rep(-8,2),upper=rep(8,2))
    sols[k,] = opt$par
    mins[k] = opt$value
  }
  best.sol = sols[which.min(mins),]
  
  Z = matrix(best.sol,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Cos.MDS = function(X.matrix,nstart=10)
{
  n = nrow(X.matrix)
  HD.dist = Cos.Dist(X.matrix)
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  if(nstart>1)
  {
    starts = rbind(8*randomLHS(n=nstart-1,k=2*n)-4,c(cmdscale(dist(X.matrix),k=2)))
  }else
  {
    starts = matrix(c(cmdscale(dist(X.matrix),k=2)),nrow=1)
  }
  sols = matrix(NA,nrow=nrow(starts),ncol=2*n)
  mins = rep(NA,nrow(starts))
  
  for(k in 1:nrow(starts))
  {
    opt = optim(par=starts[k,],fn=my.fun,method="L-BFGS-B",
                lower=rep(-4,2),upper=rep(4,2))
    sols[k,] = opt$par
    mins[k] = opt$value
    print(k)
  }
  best.sol = sols[which.min(mins),]
  
  Z = matrix(best.sol,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}

Cos.MDS.SA = function(X.matrix)
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
  
  temp = GenSA(lower = rep(-8,2*n), upper = rep(8,2*n), fn = my.fun,
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

# Inner.Product.Stress = function(x,X,Z,a)
# {
#   a = matrix(a,ncol=1)
#   x = matrix(x,ncol=1)
#   HD.dist = PCA.Dist2(X,a)
#   LD.dist = Euclidean.Dist2(Z,x)
#   stress = sum((HD.dist-LD.dist)^2)
#   return(stress)
# }

Cos.Stress = function(x,X,Z,a)
{
  a = matrix(a,ncol=1)
  x = matrix(x,ncol=1)
  HD.dist = Cos.Dist2(X,a)
  LD.dist = Euclidean.Dist2(Z,x)
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}

Manhattan.Stress = function(x,X,Z,a)
{
  a = matrix(a,nrow=1)
  x = matrix(x,nrow=1)
  HD.dist = dist2(X,a,method="manhattan")
  LD.dist = dist2(Z,x,method="euclidean")
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}

Minkowski.Stress = function(x,X,Z,a,my.power)
{
  a = matrix(a,nrow=1)
  x = matrix(x,nrow=1)
  HD.dist = dist2(X,a,method="minkowski",p=my.power)
  LD.dist = dist2(Z,x,method="euclidean")
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}

WE.Stress = function(x,X,Z,a,w)
{
  a = matrix(a,nrow=1)
  x = matrix(x,nrow=1)
  HD.dist = W.Euclidean2(w=w,X.mat=X,y.vec=a)
  LD.dist = dist2(Z,x,method="euclidean")
  stress = sum((HD.dist-LD.dist)^2)
  return(stress)
}

Sqrt.Manhattan.Stress = function(x,X,Z,a)
{
  a = matrix(a,nrow=1)
  x = matrix(x,nrow=1)
  HD.dist = sqrt(dist2(X,a,method="manhattan"))
  LD.dist = dist2(Z,x,method="euclidean")
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

Rotation = function(Z)
{
  Z = as.matrix(Z)
  Pshift=matrix(NA,nrow=3,ncol=2)
  Pshift[1,]=0
  Pshift[2,]=Z[2,]-Z[1,]
  Pshift[3,]=Z[3,]-Z[1,]
  
  theta=atan2(Pshift[2,2],Pshift[2,1])
  B=pi-theta
  
  t=matrix(c(cos(B),sin(B),-sin(B),cos(B)),ncol=2,byrow=TRUE)
  Z.star=Z%*%t
  s1=sign(Z.star[1,2]-Z.star[3,2])
  s2=sign(Z.star[2,1]-Z.star[1,1])
  r=matrix(c(s2*1,0,0,s1*1),ncol=2,byrow=TRUE)
  Z.star=Z.star%*%r
  
  return(Z.star)
}

MDS = function(X,scale=TRUE,HD.distance=c("euclidean","manhattan","sqrt.manhattan","cosine",
                               "innerproduct","minkowski"),m.power=2)
{
  if(scale==TRUE){X=scale(X)}
  
  if(HD.distance=="euclidean")
  {
    Z = Euclidean.MDS.SA(X)
  }
  
  if(HD.distance=="manhattan")
  {
    # Z= Minkowski.MDS.SA(X,my.power=1)
    Z = Manhattan.MDS.SA(X)
    # Z = Manhattan.MDS(X,nstart=1000)
  }
  
  if(HD.distance=="minkowski")
  {
    Z = Minkowski.MDS.SA(X,m.power)
  }
  
  if(HD.distance=="cosine")
  {
    Z = Cos.MDS.SA(X)
    # Z = Cos.MDS(X,nstart=1)
  }
  
  if(HD.distance=="sqrt.manhattan")
  {
    Z = Sqrt.Manhattan.MDS.SA(X)
  }
  
  if(HD.distance=="innerproduct")
  {
    V = eigen(t(X)%*%X)$vectors[,1:2]
    Z = X%*%V
  }
  
  Z = scale(Rotation(Z),scale=FALSE)
  
  return(Z)
}

MDS.G.Biplot = function(X,Z,HD.distance=c("euclidean","manhattan","sqrt.manhattan","cosine",
                                                "innerproduct","minkowski"),m.power=2,plot.seq)
{
  n = nrow(X)
  p = ncol(X)

  if(HD.distance=="euclidean")
  {
    m.power = 2
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    
    for(i in 1:p)
    {
      startx = c(0,0)
      for(j in 1:length(plot.seq))
      {
        starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
        sols = matrix(NA,nrow=nrow(starts),ncol=2)
        mins = rep(NA,nrow(starts))
        
        for(k in 1:nrow(starts))
        {
          a = rep(0,p)
          a[i] = plot.seq[j]
          opt = optim(par=starts[k,],fn=Minkowski.Stress,method="L-BFGS-B",
                      lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z,a=a,my.power=m.power)
          sols[k,] = opt$par
          mins[k] = opt$value
        }
        LD.axis.array[j,,i] = sols[which.min(mins),]
        startx = sols[which.min(mins),]
        stress.mat[i,j] = mins[which.min(mins)]
      }
    }
  }
  
  if(HD.distance=="manhattan")
  {
    m.power = 1
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    
    for(i in 1:p)
    {
      startx = c(0,0)
      for(j in 1:length(plot.seq))
      {
        starts = rbind(20*randomLHS(n=10,k=2)-10,startx)
        sols = matrix(NA,nrow=nrow(starts),ncol=2)
        mins = rep(NA,nrow(starts))
        
        for(k in 1:nrow(starts))
        {
          a = rep(0,p)
          a[i] = plot.seq[j]
          opt = optim(par=starts[k,],fn=Minkowski.Stress,method="L-BFGS-B",
                      lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z,a=a,my.power=m.power)
          sols[k,] = opt$par
          mins[k] = opt$value
        }
        LD.axis.array[j,,i] = sols[which.min(mins),]
        startx = sols[which.min(mins),]
        stress.mat[i,j] = mins[which.min(mins)]
      }
    }
  }
  
  if(HD.distance=="sqrt.manhattan")
  {
    m.power = 1
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    
    for(i in 1:p)
    {
      startx = c(0,0)
      for(j in 1:length(plot.seq))
      {
        starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
        sols = matrix(NA,nrow=nrow(starts),ncol=2)
        mins = rep(NA,nrow(starts))
        
        for(k in 1:nrow(starts))
        {
          a = rep(0,p)
          a[i] = plot.seq[j]
          opt = optim(par=starts[k,],fn=Sqrt.Manhattan.Stress,method="L-BFGS-B",
                      lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z,a=a)
          sols[k,] = opt$par
          mins[k] = opt$value
        }
        LD.axis.array[j,,i] = sols[which.min(mins),]
        startx = sols[which.min(mins),]
        stress.mat[i,j] = mins[which.min(mins)]
      }
    }
  }
  
  if(HD.distance=="cosine")
  {
    LD.axis.array = array(NA,c(1,2,p))
    stress.mat = matrix(NA,nrow=p,ncol=1)
    
    for(i in 1:p)
    {
      startx = c(0,0)
      starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
      sols = matrix(NA,nrow=nrow(starts),ncol=2)
      mins = rep(NA,nrow(starts))
      
      for(k in 1:nrow(starts))
      {
        a = rep(0,p)
        a[i] = 1
        opt = optim(par=starts[k,],fn=Cos.Stress,method="L-BFGS-B",
                    lower=rep(-8,2),upper=rep(8,2),X=X,Z=Z,a=a)
        sols[k,] = opt$par
        mins[k] = opt$value
      }
      LD.axis.array[,,i] = sols[which.min(mins),]
      stress.mat[i,] = Cos.Stress(x=sols[which.min(mins),],X=X,Z=Z,a=a)
    }
  }
  
  if(HD.distance=="innerproduct")
  {
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    V = eigen(t(X)%*%X)$vectors[,1:2]
    
    Z = X%*%V
    Pshift=matrix(NA,nrow=3,ncol=2)
    Pshift[1,]=0
    Pshift[2,]=Z[2,]-Z[1,]
    Pshift[3,]=Z[3,]-Z[1,]
    
    theta=atan2(Pshift[2,2],Pshift[2,1])
    B=pi-theta
    
    t=matrix(c(cos(B),sin(B),-sin(B),cos(B)),ncol=2,byrow=TRUE)
    Z.star=Z%*%t
    s1=sign(Z.star[1,2]-Z.star[3,2])
    s2=sign(Z.star[2,1]-Z.star[1,1])
    r=matrix(c(s2*1,0,0,s1*1),ncol=2,byrow=TRUE)
    
    for(i in 1:p)
    {
      for(j in 1:length(plot.seq))
      {
        a = rep(0,p)
        a[i] = plot.seq[j]

        LD.axis.array[j,,i] = a%*%V
        stress.mat[i,j] = sum((X%*%a-Z%*%matrix(a%*%V,ncol=1))^2)
      }
      LD.axis.array[,,i] = LD.axis.array[,,i]%*%t%*%r
      
    }
  }
  
  if(HD.distance=="minkowski")
  {
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    
    for(i in 1:p)
    {
      startx = c(0,0)
      for(j in 1:length(plot.seq))
      {
        starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
        sols = matrix(NA,nrow=nrow(starts),ncol=2)
        mins = rep(NA,nrow(starts))
        
        for(k in 1:nrow(starts))
        {
          a = rep(0,p)
          a[i] = plot.seq[j]
          opt = optim(par=starts[k,],fn=Minkowski.Stress,method="L-BFGS-B",
                      lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z,a=a,my.power=m.power)
          sols[k,] = opt$par
          mins[k] = opt$value
        }
        LD.axis.array[j,,i] = sols[which.min(mins),]
        startx = sols[which.min(mins),]
        stress.mat[i,j] = mins[which.min(mins)]
      }
    }
  }
  
  return(list(LD.axis.array=LD.axis.array,stress.mat=stress.mat))
}

MDS.G.Biplot.WE = function(X,Z,w,plot.seq)
{
  n = nrow(X)
  p = ncol(X)
  
    LD.axis.array = array(NA,c(length(plot.seq),2,p))
    stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))
    
    for(i in 1:p)
    {
      startx = c(0,0)
      for(j in 1:length(plot.seq))
      {
        starts = rbind(12*randomLHS(n=10,k=2)-6,startx)
        sols = matrix(NA,nrow=nrow(starts),ncol=2)
        mins = rep(NA,nrow(starts))
        for(k in 1:nrow(starts))
        {
          a = rep(0,p)
          a[i] = plot.seq[j]
          opt = optim(par=starts[k,],fn=WE.Stress,method="L-BFGS-B",
                      lower=rep(-10,2),upper=rep(10,2),X=X,Z=Z,a=a,w=w)
          sols[k,] = opt$par
          mins[k] = opt$value
        }
        LD.axis.array[j,,i] = sols[which.min(mins),]
        startx = sols[which.min(mins),]
        stress.mat[i,j] = mins[which.min(mins)]
      }
    }
  
  return(list(LD.axis.array=LD.axis.array,stress.mat=stress.mat))
}

# Functions Matt wrote
Manhattan.Dist <- function(X.mat){ #R
  rows <- nrow(X.mat)
  dist <- matrix(0,nrow=rows,ncol=rows)
  for (ii in 1:(rows-1)){
    for (jj in (ii+1):rows){
      dist[jj,ii] <- dist[ii,jj] <- sum(abs(X.mat[ii,]-X.mat[jj,]))
    }
  }
  return(dist)
}

Manhattan.Dist2 <- function(X.mat,x.vec){ #R
  rows <- nrow(X.mat)
  dist <- matrix(0,nrow=rows,ncol=1)
  for (ii in 1:rows){
    dist[ii,1] <- sum(abs(X.mat[ii,]-x.vec))
  }
  return(dist)
}

gower_nl_fnc <- function(filepath.and.name, 
                         raw.var.names,
                         plot.var.names,
                         high.D.dist.metric,
                         plot.seq
){
  #Read in Raw Data
  raw.data = read.csv(filepath.and.name,header=TRUE)
  
  #Get Row Names and Create Data Matrix (DM)
  observations = raw.data[,1]
  attributes.all = plot.var.names
  DM = raw.data[,colnames(raw.data) %in% raw.var.names]
  colnames(DM) <- attributes.all
  
  #Attributes of DM
  n = nrow(DM)
  p = ncol(DM)
  
  #Set Up Storage
  axis.results = array(NA, dim=c(length(plot.seq),2,p), dimnames = list(1:length(plot.seq),c("Z1","Z2"),plot.var.names))
  pt.results = array(NA,dim=c(n,2), dimnames = list(1:n,c("Z1","Z2")))
  
  #########Gower's NL Method
  #Scale DM so that all entries are mean 0 and variance 1
  DM_scaled <- apply(DM,2,function(x) (x-mean(x))/sd(x))
  
  #Convert DM_scaled to matrix for use in distance function
  DM_scaled <- as.matrix(DM_scaled)
  
  #Compute Dstar Matrix (which holds distances between rows of DM_scaled)
  if (high.D.dist.metric == "Euclidean"){
    Dstar <- Euclidean.Dist(DM_scaled) #This function took the square root
  } else{
    if (high.D.dist.metric == "Manhattan"){
      Dstar <- sqrt(Manhattan.Dist(DM_scaled)) #Needs the sqrt rt b/c of Euclidean embeddable stuff - see book
    } else {print("Incorrect Distance Name")}
  }
  
  #Perform Double Centering Technique
  D <- -0.5*Dstar^2
  J <- matrix(1,n,n)
  B <- (diag(1,n) - (1/n)*J) %*% D %*% (diag(1,n) - (1/n)*J)
  
  #Find Low-D Coordinates
  eigvals <- round(sort(eigen(B)$values,decreasing = TRUE),digits=8)
  Delta <- diag(eigvals)
  V <- eigen(B)$vectors #Columns are normalized eigenvectos
  Y_star <- V%*%Delta^.5
  Y_star <- Rotation(Y_star[,1:2])
  
  #LD Points
  pt.results = Y_star
  
  #Axis Projections
  get_trajectory <- function(axis,axis_seq,dist.met){
    axis_points_projected <- NULL
    for (i in axis_seq){
      #Create new point
      new_pt <- c(rep(0,p))
      new_pt[axis] <- i #Making one element non-zero
      
      #Compute distance between new point and original points 
      if (dist.met=="Euclidean"){
        dist_sq_new_old <- Euclidean.Dist2(DM_scaled,new_pt)^2
      } else {
        if (dist.met=="Manhattan"){
          dist_sq_new_old <- Manhattan.Dist2(DM_scaled,new_pt)
        }
      }
      
      #Compute vector of d_{ii}^2 Using Eqn (3) for i=1,2,...,n
      d2_ii <- matrix(apply(Dstar,1,function(x) sum(x^2))/n - sum(Dstar^2)/(2*n^2),ncol=1) 
      
      #Compute d=d_{ii}^2 - d_{i,n+1}^2 as defined in text before Eqn (3) but with formula corrected from other paper 
      d <- d2_ii - dist_sq_new_old
      
      #Compute y using Eqn (4) in NLBiplots Paper
      y <- (1/2)*solve(t(Y_star)%*%Y_star)%*%t(Y_star)%*%d
      axis_points_projected <- rbind(axis_points_projected,t(y))
    }
    return(axis_points_projected[,1:2]) #Return just the columns needed for plotting in 2-d
  }
  
  for (ii in 1:p){
    axis.results[,,ii] = get_trajectory(axis=ii,axis_seq=plot.seq,dist.met = high.D.dist.metric)
  }
  
  results = list(LD.pts=pt.results, axis.pts=axis.results)
}

gower_nl_fnc_we <- function(filepath.and.name, raw.var.names,plot.var.names,high.D.dist.metric,plot.seq)
{
  #Read in Raw Data
  raw.data = read.csv(filepath.and.name,header=TRUE)
  
  #Get Row Names and Create Data Matrix (DM)
  observations = raw.data[,1]
  attributes.all = plot.var.names
  DM = raw.data[,colnames(raw.data) %in% raw.var.names]
  colnames(DM) <- attributes.all
  
  #Attributes of DM
  n = nrow(DM)
  p = ncol(DM)
  
  #Set Up Storage
  axis.results = array(NA, dim=c(length(plot.seq),2,p), dimnames = list(1:length(plot.seq),c("Z1","Z2"),plot.var.names))
  pt.results = array(NA,dim=c(n,2), dimnames = list(1:n,c("Z1","Z2")))
  
  #########Gower's NL Method
  #Scale DM so that all entries are mean 0 and variance 1
  DM_scaled <- apply(DM,2,function(x) (x-mean(x))/sd(x))
  
  #Convert DM_scaled to matrix for use in distance function
  DM_scaled <- as.matrix(DM_scaled)
  
  #Compute Dstar Matrix (which holds distances between rows of DM_scaled)
  if (high.D.dist.metric == "Euclidean")
  {
    HD.dist = as.matrix(dist(DM_scaled,method="euclidean"))
    out = optim(par=rep(1,p),fn=W.Euclidean.stress,X.mat=DM_scaled,HD.dist=HD.dist,method="L-BFGS-B",
                lower=rep(0,p),upper=rep(10,p),control=list(maxit=500))
    w.opt = out$par
    Dstar = W.Euclidean(w.opt,DM_scaled)
  }
  if (high.D.dist.metric == "Manhattan")
  {
    HD.dist = as.matrix(dist(DM_scaled,method="manhattan"))
    out = optim(par=rep(1,p),fn=W.Euclidean.stress,X.mat=DM_scaled,HD.dist=HD.dist,method="L-BFGS-B",
                lower=rep(0,p),upper=rep(10,p),control=list(maxit=500))
    w.opt = out$par
    Dstar = W.Euclidean(w.opt,DM_scaled)
  }
  if(high.D.dist.metric == "Cosine")
  {
    HD.dist = Cos.Dist(DM_scaled)
    out = optim(par=rep(1,p),fn=W.Euclidean.stress,X.mat=DM_scaled,HD.dist=HD.dist,method="L-BFGS-B",
                lower=rep(0,p),upper=rep(10,p),control=list(maxit=500))
    w.opt = out$par
    Dstar = W.Euclidean(w.opt,DM_scaled)
  }
  
  #Perform Double Centering Technique
  D <- -0.5*Dstar^2
  J <- matrix(1,n,n)
  B <- (diag(1,n) - (1/n)*J) %*% D %*% (diag(1,n) - (1/n)*J)
  
  #Find Low-D Coordinates
  eigvals <- round(sort(eigen(B)$values,decreasing = TRUE),digits=8)
  Delta <- diag(eigvals)
  V <- eigen(B)$vectors #Columns are normalized eigenvectos
  Y_star <- V%*%Delta^.5
  Y_star <- Rotation(Y_star[,1:2])
  
  #LD Points
  pt.results = Y_star
  
  #Axis Projections
  get_trajectory <- function(axis,axis_seq)
  {
    axis_points_projected <- NULL
    for (i in axis_seq){
      #Create new point
      new_pt <- c(rep(0,p))
      new_pt[axis] <- i #Making one element non-zero
      
      #Compute distance between new point and original points 

      dist_sq_new_old <- W.Euclidean2(w=w.opt,X.mat=DM_scaled,y.vec=new_pt)^2

      
      #Compute vector of d_{ii}^2 Using Eqn (3) for i=1,2,...,n
      d2_ii <- matrix(apply(Dstar,1,function(x) sum(x^2))/n - sum(Dstar^2)/(2*n^2),ncol=1) 
      
      #Compute d=d_{ii}^2 - d_{i,n+1}^2 as defined in text before Eqn (3) but with formula corrected from other paper 
      d <- d2_ii - dist_sq_new_old
      
      #Compute y using Eqn (4) in NLBiplots Paper
      y <- (1/2)*solve(t(Y_star)%*%Y_star)%*%t(Y_star)%*%d
      axis_points_projected <- rbind(axis_points_projected,t(y))
    }
    return(axis_points_projected[,1:2]) #Return just the columns needed for plotting in 2-d
  }
  
  for (ii in 1:p){
    axis.results[,,ii] = get_trajectory(axis=ii,axis_seq=plot.seq)
  }
  
  results = list(LD.pts=pt.results, axis.pts=axis.results)
}

mueller_fnc <- function(filepath.and.name, 
                        raw.var.names,
                        plot.var.names,
                        high.D.obs.dist.metric,
                        mds.method
){
  #Read in Raw Data
  raw.data = read.csv(filepath.and.name,header=TRUE)
  
  #Get Row Names and Create Data Matrix (DM)
  observations = raw.data[,1]
  attributes.all = plot.var.names
  DM = raw.data[,colnames(raw.data) %in% raw.var.names]
  colnames(DM) <- attributes.all
  
  #Attributes of DM
  n = nrow(DM)
  p = ncol(DM)
  
  #Set Up Storage
  var.results = array(NA, dim=c(1,2,p), dimnames = list(1,c("Z1","Z2"),plot.var.names))
  pt.results = array(NA,dim=c(n,2), dimnames = list(1:n,c("Z1","Z2")))
  
  #########Mueller's Method
  #Scale DM so that all entries are between 0 and 1 - Section 3.1 p. 123
  DM_range <- apply(DM, 2, range)
  for (i in 1:p){
    DM[,i] <- (DM[,i]-DM_range[1,i])/(DM_range[2,i] - DM_range[1,i])
  }
  
  #Convert DM to matrix for use in distance function
  DM <- as.matrix(DM)
  
  #Compute DD Matrix (pairwise distances between observation/samples/rows of DM)
  #Compute Dstar Matrix (which holds distances between rows of DM_scaled)
  if (high.D.obs.dist.metric == "Euclidean"){
    DD <- Euclidean.Dist(DM)
  } 
  if (high.D.obs.dist.metric == "Manhattan"){
      DD <- Manhattan.Dist(DM) }
  if (high.D.obs.dist.metric == "SqrtManhattan"){
    DD <- sqrt(Manhattan.Dist(DM)) }
  if (high.D.obs.dist.metric == "Cosine"){
    DD <- Cos.Dist(DM) }
  
  #Create V Matrix (each row consists of the values observed for that variable)
  V <- as.matrix(t(DM))
  
  #Compute VV Matrix (distances between rows of V Matrix - Using Euclidean for now)
  VV <- 1 - cor(t(V))
  
  #Compute VD (distances between variables and data) which is nxm (8x15) using (1-value distance)
  VD <- as.matrix(matrix(1,nrow=ncol(DM),ncol=nrow(DM)) - V)
  
  #Compute DV
  DV <- as.matrix(matrix(1,nrow=nrow(DM),ncol=ncol(DM)) - DM)
  
  #Scale submatrices so that all four submatrices have the same mean
  DD_mean <- mean(DD)
  VD_mean <- mean(VD)
  DV_mean <- mean(DV)
  VV_mean <- mean(VV)
  M_max <- max(DD_mean,VD_mean,DV_mean,VV_mean)
  
  W_DD <- M_max/DD_mean
  W_VD <- M_max/VD_mean
  W_DV <- M_max/DV_mean
  W_VV <- M_max/VV_mean
  
  DD_theta <- DD*W_DD
  VD_theta <- VD*W_VD
  DV_theta <- DV*W_DV
  VV_theta <- VV*W_VV
  
  #Form Composite Matrix CM
  CM <- rbind(cbind(DD_theta,DV_theta),cbind(VD_theta,VV_theta))
  
  # CM_mds_results <- Rotation(Euclidean.MDS.SA(as.matrix(CM)))
  CM_mds_results <- Rotation(MDS.for.DCM(as.matrix(CM)))
  #Perform MDS on CM - Pick One
  # if (mds.method == "Euclidean")
  #   {
  #   CM_mds_results <- Rotation(Euclidean.MDS.SA(as.matrix(CM)))
  #   }
  # if (mds.method == "Manhattan")
  #   {
  #     CM_mds_results <- Rotation(Manhattan.MDS.SA(as.matrix(CM)))
  #   }
  # if (mds.method == "SqrtManhattan")
  # {
  #   CM_mds_results <- Rotation(Sqrt.Manhattan.MDS.SA(as.matrix(CM)))
  # }
  # if (mds.method == "Cosine")
  # {
  #   CM_mds_results <- Rotation(Cosine.MDS.SA(as.matrix(CM)))
  # } 
  
  #LD Points for observations
  pt.results = as.data.frame(CM_mds_results[1:n,])
  colnames(pt.results) = c("Z1","Z2")
  
  #LD Points for variables
  for (i in 1:p){
    var.results[,,i] = CM_mds_results[n+i,]
    #colnames(var.results[,,i]) = c("Z1","Z2")
  }
  
  results = list(LD.pts=pt.results, var.pts=var.results)
}


MDS.for.DCM = function(HD.dist)
{
  n = nrow(HD.dist)
  x.start = c(rnorm(2*n))
  
  my.fun = function(z.vec)
  {
    Z.mat = matrix(z.vec,nrow=n,ncol=2,byrow=TRUE)
    LD.dist = Euclidean.Dist(Z.mat)
    stress = sum((HD.dist-LD.dist)^2)
    return(stress)
  }
  
  temp = GenSA(lower = rep(-6,2*n), upper = rep(6,2*n), fn = my.fun,
               control=list(verbose=FALSE))
  Z = matrix(temp$par,nrow=n,ncol=2,byrow=TRUE)
  return(Z)
}