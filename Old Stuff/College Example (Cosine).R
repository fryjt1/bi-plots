library(ggplot2)
library(ggrepel)
library(lhs)

raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/College_Data.csv")
data = raw.data[,colnames(raw.data)%in%c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                                         "ACT_75","Percent_Admit",
                                         "Graduation_Rate","Percent_Male","Avg_Net_Cost")]
observations = raw.data[,1]
attributes.all = c("Student-to-Faculty","Enrollment","Percent Graduate Students",
                   "ACT Score","Percent Admitted","Graduation Rate","Percent Male",
                   "Average Cost")
X = scale(data)
n = nrow(X)
p = ncol(X)

Z.mds.Cos = Cos.MDS.SA(X)

stress.j = rep(NA,p)
for(j in 1:p)
{
  X.temp = X[,-j] 
  stress.j[j] = mean( (as.matrix(dist(X.temp,method="manhattan")) - 
                         as.matrix(dist(Z.mds.Cos,method="manhattan")) )^2 )
}

Z.mds.Cos = scale(Rotation(Z.mds.Cos),scale=FALSE)
LD.axis.all = matrix(NA,nrow=p,ncol=2)
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
                lower=rep(-8,2),upper=rep(8,2),X=X,Z=Z.mds.Cos,a=a)
    sols[k,] = opt$par
    mins[k] = opt$value
  }
  LD.axis.all[i,] = sols[which.min(mins),]
  stress.mat[i,1] = mins[which.min(mins)]

  print(i)
}

K = 8
LD.axis = LD.axis.all[which(rank(apply(stress.mat,1,mean))<(K+1)),]
attributes = attributes.all[which(rank(apply(stress.mat,1,mean))<(K+1))]

Z.mds.Cos.df = as.data.frame(Z.mds.Cos); colnames(Z.mds.Cos.df) = c("Z1","Z2")
LD.axis = as.data.frame(LD.axis)

mdsplot <- ggplot() + labs(title="MDS (Cosine) Projection",x="",y="") +
  geom_point(data=Z.mds.Cos.df,aes(x=Z1,y=Z2),color="black") +
  geom_text_repel(data=Z.mds.Cos.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_point(data=LD.axis,aes(x=V1,y=V2),color="blue") +
  geom_text_repel(data=LD.axis,aes(x=V1,y=V2,label=attributes),color="blue",cex=3) + 
  scale_x_continuous(expand = c(.2, .2)) +
  theme(legend.position="none")

print(mdsplot)

