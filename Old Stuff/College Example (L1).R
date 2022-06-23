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

Z.mds.L1 = Manhattan.MDS.SA(X)

stress.j = rep(NA,p)
for(j in 1:p)
{
  X.temp = X[,-j] 
  stress.j[j] = mean( (as.matrix(dist(X.temp,method="manhattan")) - 
                       as.matrix(dist(Z.mds.L1,method="manhattan")) )^2 )
}

Z.mds.L1 = scale(Rotation(Z.mds.L1),scale=FALSE)
plot.seq = seq(-4,4,by=.1)
LD.axis.array.all = array(NA,c(length(plot.seq),2,p))
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
      opt = optim(par=starts[k,],fn=Manhattan.Stress,method="L-BFGS-B",
                  lower=rep(-8,2),upper=rep(8,2),X=X,Z=Z.mds.L1,a=a)
      sols[k,] = opt$par
      mins[k] = opt$value
    }
    LD.axis.array.all[j,,i] = sols[which.min(mins),]
    startx = sols[which.min(mins),]
    stress.mat[i,j] = mins[which.min(mins)]
  }
  print(i)
}

K = 4
LD.axis.array = LD.axis.array.all[,,which(rank(apply(stress.mat,1,mean))<(K+1))]
attributes = attributes.all[which(rank(apply(stress.mat,1,mean))<(K+1))]
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = rep(attributes,length(plot.seq))[order(rep(1:K,length(plot.seq)))]

Z.mds.L1.df = as.data.frame(Z.mds.L1); colnames(Z.mds.L1.df) = c("Z1","Z2")
LD.axis = as.data.frame(LD.axis)
LD.axis$attribute = attribute.vec; colnames(LD.axis) = c("b1","b2","attribute")

endpoint.ind = seq(length(plot.seq),K*length(plot.seq),by=length(plot.seq))
endpoint.df = LD.axis[endpoint.ind,]


mdsplot <- ggplot() + labs(title="MDS (Manhattan) Projection",x="",y="") +
  geom_point(data=Z.mds.L1.df,aes(x=Z1,y=Z2),color="black") +
  geom_text_repel(data=Z.mds.L1.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_path(data=LD.axis,aes(x=b1,y=b2,color=factor(attribute)),linetype="dashed") +
  geom_text_repel(data=endpoint.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(legend.position="none")


 print(mdsplot)
