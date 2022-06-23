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
V = -eigen(t(X)%*%X)$vectors[,1:2]

Z.mds.IP = X%*%V

Z.mds.IP = scale(Rotation(Z.mds.IP),scale=FALSE)
plot.seq = seq(-10,10,by=1)
LD.axis.array.all = array(NA,c(length(plot.seq),2,p))
stress.mat = matrix(NA,nrow=p,ncol=length(plot.seq))

for(i in 1:p)
{
  for(j in 1:length(plot.seq))
  {

    a = rep(0,p)
    a[i] = plot.seq[j]
    
    LD.axis.array.all[j,,i] = a%*%V
    stress.mat[i,j] = PCA.Stress(a%*%V,X,Z.mds.IP,a)
  }
}

K = 8
LD.axis.array = LD.axis.array.all[,,which(rank(apply(stress.mat,1,mean))<(K+1))]
attributes = attributes.all[which(rank(apply(stress.mat,1,mean))<(K+1))]
LD.axis = cbind(c(apply(LD.axis.array,3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = rep(attributes,length(plot.seq))[order(rep(1:K,length(plot.seq)))]

Z.mds.IP.df = as.data.frame(Z.mds.IP); colnames(Z.mds.IP.df) = c("Z1","Z2")
LD.axis = as.data.frame(LD.axis)
LD.axis$attribute = attribute.vec; colnames(LD.axis) = c("b1","b2","attribute")
tick.seq = seq(min(plot.seq),max(plot.seq),by=1); tick.seq = tick.seq[tick.seq!=0]
LD.ticks = cbind(c(apply(LD.axis.array[which(plot.seq%in%tick.seq),,],3,function(x)x[,1])),
                           c(apply(LD.axis.array[which(plot.seq%in%tick.seq),,],3,function(x)x[,2])))
LD.ticks = as.data.frame(LD.ticks)
LD.ticks$attribute = rep(attributes,length(tick.seq))[order(rep(1:K,length(tick.seq)))]
LD.ticks$mark = rep(tick.seq,K)
endpoint.ind = seq(length(plot.seq),K*length(plot.seq),by=length(plot.seq))
endpoint.df = LD.axis[endpoint.ind,]

mdsplot <- ggplot() + labs(title="MDS (Inner-Product) Projection",x="",y="") +
  geom_point(data=Z.mds.IP.df,aes(x=Z1,y=Z2),color="black") +
  geom_text_repel(data=Z.mds.IP.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_path(data=LD.axis,aes(x=b1,y=b2,color=factor(attribute)),linetype="dotted") +
  geom_text_repel(data=endpoint.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3) +
  geom_text(data=LD.ticks,aes(x=V1,y=V2,label=mark,color=factor(attribute)),cex=2.5) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(legend.position="none")


print(mdsplot)
