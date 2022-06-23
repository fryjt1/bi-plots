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
S = t(X)%*%X/n

b = 5
SVD = svd(X)
U = SVD$u[,1:2]; V = SVD$v[,1:2]; L = diag(SVD$d[1:2])
Z.mds = U%*%L
LD.axis = V*b

Z = as.matrix(Z.mds)
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

Z.mds = Z.mds%*%t%*%r
LD.axis = LD.axis%*%t%*%r
attribute.vec = attributes.all

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec
colnames(LD.axis.df) = c("b1","b2","attribute")

mdsplot <- ggplot() + labs(title=NULL,x="",y="") +
  geom_segment(data=LD.axis.df,aes(x=0,y=0,xend=b1,yend=b2,color=factor(attribute)),
               arrow = arrow(length = unit(0.3,"cm"))) +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2),color="black",cex=2.5) +
  geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=3) + 
  geom_text_repel(data=LD.axis.df,aes(x=b1,y=b2,color=factor(attribute),label=attribute),cex=3.5) +
  scale_x_continuous(limits = c(-4, 4)) + 
  scale_y_continuous(limits = c(-3, 3)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
        axis.ticks=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), 
       axis.ticks=element_blank()) +
  theme(legend.position="none") 

path = "~/Dropbox/Leman Research/General Dynamics/Bi-Plots/JSM Materials/Short Presentation/"
file.name = "PCA_arrow.png"
# png(file = paste(path,file.name,sep=""), width = 1024, height = 768, units = "px")
print(mdsplot)
# dev.off()
