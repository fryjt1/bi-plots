library(ggplot2)
library(ggrepel)

raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/College_Data.csv")
data = raw.data[,colnames(raw.data)%in%c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                                 "ACT_75","Percent_Admit",
                                 "Graduation_Rate","Percent_Male","Avg_Net_Cost")]
observations = raw.data[,1]
attributes = c("Student-to-Faculty","Enrollment","Percent Graduate Students",
               "ACT Score","Percent Admitted","Graduation Rate","Percent Male",
               "Average Cost")
X = scale(data)
n = nrow(X)
p = ncol(X)
V = -eigen(t(X)%*%X)$vectors[,1:2]
Z = Rotation(X%*%V)
b = sqrt(n)

Z.pca.df = as.data.frame(Z/b)
V.pca.df = as.data.frame(V*b)
colnames(V.pca.df) = c("V1","V2")
colnames(Z.pca.df) = c("Z1","Z2")
pcaplot <- ggplot() + labs(title="PCA Projection",x="",y="") +
      geom_point(data=Z.pca.df,aes(x=Z1,y=Z2),color="black") + 
      geom_segment(aes(x=0,y=0,xend=V1,yend=V2,linetype="dotted"),linetype=2,data=V.pca.df,color="black") +
      geom_text_repel(data=Z.pca.df,aes(x=Z1,y=Z2,label=observations),color="blue",cex=3) + 
      geom_text(aes(x=V1,y=V2,label=attributes),data=V.pca.df,color="black",cex=3) +
      scale_x_continuous(expand = c(.1, .1))
print(pcaplot)

