library(ggplot2)
library(ggrepel)

raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/College_Data.csv")
data = raw.data[,colnames(data)%in%c("Fac_to_Student","Tot_Enrollment","Percent_Grad_Student",
                                 "Tuition_In","Tuition_Out","Avg_Scholarship","Percent_Admit",
                                 "Graduation_Rate","Percent_Male")]
observations = raw.data[,1]
attributes = c("Faculty-to-Student","Enrollment","Percent Graduate Students","In-State Tuition",
               "Out-of-State Tuition","Average Scholarship","Percent Admitted","Graduation Rate","Percent Male")
X = scale(data)
n = nrow(X)
p = ncol(X)
V = eigen(t(X)%*%X)$vectors[,1:2]
Z = X%*%V
b = sqrt(n)

library(ggplot2)

Z = as.data.frame(Z/b)
V = as.data.frame(V*b)
colnames(V) = c("V1","V2")
colnames(Z) = c("Z1","Z2")
pcaplot <- ggplot() + labs(title="PCA Projection",x="",y="") +
      geom_point(data=Z,aes(x=Z1,y=Z2),color="blue") + 
      geom_segment(aes(x=0,y=0,xend=V1,yend=V2,linetype="dotted"),linetype=2,data=V,color="black") +
      geom_text_repel(aes(x=Z1,y=Z2,label=observations),data=Z,color="blue",cex=3) + 
      geom_text(aes(x=V1,y=V2,label=attributes),data=V,color="black",cex=3) +
      scale_x_continuous(expand = c(.1, .1))
print(pcaplot)

