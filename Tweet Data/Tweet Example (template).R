library(ggplot2)
library(ggrepel)
library(lhs)

my.distance = "cosine"
my.power = 1

data = tfidf
observations = row.names(data)
attributes.all = colnames(data)

X = scale(data,scale=TRUE)
n = nrow(X)
p = ncol(X)

Z.mds = MDS(X,HD.distance=my.distance,m.power=my.power)
my.plot.seq = seq(-5,5,by=1)
fit = MDS.G.Biplot(X=X,Z=Z.mds,HD.distance=my.distance,m.power=my.power,plot.seq=my.plot.seq)
LD.axis.array = fit$LD.axis.array
stress = fit$stress

K = 50 # Number of attributes to plot
# my.ranks = rev(rank(stress)); ind = which(my.ranks<(K+1))
my.ranks = rank(colMeans(as.matrix(tfidf))); ind = which(my.ranks<(K+1))
# my.ranks = rev(rank(idf)); ind = which(my.ranks<(K+1))
attributes = attributes.all[ind]
LD.axis = t(LD.axis.array[,,ind])
# LD.axis = cbind(c(apply(LD.axis.array[,,ind],3,function(x)x[,1])),c(apply(LD.axis.array,3,function(x)x[,2])))
attribute.vec = attributes

Z.mds.df = as.data.frame(Z.mds); colnames(Z.mds.df) = c("Z1","Z2")
Z.mds.df$user = user
LD.axis.df = as.data.frame(LD.axis)
LD.axis.df$attribute = attribute.vec; colnames(LD.axis.df) = c("b1","b2","attribute")

mdsplot <- ggplot() + labs(title="MDS Projection",x="",y="") +
  geom_point(data=Z.mds.df,aes(x=Z1,y=Z2,color=factor(user))) +
  scale_shape_manual(values=c("HillaryClinton","realDonaldTrump"))+
  scale_colour_manual(values=c("blue", "red")) +
  geom_point(data=LD.axis.df,aes(x=b1,y=b2),color="black") +
  # geom_text_repel(data=Z.mds.df,aes(x=Z1,y=Z2,label=observations),color="black",cex=2.5) + 
  geom_text_repel(data=LD.axis.df,aes(x=b1,y=b2,label=attribute),color="black",cex=3) +
  scale_x_continuous(expand = c(.2, .2)) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), 
        axis.ticks=element_blank()) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), 
        axis.ticks=element_blank()) +
  theme(legend.position="none") 
print(mdsplot)
