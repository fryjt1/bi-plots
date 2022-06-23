library(tm)
library(SnowballC)
library(topicmodels)

n = 200
raw.data = read.csv("~/Dropbox/Leman Research/Bi-plots/Tweet Data/tweets.csv")
row.sample = sample(1:nrow(raw.data),size=n,replace=FALSE)
tweets = as.data.frame(raw.data[row.sample,colnames(raw.data)%in%c("text")])
user = as.data.frame(raw.data[row.sample,colnames(raw.data)%in%c("handle")])
names(user) = "user"

corpus = Corpus(DataframeSource(tweets))
corp.copy = corpus
for(j in seq(corpus))   
{   
  corpus[[j]] <- gsub('http\\S+\\s*',"", corpus[[j]])
  corpus[[j]] <- gsub("@", "", corpus[[j]])   
  corpus[[j]] <- gsub("\\|", "", corpus[[j]])
  
  corpus[[j]] = gsub("&amp", "", corpus[[j]])
  corpus[[j]] = gsub("(RT|via)((?:\\b\\W*@\\w+)+)", "", corpus[[j]])
  corpus[[j]] = gsub("@\\w+", "", corpus[[j]])
  corpus[[j]] = gsub("[[:punct:]]", "", corpus[[j]])
  corpus[[j]] = gsub("[[:digit:]]", "", corpus[[j]])
  corpus[[j]] = gsub("http\\w+", "", corpus[[j]])
  corpus[[j]] = gsub("[ \t]{2,}", "", corpus[[j]])
  corpus[[j]] = gsub("^\\s+|\\s+$", "", corpus[[j]]) 
}   

corpus = tm_map(corpus, PlainTextDocument)
corpus = tm_map(corpus, content_transformer(tolower))
corpus = tm_map(corpus, removePunctuation)
corpus = tm_map(corpus, removeNumbers)
corpus = tm_map(corpus, removeWords, stopwords("english"))
corpus = tm_map(corpus, stemDocument, language = "english") 
corpus = tm_map(corpus, stripWhitespace)
dtm = DocumentTermMatrix(corpus)
# ,control=list(weighing=weightTfIdf)

n.topics = 5
fit = LDA(x=dtm,k=n.topics,method="VEM")
terms(fit,k=3)
temp = posterior(fit)
X = temp$topics

attributes.all = c(paste("Topic",paste(1:n.topics)))
observations = c(1:n)
X = scale(X)
p = ncol(X)

Z.mds.Cos = Cos.MDS(X,nstart=1)

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

K = n.topics
LD.axis = LD.axis.all[which(rank(apply(stress.mat,1,mean))<(K+1)),]
attributes = attributes.all[which(rank(apply(stress.mat,1,mean))<(K+1))]

Z.mds.Cos.df = as.data.frame(Z.mds.Cos); colnames(Z.mds.Cos.df) = c("Z1","Z2")
Z.mds.Cos.df$user = unlist(user)
LD.axis = as.data.frame(LD.axis)
LD.axis$attributes = attributes

mdsplot <- ggplot() + labs(title="MDS (Cosine) Projection",x="",y="") +
  geom_point(data=Z.mds.Cos.df,aes(x=Z1,y=Z2,color=user)) +
  scale_colour_manual(values = c("HillaryClinton"="blue","realDonaldTrump"="red")) +
  geom_text_repel(data=Z.mds.Cos.df,aes(x=Z1,y=Z2,label=observations,color=user),cex=2.5) + 
  geom_point(data=LD.axis,aes(x=V1,y=V2),color="black") +
  geom_text_repel(data=LD.axis,aes(x=V1,y=V2,label=attributes,fontface="bold"),color="black",cex=5) + 
  scale_x_continuous(expand = c(.2, .2)) +
  theme(legend.position="none")

print(mdsplot)

terms(fit,k=3)
