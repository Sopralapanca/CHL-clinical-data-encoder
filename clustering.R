#study proteomics data
prot=read.csv("proteomics.csv")
labels=prot[,1]
prot=prot[,-1]
prot=log(prot)

#visualise using PCA and t-sne
##visualiose patients with PCA
pc=prcomp(prot)
summary(pc)
protPC=predict(pc,prot)
plot(protPC[,1],protPC[,2],cex=1, pch=labels+10, col=labels+1)

#tsne visualisation of patients
library(Rtsne)
tsne=Rtsne(prot,perplexity=10)
plot(tsne$Y,pch=labels+15, col=labels+1 )

#visualise proteins with tsne
tsne=Rtsne(t(prot),perplexity=30)
plot(tsne$Y,pch='.')

#cluster patients using k-means
#TODO homework

#cluster proteins using k-means
km=kmeans(t(prot),centers=2) 
plot(tsne$Y,pch='.',col=km$cluster)

#find best number of clusters
protDist=dist(t(prot))
ks=c(2:100)
wss=c()
silh=c()
for (k in ks){
  km=kmeans(t(prot),centers=k)
  wss=c(wss,km$tot.withinss)
  silh=c(silh,mean(silhouette(km$cluster,protDist)[,3]))
}

plot(ks,wss)

##use factoextra
library(factoextra)
f=fviz_nbclust(t(prot),kmeans,"silhouette",k.max=10)
f

##use hierarchical clustering
patientDist=dist(prot)
hc=hclust(patientDist,method="average")
plot(hc)

clst=cutree(hc,k=10)
tsne=Rtsne(prot,perplexity=10)
plot(tsne$Y,pch=labels+15, col=clst)

#study gene expression dataset GSE183356
library(GEOquery)
covid=getGEO("GSE183356")

log(exprs(covid[[1]]))->exprMatrix
as.numeric(substring(patient,9,9))->patient
organ=pData(covid[[1]])[,"source_name_ch1"]


tsne=Rtsne(t(exprMatrix))
plot(tsne$Y, col=patient,pch=5)



