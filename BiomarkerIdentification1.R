
a=sample(1:1000,10)
b=sample(1:1000,10)
a
a*2
b
a+b

a[a>500]=b[a>500]
a



data=read.csv("proteomics.csv")
data=sapply(data,as.numeric)
pdf(file="heatmap.pdf",width=5,height=5)
heatmap(data,Rowv=NA,Colv=NA, scale="column")
dev.off()

pdf(file="genehist.pdf",width=5,height=5)
hist(data[,3],breaks=20,xlab=colnames(data)[3],main="")
dev.off()


pdf(file="genebox.pdf",width=3,height=5)
boxplot(data[,3],xlab=colnames(data)[3],main="",ylab="Expression level")
dev.off()



d0_means=apply(data[data[,1]==0,-1],2,mean)
d1_means=apply(data[data[,1]==1,-1],2,mean)
LF=log(d0_means/d1_means)
pdf(file="LF.pdf",width=5,height=5)
hist(LF,breaks=20,col="green")
dev.off()


pdf(file="pmf.pdf",width=5,height=5)
hist(data[,1],xlab="Disease status", main="")
dev.off()

x=seq(-4,4,by=0.1)
y=dnorm(x)
z=pnorm(x)
pdf(file="pdf-cdf.pdf",width=6,height=3)
par(mfrow=c(1,2))
plot(x,y,type='l',col='red',ylab="PDF")
plot(x,z,type='l',col='green',ylab="CDF")
dev.off()


data0=data[data[,1]==0,]
data1=data[data[,1]==1,]
pvals=c()
for (i in 2:ncol(data0)){
	pvals=c(pvals,t.test(data0[,i],data1[,i])$p.value)
}

pvals_adj=pvals*(ncol(data0)-1)
for (i in 1:length(pvals_adj)){
	pvals_adj[i]=min(1,pvals_adj[i])
}

pdf(file="pvals.pdf",width=6,height=3)
par(mfrow=c(1,2))
hist(log(pvals,10), breaks=20,col='red',xlab='log10(pval)',main="",xlim=c(-5,0))
hist(log(pvals_adj,10), breaks=20,col='red',xlab='log10(pval adj)',main="",xlim=c(-5,0))
dev.off()

pdf(file="volcano.pdf",width=10,height=6)
par(mfrow=c(1,2))
plot(LF,-log(pvals,10),col='red',ylab="-log10(pval)",xlab="LF",main="",pch=20)
plot(LF,-log(pvals_adj,10),col='red',ylab="-log10(pval adj)",xlab="LF",main="",pch=20)
dev.off()


library(GEOquery)
data<-getGEO("GSE144459", GSEMatrix = TRUE)

head(exprs(data[[1]]))

head(pData(data[[1]]))

head(fData(data[[1]]))


pData(data[[1]])[['characteristics_ch1']]
pData(data[[1]])[['characteristics_ch1.1']]



library(limma)

ctl=pData(data[[1]])[['characteristics_ch1.1']]=="tissue: blood" & pData(data[[1]])[['characteristics_ch1']]=="strain: B6129SF2 WT"

ad=pData(data[[1]])[['characteristics_ch1.1']]=="tissue: blood" & pData(data[[1]])[['characteristics_ch1']]=="strain: 3xTg-AD"

ctl=ctl*1
ad=ad*1

design=cbind(CTL=ctl,AD=ad)
fit=lmFit(data[[1]],design)

cont.matrix <- makeContrasts(AD-CTL, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

de=topTable(fit2,adjust='BH',number=100)



data<-getGEO("GSE183356", GSEMatrix = TRUE) 
p2=startsWith(pData(data[[1]])$title,"Patient 2 Lung")
p1=startsWith(pData(data[[1]])$title,"Patient 1 Lung")
design=cbind(p1=p1,p2=p2)
cm=makeContrasts(p1-p2, levels=design)
fit=lmFit(data[[1]],design)
fit2=contrasts.fit(fit, cm)
fit2=eBayes(fit2)
degenes=topTable(fit2,adjust='BH')
decideTests(de)


datap1=exprs(data[[1]])[,design$p1]




rnaseq=read.table("GSE189990.txt",header=TRUE)
genes=rnaseq[,1]
rnaseq=log(rnaseq[,-1],2)
edata=ExpressionSet(as.matrix(rnaseq))
featureNames(edata)=genes

data<-getGEO("GSE189990", GSEMatrix = TRUE)
(pData(data[[1]])[["mortality:ch1"]]=="ok")*1->ok
(pData(data[[1]])[["mortality:ch1"]]=="fatal")*1->ko
design=cbind(fatal=ko,ok=ok)
cm=make.contrasts(fatal-ok,levels=design)

fit1=lmFit(edata,design)
fit2=contrasts.fit(fit1,cm)
fit2=eBayes(fit2)

deGenes=topTable(fit2,p.value=0.05,adjust.method="none",number=1000)

plot(deGenes$logFC,-log(deGenes$P.Value,10),xlim=c(-100,100))


