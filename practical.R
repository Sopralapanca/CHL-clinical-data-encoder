
##DE
library(GEOquery)
data<-getGEO("GSE15471", GSEMatrix = TRUE)

# explore the dataset
head(exprs(data[[1]]))

head(pData(data[[1]]))

head(fData(data[[1]]))



pData(data[[1]])[['sample:ch1']]


hist(exprs(data[[1]]))
hist(log(exprs(data[[1]])))


expression=log(exprs(data[[1]]))
rownames(expression)=fData(data[[1]])[["Gene Symbol"]]
head(expression)


# differential expression
library(limma)

ctl=pData(data[[1]])[['sample:ch1']]=="normal"

pc=pData(data[[1]])[['sample:ch1']]=="tumor"

ctl=ctl*1
pc=pc*1

design=cbind(CTL=ctl,PC=pc)
fit=lmFit(expression,design)

cont.matrix <- makeContrasts(PC-CTL, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

de=topTable(fit2,adjust='BH',number=nrow(fData(data[[1]])),p.value=0.01, lfc=0.3)

genes=unique(de[['ID']])
genes=genes[genes!=""]

write.table(genes,"degenes.txt",row.names=FALSE)

# perform enrichment analysis using DAVID, link on the slide on the genes

##random forest
deData=expression[genes,]
labels=pc


#create co-expression network
library(igraph)

#transpose the matrix and comput correlation
r <- cor(t(deData))
hist(r)
# 0.8 is the correlation threshold
adj1 <- (r>0.8)*1
adj2 <- (r< -0.8)*1
adj <- adj1+adj2

graph <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=NULL, diag=FALSE)

plot(graph, layout=layout_with_fr, vertex.size=1,vertex.label=NA)

#centrality of nodes - select the most important genes
# the more a path pass through a node, the more important it is
bet=estimate_betweenness(graph,directed=FALSE,cutoff=-1)
bet=sort(bet,decreasing=TRUE)

# another measure of centrality, another way to compute the important genes
cl=closeness(graph)
cl=sort(cl,decreasing=TRUE) # the genes are different from bet

# another measure of centrality, the easiest way to compute the important genes
dg=degree(graph)
dg=sort(dg,decreasing=TRUE) #again different genes from the other one

# we can  names(head(dg, n=10)) and search for them on david

# select top five genes from closeness metric and compute pubmed analysis
centralGenes=names(cl)[c(1:5)]

#pubmed - download infromation about the central genes
library(easyPubMed)

ids<- get_pubmed_ids(centralGenes[1]) # get the ids of the articles of where the genes appears, return the first 20 ids

pubmeddata=fetch_pubmed_data(ids) # download the abstracts and titles of the articles in xml format

titles= custom_grep(pubmeddata, "ArticleTitle", "char") # extract the titles

abstracts= custom_grep(pubmeddata, "Abstract", "char") # extract the abstracts

#named entity recognition using pubmed http://bern2.korea.ac.kr/pubmed/ and restful api
library(httr)

httr::set_config(config(ssl_verifypeer = 0L))

allAnnotations=list()
for (i in 1:length(ids$IdList)){
    # take one id at a time and make a request to the server
    resp=GET(paste("http://bern2.korea.ac.kr/pubmed/",ids$IdList[[i]],sep=''))
    # parse the response as json, handle exceptions with try statement
    annotations=try(content(resp,as="parsed",type="application/json"))
    # take all the annotation and build a graph where the links are the co-occurence of the entities
    if(length(annotations[[1]])>1){
      allAnnotations[[i]]=annotations[[1]]$annotations
    }
}

# create a list of entities for each article and build a graph
edges=c()
E=0
#this is not efficient!!!!!!!
for (articleAnnotation in allAnnotations){
  if (!is.null(articleAnnotation)){
  entities=c()
  for (entity in articleAnnotation){
    e2=entity$id[[1]]
    for(e1 in entities){
      if(e1!=e2){
        newEdge=TRUE
        if(!is.null(edges)){
          for (i in 1:nrow(edges)){
            if((edges[i,][1]==e1 && edges[i,][2]==e2 )||(edges[i,][1]==e2 && edges[i,][2]==e1 ) )
              newEdge=FALSE
          }}
        if(newEdge){
          edges=rbind(edges,c(e2,e1))
          E=E+1
        }
      }
    }
    entities=c(entities, e2)
  }
  }
}


graph=graph_from_edgelist(edges,directed=FALSE)
plot(graph, layout=layout_with_fr,vertex.size=1, vertex.label=NA)



library(MeSH.db)
select(MeSH.db, substring(bet[4], 6), columns=c("MESHTERM"), keytype="MESHID")[["MESHTERM"]]