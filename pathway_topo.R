library(KEGGgraph)
all_infile=list.files('D:/pathway_kgml')
library(RBGL)
rel<-c("activation","compound","binding/association","expression","inhibition",
       "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
       "inhibition_dephosphorylation","dissociation","dephosphorylation",
       "activation_dephosphorylation","state change","activation_indirect effect",
       "inhibition_ubiquination","ubiquination", "expression_indirect effect",
       "inhibition_indirect effect","repression","dissociation_phosphorylation",
       "indirect effect_phosphorylation","activation_binding/association",
       "indirect effect","activation_compound","activation_ubiquination")
beta<-c(1,1,1,1,-1,2,1,-2,-2,1,1,2,1,1,-1,1,1,-1,-1,1,1,1,1,1,1)
for (file_name in all_infile){
mapkKGML <- system.file(paste("extdata/",file_name,sep=''),package="KEGGgraph")
mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
mapkpathway <- parseKGML(mapkKGML)
mapkGedgedata <- getKEGGedgeData(mapkG)
g<-mapkG
numer_for_gene<-list()
all_gene<-list()
for (i in 1:length(mapkGedgedata)){
  n<-grep('name',mapkGedgedata[i],value=TRUE)
  if (length(n)>0){
  d<-strsplit(n,split=',')
  entry<-d[[1]][2]
  entryid<-substring(entry,14,nchar(entry)-1)
  name<-d[[1]][6]
  sign_gene<-substring(name,10,nchar(name)-1)
  index=which(rel==sign_gene)
  numer_for_gene<-c(numer_for_gene,list(beta[index]))
  all_gene[[entryid]]<-beta[index]
  }
}
{ 
  
  nv <- length(nodes(g))
  em <- edgeMatrix(g)
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))
  
  if (any(eW <= 0, na.rm = TRUE)) 
    stop("'brandes.betweenness.centrality' requires that all edge weights are positive")
  ans <- .Call("BGL_brandes_betweenness_centrality", as.integer(nv), 
               as.integer(ne), as.integer(em - 1), as.double(eW), PACKAGE = "RBGL")
  colnames(ans[[1]]) <- nodes(g)
  colnames(ans[[3]]) <- nodes(g)
  rownames(ans[[2]]) <- c("centrality")
  rownames(ans[[5]]) <- c("from", "to")
  ans[[5]] <- apply(ans[[5]], 2, function(x, y) y[x + 1], nodes(g))
  list(betweenness.centrality.vertices = ans[[1]], edges = ans[[5]], 
       betweenness.centrality.edges = ans[[2]], relative.betweenness.centrality.vertices = ans[[3]], 
       dominance = ans[[4]])
}
for (i in 1:length(all_gene)){
  gene_index<-which(nodes(g)==names(all_gene[i]))
  ans[[3]][gene_index]<-as.numeric(ans[[3]][gene_index])*as.numeric(all_gene[i])
  
}
for (i in 1:length(ans[[3]])){
     Gene_id<-translateKEGGID2GeneID(nodes(g))[i]
     write.table((paste(Gene_id,ans[[3]][i],sep=' ')),paste('D:/Topo_kgml/',mapkpathway@pathwayInfo@title,sep=''),append=TRUE)
}
}

