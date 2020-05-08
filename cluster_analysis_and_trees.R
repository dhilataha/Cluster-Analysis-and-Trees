#Cluster Analysis and Trees
# Ex 3
#To compute the Euclidian distance between two vectors
a <- c(1,1); b <- c(4,5) 
sqrt(sum((a-b)^2))

#ex 4
# Distances between Cyclin gene expressions
# BiocManager::install("multtest")
library(multtest); data(golub)
index <- grep("Cyclin",golub.gnames[,2])
golub.gnames[index,2]
dist.cyclin <- dist(golub[index,],method="euclidian")
diam <- as.matrix(dist.cyclin)
rownames(diam) <- colnames(diam) <- golub.gnames[index,3]
diam[1:5,1:5]

#Example 5
# Finding the ten closest genes to a given one
# BiocManager::install("genefilter")
# BiocManager::install("Rcpp")
#library("genefilter"); library("ALL"); data(ALL)
closeto1389_at <- genefinder(ALL, "1389_at", 10, method = "euc")
closeto1389_at[[1]]$indices
round(closeto1389_at[[1]]$dists,1)
featureNames(ALL)[closeto1389_at[[1]]$indices]

#7.2 Cluster Analysis
#Single Linkage
# computes the the distances between the genes and performs a single linkage cluster analysis
#Ex 1
names <- list(c("g1","g2","g3","g4","g5"),c("p1","p2"))
sl.clus.dat <- matrix(c(1,1,1,1.1,3,2,3,2.3,5,5),ncol = 2,
          byrow = TRUE,dimnames = names) 
# pict 7.1
plot(sl.clus.dat,type="n", xlim=c(0,6), ylim=c(0,6))
text(sl.clus.dat,labels=row.names(sl.clus.dat)) 
print(dist(sl.clus.dat,method="euclidian"),digits=3)
# pict 7.2 cluster dendogram
sl.out<-hclust(dist(sl.clus.dat,method="euclidian"),method="single")
plot(sl.out)

#Example 2
sl.out<-hclust(dist(rnorm(20,0,1),method="euclidian"),method="single") 
#pict 7.3
plot(sl.out)

x <- c(rnorm(10,0,0.1),rnorm(10,3,0.5),rnorm(10,10,1.0)) 
#pict 7.4
plot(hclust(dist(x,method="euclidian"),method="single")) 
#pict 7.5
plot(sl.out)

#Example 3
# Recall that the ﬁrst twenty seven patients belong to ALL and the remaining eleven to AML and that we found earlier that the expression values of the genes 
data(golub, package="multtest") 
clusdata <- data.frame(golub[1042,],golub[2124,])
colnames(clusdata)<-c("CCND3 Cyclin D3","Zyxin") 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
#pict 7.5
plot(clusdata, pch=as.numeric(gol.fac)) 
legend("topright",legend=c("ALL","AML"),pch=1:2) 
#pict 7.6
# tree from single linkage cluster analysis
plot(hclust(dist(clusdata,method="euclidian"),method="single"))
                       
#K-means cluster analysis
#ex 1
# Relating a data generation process to k-means cluster analysis
data <- rbind(matrix(rnorm(100,0,0.5), ncol = 2), 
    + matrix(rnorm(100,2,0.5), ncol = 2)) 
cl <- kmeans(data, 2)
# Specify the color of each data 
# cluster membership
#pict 7.7
plot(data, col = cl$cluster) 
points(cl$centers, col = 1:2, pch = 8, cex=2)

# the sample means of potential clusters or the hypothesized population means
initial <- matrix(c(0,0,2,2), nrow = 2, ncol=2, byrow=TRUE)
cl<- kmeans(data, initial, nstart = 10) 
#Calculate quantile
n <- 100; nboot<-1000
boot.cl <- matrix(0,nrow=nboot,ncol = 4) 
for (i in 1:nboot){ 
  dat.star <- data[sample(1:n,replace=TRUE),]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
}
quantile(boot.cl[,1],c(0.025,0.975)) 
quantile(boot.cl[,2],c(0.025,0.975)) 
quantile(boot.cl[,3],c(0.025,0.975)) 
quantile(boot.cl[,4],c(0.025,0.975))

#Ex 2
data <- data.frame(golub[1042,],golub[2124,]) 
colnames(data)<-c("CCND3 Cyclin D3","Zyxin")
cl <- kmeans(data, 2,nstart = 10) 
#K-means clustering with 2 clusters of sizes 11, 27
cl

mean(data.frame(boot.cl))
quantile(boot.cl[,1],c(0.025,0.975)) 
quantile(boot.cl[,2],c(0.025,0.975)) 
quantile(boot.cl[,3],c(0.025,0.975)) 
quantile(boot.cl[,4],c(0.025,0.975)) 

#Correlation Coefficient
#ex 1
#library(TeachingDemos)
# run.cor.examp(1000)
#ex 2
# put.points.demo()
#ex 3
library(multtest); data(golub)
x <- golub[2289,]; y <- golub[2430,]                   
cor(x,y)                       
cor.test(x,y)

#Ex 4
nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1) 
data <- matrix(c(x,y),ncol=2,byrow=FALSE) 
for (i in 1:nboot){ 
  dat.star <- data[sample(1:nrow(data),replace=TRUE),] 
  boot.cor[i,] <- cor(dat.star)[2,1]}
mean(boot.cor) 
quantile(boot.cor[,1],c(0.025,0.975)) 

#ex 5
library(multtest); data(golub) 
corgol<- apply(golub, 1, function(x) cor(x,golub.cl)) 
o <- order(corgol)

# Principal Components Analysis
#ex 1 - Compute correlation matrix
Z <- matrix(c( 1.63, 1.22, -0.40, 0.79, 0.93, 0.97, -1.38,
      -1.08, -0.17, -0.96, -0.61, -0.93), nrow=6, byrow=TRUE)      
K <- eigen(cor(Z))

#Output skor dalam 2 digit
print(K, digits = 2)
print(Z %*% K$vec, digits=2)
# perform principal components analysis - Primpcomp
pca <- princomp(Z, center = TRUE, cor=TRUE, scores=TRUE)
pca$scores

#ex 2 -  The ﬁrst ﬁve eigenvalues from the correlation matrix of golub 
eigen(cor(golub))$values[1:5]
data <- golub; p <- ncol(data); n <- nrow(data) ; nboot<-1000
eigenvalues <- array(dim=c(nboot,p)) 

for (i in 1:nboot){dat.star <- data[sample(1:n,replace=TRUE),]
eigenvalues[i,] <- eigen(cor(dat.star))$values}
for (j in 1:p) print(quantile(eigenvalues[,j],c(0.025,0.975)))
for (j in 1:5) cat(j,as.numeric(quantile(eigenvalues[,j],              
       c(0.025,0.975))),"\n" )   

# The percentages of variance 
sum(eigen(cor(golub))$values[1:2])/38*100  

# the weights
-eigen(cor(golub))$vec[,1:2]   

# The ﬁrst and the last 10 gene names on the second component 
pca <- princomp(golub, center = TRUE, cor=TRUE, scores=TRUE) 
o <- order(pca$scores[,2]) 
golub.gnames[o[1:10],2]
golub.gnames[o[3041:3051],2]

#ex 3 - biplot which is based on a two-dimensional approximation of the data 
#Biplot of selected genes from the golub data
biplot(princomp(data,cor=TRUE),pc.biplot=TRUE,cex=0.5,expand=0.8)

#ex 4
data(golub, package = "multtest") 
factor <- factor(golub.cl)
o1 <- grep("CD",golub.gnames[,2])
o2 <- grep("Op",golub.gnames[,2])
o3 <- grep("MCM",golub.gnames[,2]) 
o <- c(o1,o2,o3)

#use a two-sample t-test
pt <- apply(golub, 1, function(x) t.test(x ~ gol.fac)$p.value)
oo <- o[pt[o]<0.01]

# to identify genes in directions of large variation
# using scores on the first two principal components
Z <- as.matrix(scale(golub, center = TRUE, scale = TRUE))
K <- eigen(cor(Z)) 
P <- Z %*% -K$vec[,1:2] 
leu <- data.frame(P[oo,], row.names= oo) 
attach(leu)

#Analysis Cluster Hierarchical
cl <- hclust(dist(leu,method="euclidian"),method="single") 
plot(cl)
#Gene Order
a <- as.integer(rownames(leu)[cl$order])
for (i in 1:length(a)) cat(a[i],golub.gnames[a[i],2],"\n")

