source("https://bioconductor.org/biocLite.R")
biocLite("multtest")

library(multtest); data(golub)
library(MASS)
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 

data <- golub; p <- ncol(data); n <- nrow(data) ; nboot<-1000
eigenvalues <- array(dim=c(nboot,p))
for (i in 1:nboot){dat.star <- data[sample(1:n,replace=TRUE),]
eigenvalues[i,] <- eigen(cor(dat.star))$values}
for (j in 1:5) cat(j,as.numeric(quantile(eigenvalues[,j],
                                         + c(0.025,0.975))),"\n" )
sum(eigen(cor(golub))$values[1:2])/38*100

-eigen(cor(golub))$vec[,1:2]
pca <- princomp(golub, center = TRUE, cor=TRUE, scores=TRUE)
o <- order(pca$scores[,2])
golub.gnames[o[1:10],2]
Golub.gnames[o[3041:3051],2]
biplot(princomp(data,cor=TRUE),pc.biplot=TRUE,cex=0.5,expand=0.8)

#Example code for Hierarchical clustering

grep("MPO Myeloperoxidase", golub.gnames[,2])
grep("CST3 Cystatin C", golub.gnames[,2])
grep("INTERLEUKIN-8 PRECURSOR", golub.gnames[,2])
grep("Interleukin 8", golub.gnames[,2])
grep("TCL1 gene", golub.gnames[,2])

data= data.frame(golub[2664,])
zyxinclus <-hclust(dist(data, method="euclidian"), method="single")
plot(zyxinclus, labels= gol.fac)

data= data.frame(golub[2734,])
zyxinclus <-hclust(dist(data, method="euclidian"), method="single")
plot(zyxinclus, labels= gol.fac)


#K-mean
#Assessing the K number of cluster
mydata <- golub
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=1)

#K-Means on the complete golub dataset
fit <- kmeans(mydata,2)
plot(mydata,col=fit$cluster,pch=16:17,xlab='Component 1',ylab='Component 2',bty='L')
points(fit$centers,col='51',pch=8)
legend(2.2,-0.2,legend="Centroid",pch=8,col ="green",cex=0.8)

#K-Means on the gene CST3 Cystatin C
data1= data.frame(golub[829,])
initial <-as.matrix(tapply(golub[829,],gol.fac,mean),nrow =2, ncol=1)
cl <- kmeans(data1,initial, nstart=10)
cl
library("cluster")
clusplot(data1, cl$cluster,xlab="Expression levels",main=golub.gnames[829,2], color=TRUE, shade=TRUE, labels=4, lines=0,bty='L')
legend(2,0.8,legend=c("ALL","AML"),pch=49:50,col =c("blue","red"),cex=0.8)

#K-Means on the genes CST3 Cystatin C, Interleukin 8, DF D
data2=data.frame(golub[1009,])
initial <-as.matrix(tapply(golub[1009,],gol.fac,mean),nrow =2, ncol=1)
dfd <- kmeans(data2,initial, nstart=10)
dfd
data3= data.frame(golub[2663,])
initial <-as.matrix(tapply(golub[2663,],gol.fac,mean),nrow =2, ncol=1)
il8 <- kmeans(data3,initial, nstart=10)
il8
data4<-data.frame(golub[829,],golub[2663,],golub[1009,])
km.res <- kmeans(data4, 2, nstart = 25)
#instal.packages("factoextra")
library("factoextra")
fviz_cluster(km.res, data = data4, frame.type = "convex")+
  theme_minimal()

#zyxin n ccnd3 plot code

data_z_c <- data.frame(golub[1042,],golub[2124,])
zyx_ccnd3 <- kmeans(data_z_c, 2,nstart = 10)
zyx_ccnd3

plot(data_z_c,col=zyx_ccnd3$cluster+1,pch = c(16, 17)[as.numeric(gol.fac)])
legend("topright",legend=c("ALL","AML"),pch=16:17,col = c('51','red'))
points(cl$centers, pch=8,col=c('red','51'))
clusplot(data_z_c, cl$cluster, main="Clustplot of Zyxin and CCND3 Cyclin D3" ,color=TRUE, shade=TRUE, labels=2, lines=0)


#bootstrap
data= data.frame(golub[829,])
table(cl$cluster,gol.fac)
n <- nrow(data); nboot<-1000
boot.cl <- matrix(0,nrow=nboot,ncol = 2)


for (i in 1:nboot) {
  dat.star <- data[sample(1:n,replace=TRUE),]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
  
  
}

hist(boot.cl[,1],main="Histogram of cluster 1",xlab="Gene Expression values")
hist(boot.cl[,2],main="Histogram of cluster 2",xlab="Gene Expression values")

quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))
mean(boot.cl[,1])
mean(boot.cl[,2])

