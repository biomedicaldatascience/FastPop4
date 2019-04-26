########################################################################
#   PCA Scoring Algorithm using SNP Weights   - updated 4/24/2019      #
#   created by Jinyoung Byun, Younghun Han, and Chris Amos             #
########################################################################
PCAScore <- function(infile1,infile2,discovery.PCs,hapmap.mean,discovery.adj,outfile){
 #read SNP weights and data file coded with additive components
  SNP.weights <- read.table(infile1, header =FALSE)
  plink.raw <- read.table(infile2, header =TRUE)
  pca.discovery <- read.table(discovery.PCs,header=FALSE)
  hapmap.mean <- read.table(hapmap.mean,header=FALSE)
  adj1 <- read.table(discovery.adj,header=FALSE)
  k=2318
  U <- plink.raw[,7:length(plink.raw)]
  # Replace the missing additive components with each original SNP mean from the discovery set.
  for(j in 1:k){ replace(U[,j], is.na(U[,j]) , adj1[j,1]) -> U[,j] }

  n2<-nrow(plink.raw)
  vec1<- as.vector(rep(1,n2))

  U <- as.matrix(U)
  SNP.weights <- as.matrix(SNP.weights)
  U_std <- (U - vec1 %*% t(adj1[,1]))/ vec1%*% t(adj1[,2])
  pca.new <- U_std%*%SNP.weights
  
  png(outfile)
  par(mfrow=c(3,1))
  plot(pca.discovery[,1],pca.discovery[,2],type="p",col="cyan",pch=1,xlab="PC1", ylab="PC2")
  title(main="Plot on discovery samples and Hapmap samples", sub="on 2318 markers",col.main="black", font.main=1)
  plot(pca.new[,1],pca.new[,2],col="pink",pch="*",xlab="PC1", ylab="PC2",type="p")
  title(main="Plot on new samples using SNP Weights of discovery samples", sub="on 2318 markers",col.main="black", font.main=1)
  
  plot(pca.discovery[,1],pca.discovery[,2],type="p",col="grey",pch=1,xlab="PC1", ylab="PC2")
  points(pca.new[,1],pca.new[,2],col="pink",pch="*")
  points(pca.discovery[1:165,1],pca.discovery[1:165,2],col="red",pch=20);#CEU
  points(pca.discovery[166:302,1],pca.discovery[166:302,2],col="green",pch=20);#CHB_Asian
  points(pca.discovery[303:505,1],pca.discovery[303:505,2],col="cyan",pch=20);#YRI_African
  points(pca.discovery[506:548,1],pca.discovery[506:548,2],col="magenta",pch=20);#Native American
  points(pca.discovery[549:601,1],pca.discovery[549:601,2],col="seagreen",pch=20);#Mexican
  points(hapmap.mean[,1],hapmap.mean[,2],col="black",pch=16)
  title(main="Prediction of new samples using SNP Weights from discovery set", sub="num.discovery=46,875 on 2318 markers",col.main="black", font.main=1)
  dev.off()
  
  score1 <- cbind(plink.raw[,1:2],pca.new)
  final=list(score1=score1)
  return(final)
  }
  
####### Arguments  :
#         infile1  : 2318 SNP Weights as input file
#         infile2  : New samples as input file
#         discovery.PCs : Scores from discovery 46,875 samples
#         hapmap.mean : mean matrix from CEU, CHB, YRI, and NA
#         discovery.adj : standardizing the data
#         outfile : output file

y <- PCAScore(infile1="fastpop4_snpweight.txt",infile2="validation.30p.A.raw",discovery.PCs="fastpop4_pca.discovery.txt",hapmap.mean="fastpop4_popmean.txt",discovery.adj="fastpop4_adj.txt",outfile="PCAplot.png")

#The scores will be saved for calculating probability of each continental ethnicity definition;
write.table(y$score1,"new.pca.txt",row.names=FALSE,col.names=FALSE)





