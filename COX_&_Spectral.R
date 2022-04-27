### Loading Data
gene = read.csv(file="GBM/GLIO_Gene_Expression.txt", sep=",", header = T, row.names = 1)
mirna = read.csv(file="GBM/GLIO_Mirna_Expression.txt",sep=",", header = T, row.names = 1)
met = read.csv(file="GBM/GLIO_Methy_Expression_1.txt",sep=",", header = T, row.names = 1)
clinical = read.csv(file="GBM/GLIO_Survival.txt", sep="\t", header = T)


library(CancerSubtypes)

#### STEP 1: Feature selection using COX-model
data1 = FSbyCox(as.matrix(gene_colon), cli_colon$Survival, cli_colon$Death, cutoff = 0.05)
data2 = FSbyCox(as.matrix(mirna_colon), cli_colon$Survival, cli_colon$Death, cutoff = 0.05)
data3 = FSbyCox(as.matrix(met_colon), cli_colon$Survival, cli_colon$Death, cutoff = 0.05)
cat = rbind(data1,data2,data3)

write.csv(cat, file="GBM/cox_GBM.csv")

#### STEP 2:  The selected features are then passed to the Autoencoder and 
#### the compressed BL layers at different numbers of neurons are stored.
### In the next step spectral clustering is performed on each of the BLs obtained.
### Use the script Autoencoder_cancer for performing this step.

#### STEP 3: Spectral clustering using knn graphs

spectral_clustering <- function(X, # matrix of data points
                                nn=10, # the k nearest neighbors to consider
                                n_eig=3) # n number of eigen-vectors to keep
{
  mutual_knn_graph <- function(X, nn=10)
  {
    D <- as.matrix(dist(X) ) # matrix of euclidean distances between data points in X
    
    # initialize the knn matrix
    knn_mat <- matrix(0,
                      nrow = nrow(X),
                      ncol = nrow(X))
    
    # find the 10 nearest neighbors for each point
    for (i in 1: nrow(X)) {
      neighbor_index <- order(D[i,])[2:(nn + 1)]
      knn_mat[i,][neighbor_index] <- 1 
    }
    
    # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1 
    knn_mat <- knn_mat + t(knn_mat) # find mutual knn
    
    knn_mat[ knn_mat == 2 ] = 1
    
    return(knn_mat)
  }
  
  graph_laplacian <- function(W, normalized = TRUE)
  {
    stopifnot(nrow(W) == ncol(W)) 
    
    g = colSums(W) # degrees of vertices
    n = nrow(W)
    
    if(normalized)
    {
      D_half = diag(1 / sqrt(g) )
      return( diag(n) - D_half %*% W %*% D_half )
    }
    else
    {
      return( diag(g) - W )
    }
  }
  
  W = mutual_knn_graph(X) # 1. matrix of similarities
  L = graph_laplacian(W) # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  n = nrow(L)
  return(ei$vectors[,(n - n_eig):(n - 1)]) # return the eigenvectors of the n_eig smallest eigenvalues
  
}

## parameter selection for spectral clustering        
n_eig = k = 3 #number of clusters 
nn = 10 #the k nearest neighbors to consider

### Spectral clustering is performed on compressed data obtained from the BL layer of the AE
#### Silhouette score and p-value for log rank test (survival analysis) is recorded. The BL
### layer that gives maximum silhouette score is said to be the best compression and the obtained
### clusters are the resultant cancer subtypes.

en_dim = NULL
for(i in seq(10,90,10)){
en_dim = append(en_dim,round((i/100)*nrow(cat)))# varying from 10% to 90% of the input features (biomarkers) 
}
deco = NULL
for(i in seq(50,90,10)){
  deco = append(deco,round((i/100)*nrow(cat)))# varying from 50% to 90% of the input features (biomarkers)
}

### For Undercomplete AE
encod = NULL
decod = NULL
sil = NULL
p_value = NULL
days = as.numeric(clinical$Survival)
vital = as.numeric(paste(clinical$Death))
for(i in en_dim){
  for(j in deco){
    data = read.csv(paste0(file="Results/undercomplete_gbm_",i,"_",j,".csv"))
    set.seed(1993)
    X_sc <- spectral_clustering(data)
    X_sc_kmeans <- kmeans(X_sc, k, iter.max = 100)
    sil = append(sil, mean(cluster::silhouette(X_sc_kmeans$cluster,dist(X_sc))[,3]))
    p_value = append(p_value, survAnalysis(mainTitle = "Survival Analysis", days, vital, X_sc_kmeans$cluster))
    encod = append(encod,i)
    decod = append(encod,j)
  }
}
print(as.data.frame(c(encod,decod,sil,p_value)))

### For Sparse AE
regu = c(0.01,0.001,0.0001,0.00001,0.000001)

encod = NULL
decod = NULL
lambda = NULL
sil = NULL
p_value = NULL
days = as.numeric(clinical$Survival)
vital = as.numeric(paste(clinical$Death))
for(i in en_dim){
  for(j in deco){
    for(k in regu){
      data = read.csv(paste0(file="Results/sparse_gbm_",i,"_",j,"_",k,".csv"))
      set.seed(1993)
      X_sc <- spectral_clustering(data)
      X_sc_kmeans <- kmeans(X_sc, k, iter.max = 100)
      sil = append(sil, mean(cluster::silhouette(X_sc_kmeans$cluster,dist(X_sc))[,3]))
      p_value = append(p_value, survAnalysis(mainTitle = "Survival Analysis", days, vital, X_sc_kmeans$cluster))
      encod = append(encod,i)
      decod = append(encod,j)
      lambda = append(lambda,k)
    }
  }
}
print(as.data.frame(c(encod,decod,lambda,sil,p_value)))











