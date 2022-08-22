######################################
# Breast Cancer Diagnostic
######################################
######################################
# SECTION 1 : DATA SET
######################################
########################
## Read data 
########################

library(readr)
data <- read.csv("breastcancer.csv")
dim(data)
head(data)

########################
# Create matrix X
########################

#removing unused columns -> ID and gender
X <- data[, -c(1, 2)]
head(X)
dim(X)

# number of observations (rows)
n <- dim(X)[1]
n

# center and scaled data
X.cs <- scale(X, scale = T)
head(X.cs)
dim(X.cs)

########################
# Removing outliers
########################

# removing outliers
# from the plots below, I observe there are only points 
# beyond -3 or +3 so, I will be classifying those as outliers
# and removing them 
# mean of characteristics
pairs(X.cs[,1:10], lower.panel = NULL,
      col = c("red", "blue")[factor(data$diagnosis)],
      main = "Mean of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

# se of characteristics
pairs(X.cs[,11:20], lower.panel = NULL,
      col = c("red", "blue")[factor(data$diagnosis)],
      main = "SE of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

# "worst" of characteristics
pairs(X.cs[,21:30], lower.panel = NULL,
      col = c("red", "blue")[factor(data$diagnosis)],
      main = "Worst of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

# replacing outliers with NA values
X.cs <- apply(X.cs, 2, function(x){
  x <- ifelse(x > 3 | x < -3, NA, x)
})

# rows that contain variables that have outliers 
row.outliers <- which(apply(X.cs, 1, function(x) {sum(is.na(x)) > 0}) == TRUE )
row.outliers

# total number of outliers
length(row.outliers)

# removing outliers
X.cs <- na.omit(X.cs)      
dim(X.cs)

# checking that outliers are removed
sum(apply(X.cs, 2, is.na)) == 0 
sum (apply(X.cs, 2, function(x) {
  any(x < -3 | x > 3)
})) == 0

# new original data set 
data.new <- data[-row.outliers,]

# pairs plot after removing outliers
pairs(X.cs[,1:10], lower.panel = NULL,
      col = c("red", "blue")[factor(data.new$diagnosis)],
      main = "Mean of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

pairs(X.cs[,11:20], lower.panel = NULL,
      col = c("red", "blue")[factor(data.new$diagnosis)],
      main = "SE of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

pairs(X.cs[,21:30], lower.panel = NULL,
      col = c("red", "blue")[factor(data.new$diagnosis)],
      main = "Worst of Characteristics Explanatory Pairs Plot")
legend("bottomleft",c("benign", "malignant"),pch=c("R","B"),
       col=c("red","blue"))

########################
# Final data set 
########################

# After removing 74 observations, we obtain the final data set
dim(data.new)
n <- dim(data.new)[1]
n
head(data.new)

# Final centered and scaled matrix after removing outlier
X <- data.new[, -(1:2)]
X.cs <- scale(X, scale = TRUE)
head(X.cs)
dim(X.cs)

######################################
# SECTION 2 : CLUSTER ANALYSIS
######################################
########################
# k-means algorithm
########################
# Determine initial means using distance matrix, D
# choose the two most distant observations based on the distance matrix. 
D <- dist(X.cs)
class(D)

## tells us largest distance between rows. 
max(D) 

## find which are the rows of the data at that distance
D1 <- as.matrix(D)
rows.most.distant <- matrix(c(0, 0), ncol = 1) 
rows.most.distant

rows <- which(apply(D1, 1, function(x){sum(x == max(D1)) == 1}) == TRUE)

## Initial values: use row 85 and row 492
rows.most.distant <- matrix(c(rows[1], rows[2]), ncol = 1) 
rows.most.distant

## view initial vector for cluster 1 
X.cs[85, ] 
## view initial vector for cluster 2
X.cs[492, ] 

## Initial means, mean^(0), for variables in cluster 1
## mean^(0) for cluster 1 must be numeric
c1 <- as.numeric(as.vector(X.cs[85,])) 
c1

# Initial means, mean^(0), for the two variables in cluster 2 
# mean^(0)  for cluster 2 must be numeric
c2 <- as.numeric(as.vector(X.cs[492,])) 
c2

## Two index vectors:
## one will tell us clusters assigned in last iteration
## another index vector will tell cluster in the new iteration
## If z_(t-1) != z_(t), must continue allocating. 

## z_(t-1)
## Initial value for z
pastIndicator <- n:1 
## z_(t)
## Past indicator will be compared with new indicator
indicator <- 1:n   

while(sum(pastIndicator!=indicator)!=0) {
  pastIndicator <- indicator
  
  ## Distance to current cluster centers of each row 
  dc1 <- colSums((t(X.cs)-c1)^2) # distance of each from from mean cluster 1
  dc2 <- colSums((t(X.cs)-c2)^2) # distance of each from from mean cluster 2
  dMat <- matrix(c(dc1,dc2), ncol=2)
  
  ## Decide which cluster each point (each row) belongs to 
  indicator <- max.col(-dMat)
  
  ## Update the cluster centers
  c1 <- colMeans(X.cs[indicator==1,])
  c2 <- colMeans(X.cs[indicator==2,])
}

## Double check that z_(t-1) = z(t) (convergence)
indicator
all(indicator == pastIndicator)

########################
# Summary statitsics
########################
## Table of which type of breast cancer in each cluster 
cluster <- indicator
table <- table(cluster, data.new$diagnosis)
table

## proportion of type of breast cancer in each cluster (3dp)
round(table / c(sum(table[1,]), sum(table[2,])), 3)

## The final mean vector 
## cluster 1
c1   
round(c1,3)
## cluster 2
c2  
round(c2,3)
## the mean vector to which we converge for each cluster
mean_vec <- rbind(c1,c2)
round(mean_vec, 3)

## Intrepretable summary 
## combine cluster and original data 
data.new <-  cbind(data.new, cluster)
dim(data.new)
head(data.new)

## mean vector for each cluster  
colMeans(data.new[cluster == 1, -c(1, 2, 33)])
colMeans(data.new[cluster == 2, -c(1, 2, 33)])

## median vector for each cluster  
apply(data.new[cluster == 1, -c(1, 2, 33)], 2, median)
apply(data.new[cluster == 2, -c(1, 2, 33)], 2, median)

## standard deviation vector for each cluster 
apply(data.new[cluster == 1, -c(1, 2, 33)], 2, sd)
apply(data.new[cluster == 2, -c(1, 2, 33)], 2, sd)

########################
# Plot
########################
names(X)
mass_type <- factor(data.new$diagnosis, labels = c("benign", "malignant"))

par(mfrow=c(1,2))
# plot of mean area and mean smoothness 
plot(X[, 4], X[, 5], 
     col = c("red", "blue")[cluster],
     pch = c(25, 21)[unclass(mass_type)], 
     main = c("(A) Area clearly separates \n but not smoothness"),
     cex = 0.5,
     xlab = "mean area",
     ylab = "mean smoothness")

legend("bottomright",
       c("cluster 1","cluster 2"),
       pch = c(25, 21),
       col = c("red", "blue"),
       cex = 0.75)

# plot of mean concavity and mean symmetry 
plot(X[, 7], X[, 9], 
     col = c("red", "blue")[cluster],
     pch = c(25, 21)[unclass(mass_type)], 
     main = c("(B) concavity clearly separates \n but not symmetry"),
     cex = 0.5,
     xlab = "mean concavity",
     ylab = "mean symmetry")

legend("bottomright",
       c("cluster 1","cluster 2"),
       pch = c(25, 21),
       col = c("red", "blue"),
       cex = 0.75)

# plot of mean area and mean concavity  
plot(X[, 4], X[, 7], 
     col = c("red", "blue")[cluster],
     pch = c(25, 21)[unclass(mass_type)], 
     main = c("(C) Cluster 1 has smaller area \n and is less concave"),
     cex = 0.5,
     xlab = "mean area",
     ylab = "mean concavity")

legend("bottomright",
       c("cluster 1","cluster 2"),
       pch = c(25, 21),
       col = c("red", "blue"),
       cex = 0.75)

# plot of mean area and mean concavity points 
plot(X[, 4], X[, 8], 
     col = c("red", "blue")[cluster],
     pch = c(25, 21)[unclass(mass_type)], 
     main = c("(D) Cluster 1 has smaller area \n and concavity points"),
     cex = 0.5,
     xlab = "mean area",
     ylab = "mean concavity points")

legend("bottomright",
       c("cluster 1","cluster 2"),
       pch = c(25, 21),
       col = c("red", "blue"),
       cex = 0.75)

######################################
# SECTION 3 : PCA
######################################
# variance covariance matrix of Xcs
Sxcs <- var(X.cs)
round(Sxcs, 3)

# pairs plot of first seven variables to illustrate 
# that there is high correlation 
# it is worth it to conduct pca
pairs(X.cs[,1:7], lower.panel = NULL,
      col=c("red","blue")[factor(data.new[,2])],
      pch=c(8,18)[factor(data.new[,2])],
      main="pairs plot of data\n (standard deviation units)")

legend("bottomleft",c("benign","malignant"),
       pch=c("R","B"),
       col=c("red","blue"))

# eigenvalues (variances of the PC) and 
# eigenvectors (the axes of the projection)
EP <- eigen(Sxcs)
lambda <- EP$values
round(lambda,3)

# sum of eigenvalues 
sum(lambda)

# proportion of variance of each PC
prop <- 100*(lambda/sum(lambda))
round(prop, 3)
round(cumsum(prop), 3) #7 PC to capture 90% variability 

# eigenvectors
V <- EP$vectors
round(V,3) 

# check that eigenvectors are orthonormal
round(crossprod(V,V))

##########################
# Obtain the principal components 
###########################
## Note: X.tilde is the PC matrix 
X.tilde <- X.cs %*% V
round(X.tilde, 2)

head(X.tilde)
colnames(X.tilde)=c("PC1","PC2","PC3","PC4", "PC5", "PC6", "PC7",
                    "PC8", "PC9", "PC10", "PC11", "PC12", "PC13",
                    "PC14", "PC15", "PC16", "PC17", "PC18", "PC19",
                    "PC20", "PC21", "PC22", "PC23", "PC24", "PC25",
                    "PC26", "PC27", "PC28", "PC29", "PC30")
head(X.tilde)
round(X.tilde, 3)

#######Check that the PC are orthogonal

round(t(X.tilde)%*%X.tilde,3)

#### Check the variance-covariance  of the PC 
round(var(X.tilde),3)

scale(X.tilde,scale=T)

############Pairwise plot of the PCs

pairs(X.tilde[,1:7], lower.panel = NULL,
      col=c("red","blue")[factor(data.new[,2])],
      pch=c(8,18)[factor(data.new[,2])],
      main="pairs plot of PC of data \n (standard deviation units)")

legend("bottomleft",c("benign","malignant"),
       pch=c("R","B"),
       col=c("red","blue"))


### Latent variables interpretation 
cor(X.tilde,X.cs)
t(cor(X.tilde,X.cs))

######################################
# SECTION 4 : Model selection
######################################
###########################
#find probabilistic model
###########################
par(mfrow = c(2,2))
# define discrete values of x over specified range
x = c(seq(0,max(data.new$area_mean),by=0.1))

# histogram of area_mean
# simulate Log normal distributions 
# hist 1 
hist(data.new$area_mean, prob=T, ylim = c(0, 0.003), 
     xlab = "mean area", main = "histogram of area mean \n meanlog = 6, sdlog = 0.8")
points(x,dlnorm(x, meanlog=6, sdlog=0.8, log=FALSE), 
       col="blue", type="o", pch=24, bg="blue", cex = 0.25)

# hist 2
hist(data.new$area_mean, prob=T, ylim = c(0, 0.003), 
     xlab = "mean area", main = "histogram of area mean \n meanlog = 6.5, sdlog = 0.5")
points(x,dlnorm(x, meanlog=6.5, sdlog=0.5, log=FALSE), 
       col="purple", type="o", pch=24, bg="purple", cex = 0.25)

# hist 3
hist(data.new$area_mean, prob=T, ylim = c(0, 0.003), 
     xlab = "mean area", main = "histogram of area mean \n meanlog = 6.4, sdlog = 0.3")
points(x,dlnorm(x, meanlog=6.4, sdlog=0.3, log=FALSE), 
       col="red", type="o", pch=24, bg="red", cex = 0.25)

# hist 4 (best fit - meanlog = 6.3, sdlog = 0.3)
hist(data.new$area_mean, prob=T, ylim = c(0, 0.003), 
     xlab = "mean area", main = "histogram of area mean \n meanlog = 6.3, sdlog = 0.3")
points(x,dlnorm(x, meanlog=6.3, sdlog=0.3, log=FALSE), 
       col="green", type="o",pch=21, bg="green", cex = 0.25)

###########################
# Newton's algorithm
###########################
y <- data.new$area_mean
n <- length(y)
n

# Function 
# Let p[1] denote mu and p[2] denote sigma
fn <- function(p, y, n){
  -sum(log(y)) -n*log(p[2])-(1/(2*(p[2])^2 ))*sum((log(y)-p[1])^2)
}

##################################################
# (a) Finding initial values 
##################################################
# using meanlog = 6.3, sdlog = 0.3 based on histogram 

##################################################
# (b) Use Newton's method to find critical values
##################################################
## Initial values for xt
xt <- c(0,0)
## tolerance for (xtp1 - xt)
tol <- 0.0000000001
## Initial values for x1 and x2 
xtp1 <- c(6.3, 0.3)
## Saves history of xt
## : row 1 and 2 are values of x1 and x2 respectively
xHist <- matrix(xtp1, nrow = 2, ncol = 1)
## Objective function 
f <- fn(xtp1, y, n) 
## History of objective function 
fHist <- f
## History of gradient 
gHist <- c()

iter <- 0 
while (sum((xtp1 - xt)^2) > tol) {
  xt <- xtp1
  gradient <- as.vector(c(sum(log(y) - xt[1])/xt[2]^2,
                          -n/xt[2] + sum((log(y)-xt[1])^2)/xt[2]^3))
  hessian <- matrix(c(-n/xt[2]^2, -2*sum(log(y)-xt[1])/xt[2]^3,
                      -2*sum(log(y)-xt[1])/xt[2]^3, n/xt[2]^2 - 3*sum((log(y)-xt[1])^2)/xt[2]^4),
                    nrow = 2, ncol = 2)
  xtp1 <- xt - solve(hessian)%*%gradient
  xHist <- matrix(c(xHist, xtp1), nrow = 2)
  f <- fn(xtp1, y, n)
  fHist <- c(fHist, f)
  gHist <- cbind(gHist, gradient)
  iter <- iter + 1
}

## History of x1 and x2 
xHist

## History of objective function
fHist

## History of gradient
gHist

## Hessian matrix at convergence
hessian

##################################################
# (c) Number of iterations and gradient
##################################################
## The number of iterations is 6
iter 

## Value of x1 and x2 at convergence
xHist[,7]

## Value of gradient 
gradient

##################################################
# (d) Finding eigenvalues
##################################################
eigen(hessian)$values

## Since both the eigenvalues are negative, the critical
## point (6.3382741, 0.4329019) is a maximum point

##################################################
# (e) 95% confidence interval 
##################################################
pt_estimate <- c(xHist[,7])
se <- sqrt(diag(solve(-hessian)))
# CI for mu 
CI_mu <- pt_estimate[1] + c(-1,1)*qnorm(0.975)*se[1]
CI_mu
# CI for sigma 
CI_sigma <- pt_estimate[2] + c(-1,1)*qnorm(0.975)*se[2]
CI_sigma

##################################################
# (f) Model with MLE estimates 
##################################################
dev.off()
hist(data.new$area_mean, prob=T,main="Log normal with MLE estimates",
     xlab="area mean ")
points(x,dlnorm(x, meanlog=pt_estimate[1],
                sdlog=pt_estimate[2], log=FALSE), col="red", type="o", pch=21, bg="red")

