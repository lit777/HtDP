

## check whether packages are installed, then load
pacotes = c("ggplot2", "gridExtra", "invgamma", "dplyr", "MASS", "LaplacesDemon",
            "armspp", "CholWishart", "mvtnorm", "mnormt", "progress",
            "Rcpp", "mniw", "plyr", "mcclust", "cluster")

package.check <- lapply(pacotes, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
})

## set working directory
CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

## load functions and datasets
load("code_functions.RData")
sourceCpp('HtDP_update.cpp')
sourceCpp('CDP_update.cpp')
## If you get error, please put the path directly on your own.
#sourceCpp('.../HtDP_update.cpp')
#sourceCpp('.../CDP_update.cpp')
sourceCpp('C:/rcpp/HtDP_update.cpp')
sourceCpp('C:/rcpp/CDP_update.cpp')

geyser <- read.table('old_faithful.txt', header=TRUE)
syn_geyser <- read.table('syn_old_faithful.txt', header=TRUE)
toy1 <- read.table('toy1.txt', header=TRUE)
toy2 <- read.table('toy2.txt', header=TRUE)
toy3 <- read.table('toy3.txt', header=TRUE)


#-----------------------------------------------------------

## Old Faithful dataset modeling
iter <- 5000
xx <- 1:iter

CDP_geyser <- CDP(geyser, iter)
HtDP_geyser <- HtDP(geyser, iter)


## PAM
cdp.mat <- matrix(CDP_geyser$indicator[[iter-250+1]], nrow=1, ncol=nrow(geyser))
for (i in 2:250) {
  cdp.mat <- rbind(cdp.mat, CDP_geyser$indicator[[iter-250+i]])
}
cdp.psm <- comp.psm(cdp.mat)
htdp.mat <- matrix(HtDP_geyser$indicator[[iter-250+1]], nrow=1, ncol=nrow(geyser))
for (i in 2:250) {
  htdp.mat <- rbind(htdp.mat, HtDP_geyser$indicator[[iter-250+i]])
}
htdp.psm<-comp.psm(htdp.mat)

cdp.best.sim<-NULL  
for(i in 1:9){
  cdp.best.sim[i] <- pam(1-cdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
cdp.num.clus <- which.max(cdp.best.sim)
cdp.best <- pam(1-cdp.psm, cdp.num.clus, diss=TRUE)
htdp.best.sim<-NULL  
for(i in 1:9){
  htdp.best.sim[i] <- pam(1-htdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
htdp.num.clus <- which.max(htdp.best.sim)
htdp.best <- pam(1-htdp.psm, htdp.num.clus, diss=TRUE)

geyser_df <- cbind(geyser, cdp.best$clustering, htdp.best$clustering)
geyser_df <- as.data.frame(geyser_df)
colnames(geyser_df)[3:4] <- c('CDP', 'HtDP')
geyser_df$CDP <- as.factor(geyser_df$CDP)
geyser_df$HtDP <- as.factor(geyser_df$HtDP)
cdp_geyser <- ggplot(geyser_df, aes(x=Previous_Duration, y=Current_Duration, fill=CDP, colour=CDP)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
htdp_geyser <- ggplot(geyser_df, aes(x=Previous_Duration, y=Current_Duration, fill=HtDP, colour=HtDP)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
grid.arrange(cdp_geyser, htdp_geyser, ncol=2)


#-----------------------------------------------------------

## Synthetic Old Faithful dataset modeling
CDP_syn.geyser <- CDP(syn_geyser, iter)
HtDP_syn.geyser <- HtDP(syn_geyser, iter)


## PAM
cdp.mat <- matrix(CDP_syn.geyser$indicator[[iter-250+1]], nrow=1, ncol=nrow(syn_geyser))
for (i in 2:250) {
  cdp.mat <- rbind(cdp.mat, CDP_syn.geyser$indicator[[iter-250+i]])
}
cdp.psm <- comp.psm(cdp.mat)
htdp.mat <- matrix(HtDP_syn.geyser$indicator[[iter-250+1]], nrow=1, ncol=nrow(syn_geyser))
for (i in 2:250) {
  htdp.mat <- rbind(htdp.mat, HtDP_syn.geyser$indicator[[iter-250+i]])
}
htdp.psm<-comp.psm(htdp.mat)

cdp.best.sim<-NULL  
for(i in 1:9){
  cdp.best.sim[i] <- pam(1-cdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
cdp.num.clus <- which.max(cdp.best.sim)
cdp.best <- pam(1-cdp.psm, cdp.num.clus, diss=TRUE)
htdp.best.sim<-NULL  
for(i in 1:9){
  htdp.best.sim[i] <- pam(1-htdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
htdp.num.clus <- which.max(htdp.best.sim)
htdp.best <- pam(1-htdp.psm, htdp.num.clus, diss=TRUE)

df <- cbind(syn_geyser, cdp.best$clustering, htdp.best$clustering)
df <- as.data.frame(df)
colnames(df)[3:4] <- c('CDP', 'HtDP')
df$CDP <- as.factor(df$CDP)
df$HtDP <- as.factor(df$HtDP)
cdp_syn.geyser <- ggplot(df, aes(x=Previous_Duration, y=Current_Duration, fill=CDP, colour=CDP)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
htdp_syn.geyser <- ggplot(df, aes(x=Previous_Duration, y=Current_Duration, fill=HtDP, colour=HtDP)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
grid.arrange(cdp_syn.geyser, htdp_syn.geyser, ncol=2)


#-----------------------------------------------------------

## toy 1 modeling
iter <- 5000
xx <- 1:iter
CDP_toy1 <- CDP(toy1[,-3], iter)
HtDP_toy1 <- HtDP(toy1[,-3], iter)

## PAM
cdp.mat <- matrix(CDP_toy1$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy1))
for (i in 2:250) {
  cdp.mat <- rbind(cdp.mat, CDP_toy1$indicator[[iter-250+i]])
}
htdp.mat <- matrix(HtDP_toy1$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy1))
for (i in 2:250) {
  htdp.mat <- rbind(htdp.mat, HtDP_toy1$indicator[[iter-250+i]])
}
cdp.psm <- comp.psm(cdp.mat)
htdp.psm <- comp.psm(htdp.mat)

cdp.best.sim<-NULL; htdp.best.sim<-NULL    
for(i in 1:9){
  cdp.best.sim[i] <- pam(1-cdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
for(i in 1:9){
  htdp.best.sim[i] <- pam(1-htdp.psm, i, diss=TRUE)$silinfo$avg.width
}  

cdp.num.clus <- which.max(cdp.best.sim)
cdp.best <- pam(1-cdp.psm, cdp.num.clus, diss=TRUE)
htdp.num.clus <- which.max(htdp.best.sim)
htdp.best <- pam(1-htdp.psm, htdp.num.clus, diss=TRUE)

df <- cbind(toy1, cdp.best$clustering, htdp.best$clustering)
df <- as.data.frame(df)
colnames(df)[4:5] <- c('cdp_cluster', 'htdp_cluster')
df$toy1_true <- as.factor(df$toy1_true)
df$cdp_cluster <- as.factor(df$cdp_cluster)
df$htdp_cluster <- as.factor(df$htdp_cluster)


true <- ggplot(df, aes(x=X1, y=X2, fill=toy1_true, colour=toy1_true)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
cdp <- ggplot(df, aes(x=X1, y=X2, fill=cdp_cluster, colour=cdp_cluster)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
htdp <- ggplot(df, aes(x=X1, y=X2, fill=htdp_cluster, colour=htdp_cluster)) + 
  geom_point(cex=2) + 
  coord_fixed() + 
  theme(legend.position = "none")
grid.arrange(true, cdp, htdp, ncol=3)


#-----------------------------------------------------------

## toy 2 modeling
iter <- 5000
xx <- 1:iter
CDP_toy2 <- CDP(toy2[,-3], iter)
HtDP_toy2 <- HtDP(toy2[,-3], iter)


## PAM
cdp.mat <- matrix(CDP_toy2$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy2))
for (i in 2:250) {
  cdp.mat <- rbind(cdp.mat, CDP_toy2$indicator[[i+iter-250]])
}
htdp.mat <- matrix(HtDP_toy2$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy2))
for (i in 2:250) {
  htdp.mat <- rbind(htdp.mat, HtDP_toy2$indicator[[i+iter-250]])
}
cdp.psm <- comp.psm(cdp.mat)
htdp.psm <- comp.psm(htdp.mat)

cdp.best.sim<-NULL; htdp.best.sim<-NULL 
for(i in 1:9){
  cdp.best.sim[i] <- pam(1-cdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
for(i in 1:9){
  htdp.best.sim[i] <- pam(1-htdp.psm, i, diss=TRUE)$silinfo$avg.width
}  

cdp.num.clus <- which.max(cdp.best.sim)
cdp.best <- pam(1-cdp.psm, cdp.num.clus, diss=TRUE)
htdp.num.clus <- which.max(htdp.best.sim)
htdp.best <- pam(1-htdp.psm, htdp.num.clus, diss=TRUE)

df <- cbind(toy2, cdp.best$clustering, htdp.best$clustering)
df <- as.data.frame(df)
colnames(df)[4:5] <- c('cdp_cluster', 'htdp_cluster')
df$toy2_true <- as.factor(df$toy2_true)
df$cdp_cluster <- as.factor(df$cdp_cluster)
df$htdp_cluster <- as.factor(df$htdp_cluster)
  
  
true <- ggplot(df, aes(x=X1, y=X2, fill=toy2_true, colour=toy2_true)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
cdp <- ggplot(df, aes(x=X1, y=X2, fill=cdp_cluster, colour=cdp_cluster)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
htdp <- ggplot(df, aes(x=X1, y=X2, fill=htdp_cluster, colour=htdp_cluster)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
grid.arrange(true, cdp, htdp, ncol=3)


#-----------------------------------------------------------

## toy 3 modeling

CDP_toy3 <- CDP(toy3[,-3], iter)
HtDP_toy3 <- HtDP(toy3[,-3], iter)


# PAM
cdp.mat <- matrix(CDP_toy3$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy3))
for (i in 2:250) {
  cdp.mat <- rbind(cdp.mat, CDP_toy3$indicator[[i+iter-250]])
}
htdp.mat <- matrix(HtDP_toy3$indicator[[iter-250+1]], nrow=1, ncol=nrow(toy3))
for (i in 2:250) {
  htdp.mat <- rbind(htdp.mat, HtDP_toy3$indicator[[i+iter-250]])
}
cdp.psm <- comp.psm(cdp.mat)
htdp.psm <- comp.psm(htdp.mat)

cdp.best.sim<-NULL ; htdp.best.sim<-NULL 
for(i in 1:9){
  cdp.best.sim[i] <- pam(1-cdp.psm, i, diss=TRUE)$silinfo$avg.width
}  
for(i in 1:9){
  htdp.best.sim[i] <- pam(1-htdp.psm, i, diss=TRUE)$silinfo$avg.width
}  

cdp.num.clus <- which.max(cdp.best.sim)
cdp.best <- pam(1-cdp.psm, cdp.num.clus, diss=TRUE)
htdp.num.clus <- which.max(htdp.best.sim)
htdp.best <- pam(1-htdp.psm, htdp.num.clus, diss=TRUE)

df <- cbind(toy3, cdp.best$clustering, htdp.best$clustering)
df <- as.data.frame(df)
colnames(df)[4:5] <- c('cdp_cluster', 'htdp_cluster')
df$toy3_true <- as.factor(df$toy3_true)
df$cdp_cluster <- as.factor(df$cdp_cluster)
df$htdp_cluster <- as.factor(df$htdp_cluster)
  
  
true <- ggplot(df, aes(x=X1, y=X2, fill=toy3_true, colour=toy3_true)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
cdp <- ggplot(df, aes(x=X1, y=X2, fill=cdp_cluster, colour=cdp_cluster)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
htdp <- ggplot(df, aes(x=X1, y=X2, fill=htdp_cluster, colour=htdp_cluster)) + 
    geom_point(cex=2) + 
    coord_fixed() + 
    theme(legend.position = "none")
grid.arrange(true, cdp, htdp, ncol=3)
























