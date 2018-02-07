# setwd("C:/Users/jyotishka/Documents/R/modelselect")
# setwd("X:/R/modelselect")

setwd("C:/Users/Jyotishka Datta/Dropbox/hs_review/R")
rm(list=ls())
library(glmnet)
# library(BayesVarSel)
library(monomvn)
library(bayesm)
library(plyr)
library(lars)
library(horseshoe)
library(data.table)
library("parallel")
library("foreach")
library("doParallel")

no_cores <- detectCores()-1
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(no_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 13)


niter = 100
n_inf = rep(0,niter)
prop_true_hs = rep(0,niter)
prop_true_lasso = rep(0,niter)

# pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)

n=100;p=32
p1 = 5
beta <- c(c(7,4,2,1,1),rep(0,p-p1))
index <- c(rep(1,p1),rep(0,p-p1))

res<- foreach(i = 1:niter, .packages = c("bayesm","horseshoe","glmnet","data.table"))%dopar%{
  
    # sigma <- rwishart(p,diag(p))
    rho = 0.7
    sigma = rho^{abs((1:p)%*%t(rep(1,p)) - t((1:p)%*%t(rep(1,p))))}
    #different X matrix
    x=matrix(rnorm(n*p),n,p)
    #R <- chol(sigma$W)
    R <- chol(sigma)
    x <- x%*%R
    fx <- x%*% beta
    A <- t(x[,-(1:p1)])%*%x[,(1:p1)]%*%solve((t(x[,(1:p1)])%*%x[,(1:p1)]))%*%rep(1,p1)
    n_inf <- 1-norm(A,'i')
    mc <- setDT(melt(cor(x)))[Var1 != Var2, .SD[which.max(value)], keyby=Var2][1]$value
    reps = 100
    mp_hs = rep(0,reps)
    mp_lasso = rep(0, reps)
    
    # foreach(j = 1:reps)%dopar%{
    for (j in 1:reps)
    { 
      eps=sqrt(0.1)*rnorm(n)
      y=drop(fx+eps)
      data = data.frame(y,x)
      ## Using horseshoe
      res <- horseshoe(y, x, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys", burn = 1000, nmc = 1000)
      hs_ind <- HS.var.select(res, y, method = "intervals") #selected betas
      
      mp_hs[j] <- sum((hs_ind!=index))
      
      lasso.mod = glmnet(x,y,alpha = 1)
      cv.out = cv.glmnet(x,y,alpha = 1, nfolds=10)
      bestlam = cv.out$lambda.min
      lasso.coef = predict(lasso.mod,type ="coefficients",s=bestlam)
      lasso_cv_est = lasso.coef[-1]
      lasso_ind <- ((lasso_cv_est!=0))
      mp_lasso[j] <- sum((lasso_ind!=index))
    }
    # stopCluster(cl)
    prop_true_hs <- sum((mp_hs==0))/reps
    prop_true_lasso <- sum((mp_lasso==0))/reps
    
    # info <- sprintf("%d%% done", round((i/niter)*niter))
    # setWinProgressBar(pb, i/(niter)*niter, label=info)
    
    cbind(mc,n_inf,prop_true_hs,prop_true_lasso)
}

# close(pb)
stopCluster(cl)
registerDoSEQ()


## prop_true <- cbind(n_inf,prop_true_lasso,prop_true_hs)
reps = 100

prop_true <- matrix(unlist(res), nrow = reps, ncol = 4, byrow = TRUE)
write.table(prop_true, "prop_true_100_100_parallel_with_mc.csv", sep=",")

####

eta.data <- rbind(data.frame(prop = prop_true[,3], type = "Horseshoe"),
                  data.frame(prop = prop_true[,4], type = "Lasso"))
eta.data <- cbind(eta.data, eta = prop_true[,2])

library(ggplot2)
eta.plot = ggplot(eta.data, aes(x=eta, y = prop)) + ylim(0,1)+
  geom_point(size=1) + facet_grid(type ~ ., scale="free")+ylab("Percentage of Correct Selection")+
  xlab(expression(eta[infty]))+theme_bw()

eta.plot <- eta.plot+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
eta.plot <- eta.plot + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.2)))

eta.plot<- eta.plot+ theme(axis.text = element_text(size = rel(1.2)))
eta.plot <- eta.plot+theme(strip.text.x = element_text(size=13.5, face="bold"),
                           strip.text.y = element_text(size=13.5, face="bold"))+theme(legend.position="none")
#   theme(legend.position = c(0.5,0.5),legend.key = element_rect(colour = "black"),
#         legend.text=element_text(size=rel(2)),legend.title=element_blank())
print(eta.plot)

ggsave("irrepresentability_effect_100_100.eps", eta.plot, width=7, height=5)
ggsave("irrepresentability_effect_100_100.pdf", eta.plot, width=7, height=5)

#### MC Plot

mc.data <- rbind(data.frame(prop = prop_true[,3], type = "Horseshoe"),
                  data.frame(prop = prop_true[,4], type = "Lasso"))
mc.data <- cbind(mc.data, mc = prop_true[,1])

library(ggplot2)
mc.plot = ggplot(mc.data, aes(x=mc, y = prop)) + ylim(0,1)+xlim(0,1)+
  geom_point(size=1) + facet_grid(type ~ ., scale="free")+ylab("Percentage of Correct Selection")+
  xlab("mutual coherence")+theme_bw()

mc.plot <- mc.plot+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
mc.plot <- mc.plot + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.2)))

mc.plot<- mc.plot+ theme(axis.text = element_text(size = rel(1.2)))
mc.plot <- mc.plot+theme(strip.text.x = element_text(size=13.5, face="bold"),
                           strip.text.y = element_text(size=13.5, face="bold"))+theme(legend.position="none")
#   theme(legend.position = c(0.5,0.5),legend.key = element_rect(colour = "black"),
#         legend.text=element_text(size=rel(2)),legend.title=element_blank())
print(mc.plot)

ggsave("mc_effect_100_100.eps", mc.plot, width=7, height=5)
ggsave("mc_effect_100_100.pdf", mc.plot, width=7, height=5)


library(cowplot)
comb.plot<-plot_grid(eta.plot,mc.plot, labels = "AUTO")
ggsave("irrep_mc_effect_100_100.eps", comb.plot, width=7, height=5)
ggsave("irrep_mc_effect_100_100.pdf", comb.plot, width=7, height=5)


# library(data.table)
# setDT(melt(cor_mat))[Var1 != Var2, .SD[which.max(value)], keyby=Var2]
