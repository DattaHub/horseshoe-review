## Horseshoe for logistic regression 
## Stan code from the paper: https://projecteuclid.org/download/pdfview_1/euclid.ejs/1513306866
## "Sparsity information and regularization in the horseshoe and other shrinkage priors"
## Authors: Juho Piironen and Aki Vehtari


rm(list=ls())
library(rstan)
# set_cppo("fast")
library(ggplot2)
library(plyr)
library(reshape2)


#//// Horseshoe Logistic Regression ////

hs.logis.code = "
data {
int<lower=0> n;         // number of observations
int<lower=0> d;         // number of predictors
int<lower=0,upper=1> y[n]; // outputs
matrix[n,d] x;          // inputs
real<lower=0> scale_icept ;  // prior std for the intercept
real<lower=0> scale_global ; // scale for the half -t prior for tau
real<lower=1> nu_global ; // degrees of freedom for the half -t prior for tau
real<lower=1> nu_local ; // degrees of freedom for the half - t priors for lambdas
real<lower=0> slab_scale ; // slab scale for the regularized horseshoe
real<lower=0> slab_df ; // slab degrees of freedom for the regularized HS
}
parameters {
real beta0 ;
vector [d] z;
real<lower=0> aux1_global ;
real<lower=0> aux2_global ;
vector<lower=0>[d] aux1_local ;
vector<lower=0>[d] aux2_local ;
real<lower=0> caux ;
}
transformed parameters {
real<lower=0> tau ;         // global shrinkage parameter
vector<lower=0>[d] lambda ; // local shrinkage parameter
vector<lower=0>[d] lambda_tilde ; // ’ truncated ’ local shrinkage parameter
real<lower=0> c; // slab scale
vector[d] beta ; // regression coefficients
vector[n] f; // latent function values

lambda = aux1_local.*sqrt(aux2_local);
tau = aux1_global*sqrt(aux2_global)*scale_global ;
c = slab_scale*sqrt(caux);
lambda_tilde = sqrt (c^2*square(lambda)./(c^2 + tau^2*square(lambda)));
beta = z.*lambda_tilde*tau ;
f = beta0+x*beta ;
}
model {
// half -t priors for lambdas and tau , and inverse - gamma for c ^2
z ∼ normal(0,1);
aux1_local ∼ normal(0,1);
aux2_local ∼ inv_gamma (0.5*nu_local , 0.5*nu_local );
aux1_global ∼ normal(0,1);
aux2_global ∼ inv_gamma (0.5* nu_global , 0.5*nu_global );
caux ∼ inv_gamma (0.5* slab_df , 0.5* slab_df );
beta0 ∼ normal(0, scale_icept );
y ∼ bernoulli_logit(f);
}
"

hs.logis.fit = stan_model(model_code=hs.logis.code , model_name="HS logit")

set.seed(25) 

n = 30
d = 100
y = sample(c(0,1),30, replace = T)
x1 = ifelse(y == 1, rnorm(1,1,0.5), rnorm(1,-1,0.5))
x2 = ifelse(y == 1, rnorm(1,1,0.5), rnorm(1,-1,0.5))
x = cbind(x1,x2, matrix(rnorm(n*(d-2)),n,(d-2)))

scale_icept = 5;
scale_global = 2.5
nu_global = 3
nu_local = 3
slab_scale = 2
slab_df = 4

test.data = list('n'=n,'d'=d,'y'=y,'x'=x,
                 'scale_icept' = scale_icept,
                 'scale_global'=scale_global,
                 'nu_global'=nu_global,
                 'nu_local' = nu_local,
                 'slab_scale' = slab_scale,
                 'slab_df' = slab_df)

n.iters = 5000
n.chains = 2
rng.seed = 234
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


## Horseshoe Sampling ##

smpls.hslogit.res = sampling(hs.logis.fit, 
                            data = test.data, 
                            iter = n.iters,
                            warmup = floor(n.iters/2),
                            # init = 0,
                            seed = rng.seed, 
                            chains = n.chains)

beta.smpls.hslogit = extract(smpls.hslogit.res, pars=c("beta"), permuted=TRUE)[[1]]
beta.mean = apply(beta.smpls.hslogit, 2, mean)

# par(mfrow=c(2,1))
# plot(hs.mean.beta,col=1+y)
# plot(beta.smpls.hslogit[,1],beta.smpls.hslogit[,2],type="p")

beta0.smpls.hslogit = extract(smpls.hslogit.res, pars=c("beta0"), permuted=TRUE)[[1]]
beta0.mean = mean(beta0.smpls.hslogit)

f.all = rep(beta0.smpls.hslogit,30)+tcrossprod(beta.smpls.hslogit,x)

f = beta0.mean+x%*%beta.mean;
pred.prob <- exp(f)/(1+exp(f))

par(mfrow=c(1,1))
plot(pred.prob,col=1+y)


####### Supervised Learning ####### 
#### Compare AUC #####


