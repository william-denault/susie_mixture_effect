#simulation_par
library(gtools)
library(MixSER)
#simulation per scenario
n_sim =300
n=100
MAF=  .3
sds=.25
beta= 1
generate_data= function (n,
                         mixture_prop= c(1/3,
                                         1/3,
                                         1/3),
                         beta=1,
                         sds=.25,
                         MAF=.3){


  if(length(beta==1)){
    beta= rep(beta,3)
  }
  g  <- rbinom(n, 2, 0.3); while (length(unique(g)) < 3) g <- rbinom(n, 2, 0.45)
  rc <- recode_snp(g)
  z= sample(c("A","D","R"),
            prob= mixture_prop,
            size=n,
            replace=TRUE)
  y=   rnorm(n,sd  = sds)
  y_true=   0*y

  for ( i in 1:n){

    if( z[i]=="A"){
      y_true[i]= rc$additive[i]*beta[1]
      y[i]=y[i]+rc$additive[i]*beta[1]
    }
    if( z[i]=="D"){
      y_true[i]= rc$dominant[i]*beta[1]
      y[i]=y[i]+rc$dominant[i]*beta[2]
    }
    if( z[i]=="R"){
      y_true[i]=  rc$recessive[i]*beta[3]
      y[i]=y[i]+ rc$recessive[i]*beta[3]
    }
  }

  return(list(y=y,
              y_true=y_true,
              z=z,
              rc= rc,
              beta=beta,
              mixture_prop=mixture_prop,
              sds=sds,
              MAF=MAF))

}

## A genuine per-individual mixture at ONE SNP: a fraction pi_R of individuals
## respond recessively, the rest additively.
res=list()
rmse= function(x,y){ sqrt(mean((x-y)^2))}

m=1

for ( o in 1: n_sim ){
    mix_pi <- c(rdirichlet(n = 1, alpha = rep(1, 3)))

    dat=  generate_data  (n=n,
                          mixture_prop= mix_pi,
                          beta=  beta,
                          sds=sds,
                          MAF=  MAF)



    fit <- ser_pim_em(
      codings = list(additive = dat$rc$additive,
                     dominant  = dat$rc$dominant,
                     recessive = dat$rc$recessive),
      y = dat$y, residual_variance = dat$sds^2,
      n_start = 10,
      pi_init = c(.9,.05,.05))

    susie_res= susie_res <- susieR::susie(
      X = matrix(as.double(dat$rc$additive), ncol = 1),
      y = dat$y
    )




    res[[m]]=list(fit=fit,
                  dat=dat,
                  rmse_sample= rmse(dat$y, dat$ y_true),
                  rmse_mix=rmse(fit$fitted, dat$y),
                  rmse_susie= rmse(susie_res$fitted, dat$y)
    )


    m=m+1
    print(m)
}




save(res,file =paste0(getwd(),"/sim/mixture_model_effect.RData"))



