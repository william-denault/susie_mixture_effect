#simulation_par

library(MixSER)
#simulation per scenario
n_sim =100
n=100
MAF=  .3
sds=.25
beta=  beta
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

  for ( i in 1:n){

    if( z[i]=="A"){
      y[i]=y[i]+rc$additive[i]*beta[1]
    }
    if( z[i]=="D"){
      y[i]=y[i]+rc$dominant[i]*beta[2]
    }
    if( z[i]=="R"){
      y[i]=y[i]+rc$dominant[i]*beta[3]
    }
  }

  return(list(y=y,
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

m=1
for ( k in 1:3){

  if (k ==1){
    mix_pi=c(1,0,0)
  }
  if (k ==2){
    mix_pi=c(0,1,0)
  }
  if (k ==3){
    mix_pi=c(0,0,1)
  }
  mix_pi
 for ( o in 1: n_sim ){


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
    pi_init = c(.8,.1,.1))



  res[[m]]=list(fit=fit,
                dat=dat
                )
  m=m+1
  print(m)
 }


}


save(res,file =paste0(getwd(),"/sim/single_model_effect.RData"))



