load("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mixture_effect/sim/mixture_model_effect.RData")



pi_est= do.call(rbind,
                lapply(1: length(res), function(k) res[[k]]$fit$pi)
)


true_pi=  do.call(rbind,
                  lapply(1: length(res), function(k) res[[k]]$dat$mixture_prop)
)

ind_est= do.call(rbind,
                 lapply(1: length(res), function(k) res[[k]]$fit$r)
)


true_ind=  do.call(c,
                   lapply(1: length(res), function(k) res[[k]]$dat$z)
)

rmse_res=  do.call(rbind,
                   lapply(1: length(res), function(k)  data.frame(rmse= c(res[[k]]$rmse_sample,
                                                                          res[[k]]$rmse_mix,
                                                                          res[[k]]$rmse_susie),
                                                                  type=c("in sample",
                                                                         "mix SER",
                                                                         "SER"),
                                                                  reg_type= rep(res[[k]]$dat$mixture_prop[1],3)
                   )
                   ))

ggplot(rmse_res, aes(y =rmse, x=type,))+geom_boxplot()
ggplot(rmse_res, aes(y =rmse,x = reg_type,
                     color = type,
                     shape=type))+geom_point()+
  xlab("proportion of additive individual")



true_ind_mat=0*ind_est
true_ind_mat[which(true_ind=="A"),1]=1
true_ind_mat[which(true_ind=="D"),2]=1
true_ind_mat[which(true_ind=="R"),3]=1



plot(c(pi_est[,1]),c(true_pi[,1]), xlab="proportion ind additive regulatio" )

plot(c(pi_est[,2]),c(true_pi[,2]) )

plot(c(pi_est[,3]),c(true_pi[,3]) , main="Estiat")



plot(c(pi_est[,2]+pi_est[,3]),c(true_pi[,2]+true_pi[,3]) )

library(pROC)
roc_curve <- roc(c(true_pi[,1]), c(pi_est[,1]))
plot(roc_curve)






roc_curve <- roc( (true_ind_mat[,1]),  (ind_est[,1]))
plot(roc_curve, main="classification of individual being under additive regulation")
roc_curve <- roc( (true_ind_mat[,2])+(true_ind_mat[,3]),  (ind_est[,2]+ind_est[,3]) )
plot(roc_curve, main="classification of individual being under non additive regulation")
roc_curve <- roc( (true_ind_mat[,3]),  (ind_est[,3]) )
plot(roc_curve, main="classification of individual being under recessive regulation")
roc_curve <- roc( (true_ind_mat[,2]),  (ind_est[,2]) )
plot(roc_curve, main="classification of individual being under dominant regulation")

