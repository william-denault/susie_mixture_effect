load("C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mixture_effect/sim/single_model_effect.RData")



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

true_ind_mat=0*ind_est
true_ind_mat[which(true_ind=="A"),1]=1
true_ind_mat[which(true_ind=="D"),2]=1
true_ind_mat[which(true_ind=="R"),3]=1



plot(c(pi_est[,1]),c(true_pi[,1]) )



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
