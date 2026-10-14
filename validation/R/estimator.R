kernel_config <- function(level) {
 k<-list(n_states=21L,per_state=30000L,substeps=32L,nodes=31L,rho_grid=seq(.7,2.3,.1),seed=5000L)
 if(level=='refined')k$substeps<-64L
 if(level=='seed2')k$seed<-5001L
 if(level=='seed3')k$seed<-5002L
 k
}
get_kernel <- function(level,out) {
 k<-kernel_config(level)
 if(level=='standard') {
  c<-readRDS('original_source/outputs/cache/conditional_cf_endpoint_kernel_production.rds')
  stopifnot(c$config$n_states==21,c$config$substeps==32,c$config$max_path_nodes==31,
   c$config$transition_sim_blocks==630000,c$config$kernel_seed==5000,
   abs(c$config$block_horizon-.05)<1e-12,isTRUE(all.equal(c$config$rho_grid,k$rho_grid)))
  return(c$models)
 }
 p<-file.path(out,'cache',paste0('kernel_',level,'.rds'))
 if(file.exists(p)){x<-readRDS(p);stopifnot(identical(x$config,k));return(x$models)}
 state<-make_q_endpoint_breaks(k$n_states)
 models<-lapply(k$rho_grid,function(r)simulate_q_endpoint_cf_model(r,.05,state,k$per_state*k$n_states,
  substeps_per_block=k$substeps,max_path_nodes_per_pair=k$nodes,state_model='conditional',seed=k$seed))
 atomic_rds(list(config=k,models=models),p);models
}
get_transfer <- function(row,out) {
 q<-1/log(exp(1)+row$n);h<-floor(sqrt(row$n))*row$dt
 label<-sprintf('transfer_n%d_T%s_%s',row$n,row$T,row$transfer);p<-file.path(out,'cache',paste0(label,'.rds'))
 if(file.exists(p)) {x<-readRDS(p);stopifnot(all(abs(x$q-q)<1e-12),all(abs(x$h-h)<1e-12));return(x)}
 if(row$transfer=='standard') {
  candidates<-c('original_source/finite_window_correction_numerical_profiles.csv',list.files('original_source/overarching_profiles',full.names=TRUE))
  for(f in candidates) {
   x<-read.csv(f);if('replication'%in%names(x))x<-x[x$replication=='primary',]
   if(nrow(x)&&all(abs(x$q-q)<1e-12)&&all(abs(x$h-h)<1e-12)){atomic_rds(x,p);return(x)}
  }
 }
 refined<-row$transfer=='refined';finegrid<-row$transfer=='finegrid'
 x<-rho_fw_build_profile(q,h,rho_grid=seq(.4,3.5,if(finegrid).05 else .1),
  pair_count=if(refined)100000L else 50000L,fine_steps=if(refined)512L else 256L,
  resolution_divisors=c(4L,2L,1L),seed=if(row$transfer=='independent')941002L else 941001L,chunk_size=1000L)
 # Invalid/nonmonotone numerical maps are retained and reported, never repaired.
 atomic_rds(x,p);x
}
fit_fixed <- function(obs,row,models,profile,out,reference=NULL) {
 estimate<-setNames(rep(NA_real_,5),c('alpha','sigma2','beta','rho','sigma1'))
 stages<-list();stage<-'tail';warnings<-character();start<-proc.time()[3]
 status<-'success';error<-''
 tryCatch(withCallingHandlers({
  if(row$variant %in% c('oracle_tail','oracle_tail_beta')) {
   tail<-list(alpha_hat=row$alpha,sigma2_hat=row$sigma2,information='oracle')
  } else tail<-estimate_alpha_sigma2_hill(obs$R,rep(obs$dt,length(obs$R)))
  stages$tail<-tail
  if(row$variant=='oracle_alpha') {
   tail$alpha_hat<-row$alpha
   tail$sigma2_hat<-estimate_sigma2_tail_robust(row$alpha,sort(abs(obs$R),decreasing=TRUE),tail$k_hat,rep(obs$dt,length(obs$R)))$sigma2_hat
  }
  estimate[c('alpha','sigma2')]<-c(tail$alpha_hat,tail$sigma2_hat)
  stage<-'beta'
  if(row$variant=='oracle_tail_beta') estimate['beta']<-row$sigma1*row$rho else if(row$component=='transfer') {
   if(is.null(reference))stop('Transfer diagnostics require completed ordinary fit')
   stages$endpoint<-reference$fit$stages$endpoint
   estimate['beta']<-reference$fit$estimates['beta']
  } else {
   k<-kernel_config(row$kernel)
   cp<-file.path(out,'checkpoints',row$job_id);ensure(cp)
   endpoint<-estimate_exact_cf_endpoint_hmm(obs$R,rep(obs$dt,length(obs$R)),tail$alpha_hat,tail$sigma2_hat,
    rho_grid=k$rho_grid,kernel_models=models,covariance_model='conditional_diag',
    n_states=k$n_states,transition_sim_blocks=k$n_states*k$per_state,substeps_per_block=k$substeps,
    max_path_nodes_per_pair=k$nodes,kernel_seed=k$seed,beta_profile_grid_size=25L,
    checkpoint_dir=cp,checkpoint_signature=paste(row$job_id,sha(file.path(out,'design.csv')),row$kernel,row$variant,sep=':'),verbose=FALSE)
   stages$endpoint<-endpoint;estimate['beta']<-endpoint$beta_hat
  }
  stage<-'rho'
  u<-1/(estimate['beta']*sqrt(log(exp(1)+length(obs$R))))
  local<-rho_level4_estimate(obs$R,obs$dt,tail$alpha_hat,tail$sigma2_hat,unname(estimate['beta']),u=unname(u))
  stages$local<-local
  corrected<-rho_fw_correct_level4(local,profile);stages$corrected<-corrected
  estimate[c('rho','sigma1')]<-c(corrected$rho_hat,corrected$sigma1_hat)
  if(corrected$status!='success')status<-paste0('rho:',corrected$status)
  if(any(!is.finite(estimate))&&status=='success')status<-'nonfinite'
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){status<<-paste0(stage,':error');error<<-conditionMessage(e)})
 list(estimates=estimate,status=status,error=error,warnings=warnings,stages=stages,seconds=unname(proc.time()[3]-start))
}
