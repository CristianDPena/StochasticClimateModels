options(stringsAsFactors=FALSE)
validation_root <- normalizePath('.', mustWork=TRUE)
if (!dir.exists('original_source')) stop('Run from repository root')
lib <- file.path(validation_root,'validation','.library')
if(dir.exists(lib)) .libPaths(c(lib,.libPaths()))
args <- commandArgs(TRUE)
arg <- function(key, default=NULL) {
  a <- grep(paste0('^--',key,'='), args, value=TRUE)
  if(length(a)) sub(paste0('^--',key,'='),'',a[1]) else default
}
ensure <- function(path) dir.create(path,recursive=TRUE,showWarnings=FALSE)
atomic_rds <- function(x,path) {
  ensure(dirname(path)); tmp<-paste0(path,'.',Sys.getpid(),'.tmp')
  saveRDS(x,tmp,compress=FALSE); if(!file.rename(tmp,path))stop('Atomic write failed: ',path)
}
sha <- function(path) unname(tools::md5sum(path)) # integrity, not cryptographic authentication
bind_rows <- function(xs) {
  xs<-Filter(function(x)!is.null(x)&&nrow(x)>0,xs);if(!length(xs))return(data.frame())
  cols<-unique(unlist(lapply(xs,names)));do.call(rbind,lapply(xs,function(x){for(k in setdiff(cols,names(x)))x[[k]]<-NA;x[cols]}))
}
load_frozen <- function() {
  old<-getwd();on.exit(setwd(old));setwd(file.path(validation_root,'original_source'))
  source('work/diagnose_exact_cf_endpoint_hmm.R',local=.GlobalEnv)
  assign('unchunked_forward_reference',conditional_exact_cf_forward_loglik,.GlobalEnv)
  source('work/chunked_exact_endpoint_likelihood.R',local=.GlobalEnv)
  source('rho_level4_oracle_residual_estimator.R',local=.GlobalEnv)
  source('finite_window_corrected_rho_estimator.R',local=.GlobalEnv)
  assign('conditional_exact_cf_forward_loglik',function(y_mat,model,beta,freqs,alpha,sigma2,dt_reference,block_size,full_covariance=FALSE,return_filter=FALSE){
    if(full_covariance)stop('Frozen covariance is diagonal')
    conditional_exact_cf_forward_loglik_chunked_library(y_mat,model,beta,freqs,alpha,sigma2,dt_reference,block_size,groups_per_chunk=16L,return_filter=return_filter)
  },.GlobalEnv)
}
verify_frozen <- function() {
  m<-read.csv('validation/config/frozen_sources.csv'); stopifnot(all(file.exists(m$path)))
  if(!identical(as.character(m$md5),as.character(unname(tools::md5sum(m$path)))))stop('Frozen source/cache changed')
}
read_design <- function(out) read.csv(file.path(out,'design.csv'))
