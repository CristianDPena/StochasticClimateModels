compile_simulator <- function(out) {
 ensure(out);out<-normalizePath(out);dll<-file.path(out,paste0('simulator',.Platform$dynlib.ext))
 if(!file.exists(dll)) {
  file.copy('validation/src/simulator.c',file.path(out,'simulator.c'),overwrite=TRUE)
  old<-getwd();on.exit(setwd(old));setwd(out)
  status<-system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB','simulator.c','-o',paste0('simulator',.Platform$dynlib.ext)))
  if(status!=0)stop('C simulator compilation failed')
 }
 dyn.load(dll);invisible(dll)
}
generate_trajectory <- function(row,d,out) {
 path<-file.path(out,'raw',paste0(row$trajectory_id,'.rds'))
 if(file.exists(path)){
  if(!file.exists(paste0(path,'.md5'))||sha(path)!=readLines(paste0(path,'.md5')))stop('Raw-data checksum mismatch')
  return(path)
 }
 lock<-file.path(out,'locks',paste0('trajectory_',row$trajectory_id,'.lock'))
 if(!dir.create(lock,showWarnings=FALSE))stop('Trajectory locked: ',row$trajectory_id)
 writeLines(as.character(Sys.getpid()),file.path(lock,'pid'));on.exit(unlink(lock,recursive=TRUE))
 peers<-d[d$trajectory_id==row$trajectory_id,];delta<-min(peers$dt);T<-max(peers$T);n<-as.integer(round(T/delta))
 levels<-as.integer(sort(unique(peers$sim_level)))
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(row$seed_v)
 V<-.Call('coupled_latent',n,delta,row$rho,levels)
 if(any(!is.finite(V)))stop('Simulator produced nonfinite latent path; no resampling allowed')
 set.seed(row$seed_l);L<-row$sigma2*delta^(1/row$alpha)*symmetric_stable_r(n,row$alpha)
 residuals<-setNames(lapply(seq_along(levels),function(l)row$sigma1*diff(V[,l])+L),as.character(levels))
 horizons<-sort(unique(peers$T))
 diagnostics<-do.call(rbind,lapply(seq_along(levels),function(l)data.frame(level=levels[l],T=horizons,
  max_abs_V=vapply(horizons,function(t)max(abs(V[seq_len(round(t/delta)+1L),l])),numeric(1)))))
 atomic_rds(list(residuals=residuals,diagnostics=diagnostics,delta=delta,T=T,truth=row[c('alpha','sigma2','sigma1','rho')],
  seeds=row[c('seed_v','seed_l')],generator='Euler coupled 16 Brownian increments per observation; V0=0',
  source_md5=sha('validation/src/simulator.c')),path)
 writeLines(sha(path),paste0(path,'.md5'));path
}
observe <- function(raw,row) {
 stopifnot(isTRUE(all.equal(unname(unlist(raw$truth)),unname(unlist(row[c('alpha','sigma2','sigma1','rho')])))))
 stride<-as.integer(round(row$dt/raw$delta));stopifnot(abs(stride*raw$delta-row$dt)<1e-12)
 n<-as.integer(row$n);R<-raw$residuals[[as.character(row$sim_level)]]
 stopifnot(!is.null(R),n*stride<=length(R));R<-R[seq_len(n*stride)]
 if(stride>1L)R<-colSums(matrix(R,nrow=stride))
 diag<-raw$diagnostics[raw$diagnostics$level==row$sim_level & abs(raw$diagnostics$T-row$T)<1e-10,]
 stopifnot(nrow(diag)==1L)
 list(R=as.numeric(R),dt=row$dt,max_abs_V=diag$max_abs_V)
}
