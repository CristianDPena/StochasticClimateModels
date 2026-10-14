# Both benchmarks use the SAME residual data, no latent path or true parameter.
# Neither estimates beta/rho/sigma1. No iid model-based CIs are used.
benchmark_tail <- function(R,dt,method) {
 start<-proc.time()[3];fit<-NULL;est<-c(alpha=NA_real_,sigma2=NA_real_);warnings<-character();status<-'success';error<-''
 tryCatch(withCallingHandlers({
  A<-sort(abs(R[R!=0&is.finite(R)]),decreasing=TRUE);n<-length(A)
  if(method=='hill_sqrt') {
   k<-floor(sqrt(n));alpha<-1/mean(log(A[seq_len(k)]/A[k+1L]));logconstant<-log(k/n)+alpha*log(A[k+1L])
   fit<-list(k=k,threshold=A[k+1L],rule='floor(sqrt(n)); no threshold search')
  } else if(method=='gpd_ml_99') {
   if(!requireNamespace('evir',quietly=TRUE))stop('Required benchmark dependency evir is missing')
   threshold<-unname(quantile(A,.99));normalizer<-median(A);fit<-evir::gpd(A/normalizer,threshold=threshold/normalizer,method='ml')
   if(!is.null(fit$converged)&&fit$converged!=0)stop('GPD optimizer did not converge')
   xi<-unname(fit$par.ests['xi']);scale<-unname(fit$par.ests['beta'])*normalizer
   alpha<-1/xi
   logconstant<-log(mean(A>threshold))+alpha*log(scale/xi)
  } else stop('Unknown benchmark')
  est['alpha']<-alpha
  if(!is.finite(alpha)||alpha<=0||alpha>=2) {status<-'outside_stable_domain'} else {
   constant<-gamma(alpha)*sin(pi*alpha/2)/pi
   est['sigma2']<-exp((logconstant-log(2*constant*dt))/alpha)
   if(!is.finite(est['sigma2'])||est['sigma2']<=0)status<-'invalid_scale'
  }
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){status<<-'error';error<<-conditionMessage(e)})
 list(estimates=est,status=status,error=error,warnings=warnings,fit=fit,seconds=unname(proc.time()[3]-start))
}
