metric_row <- function(x) {
 good<-is.finite(x$estimate);e<-x$estimate[good]-x$truth[good];m<-length(e);N<-nrow(x)
 relative<-if(x$parameter[1]=='alpha')e else e/x$truth[good]
 success<-good & x$status=='success' & if(x$parameter[1]=='alpha')abs(x$estimate-x$truth)<.1 else abs(x$estimate/x$truth-1)<.1
 p<-mean(success);z<-qnorm(.975);den<-1+z*z/N;mid<-(p+z*z/(2*N))/den;half<-z*sqrt(p*(1-p)/N+z*z/(4*N*N))/den
 data.frame(N=N,finite=m,bias=if(m)mean(e) else NA,sd=if(m>1)sd(e) else NA,
  rmse=if(m)sqrt(mean(e^2)) else NA,mae=if(m)mean(abs(e)) else NA,median_abs=if(m)median(abs(e)) else NA,
  q90_abs=if(m)unname(quantile(abs(e),.9)) else NA,relative_rmse=if(m)sqrt(mean(relative^2)) else NA,
  bias_mcse=if(m>1)sd(e)/sqrt(m) else NA,
  rmse_mcse=if(m>1&&sum(e^2)>0)sd(e^2)/(2*sqrt(mean(e^2))*sqrt(m)) else NA,
  tolerance_fraction=p,tolerance_low=mid-half,tolerance_high=mid+half)
}
summarize_study <- function(out) {
 d<-read_design(out);files<-list.files(file.path(out,'results'),pattern='[.]rds$',full.names=TRUE)
 objs<-lapply(files,readRDS);rows<-list();joint<-list()
 for(o in objs) {
  r<-o$design;fits<-c(list(fixed=o$fit),o$benchmarks)
  for(method in names(fits)) {
   f<-fits[[method]];truth<-c(alpha=r$alpha,sigma2=r$sigma2,beta=r$sigma1*r$rho,rho=r$rho,sigma1=r$sigma1)
   for(p in names(f$estimates))rows[[length(rows)+1L]]<-cbind(r,data.frame(method=method,parameter=p,truth=unname(truth[p]),estimate=unname(f$estimates[p]),status=f$status,seconds=f$seconds,max_abs_V=o$max_abs_V))
  }
  e<-o$fit$estimates;finite<-all(is.finite(e));ok<-finite&&o$fit$status=='success'
  joint[[length(joint)+1L]]<-cbind(r,data.frame(status=o$fit$status,error=o$fit$error,finite_complete=finite,
   joint_success=ok&&abs(e['alpha']-r$alpha)<.1&&all(abs(e[c('sigma2','rho','sigma1')]/c(r$sigma2,r$rho,r$sigma1)-1)<.1)))
 }
 long<-bind_rows(rows);j<-bind_rows(joint);ensure(file.path(out,'summary'));ensure(file.path(out,'figures'))
 write.csv(long,file.path(out,'summary','estimates.csv'),row.names=FALSE);write.csv(j,file.path(out,'summary','failures_and_joint_recovery.csv'),row.names=FALSE)
 done<-if(nrow(j))j$job_id else character();missing<-d[!d$job_id%in%done,];write.csv(missing,file.path(out,'summary','pending.csv'),row.names=FALSE)
 report<-c('# Validation report','',sprintf('Completed scheduled jobs: %d / %d.',sum(done%in%d$job_id),nrow(d)),
  sprintf('Principal-study jobs completed: %d / 216.',sum(done%in%d$job_id[d$component%in%c('main','stress','infill','horizon')])),
  '',if(nrow(missing))'**INCOMPLETE STUDY. These outputs do not establish robustness or comparative superiority.**' else '**All scheduled jobs have result records; inspect failures and numerical refinement before scientific claims.**',
  '', 'The fixed estimator sources are verified against the frozen checksums. Oracle and refinement results are separate from ordinary estimates. Undefined estimates remain in failure and tolerance denominators; error moments use finite outputs.',
  '', 'Benchmark scope: fixed-rank Hill and evir generalized-Pareto maximum likelihood estimate only the stable tail. Their jump-scale estimates use stable-tail conversion. They do not estimate the latent amplitude or clock. No full-pipeline superiority claim is supported by these partial comparisons.')
 if(nrow(long)) {
  keys<-c('component','regime','T','n','variant','kernel','transfer','sim_level','method','parameter')
  groups<-split(seq_len(nrow(long)),do.call(paste,c(long[keys],sep='|')))
  summary<-do.call(rbind,lapply(groups,function(ix)cbind(long[ix[1],keys,drop=FALSE],metric_row(long[ix,]))))
  write.csv(summary,file.path(out,'summary','metrics.csv'),row.names=FALSE)
  main<-summary[summary$component%in%c('main','stress','smoke')&summary$method=='fixed',]
  tex<-c('\\begin{tabular}{llrrrrr}','\\toprule','Regime & Parameter & $N$ & Finite & Bias & RMSE & Median AE \\\\','\\midrule')
  for(i in seq_len(nrow(main)))tex<-c(tex,sprintf('%s & %s & %d & %d & %.4f & %.4f & %.4f \\\\',main$regime[i],main$parameter[i],main$N[i],main$finite[i],main$bias[i],main$rmse[i],main$median_abs[i]))
  writeLines(c(tex,'\\bottomrule','\\end{tabular}'),file.path(out,'summary','recovery_table.tex'))
  # Paired contrasts against the ordinary reference of the SAME trajectory and parameter.
  ref<-long[long$method=='fixed'&long$component%in%c('main','stress','smoke'),c('trajectory_id','parameter','estimate','truth','status')]
  names(ref)[3:5]<-c('reference','reference_truth','reference_status')
  pair<-merge(long,ref,by=c('trajectory_id','parameter'),all.x=TRUE)
  # Competitors must be paired to the fixed estimator on the identical observed dataset.
  fixed_same<-long[long$method=='fixed',c('job_id','parameter','estimate','truth','status')]
  for(ii in which(pair$method!='fixed')) {
   hit<-which(fixed_same$job_id==pair$job_id[ii]&fixed_same$parameter==pair$parameter[ii])
   if(length(hit)==1L){pair$reference[ii]<-fixed_same$estimate[hit];pair$reference_truth[ii]<-fixed_same$truth[hit];pair$reference_status[ii]<-fixed_same$status[hit]}
  }
  pair$reference_kind<-ifelse(pair$method=='fixed','same_trajectory_reference_design','same_observed_dataset')
  pair$paired_finite<-is.finite(pair$estimate)&is.finite(pair$reference)
  pair$absolute_error_difference<-abs(pair$estimate-pair$truth)-abs(pair$reference-pair$reference_truth)
  pair$squared_error_difference<-(pair$estimate-pair$truth)^2-(pair$reference-pair$reference_truth)^2
  pair$estimate_shift<-pair$estimate-pair$reference
  write.csv(pair,file.path(out,'summary','paired_contrasts.csv'),row.names=FALSE)
  pg<-split(seq_len(nrow(pair)),do.call(paste,c(pair[keys],sep='|')))
  set.seed(773001)
  ps<-bind_rows(lapply(pg,function(ix){x<-pair[ix,];v<-x$absolute_error_difference[x$paired_finite];m<-length(v)
   interval<-if(m>=2)quantile(replicate(999,mean(sample(v,m,replace=TRUE))),c(.025,.975)) else c(NA,NA)
   cbind(x[1,keys,drop=FALSE],data.frame(N=nrow(x),paired_finite=m,mean_AE_difference=if(m)mean(v) else NA,
    mcse=if(m>1)sd(v)/sqrt(m) else NA,bootstrap_low=interval[1],bootstrap_high=interval[2]))}))
  write.csv(ps,file.path(out,'summary','paired_summary.csv'),row.names=FALSE)
  params<-c('alpha','sigma2','sigma1','rho');col<-c('#0072B2','#D55E00','#009E73','#CC79A7')
  pdf(file.path(out,'figures','recovery.pdf'),width=10,height=7,useDingbats=FALSE);par(mfrow=c(2,2),mar=c(5,4,2,1))
  for(p in params) {
   x<-long[long$component%in%c('main','stress','smoke')&long$method=='fixed'&long$parameter==p&is.finite(long$estimate),]
   if(nrow(x)){err<-if(p=='alpha')x$estimate-x$truth else x$estimate/x$truth-1
    boxplot(err~x$regime,outline=FALSE,las=2,col='#DCEAF2',ylab=if(p=='alpha')'Signed error' else 'Relative signed error',main=p)
    stripchart(err~x$regime,vertical=TRUE,add=TRUE,method='jitter',pch=16,cex=.5,col='#0072B280');abline(h=0,lty=2)}else {plot.new();title(p)}
  };dev.off()
  growth<-long[long$method=='fixed'&long$regime%in%c('B','H')&long$replication<=12&long$component%in%c('main','infill','horizon'),]
  gg<-split(seq_len(nrow(growth)),do.call(paste,c(growth[c('regime','parameter','T','n','component')],sep='|')))
  gm<-bind_rows(lapply(gg,function(ix)cbind(growth[ix[1],c('regime','parameter','T','n','component')],metric_row(growth[ix,]))))
  write.csv(gm,file.path(out,'summary','information_growth.csv'),row.names=FALSE)
  pdf(file.path(out,'figures','information_growth.pdf'),width=10,height=11,useDingbats=FALSE);par(mfrow=c(4,2),mar=c(4,4,2,1))
  for(p in params)for(comp in c('infill','horizon')) {
   x<-gm[gm$parameter==p&gm$component%in%c('main',comp)&is.finite(gm$rmse),]
   axisvar<-if(comp=='infill')'n' else 'T'
   if(nrow(x)) {
    xx<-x[[axisvar]]/if(comp=='infill')1e6 else 1
    plot(range(xx)+c(-.01,.01),range(c(0,x$rmse+1.96*x$rmse_mcse),na.rm=TRUE)+c(0,1e-6),type='n',
     xlab=if(comp=='infill')'n (millions), T=50' else 'T, spacing=0.000005',ylab='RMSE',main=p)
    for(a in c('B','H')){y<-x[x$regime==a,];y<-y[order(y[[axisvar]]),];if(nrow(y)) {
     xx<-y[[axisvar]]/if(comp=='infill')1e6 else 1;color<-col[match(a,c('B','H'))]
     lines(xx,y$rmse,type='b',pch=16,col=color)
     good<-is.finite(y$rmse_mcse);arrows(xx[good],pmax(0,y$rmse[good]-1.96*y$rmse_mcse[good]),xx[good],y$rmse[good]+1.96*y$rmse_mcse[good],angle=90,code=3,length=.03,col=color)
    }}
    legend('topright',legend=c('B','H'),col=col[1:2],pch=16,bty='n',cex=.8)
   }else{plot.new();title(p)}
  };dev.off()
  pdf(file.path(out,'figures','benchmarks.pdf'),width=10,height=5,useDingbats=FALSE);par(mfrow=c(1,2),mar=c(8,4,2,1))
  for(p in c('alpha','sigma2')){x<-long[long$component%in%c('main','stress','smoke')&long$parameter==p&is.finite(long$estimate),]
   if(nrow(x))boxplot(abs(x$estimate-x$truth)~interaction(x$regime,x$method),las=2,cex.axis=.55,col='#DCEAF2',main=p,ylab='Absolute error') else plot.new()};dev.off()
  if(nrow(j))report<-c(report,'',sprintf('Recorded estimator statuses: %s.',paste(names(table(j$status)),as.integer(table(j$status)),sep='=',collapse='; ')))
  report<-c(report,'','## Manuscript selection','Use ordinary main-regime recovery distributions and paired information-growth results as primary evidence. Include numerical-stability results only after both refinements complete. Keep oracle comparisons and partial tail benchmarks explicitly labeled; report materially adverse outcomes in the main text. Do not pool stress regimes into a single headline success rate.','', '## Reproducibility','Raw trajectory RDS files, integrity sidecars, per-job RDS objects (including profiles and warnings), the immutable design, session information, and source manifests are retained under this run directory. CSV tables and PDF figures are generated solely from these records.')
 }
 if(any(d$component=='smoke'))report<-c('# COMPUTATIONAL SMOKE TEST — NOT SCIENTIFIC VALIDATION','',report)
 writeLines(report,file.path(out,'summary','validation_report.md'))
}
