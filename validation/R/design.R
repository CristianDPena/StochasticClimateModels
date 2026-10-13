regimes <- function() data.frame(regime=c('B','Aminus','Aplus','Jminus','Cminus','Sminus','H'),
 alpha=c(1.3,1.1,1.6,1.3,1.3,1.3,1.55),sigma2=c(1,1,1,.5,1,1,1.4),
 sigma1=c(1,1,1,1,1.875,.5,.7),rho=c(1.5,1.5,1.5,1.5,.8,1.5,2.2),reps=rep(24L,7))
make_design <- function(mode='study') {
 rows<-list();i<-0L
 add<-function(r,rep,T,n,component,batch,variant='ordinary',kernel='standard',transfer='standard',sim_level=4L){
  i<<-i+1L;rows[[i]]<<-data.frame(job_id=sprintf('v%04d',i),regime=r$regime,replication=rep,
    trajectory_id=sprintf('%s_%03d',r$regime,rep),seed_v=100000L+match(r$regime,regimes()$regime)*1000L+rep,
    seed_l=200000L+match(r$regime,regimes()$regime)*1000L+rep,alpha=r$alpha,sigma2=r$sigma2,sigma1=r$sigma1,rho=r$rho,
    T=T,n=as.integer(n),dt=T/n,component=component,batch=batch,variant=variant,kernel=kernel,transfer=transfer,sim_level=sim_level)
 }
 rr<-regimes()
 for(a in seq_len(nrow(rr)))for(b in seq_len(rr$reps[a])) {
  r<-rr[a,];add(r,b,50,1e7,'main',paste0('recovery',ceiling(b/6)))
  if(r$regime %in% c('B','H') && b<=12) {
   batch<-paste0('information',ceiling(b/6))
   add(r,b,50,2.5e6,'infill',batch);add(r,b,100,2e7,'horizon',batch)
  }
  if(r$regime %in% c('B','Jminus','H') && b<=4)
   for(v in c('oracle_tail','oracle_tail_beta'))add(r,b,50,1e7,'oracle','diagnostics',v)
  if(r$regime %in% c('B','H') && b<=2) {
   add(r,b,50,1e7,'numerical','diagnostics',kernel='refined')
   add(r,b,50,1e7,'simulator','diagnostics',sim_level=8L)
   for(t in c('refined','independent'))add(r,b,50,1e7,'transfer','diagnostics',transfer=t)
  }
 }
 d<-do.call(rbind,rows)
 stopifnot(sum(d$component=='main')==168L,sum(d$component%in%c('infill','horizon'))==48L,nrow(d)==256L)
 if(mode=='smoke') {
  d<-d[1,,drop=FALSE];d$job_id<-'smoke';d$trajectory_id<-'smoke';d$T<-1;d$n<-10000L;d$dt<-d$T/d$n;d$component<-'smoke';d$batch<-'smoke'
  oracle<-d;oracle$job_id<-'smoke-oracle';oracle$component<-'smoke_oracle';oracle$variant<-'oracle_tail';d<-rbind(d,oracle)
 }
 d
}
framework_files <- function() c(list.files('validation/R',full.names=TRUE),'validation/src/simulator.c','validation/run.R','validation/config/frozen_sources.csv')
initialize <- function(out,mode='study') {
 ensure(out);d<-make_design(mode)
 if(file.exists(file.path(out,'design.csv'))) {
  if(!isTRUE(all.equal(read_design(out),d,check.attributes=FALSE)))stop('Output design differs; use a new output directory')
  p<-readRDS(file.path(out,'provenance.rds'))
  if(!identical(p$framework_hashes,tools::md5sum(framework_files())))stop('Framework changed since initialization; use a new output directory')
  if(sha(file.path(out,'design.csv'))!=p$design_md5)stop('Design integrity mismatch')
  return(invisible(d))
 }
 write.csv(d,file.path(out,'design.csv'),row.names=FALSE)
 write.csv(regimes(),file.path(out,'regimes.csv'),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
 writeLines(system('git rev-parse HEAD',intern=TRUE),file.path(out,'git_commit.txt'))
 for(dir in c('raw','results','cache','logs','summary','figures','locks','checkpoints'))ensure(file.path(out,dir))
 file.copy('validation/config/frozen_sources.csv',file.path(out,'frozen_sources.csv'))
 atomic_rds(list(mode=mode,design_md5=sha(file.path(out,'design.csv')),created=format(Sys.time(),tz='UTC'),
  framework_hashes=tools::md5sum(framework_files()),R=R.version.string,
  evir=if(requireNamespace('evir',quietly=TRUE))as.character(packageVersion('evir')) else NA_character_),file.path(out,'provenance.rds'))
 invisible(d)
}
