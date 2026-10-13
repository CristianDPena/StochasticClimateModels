#!/usr/bin/env Rscript
source('validation/R/common.R')
for(f in c('design','simulate','estimator','benchmarks'))source(paste0('validation/R/',f,'.R'))
command<-arg('command','help');mode<-arg('mode',if(command=='smoke')'smoke' else 'study')
out<-arg('out','validation/results/laptop-v1');ensure(out);out<-normalizePath(out)
if(command=='help') {cat('Commands: init, prepare, run, summarize, smoke. Options: --out=PATH --job=v0001 --batch=recovery1 --mode=study|smoke\n');quit()}
verify_frozen();initialize(out,mode)
if(command=='init')quit()
if(command=='summarize'){source('validation/R/summarize.R');summarize_study(out);quit()}
load_frozen();compile_simulator(file.path(out,'cache'));d<-read_design(out)
selected<-d
if(!is.null(arg('job')))selected<-selected[selected$job_id%in%strsplit(arg('job'),',')[[1]],]
if(!is.null(arg('component')))selected<-selected[selected$component%in%strsplit(arg('component'),',')[[1]],]
if(!is.null(arg('batch')))selected<-selected[selected$batch%in%strsplit(arg('batch'),',')[[1]],]
if(!is.null(arg('limit')))selected<-head(selected,as.integer(arg('limit')))
if(!nrow(selected))stop('No jobs selected')
if(command=='prepare') {
 for(k in unique(selected$kernel))get_kernel(k,out)
 combos<-unique(selected[c('n','T','dt','transfer')])
 for(i in seq_len(nrow(combos)))get_transfer(combos[i,],out)
 quit()
}
if(!command%in%c('run','smoke'))stop('Unknown command')
failed<-FALSE
for(i in seq_len(nrow(selected))) {
 row<-selected[i,,drop=FALSE];target<-file.path(out,'results',paste0(row$job_id,'.rds'))
 if(file.exists(target))next
 lock<-file.path(out,'locks',paste0('job_',row$job_id,'.lock'))
 if(!dir.create(lock,showWarnings=FALSE))stop('Locked job: ',row$job_id)
 writeLines(as.character(Sys.getpid()),file.path(lock,'pid'))
 tryCatch({
  message('START ',row$job_id,' ',row$regime,' ',row$component)
  p<-generate_trajectory(row,d,out);raw<-readRDS(p);obs<-observe(raw,row);rm(raw);gc(FALSE)
  reference<-NULL
  if(row$component=='transfer') {
   ref<-d[d$trajectory_id==row$trajectory_id&d$component=='main',]
   rp<-file.path(out,'results',paste0(ref$job_id,'.rds'))
   if(length(rp)!=1||!file.exists(rp))stop('Complete ordinary reference first')
   reference<-readRDS(rp)
  }
  models<-if(row$variant=='oracle_tail_beta'||row$component=='transfer')NULL else get_kernel(row$kernel,out)
  profile<-get_transfer(row,out)
  fit<-fit_fixed(obs,row,models,profile,out,reference)
  benchmarks<-if(row$component%in%c('main','infill','horizon','smoke'))setNames(lapply(c('hill_sqrt','gpd_ml_99'),function(m)benchmark_tail(obs$R,obs$dt,m)),c('hill_sqrt','gpd_ml_99')) else list()
  result<-list(design=row,fit=fit,benchmarks=benchmarks,max_abs_V=obs$max_abs_V,
    data_md5=readLines(paste0(p,'.md5')),completed=format(Sys.time(),tz='UTC'),session=capture.output(sessionInfo()))
  atomic_rds(result,target);message('DONE ',row$job_id,' ',fit$status)
 },error=function(e){
  failed<<-TRUE
  writeLines(conditionMessage(e),file.path(out,'logs',paste0(row$job_id,'.error.txt')));message('ERROR ',row$job_id,': ',conditionMessage(e))
 },finally=unlink(lock,recursive=TRUE))
}
if(failed)quit(status=1L)
