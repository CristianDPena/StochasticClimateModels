#!/usr/bin/env Rscript
# Post-run integrity audit; never modifies simulation or estimator records.
source('validation/R/common.R')
out<-arg('out','validation/results/laptop-v1');verify_frozen()
d<-read_design(out);stopifnot(nrow(d)==256L,!anyDuplicated(d$job_id))
p<-readRDS(file.path(out,'provenance.rds'))
stopifnot(p$mode=='study',sha(file.path(out,'design.csv'))==p$design_md5,
 identical(unname(p$framework_hashes),unname(tools::md5sum(names(p$framework_hashes)))))
files<-file.path(out,'results',paste0(d$job_id,'.rds'));stopifnot(all(file.exists(files)))
stopifnot(setequal(basename(files),list.files(file.path(out,'results'),pattern='[.]rds$')))
raw_ids<-unique(d$trajectory_id);raw_files<-file.path(out,'raw',paste0(raw_ids,'.rds'))
stopifnot(length(raw_ids)==168L,all(file.exists(raw_files)))
raw_hashes<-unname(tools::md5sum(raw_files));names(raw_hashes)<-raw_ids
for(i in seq_along(raw_files))stopifnot(identical(readLines(paste0(raw_files[i],'.md5')),unname(raw_hashes[i])))
records<-lapply(files,readRDS)
for(i in seq_along(records)){
 o<-records[[i]];r<-d[i,,drop=FALSE]
 stopifnot(isTRUE(all.equal(o$design,r,check.attributes=FALSE)),
 o$data_md5==raw_hashes[[r$trajectory_id]],
 setequal(names(o$fit$estimates),c('alpha','sigma2','beta','rho','sigma1')),
 length(o$fit$status)==1L,nzchar(o$fit$status))
 if(r$component %in% c('main','infill','horizon'))stopifnot(setequal(names(o$benchmarks),c('hill_sqrt','gpd_ml_99')))
}
s<-read.csv(file.path(out,'summary','failures_and_joint_recovery.csv'))
stopifnot(nrow(s)==256L,setequal(s$job_id,d$job_id),nrow(read.csv(file.path(out,'summary','pending.csv')))==0L)
for(i in seq_along(records))stopifnot(s$status[match(d$job_id[i],s$job_id)]==records[[i]]$fit$status)
ensure(file.path(out,'audit'))
manifest<-data.frame(path=c(file.path('raw',basename(raw_files)),file.path('results',basename(files))),
 md5=c(unname(raw_hashes),unname(tools::md5sum(files))))
write.csv(manifest,file.path(out,'audit','artifact_integrity.csv'),row.names=FALSE)
writeLines(c('# Completed-study integrity audit','',paste('Verified:',format(Sys.time(),tz='UTC',usetz=TRUE)),
 '256 unique scheduled result records; 168 raw trajectory files and checksum sidecars verified.',
 'Frozen estimator sources, framework hashes, design identity, per-result design and raw-data links verified.',
 'All 216 principal fits contain both declared benchmark records; summary statuses match result records.',
 'No pending jobs. Statistical failures are retained and are not audit failures.',
 'This audit checks reproducible records, not theoretical validity or scientific accuracy.'),file.path(out,'audit','README.md'))
cat('PASS: complete design, frozen framework, 168 raw checksums, 256 result records and summary statuses.\n')
