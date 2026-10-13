source('validation/R/common.R')
for(f in c('design','simulate','estimator','benchmarks','summarize'))source(paste0('validation/R/',f,'.R'))
verify_frozen();d<-make_design();stopifnot(nrow(d)==256,all(table(d$batch)[paste0('recovery',1:4)]==42),!anyDuplicated(d$job_id))
stopifnot(length(unique(d$trajectory_id))==168,sum(d$component%in%c('main','infill','horizon'))==216)
for(id in unique(d$trajectory_id)) {
 x<-d[d$trajectory_id==id,];stopifnot(length(unique(x$seed_v))==1,length(unique(x$seed_l))==1)
}
load_frozen();td<-tempfile('validation-tests-');ensure(td)
compile_simulator(td)
set.seed(234);native<-.Call('coupled_latent',10L,.01,1.5,as.integer(c(4,8,16)))
set.seed(234);z<-matrix(rnorm(160)*sqrt(.01/16),16)
ref<-matrix(0,11,3)
for(l in 1:3){m<-c(4,8,16)[l];for(i in 1:10){v<-ref[i,l];for(j in 1:m){dw<-sum(z[((j-1)*(16/m)+1):(j*(16/m)),i]);v<-v-1.5^2*v*.01/m+1.5*sqrt(1+v^2)*dw};ref[i+1,l]<-v}}
stopifnot(max(abs(native-ref))<1e-12)
row<-make_design('smoke')[1,,drop=FALSE];row$T<-.04;row$n<-4L;row$dt<-.01
raw<-list(residuals=list('4'=1:8),delta=.005,T=.04,truth=row[c('alpha','sigma2','sigma1','rho')],diagnostics=data.frame(level=4,T=.04,max_abs_V=1))
stopifnot(identical(observe(raw,row)$R,as.numeric(c(3,7,11,15))))
x<-data.frame(estimate=c(1.1,NA),truth=c(1,1),status=c('success','rho:invalid'),parameter='sigma1')
m<-metric_row(x);stopifnot(m$N==2,m$finite==1,abs(m$rmse-.1)<1e-12)
# A real data-generating law checks the tail adapters and scale units.
set.seed(4242);r<-symmetric_stable_r(100000L,1.3)*.01^(1/1.3)
h<-benchmark_tail(r,.01,'hill_sqrt');g<-benchmark_tail(r,.01,'gpd_ml_99')
stopifnot(h$status=='success',g$status=='success',abs(h$estimates['alpha']-1.3)<.3,abs(g$estimates['alpha']-1.3)<.3)
stopifnot(abs(g$estimates['sigma2']-1)<.4)
# Check no scalar-dt truncation in the existing beta interface: supply n dt entries.
model<-get_kernel('standard',td)[[1]]
obs<-block_cf_log_observations(r[1:1000],rep(.001,1000),1.3,1,50L,freq_mults=c(.25,.5,.75))
source('original_source/work/chunked_exact_endpoint_likelihood.R')
stopifnot(obs$n_blocks==20)
# Compare unchanged full and chunked forward likelihoods on exactly the same summaries.
y<-obs$y_dejumped
if(is.null(y)) {print(names(obs));stop('Need exact observation field')}
a<-unchunked_forward_reference(y,model,1.5,obs$freqs,1.3,1,.001,50L)
b<-conditional_exact_cf_forward_loglik(y,model,1.5,obs$freqs,1.3,1,.001,50L)
stopifnot(abs(a$loglik-b$loglik)<1e-8)

cat('PASS: manifest, design counts, stream identities, coupled simulator, exact aggregation, failure denominators, standard tail benchmarks and observed-block shape\n')
unlink(td,recursive=TRUE)
