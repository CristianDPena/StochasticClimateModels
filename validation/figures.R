#!/usr/bin/env Rscript
# Presentation is separate from estimation: rerendering never changes fit records.
source('validation/R/common.R')
out<-normalizePath(arg('out','validation/results/laptop-v1'),mustWork=TRUE)
x<-read.csv(file.path(out,'summary','estimates.csv'))
if(!nrow(x))quit()
ensure(file.path(out,'figures'))
parameters<-c('alpha','sigma2','sigma1','rho')
titles<-expression(alpha,sigma[2],sigma[1],rho)
ordinary<-x[x$component%in%c('main','smoke')&x$method=='fixed',]
reg_order<-c('B','Aminus','Aplus','Jminus','Cminus','Sminus','H')
reg<-reg_order[reg_order%in%ordinary$regime]
blue<-'#0072B2';colors<-c('#0072B2','#D55E00','#009E73')
set.seed(442901)
pdf(file.path(out,'figures','recovery.pdf'),width=10,height=7,useDingbats=FALSE)
par(mfrow=c(2,2),mar=c(4,5,2,1),mgp=c(3,1,0))
for(i in seq_along(parameters)) {
 y<-ordinary[ordinary$parameter==parameters[i],];y$error<-if(i==1)y$estimate-y$truth else y$estimate/y$truth-1
 finite<-is.finite(y$error);rng<-range(c(0,y$error[finite]));pad<-max(diff(rng)*.1,.005)
 plot(c(.5,length(reg)+.5),rng+c(-pad,pad),type='n',xaxt='n',xlab='Regime',
  ylab=if(i==1)'Signed error' else 'Relative signed error',main=titles[i])
 axis(1,at=seq_along(reg),labels=reg,cex.axis=.85);abline(h=0,lty=2,col='grey50')
 for(j in seq_along(reg)) {
  vals<-y$error[y$regime==reg[j]];good<-vals[is.finite(vals)]
  if(length(good)) {
   boxplot(good,at=j,add=TRUE,axes=FALSE,outline=FALSE,boxwex=.5,col='#DCEAF2')
   points(j+runif(length(good),-.12,.12),good,pch=16,cex=.65,col=adjustcolor(blue,.65))
  }
  text(j,rng[1]-pad*.65,sprintf('%d/%d',length(good),length(vals)),cex=.7,col='grey35')
 }
}
dev.off()
methods<-c('fixed','hill_sqrt','gpd_ml_99');labels<-c('Fixed method','Hill: square-root rank','GPD ML: 99% threshold')
pdf(file.path(out,'figures','benchmarks.pdf'),width=11,height=8,useDingbats=FALSE)
par(mfrow=c(2,1),mar=c(4,5,3,1),mgp=c(3,1,0))
for(i in 1:2) {
 y<-x[x$component%in%c('main','smoke')&x$parameter==parameters[i],];y$error<-abs(y$estimate-y$truth)
 ymax<-max(c(.01,y$error[is.finite(y$error)]))*1.2
 plot(c(.4,length(reg)+.6),c(-.09*ymax,ymax),type='n',xaxt='n',xlab='Regime',ylab='Absolute error',main=titles[i])
 axis(1,at=seq_along(reg),labels=reg);abline(h=0,col='grey80')
 for(j in seq_along(reg))for(m in 1:3) {
  vals<-y$error[y$regime==reg[j]&y$method==methods[m]];good<-vals[is.finite(vals)];pos<-j+(m-2)*.22
  if(length(good)) {
   boxplot(good,at=pos,add=TRUE,axes=FALSE,outline=FALSE,boxwex=.18,col=adjustcolor(colors[m],.25))
   points(pos+runif(length(good),-.035,.035),good,pch=16,cex=.5,col=adjustcolor(colors[m],.6))
  }
  text(pos,-.045*ymax,sprintf('%d/%d',length(good),length(vals)),cex=.6,col=colors[m])
 }
 legend('topright',legend=labels,col=colors,pch=15,bty='n',cex=.75)
}
dev.off()
# The information-growth figure from summarize.R already uses matched trajectory IDs.
# Annotate rather than imply information-growth evidence when only reference fits exist.
if(!any(x$component%in%c('infill','horizon'))) {
 pdf(file.path(out,'figures','information_growth.pdf'),width=8,height=4,useDingbats=FALSE)
 plot.new();text(.5,.6,'Information-growth results are pending',cex=1.2)
 text(.5,.4,'Reference fits alone do not assess improvement with n or T.',cex=.9);dev.off()
}
writeLines(c('Generated figures from saved estimates; finite/total counts appear below groups.',
 'Recovery uses signed alpha error and relative signed scale/clock errors.',
 'Benchmark alpha values outside the stable domain remain visible; undefined scale estimates are counted.',
 'Partial studies and smoke runs are not manuscript evidence.',
 paste('Plot source MD5:',sha('validation/figures.R')),
 paste('Estimates CSV MD5:',sha(file.path(out,'summary','estimates.csv')))),file.path(out,'figures','README.txt'))
