pdf('pairwise_cor_to_Mm2.pdf',15,8.5)
par(font=2,lwd=3,font.axis=2,font.main=2,font.lab=2,mfrow=c(2,4),mgp=c(2,1,0))
sps<-c('Mm','Gg','Ps','Xl','Xt','Dr','Bf','Ci','Cg')
sps<-combn(sps,2)
dir1<-'/picb/compbio.work/Xuchuan/2015_Dec_linearship/poly/corr_0.7_spline_both_age_related/'
dir2<-'/picb/compbio.work/Xuchuan/2015_Dec_linearship/poly/corr_0.7_age_model_both_age_related/'
func<-function(.x)
{
sp1<-.x[1];sp2<-.x[2]

load(paste(dir1,'real_number/',sp2,'_to_',sp1,'.Rdata',sep=''));real1<-result[-1]
per1<-c()
for(j in 1:200){load(paste(dir1,'permute_data/','cor_permutation_',sp2,'_to_',sp1,'_',j,'.Rdata',sep=''));per1<-rbind(per1,result[-1])}
per_mean1<-colMeans(per1)
#sd1<-apply(per1,2,sd);ydown1<-real1-per_mean1-sd1;yup1<-real1-per_mean1+sd1

load(paste(dir2,'real_number/',sp2,'_to_',sp1,'.Rdata',sep=''));real2<-result[-1]
per2<-c()
for(j in 1:200){load(paste(dir2,'permute_data/','cor_permutation_',sp2,'_to_',sp1,'_',j,'.Rdata',sep=''));per2<-rbind(per2,result[-1])}
per_mean2<-colMeans(per2)
#sd2<-apply(per2,2,sd);ydown2<-real2-per_mean2-sd2;yup2<-real2-per_mean2+sd2


xs<-barplot(border=NA,width=0.75,space=0.53,ylim=c(min(real1-per_mean1,real2-per_mean2)-50,max(real1-per_mean1,real2-per_mean2)+50),
		    names.arg=NA,real1-per_mean1,col='#FF0000BB',xlab=paste('stages of ',sp1,sep=''),ylab='number of genes')
xs<-as.numeric(xs)
if(sp2=='Xl')rect(lwd=0.5,xs[13]-0.375,-20,xs[13]+0.375,max((real1-per_mean1)[13],(real2-per_mean2)[13])+20,col=NA,border='blue')
if(sp2=='Xt')rect(lwd=0.5,xs[14]-0.375,-20,xs[15]+0.375,max((real1-per_mean1)[14:15],(real2-per_mean2)[14:15])+20,col=NA,border='blue')
if(sp2=='Dr')rect(lwd=0.5,xs[8]-0.375,-20,xs[10]+0.375,max((real1-per_mean1)[8:10],(real2-per_mean2)[8:10])+20,col=NA,border='blue')
if(sp2=='Bf')rect(lwd=0.5,xs[5]-0.375,-20,xs[6]+0.375,max((real1-per_mean1)[5:6],(real2-per_mean2)[5:6])+20,col=NA,border='blue')
if(sp2=='Ci')rect(lwd=0.5,xs[8]-0.375,-20,xs[9]+0.375,max((real1-per_mean1)[8:9],(real2-per_mean2)[8:9])+20,col=NA,border='blue')
if(sp2=='Cg')rect(lwd=0.5,xs[3]-0.375,-20,xs[4]+0.375,max((real1-per_mean1)[3:4],(real2-per_mean2)[3:4])+20,col=NA,border='blue')
sapply(xs,function(x) {index=which(xs==x);rect(x-0.45,0,x+0.45,(real2-per_mean2)[index],border=NA,col='#0000FF55');} )
title(main=paste(sp2,'to',sp1),line=0)
legend('topleft',c('smooth spline','polynomial'),col=c('#FF0000BB','#0000FF33'),bty='n',pch=15)
}
apply(sps[,1:8],2,func)
dev.off()
