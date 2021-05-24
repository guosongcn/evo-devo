library(plotrix)
ccs<-c('unlimited')
cols<-c('#d95f02AA','#1b9e77AA')
pdf('dNdS_6_sp_primary_all2.pdf',width=3,height=6,useDingbats=F)
par(mfrow=c(1,1),font=1,font.main=1,font.axis=1,font.lab=1,mgp=c(2,1,0))
for(cc in ccs)
{
load(paste0('KaKs_all_list_',cc,'.Rdata'))
ancient<-c(res[[1]][[1]],res[[2]][[1]],res[[3]][[1]],res[[4]][[1]],res[[5]][[1]],res[[6]][[1]],res[[7]][[1]],res[[8]][[1]])
linear<-c(res[[1]][[2]],res[[2]][[2]],res[[3]][[2]],res[[4]][[2]],res[[5]][[2]],res[[6]][[2]],res[[7]][[2]],res[[8]][[2]])
p<-wilcox.test(ancient,linear,alternative='l')$p.value
print(p)
print(head(rev(sort(c(ancient,linear)))))
if(cc=='unlimited') ancient<-ancient[ancient<0.6]
if(cc=='unlimited') linear<-linear[linear<0.6]

#boxplot(ancient,linear,names=c('PO','LO'),main='raw',ylab='Ka/Ks',cex.lab=1.3,
#		cex.main=1.5,cex.axis=1.2,col=cols[c(1,2)],boxwex=0.3)

tmp2<-boxplot(linear,main='not draw outliers',ylab='Ka/Ks',cex.lab=1.3,border=cols[c(2)],lwd=1.5,
		cex.main=1.5,cex.axis=1.2,col='#1b9e7777',boxwex=0.6,yaxt='n',outline=F,at=2,xlim=c(0.5,2.5),ylim=c(0,0.33))
tmp1<-boxplot(ancient,main='',ylab='',border=cols[c(1)],lwd=1.5,bty='n',
		cex.axis=1.2,col='#d95f0277',boxwex=0.6,outline=F,at=1,add=T,yaxt='n')
axis(2,at=c(seq(0,0.25,by=0.05),0.2753399,0.33),labels=c(seq(0,0.25,by=0.05),0.275,0.6))
axis(1,at=1:2,labels=c('PO','LO'))
print(range(tmp1$out))
print(range(tmp2$out))
#axis.break(2,breakpos=0.2753399,style='zigzag',brw=0.01)
a1<-tmp1$out;a2<-tmp2$out
a1<-sapply(a1,function(x) if(x<=0.2753399) return(x) else (x-0.2753399)/(0.6-0.2753399)*(0.33-0.2753399)+0.2753399 )
a2<-sapply(a2,function(x) (x-0.2753399)/(0.6-0.2753399)*(0.33-0.2753399)+0.2753399)
points(rep(1,length(a1)),a1,cex=0.5,col='#d95f0255')
points(rep(2,length(a2)),a2,cex=0.5,col='#1b9e7755')
segments(1,0.327,2,0.327)
text(1.5,0.332,'***',cex=1.4)

#if(cc=='unlimited') {middle<-0.3;all<-middle*1.2;max<-0.6}
#ancient<-sapply(ancient,function(x) if(x<=middle) return(x) else (x-middle)/(max-middle)*(all-middle)+middle )
#linear<-sapply(linear,function(x) if(x<=middle) return(x) else (x-middle)/(max-middle)*(all-middle)+middle )
#boxplot(ancient,linear,names=c('PO','LO'),main='rm',ylab='Ka/Ks',
#						cex.lab=1.3,cex.main=1.5,cex.axis=1.2,col=cols[c(1,2)],boxwex=0.3,yaxt='n') 
#axis(2,at=c(seq(0,0.3,by=0.05),0.36),labels=c(seq(0,0.3,by=0.05),0.6))
#axis.break(2,breakpos=0.33,style='zigzag',brw=0.15)

}
dev.off()

#plot(ecdf(ancient),cex.points=0.3,col.points='#d95f02AA',ylab='cumulative frequency',xlab='Ka/Ks')
#plot(ecdf(linear),add=T,cex.points=0.3,col.points='#1b9e77AA')
#segments(1,0.229,2,0.229)
#text(1.5,0.232,'***',cex=2)

#ancient<-sapply(ancient,function(x) if(x>0.2343748) else x)
#boxplot(ancient,other,names=c('ORP','parallel'),main='not draw outliers',ylab='Ka/Ks',cex.lab=1.3,
#		cex.main=1.5,cex.axis=1.2,col=cols[c(1,3)],boxwex=0.5,outline=F,ylim=c(0,0.2343748*2))
#segments(1,0.229,2,0.229)
#text(1.5,0.232,'***',cex=2)
#boxplot(ancient[ancient>0 & ancient<0.3],other[other>0 & other<0.3],main='',ylab='Ka/Ks',cex.lab=1.3,
#		cex.main=1.5,cex.axis=1.2,col=cols[c(1,3)])
#text(1:2,-0.05,c('ORP','other'),xpd=T,cex=1.4)
#plot(rep(1:2,c(length(ancient[ancient>=0.2]),length(other[other>=0.2]))),c(ancient[ancient>=0.2],other[other>=0.2]),
#	main='',ylab='Ka/Ks',cex.lab=1.3,cex.main=1.5,cex.axis=1.2,col=rep(cols[c(1,3)],c(length(ancient[ancient>=0.2]),length(other[other>=0.2]))))
#text(1:2,-0.05,c('ORP','other'),xpd=T,cex=1.4)

#plot(c(ancient,other),rep(1:2,c(length(ancient),length(other))),col=rep(cols[c(1,3)],c(length(ancient),length(other))),type='p')
#plot(jitter(rep(c(10,30),c(length(ancient),length(other))),amount=9),c(ancient,other),col=rep(cols[c(1,3)],c(length(ancient),length(other))),
#	 type='p',xlim=c(0,40),cex=0.4)
