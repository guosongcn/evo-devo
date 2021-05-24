load('sp2Mm_expreesion.Rdata')
source('/picb/compbio.work/Xuchuan/script/descending.cl.R')
load('ant_lne_gene.Rdata') #ancient_genes linear_genes
load('cls.Rdata')
load('species_stage_names.Rdata')
sps<-c('Xl','Xt','Dr','Bf','Ci','Cg')
pdf('./pdf/ancient_gene_pattern3.pdf',16,24)
par(mfrow=c(6,4),font=2,font.main=2,font.axis=2,font.lab=2,mgp=c(2,1,0))
for(sp in sps)
	{
	data<-get(paste(sp,'ancient_expression_norm',sep='_'))
	Mm_data<-data[,3:19];rownames(Mm_data)<-data[,'Mm'];Mm_time<-ncol(Mm_data)
	sp_data<-data[,20:ncol(data)];rownames(sp_data)<-data[,'sp'];sp_time<-ncol(sp_data)
	cl<-cls[[sp]][['cl']]
	stopifnot(length(cl)==nrow(Mm_data),length(cl)==nrow(sp_data))
	splt_Mm<-split(Mm_data,cl);splt_sp<-split(sp_data,cl)
	splt_Mm<-splt_Mm[order(-sapply(splt_Mm,nrow))]
	splt_sp<-splt_sp[order(-sapply(splt_sp,nrow))]
	FUN<-function(x)
		{
		Mm<-splt_Mm[[x]];Sp<-splt_sp[[x]]
		stopifnot(nrow(Mm)==nrow(Sp))
		Mm_mns<-colMeans(Mm);sp_mns<-colMeans(Sp)
		Mm_sds<-apply(Mm,2,sd);sp_sds<-apply(Sp,2,sd)
		Mm_up<-Mm_mns+Mm_sds;sp_up<-sp_mns+sp_sds
		Mm_down<-Mm_mns-Mm_sds;sp_down<-sp_mns-sp_sds
		if(sp=='Xl')mouse=14;if(sp=='Xt')mouse=15.5;if(sp=='Dr')mouse=10
		if(sp=='Bf')mouse=6.5;if(sp=='Ci')mouse=9.5;if(sp=='Cg')mouse=4.5
		new_time<-(1:sp_time-1)*(mouse-1)/(sp_time-1)+1;max_sp_time<-max(new_time)
		plot(1,type='n',xaxt='n',xlim=c(1,max(Mm_time,max_sp_time)),ylim=c(min(Mm_down,sp_down)-1,max(Mm_up,sp_up)+1),xlab='',ylab='scaled expression')
		title(main=paste(sp,' cluster_',x,':',nrow(Mm),sep=''),line=-1.5)
		sp_stage_names<-get(paste0(sp,'_stage_names'))
		axis(side=1,at=1:17,labels=Mm_stage_names,las=2,cex.axis=0.82,col.ticks=rgb(0,0,1,0.9),col.axis=rgb(0,0,1,0.9))
		if(sp=='Cg') axis(side=3,at=new_time,labels=c(sp_stage_names[1],'','','','',sp_stage_names[6],'','','','',sp_stage_names[11],'','','',sp_stage_names[15],'','','',sp_stage_names[sp_time]),col.axis=rgb(1,0,0,0.9),las=2,cex.axis=0.82,col.ticks=rgb(1,0,0,0.9))
		if(sp!='Cg') axis(side=3,at=new_time,labels=sp_stage_names,col.axis=rgb(1,0,0,0.9),las=2,cex.axis=0.82,col.ticks=rgb(1,0,0,0.9))
		polygon(c(1:Mm_time,Mm_time:1),c(Mm_down,rev(Mm_up)),border=NA,col=rgb(0,0,1,0.2))
		polygon(c(new_time,rev(new_time)),c(sp_down,rev(sp_up)),border=NA,col=rgb(1,0,0,0.4))
		#lines(1:Mm_time,Mm_mns,col=rgb(0,0,1,0.7),lwd=2)
		#lines(new_time,sp_mns,col=rgb(1,0,0,0.9),lwd=2)
		points(1:Mm_time,Mm_mns,col=rgb(0,0,1,0.7),pch=19)
		points(new_time,sp_mns,col=rgb(1,0,0,0.9),pch=19)		
		lines(smooth.spline(1:Mm_time,Mm_mns,df=6),col=rgb(0,0,1,0.7),lwd=2)
		lines(smooth.spline(new_time,sp_mns,df=ifelse(sp=='Bf',3,6)),col=rgb(1,0,0,0.9),lwd=2)
		legend('bottomright',legend=rev(c('Mm',sp)),lty=1,lwd=2,col=rev(c(rgb(0,0,1,0.7),rgb(1,0,0,0.9))))
		}
	X<-length(unique(cl));lapply(1:X,FUN)
	cl<-cls[[sp]][[2]]
	cl<-descending.cl(cl)
	ancient<-ancient_genes[[sp]]
	if(sp=='Bf') ancient<-sapply(strsplit(ancient,'-'),'[',2)
	data<-read.table(paste0('../KaKs/',sp,'/',sp,'_dNdS_res.txt'),header=F,as.is=T)	
	for(i in 5)
			{
				me<-quantile(data[,3],0.5)
				ancient_kaks<-data[ (as.character(data[,1]) %in% ancient[cl==1]) & data[,3]>me,i]
				linear_kaks<-data[ (as.character(data[,1]) %in% ancient[cl==2]) & data[,3]>me,i]
				other_kaks<-data[ (as.character(data[,1]) %in% ancient[cl==3]) & data[,3]>me,i]

				ancient_kaks<-na.omit(ancient_kaks)
				linear_kaks<-na.omit(linear_kaks)
				other_kaks<-na.omit(other_kaks)
				p1<-wilcox.test(ancient_kaks,linear_kaks)$p.value
				p2<-wilcox.test(ancient_kaks,other_kaks)$p.value
				p3<-wilcox.test(linear_kaks,other_kaks)$p.value
				cols<-c(grey(0.7))
				if(sp=='Xl')boxplot(ancient_kaks,linear_kaks,other_kaks,main=sp,ylab='Ka/Ks',cex.lab=1.3,cex.main=1.5,cex.axis=1.2,ylim=c(0,0.4),col=cols)
				if(sp %in% c('Xt','Dr','Bf'))boxplot(ancient_kaks,linear_kaks,other_kaks,main=sp,ylab='Ka/Ks',cex.lab=1.3,cex.main=1.5,cex.axis=1.2,ylim=c(0,0.5),col=cols)
				if(sp=='Ci')boxplot(ancient_kaks,linear_kaks,other_kaks,main=sp,ylab='Ka/Ks',cex.lab=1.3,cex.main=1.5,cex.axis=1.2,ylim=c(0,0.4),col=cols)
				if(sp=='Cg')boxplot(ancient_kaks,linear_kaks,other_kaks,main=sp,ylab='Ka/Ks',cex.lab=1.3,cex.main=1.5,cex.axis=1.2,ylim=c(0,0.7),col=cols)

				if(sp=='Xl') text(1:3,-0.05,c('cluster1','cluster2','cluster3'),xpd=T,cex=1.4)
				if(sp=='Cg') text(1:3,-0.09,c('cluster1','cluster2','cluster3'),xpd=T,cex=1.4)
				if(sp=='Ci') text(1:3,-0.05,c('cluster1','cluster2','cluster3'),xpd=T,cex=1.4)
				if(sp %in% c('Xt','Dr','Bf')) text(1:3,-0.06,c('cluster1','cluster2','cluster3'),xpd=T,cex=1.4)

				print(c(p1,p2,p3))
				p1<-round(p1,3);p2<-round(p2,3);p3<-round(p3,3)	
				if(sp=='Xt') {p2<-'3.89E-7';p3<-'1.50E-4'}
				#if(sp=='Bf') p2<-'2.36E-4'
				if(sp %in% c('Xl','Ci')) {segments(c(1,1,1.9,2.1,2.1,3,1,1,3),c(0.34,0.36,0.36,0.34,0.36,0.36,0.38,0.40,0.40),
						             c(1,1.9,1.9,2.1,3,3,1,3,3),c(0.36,0.36,0.34,0.36,0.36,0.34,0.40,0.40,0.38));
							  text(c(1.5,2.5,2),c(0.347,0.347,0.385),c(p1,p3,p2))	
								}
				if(sp %in% c('Xt','Dr','Bf')) {segments(c(1,1,1.9,2.1,2.1,3,1,1,3),c(0.34,0.36,0.36,0.34,0.36,0.36,0.38,0.40,0.40)+0.1,
						             c(1,1.9,1.9,2.1,3,3,1,3,3),c(0.36,0.36,0.34,0.36,0.36,0.34,0.40,0.40,0.38)+0.1);
							  text(c(1.5,2.5,2),c(0.347,0.347,0.385)+0.095,c(p1,p3,p2))	
							}
				if(sp %in% c('Cg')) {segments(c(1,1,1.9,2.1,2.1,3,1,1,3),c(0.34,0.35,0.35,0.34,0.35,0.35,0.38,0.40,0.40)+0.3,
						             c(1,1.9,1.9,2.1,3,3,1,3,3),c(0.35,0.35,0.34,0.35,0.35,0.34,0.40,0.40,0.38)+0.3);
							  text(c(1.5,2.5,2),c(0.335,0.335,0.385)+0.291,c(p1,p3,p2))	
							}
			}
	}
dev.off()
