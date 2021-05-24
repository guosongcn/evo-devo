sps<-c('Gg','Ps','Xl','Xt','Dr','Bf','Ci','Cg')
load('../back_age_genes.Rdata') #sp_back_age_genes
Bf_back_age_genes<-sapply(strsplit(Bf_back_age_genes,'-'),'[',2)
load('../modify_AL.Rdata') #ancient_genes linear_genes

FUN<-function(sp)
	{
	back<-get(paste0(sp,'_back_age_genes'))
	ancient<-ancient_genes[[sp]]
	if(sp=='Bf') ancient<-sapply(strsplit(ancient,'-'),'[',2)
	linear<-linear_genes[[sp]]
	if(sp=='Bf')linear<-sapply(strsplit(linear,'-'),'[',2)
	other<-setdiff(back,c(ancient,linear))	
	data<-read.table(paste0(sp,'/',sp,'_dNdS_res.txt'),header=F,as.is=T)	
	for(i in 5)
			{
				me<- -100
				ancient_kaks<-data[ (as.character(data[,1]) %in% ancient) & data[,3]>me,i]
				linear_kaks<-data[ (as.character(data[,1]) %in% linear) & data[,3]>me,i]
				other_kaks<-data[(as.character(data[,1]) %in% other) & data[,3]>me,i]	

				ancient_kaks<-na.omit(ancient_kaks)
				linear_kaks<-na.omit(linear_kaks)
				other_kaks<-na.omit(other_kaks)
			}	
	list(ancient=ancient_kaks,linear=linear_kaks,other=other_kaks)
	}
res<-lapply(sps,FUN)
names(res)<-sps
save(res,file='KaKs_all_list_unlimited.Rdata')
