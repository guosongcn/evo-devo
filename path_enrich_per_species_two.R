load('/picb/compbio.work/Xuchuan/2015_Dec_linearship/poly/corr_0.7_spline_both_age_related/real_number/mmu_pathway_name.Rdata') #t
load('/picb/compbio.work/Xuchuan/2015_Dec_linearship/poly/corr_0.7_spline_both_age_related/real_number/mmu_keggpathway_withensembl72_gene.Rdata') #mmu_path_ensembl72gene 
sps<-c('Xl','Xt','Dr','Bf','Ci','Cg')
load('KaKs/ant_lne_gene.Rdata') #ancient_genes
#load('spline_ancient_linear_gene_list.Rdata')
load('KaKs/back_age_genes.Rdata') #sp_back_age_genes
ancient_genes<-lapply(names(ancient_genes),function(sp) names(get(paste0(sp,'_back_age_genes')))[get(paste0(sp,'_back_age_genes')) %in% ancient_genes[[sp]]])
names(ancient_genes)<-sps
linear_genes<-lapply(names(linear_genes),function(sp) names(get(paste0(sp,'_back_age_genes')))[get(paste0(sp,'_back_age_genes')) %in% linear_genes[[sp]]])
names(linear_genes)<-sps

all_test_genes<-unique(c(unlist(ancient_genes),unlist(linear_genes)))
print(length(all_test_genes))
ps<-unique(as.character(mmu_path_ensembl72gene[mmu_path_ensembl72gene[,'gene'] %in% all_test_genes,'path']))
print(length(ps))

pathway_enrich_test_whole<-function(genes,universe,pathway,sS=0,name2func)
	{
	genes<-unique(genes);universe<-unique(universe)
	pathway<-as.data.frame(pathway,stringsAsFactors=F)
    genes<-data.frame(gene=genes,stringsAsFactors=F)                           
	universe<-data.frame(gene=universe,stringsAsFactors=F)
	test_dataframe<-merge(pathway,genes,by='gene',all=F)[,2:1]   ;#test_dataframe<-test_dataframe[order(test_dataframe[,'path']),]
	back_dataframe<-merge(pathway,universe,by='gene',all=F)[,2:1];#back_dataframe<-back_dataframe[order(back_dataframe[,'path']),]
	FUN<-function(p)
		{
		no_in_p<-nrow(test_dataframe[test_dataframe[,'path']==p,])
		no_in_u<-nrow(back_dataframe[back_dataframe[,'path']==p,])
		pvalue<-1-phyper(no_in_p-1,no_in_u,nrow(universe)-no_in_u,nrow(genes))
		func<-as.character(name2func[name2func[,'pathway']==p,'function'])
        c(p,func,no_in_u,no_in_p,pvalue)
		}
	result<-t(sapply(ps,FUN))
	colnames(result)<-c('Pathway','Function','Annotated','Significant','pvalue')
	katy<-data.frame(Pathway=as.character(result[,'Pathway']),Function=as.character(result[,'Function'])
			        ,Annotated=as.numeric(result[,'Annotated']),Significant=as.numeric(result[,'Significant']),
					pvalue=as.numeric(result[,'pvalue']),stringsAsFactors=F)
	katy<-katy[katy[,'Significant']>=sS,]
	#invisible(katy[order(katy[,'pvalue']),])
	katy
	}
res_ancient<-lapply(sps,function(sp) pathway_enrich_test_whole(ancient_genes[[sp]],unique(c(ancient_genes[[sp]],linear_genes[[sp]])),mmu_path_ensembl72gene,0,t))
names(res_ancient)<-sps
res_linear<-lapply(sps,function(sp) pathway_enrich_test_whole(linear_genes[[sp]],unique(c(ancient_genes[[sp]],linear_genes[[sp]])),mmu_path_ensembl72gene,0,t))
names(res_linear)<-sps

#final_res_ancient<-sapply(res_ancient,function(x) {sig<-x[,'Significant'];p<-x[,'pvalue'];index1<-sig<5;index2<-sig>=5;p[index2]<-p.adjust(p[index2],'BH');p[index1]<-1;p})
#final_res_linear<-sapply(res_linear,function(x) {sig<-x[,'Significant'];p<-x[,'pvalue'];index1<-sig<5;index2<-sig>=5;p[index2]<-p.adjust(p[index2],'BH');p[index1]<-1;p})
#final_res<-cbind(final_res_ancient,final_res_linear)

final_res_ancient<-sapply(res_ancient,function(x) {sig<-x[,'Significant'];p<-x[,'pvalue'];p})
final_res_linear<-sapply(res_linear,function(x) {sig<-x[,'Significant'];p<-x[,'pvalue'];p})
final_res<-cbind(final_res_ancient,final_res_linear)

colnames(final_res)<-paste(sps,rep(c('PO','LO'),each=6),sep='_')
rownames(final_res)<-res_ancient[[1]][,'Function']
save(final_res,file='path_enrich_per_species_two.Rdata')
