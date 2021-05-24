library(Biostrings)
library(seqinr)
#sps<-c('Xl','Xt','Dr','Bf','Ci','Cg')
sps<-c('Gg','Ps')
for(sp in sps) {assign(sp,readDNAStringSet(paste0(sp,'_stranded_gene.fa'),format='fasta',use.names=T))}
Mm<-readDNAStringSet(paste0('Mm','_stranded_gene.fa'),format='fasta',use.names=T)

load('../modify_AL.Rdata')
load('../back_age_genes.Rdata')
#if(sp=='Bf') Bf_back_age_genes<-sapply(strsplit(Bf_back_age_genes,'-'),'[',2)

for(sp in sps)
		{
		if(file.exists(paste0(sp,'.txt'))) next
		file.create(paste0(sp,'.txt'))
		print(sp)
		sp_seq<-get(sp)
		sp_seq<-sp_seq[unique(c(ancient_genes[[sp]],linear_genes[[sp]]))]
		back_age_genes<-get(paste0(sp,'_back_age_genes'))	
		#and with Mm_back_age_genes
		#file.create(paste0(sp,'/','length',length(sp_seq)))	
		long<-0
		lengths<-0	
		for(i in 1:length(sp_seq))
				{
					name<-names(sp_seq)[i]
					Mm_name<-names(back_age_genes)[back_age_genes==name]

					data<-sp_seq[[i]]
					data<-subseq(data,start=1,end=length(data)-3)
					pro<-Biostrings::translate(data,if.fuzzy.codon='solve')
					whole_pro<-toString(pro)
					if(grepl('\\*',whole_pro)) {long<-long+1;next}

					Mm_data<-Mm[[Mm_name]]
					Mm_data<-subseq(Mm_data,start=1,end=length(Mm_data)-3)
					Mm_pro<-Biostrings::translate(Mm_data,if.fuzzy.codon='solve')
					whole_Mm_pro<-toString(Mm_pro)
					if(grepl('\\*',whole_Mm_pro)) {long<-long+1;next}

					if(max(length(pro),length(Mm_pro))>=32760 ) {long<-long+1;next}
					lengths<-lengths+1

					align<-pairwiseAlignment(pro,Mm_pro,substitutionMatrix="BLOSUM50",gapOpening=0, gapExtension=8)
					pattern<-toString(pattern(align))
					subject<-toString(subject(align))

					sp_gene_vec<-s2c(toString(data))
					Mm_gene_vec<-s2c(toString(Mm_data))
					sp_pro_vec<-s2c(pattern)
					Mm_pro_vec<-s2c(subject)

					sp_index<-grep('-',sp_pro_vec)
					Mm_index<-grep('-',Mm_pro_vec)

					if(length(sp_index)==0){sp_index<-1000000000000000}
					if(length(Mm_index)==0){Mm_index<-1000000000000000} 
					clean_pattern<-c2s(sp_pro_vec[-sp_index])
					clean_l1<-nchar(clean_pattern)
					clean_subject<-c2s(Mm_pro_vec[-Mm_index])
					clean_l2<-nchar(clean_subject)

					sp_index1<-regexpr(clean_pattern,whole_pro,perl=T)
					sp_index2<-sp_index1+clean_l1-1
					Mm_index1<-regexpr(clean_subject,whole_Mm_pro,perl=T)
					Mm_index2<-Mm_index1+clean_l2-1

					cc1<-character(nchar(pattern)*3)
					cc2<-character(nchar(subject)*3)
					
					if(sp_index!=1000000000000000)
					{
					cc1[(sp_index-1)*3+1]<-'-'
					cc1[(sp_index-1)*3+2]<-'-'
					cc1[(sp_index-1)*3+3]<-'-'
					}
					cc1[cc1!='-']<-sp_gene_vec[ ((sp_index1-1)*3+1):(3*sp_index2) ]


					if(Mm_index!=1000000000000000)
					{
					cc2[(Mm_index-1)*3+1]<-'-'
					cc2[(Mm_index-1)*3+2]<-'-'
					cc2[(Mm_index-1)*3+3]<-'-'
					}
					cc2[cc2!='-']<-Mm_gene_vec[ ((Mm_index1-1)*3+1):(3*Mm_index2) ]

					length<-length(cc1)
					cc1<-c2s(cc1)
					cc2<-c2s(cc2)

					cat(paste0(2,' ',length,'\n'),file=paste0(sp,'/',sp,'.nuc'),append=T)
					cat(paste0(sp,'  ',cc1,'\n'),file=paste0(sp,'/',sp,'.nuc'),append=T)
					cat(paste0('Mm','  ',cc2,'\n'),file=paste0(sp,'/',sp,'.nuc'),append=T)
					cat(paste0(name,' ',align@score,'\n'),file=paste0(sp,'/',sp,'.gene'),append=T)
					cat(paste0(Mm_name,'\n'),file=paste0(sp,'/','Mm','.gene'),append=T)
				}
		file.create(paste0(sp,'/','length_',lengths))
		file.create(paste0(sp,'/','long_',long))
		}
