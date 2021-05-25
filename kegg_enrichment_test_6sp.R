load("../data/keggpathway_withensembl72_gene.Rdata")

KEGG_test<-function(cluster_gene,all_gene,path_mat)
{
    path_mat=path_mat[path_mat[,2] %in% all_gene,]
    all_path=unique(path_mat[,1])
    notcluster_gene=setdiff(all_gene,cluster_gene)
    path=read.table("../data/mmu_pathway_name2.txt",sep="\t")

    sign_path=matrix(ncol=2,nrow=0)
    for(eachpath in all_path)
    {
        onepath_gene=unique(path_mat[path_mat[,1] %in% eachpath,2])
        otherpath_gene=unique(path_mat[!(path_mat[,1] %in% eachpath),2])
        
        A=length(intersect(cluster_gene,onepath_gene))
        B=length(intersect(notcluster_gene,onepath_gene))
        C=length(intersect(cluster_gene,otherpath_gene))
        D=length(intersect(notcluster_gene,otherpath_gene))
        p1=phyper(A-1,A+C,B+D,A+B,lower.tail=FALSE)
        #p2=fisher.test(matrix(c(A,B,C,D),2,2),alternative="greater")$p.val
        p_BF=p1*length(all_path)
        
        if(p_BF<0.05)
        {
            eachpath_name=as.character(path[path[,1] %in% eachpath,2])
            paths=c(eachpath_name,p_BF)
            sign_path=rbind(sign_path,paths)
            
            
        }
        
    }
    return(sign_path)
}

load("needleman_fromstat_d1_sp6_cluster6_hclust.Rdata")

load("../data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")

load("../data/age_related_gene_list.Rdata")
sign_mat=orth6_mat_pre_norm[orth6_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth6_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth6_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth6_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth6_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth6_mat_pre_norm$Dr %in% age_related_gene_list$Dr,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,7:ncol(sign_mat)]



sign_mat[is.na(sign_mat)]<-0
all_gene=rownames(sign_mat)
all_path=unique(mmu_path_ensembl72gene$path)

clusters=sp6_cluster6

sign_path_cluster=list()
sign_mat=matrix(ncol=3,nrow=0)

for(i in 1:length(clusters))
{
    print(i)
    Mm=clusters[[i]]
    Mm_path=KEGG_test(Mm,all_gene,mmu_path_ensembl72gene)
    
 
    
    sign_path_cluster[[i]]=Mm_path
    ss=cbind(cluster=rep(i,nrow(Mm_path)),Mm_path)
    sign_mat=rbind(sign_mat,ss)
    
}
save(sign_path_cluster,file="sign_path_6sp.Rdata")

write.table(sign_mat,file="sign_path_6sp_BF.csv",quote=F,col.names=F,row.names=F,sep=",")


q(save="no")

