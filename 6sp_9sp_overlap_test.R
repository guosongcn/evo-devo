sign_mat=orth6_mat_pre_norm[orth6_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth6_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth6_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth6_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth6_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth6_mat_pre_norm$Dr %in% age_related_gene_list$Dr,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,7:ncol(sign_mat)]

sign_mat[is.na(sign_mat)]<-0
# all_gene=orth6_mat_pre_norm$Mm
all_gene6=rownames(sign_mat)



sign_mat=orth9_mat_pre_norm[orth9_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth9_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth9_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth9_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth9_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth9_mat_pre_norm$Dr %in% age_related_gene_list$Dr & orth9_mat_pre_norm$Bf %in% age_related_gene_list$Bf & orth9_mat_pre_norm$Ci %in% age_related_gene_list$Ci & orth9_mat_pre_norm$Cg %in% age_related_gene_list$Cg,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,10:ncol(sign_mat)]

sign_mat[is.na(sign_mat)]<-0


all_gene9=rownames(sign_mat)

overlap_check=function(cluster1,cluster2,name1,name2,all_gene1,all_gene2)
{
    
    
    overlap_cluster_fishertest=matrix(ncol=length(cluster2),nrow=0)
    
    for(i in 1:length(cluster1))
    {
        pval=c()
        for(j in 1:6)
        {
            
            gene1=cluster1[[i]]
            gene2=cluster2[[j]]
            overlap=intersect(gene1,gene2)
            not_cluster1=intersect(gene1,setdiff(all_gene2,gene2))
            not_cluster2=intersect(setdiff(all_gene1,gene1),gene2)
            non_all=intersect(setdiff(all_gene1,gene1),setdiff(all_gene2,gene2))
            fisher.p=fisher.test(matrix(c(length(overlap),length(not_cluster1),length(not_cluster2),length(non_all)),2),alternative="greater")$p.val
            pval=c(pval,fisher.p)
            
            
        }
        
        overlap_cluster_fishertest=rbind(overlap_cluster_fishertest,pval)
        
    }
    overlap_cluster_fishertest=overlap_cluster_fishertest*(length(cluster1)*length(cluster2))
    rownames(overlap_cluster_fishertest)=paste(1:length(cluster1),name1,sep="/")
    colnames(overlap_cluster_fishertest)=paste(1:length(cluster2),name2,sep="/")
    
    return(overlap_cluster_fishertest)
}


overlap_cluster_fishertest=overlap_check(sp6_cluster6,sp9_cluster6,"CL6sp","CL9sp",all_gene6,all_gene9)

plot_chord_diagram=function(mat)
{
    #library(circlize,lib="/home/song/R_lib/")
    library(circlize)
    cor_mat=-log10(mat*nrow(mat)*ncol(mat))
    # for(i in 1:36){cor_mat[i,i]<-0}
   col_mat = rand_color(length(cor_mat), transparency = 0.5)
   dim(col_mat) = dim(cor_mat) # to make sure it is a matrix
    col_mat[cor_mat < 5] = "#00000000"
    chordDiagram(cor_mat,  col = col_mat,transparency = 0.3)
    circos.clear()
}


pdf("transition_from_6sp_chordgraph.pdf")
par(cex = 1, mar = c(0, 0, 0, 0))
plot_chord_diagram(overlap_cluster_fishertest)
#plot_chord_diagram(t(overlap_cluster_fishertest))

dev.off()