load("../Data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")
load("../Data/1.mean_stage_associated_gene_matrix_rpkm2.Rdata")


sign_mat=orth6_mat_pre_norm[orth6_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth6_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth6_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth6_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth6_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth6_mat_pre_norm$Dr %in% age_related_gene_list$Dr,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,7:ncol(sign_mat)]

sign_mat[is.na(sign_mat)]<-0


dis=1-cor(t(sign_mat),method="spearman")
hc=hclust(as.dist(dis))
hc_clust=cutree(hc,k=6)

pdf("Fig3S.6sp_hclust6.pdf")

A2Rplot(hc, k = 6, boxes = FALSE, col.up = "gray50",
col.down = colorRampPalette(brewer.pal(9,"Set1"))(12),main="6species hierachical cluster")

dev.off()

sp6_cluster6=list()
sort_order=sort(table(hc_clust),decreasing=T)
for(i in 1:6)
{
    sp6_cluster6[[i]]=names(hc_clust[hc_clust==names(sort_order)[i]])
}
save(sp6_cluster6,file="needleman_fromstat_d1_sp6_cluster6_hclust.Rdata")


plot_6sp=function(Mm_cluster,sign_mat)
{
    for(ii in 1:length(Mm_cluster))
    {
        gene=Mm_cluster[[ii]]
        
        gene_mat=sign_mat[rownames(sign_mat) %in% gene,]
        gene_means=apply(gene_mat,2,mean)
        gene_sds=apply(gene_mat,2,sd)
        Mm=smooth.spline(1:33,gene_means[1:33],df=6)
        
        Gg=smooth.spline(1:33,gene_means[34:66],df=6)
        Ps=smooth.spline(1:33,gene_means[67:99],df=6)
        Xl=smooth.spline(1:33,gene_means[100:132],df=6)
        Xt=smooth.spline(1:33,gene_means[133:165],df=6)
        Dr=smooth.spline(1:33,gene_means[166:198],df=6)
        
        plot(1:33,gene_means[1:33],type="n",ylim=c(-1.5,1.5),ylab="",xlab="",,xaxt="n")
        # axis(side = 1,at=1:33,labels=FALSE ,las=2,cex.axis=0.8)
        
        
        
        
        lines(Mm,lwd=1.5,col=cols[1])
        lines(Gg,lwd=1.5,col=cols[2])
        lines(Ps,lwd=1.5,col=cols[3])
        lines(Xl,lwd=1.5,col=cols[4])
        lines(Xt,lwd=1.5,col=cols[5])
        lines(Dr,lwd=1.5,col=cols[6])
               
        
        
        mtext(nrow(gene_mat))
        
    }
    
    
    
    
}
Mm_cluster=sp6_cluster6
cols=rev(colorRampPalette(brewer.pal(9,"PuBu"))(12))
pdf("Fig3B.6sp_6cluster_plot.pdf",height=1.5,width=10)
par(mfrow=c(1,6))
par(mar=c(1,2,2,1))
#plot(hc)
plot_6sp(Mm_cluster,sign_mat)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("topright",c("Mm","Gg","Ps","Xl","Xt","Dr","Bf"),col=cols[1:6],pch=19,cex=0.8)
dev.off()
