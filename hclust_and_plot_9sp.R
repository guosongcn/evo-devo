

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")

library(RColorBrewer)

plot_9sp=function(Mm_cluster,sign_mat)
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
        Bf=smooth.spline(1:33,gene_means[199:231],df=6)
        Ci=smooth.spline(1:33,gene_means[232:264],df=6)
        Cg=smooth.spline(1:33,gene_means[265:297],df=6)
        
        plot(1:33,gene_means[1:33],type="n",ylim=c(-1.5,1.5),ylab="",xlab="",,xaxt="n")
               
        lines(Mm,lwd=1.5,col=cols[1])
        lines(Gg,lwd=1.5,col=cols[2])
        lines(Ps,lwd=1.5,col=cols[3])
        lines(Xl,lwd=1.5,col=cols[4])
        lines(Xt,lwd=1.5,col=cols[5])
        lines(Dr,lwd=1.5,col=cols[6])
        lines(Bf,lwd=1.5,col=cols[7])
        lines(Ci,lwd=1.5,col=cols[8])
        lines(Cg,lwd=1.5,col=cols[9])
        
        
        
        
        
        # mtext(nrow(gene_mat))
        
    }




}

load("../Data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")
load("../Data/1.mean_stage_associated_gene_matrix_rpkm2.Rdata")




cols<-c(rgb(238/255,94/255,89/255),rgb(200/255,130/255,33/255),rgb(127/255,159/255,46/255),rgb(39/255,166/255,58/255),rgb(31/255,173/255,136/255),rgb(28/255,170/255,220/255),rgb(92/255,129/255,194/255),rgb(162/255,100/255,166/255),rgb(221/255,176/255,150/255))



###6cluster

##age related
load("../Data/age_related_gene_list.Rdata")
###age_related
sign_mat=orth9_mat_pre_norm[orth9_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth9_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth9_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth9_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth9_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth9_mat_pre_norm$Dr %in% age_related_gene_list$Dr & orth9_mat_pre_norm$Bf %in% age_related_gene_list$Bf & orth9_mat_pre_norm$Ci %in% age_related_gene_list$Ci & orth9_mat_pre_norm$Cg %in% age_related_gene_list$Cg,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,10:ncol(sign_mat)]

sign_mat[is.na(sign_mat)]<-0
dis=1-cor(t(sign_mat),method="spearman")

hc=hclust(as.dist(dis))
hc_clust=cutree(hc,k=6)


pdf("Fig3A.9sp_hclust_tree_age_related_6cluster.pdf")

A2Rplot(hc, k = 6, boxes = FALSE, col.up = "gray50",
col.down = colorRampPalette(brewer.pal(9,"Set1"))(12),main="9species hierachical cluster")

dev.off()

sp9_cluster6=list()
sort_order=sort(table(hc_clust),decreasing=T)

for(i in 1:6)
{
	sp9_cluster6[[i]]=names(hc_clust[hc_clust==names(sort_order)[i]])
}



save(sp9_cluster6,file="needleman_fromstat_d1_sp9_cluster6_hclust_age_related_df10.Rdata")

Mm_cluster=sp9_cluster6
#cols<-c(rgb(238/255,94/255,89/255),rgb(200/255,130/255,33/255),rgb(127/255,159/255,46/255),rgb(39/255,166/255,58/255),rgb(31/255,173/255,136/255),rgb(28/255,170/255,220/255),rgb(92/255,129/255,194/255),rgb(162/255,100/255,166/255),rgb(221/255,176/255,150/255))

cols=rev(colorRampPalette(brewer.pal(9,"PuBu"))(12))
pdf("Fig3B.9sp_6cluster_plot2.pdf",height=1.5,width=10)
par(mfrow=c(1,6))
par(mar=c(1,2,2,1))
#plot(hc)
plot_9sp(Mm_cluster,sign_mat)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
#legend("topright",c("Mm","Gg","Ps","Xl","Xt","Dr","Bf","Ci","Cg"),col=cols[1:9],pch=19,cex=0.8)
dev.off()











