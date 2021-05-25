source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")

library(RColorBrewer)


load("../Data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")
load("../Data/age_related_gene_list_BH0.05.Rdata")
load("../Data/needleman_fromstat_d1_sp9_cluster6_hclust_age_related_df10.Rdata")


cor_mean=function(Mm_cluster,sign_mat)
{
    cors_table=matrix(ncol=8,nrow=0)
    
    for(ii in 1:length(Mm_cluster))
    {
        gene=Mm_cluster[[ii]]
        
        gene_mat=sign_mat[rownames(sign_mat) %in% gene,]
        cors_mat=matrix(ncol=8,nrow=0)
        for(j in 1:nrow(gene_mat))
        {
            gene_means=as.numeric(gene_mat[j,])
            Mm=gene_means[1:33]
            Gg=gene_means[34:66]
            Ps=gene_means[67:99]
            Xl=gene_means[100:132]
            Xt=gene_means[133:165]
            Dr=gene_means[166:198]
            Bf=gene_means[199:231]
            Ci=gene_means[232:264]
            Cg=gene_means[265:297]
            cors=c(cor(Mm,Gg,method="spearman"),cor(Mm,Ps,method="spearman"),cor(Mm,Xl,method="spearman"),cor(Mm,Xt,method="spearman"),cor(Mm,Dr,method="spearman"),cor(Mm,Bf,method="spearman"),cor(Mm,Ci,method="spearman"),cor(Mm,Cg,method="spearman"))
            cors_mat=rbind(cors_mat,cors)
        }
        cors_mat_mean=apply(cors_mat,2,mean)
        cors_table=rbind(cors_table,cors_mat_mean)
    }
    return(cors_table)
    
}

sign_mat=orth9_mat_pre_norm[orth9_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth9_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth9_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth9_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth9_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth9_mat_pre_norm$Dr %in% age_related_gene_list$Dr & orth9_mat_pre_norm$Bf %in% age_related_gene_list$Bf & orth9_mat_pre_norm$Ci %in% age_related_gene_list$Ci & orth9_mat_pre_norm$Cg %in% age_related_gene_list$Cg,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,10:ncol(sign_mat)]

sign_mat[is.na(sign_mat)]<-0

s=cor_mean(sp9_cluster6,sign_mat)
Mm_sep_time=c(312,312,351.8,351.8,435,684,676,797)




s=cor_mean(sp9_cluster6,sign_mat)


pdf("FigureS.cluster.dist2.pdf",height=3,width=18)
par(mfrow=c(1,6))
par(mar=c(4,3,2,1))
for(i in 1:6)
{
    plot(Mm_sep_time,s[i,],main="",pch=16,cex=2,xlab="",ylab="",col="blue",ylim=c(-0.5,1),xlim=c(300,800))
    abline(lm(s[i,]~Mm_sep_time),col="grey",lty=2,lwd=2)
       text(630,0.8,paste("rho = ",round(cor(Mm_sep_time,s[i,],method="spearman"),2),sep=""),cex=1.5)
}



dev.off()



pdf("FigureS.1cluster.dist2.pdf",height=4,width=4)
i=1   

plot(Mm_sep_time,s[i,],main="",pch=16,cex=2,xlab="",ylab="",col="darkred",ylim=c(0.2,0.7),xlim=c(300,800))
    abline(lm(s[i,]~Mm_sep_time),col="grey",lty=2,lwd=2)
       text(630,0.8,paste("rho = ",round(cor(Mm_sep_time,s[i,],method="spearman"),2),sep=""),cex=1.5)




dev.off()



