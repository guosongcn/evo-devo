
load("../data/fromstart.alignmentd-1.order.Rdata")

load("../data/9species.coding.rpkm1.meanstage.Rdata")
 
 Mm=Mm_mat_cufflink_coding_rpkm1_mean
 Bf=Bf_mat_cufflink_coding_rpkm1_mean[,Mm_Bf[,2]]
 Bf[is.na(Bf)]<-0
 Bf=Bf[apply(Bf,1,sum)>0,]


Ci=Ci_mat_cufflink_coding_rpkm1_mean[,Mm_Ci[,2]]
Ci[is.na(Ci)]<-0
 Ci=Ci[apply(Ci,1,sum)>0,]

Cg=Cg_mat_cufflink_coding_rpkm1_mean[,Mm_Cg[,2]]
Cg[is.na(Cg)]<-0
 Cg=Cg[apply(Cg,1,sum)>0,]
Ps=Ps_mat_cufflink_coding_rpkm1_mean[,Mm_Ps[,2]]
Ps[is.na(Ps)]<-0
Ps=Ps[apply(Ps,1,sum)>0,]
Gg=Gg_mat_cufflink_coding_rpkm1_mean[,Mm_Gg[,2]]
Gg[is.na(Gg)]<-0
 Gg=Gg[apply(Gg,1,sum)>0,]

Dr=Dr_mat_cufflink_coding_rpkm1_mean[,Mm_Dr[,2]]
Dr[is.na(Dr)]<-0
Dr=Dr[apply(Dr,1,sum)>0,]
Xl=Xl_mat_cufflink_coding_rpkm1_mean[,Mm_Xl[,2]]
Xl[is.na(Xl)]<-0
Xl=Xl[apply(Xl,1,sum)>0,]
Xt=Xt_mat_cufflink_coding_rpkm1_mean[,Mm_Xt[,2]]
Xt[is.na(Xt)]<-0
Xt=Xt[apply(Xt,1,sum)>0,]
norm=function(x)
{
   y=(x-mean(x))/sd(x)
   
}


predict_mat=function(mat,real,predict_point)
{
    mat_pre=t(apply(mat,1,function(x){
        y=predict(smooth.spline(real,x,df=10),predict_point)$y
        
    }))
    
    return(mat_pre)
}

 predict_mat_norm=function(mat,real,predict_point)
 {
    mat_pre=t(apply(mat,1,function(x){
        y=norm(predict(smooth.spline(real,x,df=10),predict_point)$y)
        
    }))
 
    return(mat_pre)
 }
 
 Mm_pre=predict_mat(Mm,1:17,seq(1,17,by=0.5))
 Bf_pre=predict_mat(Bf,Mm_Bf[,1],seq(1,17,by=0.5))
 Ci_pre=predict_mat(Ci,Mm_Ci[,1],seq(1,17,by=0.5))
 Cg_pre=predict_mat(Cg,Mm_Cg[,1],seq(1,17,by=0.5))
 Ps_pre=predict_mat(Ps,Mm_Ps[,1],seq(1,17,by=0.5))
 Gg_pre=predict_mat(Gg,Mm_Gg[,1],seq(1,17,by=0.5))
 Xl_pre=predict_mat(Xl,Mm_Xl[,1],seq(1,17,by=0.5))
Xt_pre=predict_mat(Xt,Mm_Xt[,1],seq(1,17,by=0.5))
Dr_pre=predict_mat(Dr,Mm_Dr[,1],seq(1,17,by=0.5))

colnames(Bf_pre)=paste("Bf_",1:33,sep="")
colnames(Ci_pre)=paste("Ci_",1:33,sep="")
colnames(Cg_pre)=paste("Cg_",1:33,sep="")
colnames(Ps_pre)=paste("Ps_",1:33,sep="")
colnames(Gg_pre)=paste("Gg_",1:33,sep="")
colnames(Xt_pre)=paste("Xt_",1:33,sep="")
colnames(Xl_pre)=paste("Xl_",1:33,sep="")
colnames(Dr_pre)=paste("Dr_",1:33,sep="")
colnames(Mm_pre)=paste("Mm_",1:33,sep="")

rpkm1=function(mat)
{
    
    mat_max=apply(mat,1,max)
    mat_0=mat[rownames(mat) %in% names(mat_max[mat_max>1]),]
    
    return(mat_0)
}


Mm_pre=rpkm1(Mm_pre)
Gg_pre=rpkm1(Gg_pre)
Ps_pre=rpkm1(Ps_pre)
Xl_pre=rpkm1(Xl_pre)
Xt_pre=rpkm1(Xt_pre)
Dr_pre=rpkm1(Dr_pre)
Bf_pre=rpkm1(Bf_pre)
Ci_pre=rpkm1(Ci_pre)
Cg_pre=rpkm1(Cg_pre)

all_species=rbind(Mm_pre,Gg_pre,Ps_pre,Xl_pre,Xt_pre,Dr_pre,Bf_pre,Ci_pre,Cg_pre)
#save(all_species,file="traceback_all_sepcies_combine_renorm.Rdata")
#save(all_species,file="fromstart_needleman_all_sepcies_combine_renorm.Rdata")
save(all_species,file="../data/predicte2Mm_all_sepcies_combine_raw_df10.Rdata")



Mm_pre_norm=predict_mat_norm(Mm,1:17,seq(1,17,by=0.5))
Bf_pre_norm=predict_mat_norm(Bf,Mm_Bf[,1],seq(1,17,by=0.5))
Ci_pre_norm=predict_mat_norm(Ci,Mm_Ci[,1],seq(1,17,by=0.5))
Cg_pre_norm=predict_mat_norm(Cg,Mm_Cg[,1],seq(1,17,by=0.5))
Ps_pre_norm=predict_mat_norm(Ps,Mm_Ps[,1],seq(1,17,by=0.5))
Gg_pre_norm=predict_mat_norm(Gg,Mm_Gg[,1],seq(1,17,by=0.5))
Xl_pre_norm=predict_mat_norm(Xl,Mm_Xl[,1],seq(1,17,by=0.5))
Xt_pre_norm=predict_mat_norm(Xt,Mm_Xt[,1],seq(1,17,by=0.5))
Dr_pre_norm=predict_mat_norm(Dr,Mm_Dr[,1],seq(1,17,by=0.5))

colnames(Bf_pre_norm)=paste("Bf_",1:33,sep="")
colnames(Ci_pre_norm)=paste("Ci_",1:33,sep="")
colnames(Cg_pre_norm)=paste("Cg_",1:33,sep="")
colnames(Ps_pre_norm)=paste("Ps_",1:33,sep="")
colnames(Gg_pre_norm)=paste("Gg_",1:33,sep="")
colnames(Xt_pre_norm)=paste("Xt_",1:33,sep="")
colnames(Xl_pre_norm)=paste("Xl_",1:33,sep="")
colnames(Dr_pre_norm)=paste("Dr_",1:33,sep="")
colnames(Mm_pre_norm)=paste("Mm_",1:33,sep="")


Mm_pre_norm=Mm_pre_norm[rownames(Mm_pre_norm) %in% rownames(Mm_pre),]
Gg_pre_norm=Gg_pre_norm[rownames(Gg_pre_norm) %in% rownames(Gg_pre),]
Ps_pre_norm=Ps_pre_norm[rownames(Ps_pre_norm) %in% rownames(Ps_pre),]
Xl_pre_norm=Xl_pre_norm[rownames(Xl_pre_norm) %in% rownames(Xl_pre),]
Xt_pre_norm=Xt_pre_norm[rownames(Xt_pre_norm) %in% rownames(Xt_pre),]
Dr_pre_norm=Dr_pre_norm[rownames(Dr_pre_norm) %in% rownames(Dr_pre),]
Bf_pre_norm=Bf_pre_norm[rownames(Bf_pre_norm) %in% rownames(Bf_pre),]
Ci_pre_norm=Ci_pre_norm[rownames(Ci_pre_norm) %in% rownames(Ci_pre),]
Cg_pre_norm=Cg_pre_norm[rownames(Cg_pre_norm) %in% rownames(Cg_pre),]

all_species_norm=rbind(Mm_pre_norm,Gg_pre_norm,Ps_pre_norm,Xl_pre_norm,Xt_pre_norm,Dr_pre_norm,Bf_pre_norm,Ci_pre_norm,Cg_pre_norm)
#save(all_species,file="traceback_all_sepcies_combine_renorm.Rdata")
#save(all_species,file="fromstart_needleman_all_sepcies_combine_renorm.Rdata")
save(all_species_norm,file="../data/predicte2Mm_all_sepcies_combine_norm_df10.Rdata")



load("/picb/compbio.work/IRIE_RNAseq/Analysis/2015_08/new_dataset_expression/orth3_9_ortholog_gene_list.Rdata")
Mm_pre=data.frame(Mm=rownames(Mm_pre),Mm_pre)
 Bf_pre=data.frame(Bf=rownames(Bf_pre),Bf_pre)
 Ci_pre=data.frame(Ci=rownames(Ci_pre),Ci_pre)
 Cg_pre=data.frame(Cg=rownames(Cg_pre),Cg_pre)
 Dr_pre=data.frame(Dr=rownames(Dr_pre),Dr_pre)
 Ps_pre=data.frame(Ps=rownames(Ps_pre),Ps_pre)
 Gg_pre=data.frame(Gg=rownames(Gg_pre),Gg_pre)
Xl_pre=data.frame(Xl=rownames(Xl_pre),Xl_pre)
Xt_pre=data.frame(Xt=rownames(Xt_pre),Xt_pre)

 
orth9_mat_pre=merge(orth9,Mm_pre,by.x="Mm",by.y="Mm")
orth9_mat_pre=merge(orth9_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth9_mat_pre=merge(orth9_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth9_mat_pre=merge(orth9_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth9_mat_pre=merge(orth9_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth9_mat_pre=merge(orth9_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth9_mat_pre=merge(orth9_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")
orth9_mat_pre=merge(orth9_mat_pre,Ci_pre,by.x="Ci",by.y="Ci")
orth9_mat_pre=merge(orth9_mat_pre,Cg_pre,by.x="Cg",by.y="Cg")
 
 orth8_mat_pre=merge(orth8,Mm_pre,by.x="Mm",by.y="Mm")
orth8_mat_pre=merge(orth8_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth8_mat_pre=merge(orth8_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth8_mat_pre=merge(orth8_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth8_mat_pre=merge(orth8_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth8_mat_pre=merge(orth8_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth8_mat_pre=merge(orth8_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")
orth8_mat_pre=merge(orth8_mat_pre,Ci_pre,by.x="Ci",by.y="Ci")


 orth7_mat_pre=merge(orth7,Mm_pre,by.x="Mm",by.y="Mm")
orth7_mat_pre=merge(orth7_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth7_mat_pre=merge(orth7_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth7_mat_pre=merge(orth7_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth7_mat_pre=merge(orth7_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth7_mat_pre=merge(orth7_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth7_mat_pre=merge(orth7_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")

 #sign_mat=orth9_mat_pre[,10:306]
 #rownames(sign_mat)=orth9_mat_pre$Mm
 #sign_mat[is.na(sign_mat)]<-0
 #sp9_cluster12=k_means(sign_mat,12)
 
 
 
orth6_mat_pre=merge(orth6,Mm_pre,by.x="Mm",by.y="Mm")
orth6_mat_pre=merge(orth6_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth6_mat_pre=merge(orth6_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth6_mat_pre=merge(orth6_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth6_mat_pre=merge(orth6_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth6_mat_pre=merge(orth6_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")


orth5_mat_pre=merge(orth5,Mm_pre,by.x="Mm",by.y="Mm")
orth5_mat_pre=merge(orth5_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth5_mat_pre=merge(orth5_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth5_mat_pre=merge(orth5_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth5_mat_pre=merge(orth5_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")



orth3_mat_pre=merge(orth3,Mm_pre,by.x="Mm",by.y="Mm")
orth3_mat_pre=merge(orth3_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth3_mat_pre=merge(orth3_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")

#save(orth3_mat_pre,orth5_mat_pre,orth6_mat_pre,orth7_mat_pre,orth8_mat_pre,orth9_mat_pre,Mm_pre,file="traceback_other_to_Mm_by_prediction_expression.Rdata")

#save(orth3_mat_pre,orth5_mat_pre,orth6_mat_pre,orth7_mat_pre,orth8_mat_pre,orth9_mat_pre,Mm_pre,file="fromstart_needleman_other_to_Mm_by_prediction_expression_renorm.Rdata")

save(orth3_mat_pre,orth5_mat_pre,orth6_mat_pre,orth7_mat_pre,orth8_mat_pre,orth9_mat_pre,Mm_pre,file="../data/orth3_9_other_to_Mm_by_prediction_expression_raw_df10.Rdata")








####norm


load("/picb/compbio.work/IRIE_RNAseq/Analysis/2015_08/new_dataset_expression/orth3_9_ortholog_gene_list.Rdata")
Mm_pre=data.frame(Mm=rownames(Mm_pre_norm),Mm_pre_norm)
Bf_pre=data.frame(Bf=rownames(Bf_pre_norm),Bf_pre_norm)
Ci_pre=data.frame(Ci=rownames(Ci_pre_norm),Ci_pre_norm)
Cg_pre=data.frame(Cg=rownames(Cg_pre_norm),Cg_pre_norm)
Dr_pre=data.frame(Dr=rownames(Dr_pre_norm),Dr_pre_norm)
Ps_pre=data.frame(Ps=rownames(Ps_pre_norm),Ps_pre_norm)
Gg_pre=data.frame(Gg=rownames(Gg_pre_norm),Gg_pre_norm)
Xl_pre=data.frame(Xl=rownames(Xl_pre_norm),Xl_pre_norm)
Xt_pre=data.frame(Xt=rownames(Xt_pre_norm),Xt_pre_norm)


orth9_mat_pre=merge(orth9,Mm_pre,by.x="Mm",by.y="Mm")
orth9_mat_pre=merge(orth9_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth9_mat_pre=merge(orth9_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth9_mat_pre=merge(orth9_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth9_mat_pre=merge(orth9_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth9_mat_pre=merge(orth9_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth9_mat_pre=merge(orth9_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")
orth9_mat_pre=merge(orth9_mat_pre,Ci_pre,by.x="Ci",by.y="Ci")
orth9_mat_pre_norm=merge(orth9_mat_pre,Cg_pre,by.x="Cg",by.y="Cg")

orth8_mat_pre=merge(orth8,Mm_pre,by.x="Mm",by.y="Mm")
orth8_mat_pre=merge(orth8_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth8_mat_pre=merge(orth8_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth8_mat_pre=merge(orth8_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth8_mat_pre=merge(orth8_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth8_mat_pre=merge(orth8_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth8_mat_pre=merge(orth8_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")
orth8_mat_pre_norm=merge(orth8_mat_pre,Ci_pre,by.x="Ci",by.y="Ci")


orth7_mat_pre=merge(orth7,Mm_pre,by.x="Mm",by.y="Mm")
orth7_mat_pre=merge(orth7_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth7_mat_pre=merge(orth7_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth7_mat_pre=merge(orth7_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth7_mat_pre=merge(orth7_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth7_mat_pre=merge(orth7_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")
orth7_mat_pre_norm=merge(orth7_mat_pre,Bf_pre,by.x="Bf",by.y="Bf")

#sign_mat=orth9_mat_pre[,10:306]
#rownames(sign_mat)=orth9_mat_pre$Mm
#sign_mat[is.na(sign_mat)]<-0
#sp9_cluster12=k_means(sign_mat,12)



orth6_mat_pre=merge(orth6,Mm_pre,by.x="Mm",by.y="Mm")
orth6_mat_pre=merge(orth6_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth6_mat_pre=merge(orth6_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth6_mat_pre=merge(orth6_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth6_mat_pre=merge(orth6_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")
orth6_mat_pre_norm=merge(orth6_mat_pre,Dr_pre,by.x="Dr",by.y="Dr")


orth5_mat_pre=merge(orth5,Mm_pre,by.x="Mm",by.y="Mm")
orth5_mat_pre=merge(orth5_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth5_mat_pre=merge(orth5_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")
orth5_mat_pre=merge(orth5_mat_pre,Xl_pre,by.x="Xl",by.y="Xl")
orth5_mat_pre_norm=merge(orth5_mat_pre,Xt_pre,by.x="Xt",by.y="Xt")



orth3_mat_pre=merge(orth3,Mm_pre,by.x="Mm",by.y="Mm")
orth3_mat_pre=merge(orth3_mat_pre,Gg_pre,by.x="Gg",by.y="Gg")
orth3_mat_pre_norm=merge(orth3_mat_pre,Ps_pre,by.x="Ps",by.y="Ps")

#save(orth3_mat_pre,orth5_mat_pre,orth6_mat_pre,orth7_mat_pre,orth8_mat_pre,orth9_mat_pre,Mm_pre,file="traceback_other_to_Mm_by_prediction_expression.Rdata")

#save(orth3_mat_pre,orth5_mat_pre,orth6_mat_pre,orth7_mat_pre,orth8_mat_pre,orth9_mat_pre,Mm_pre,file="fromstart_needleman_other_to_Mm_by_prediction_expression_renorm.Rdata")

save(orth3_mat_pre_norm,orth5_mat_pre_norm,orth6_mat_pre_norm,orth7_mat_pre_norm,orth8_mat_pre_norm,orth9_mat_pre_norm,Mm_pre_norm,file="../data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")

