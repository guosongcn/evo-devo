load("../Data/9species.cufflink.expression.unique.coding.maxrpkm1.withoutadult.Rdata")
load("../Data/sample_reads_count_from_naoki.Rdata")

##MM
stage=sapply(colnames(Mm_mat_cufflink_coding_rpkm1),function(x){substr(x,1,nchar(x)-2)})
rep=sapply(colnames(Mm_mat_cufflink_coding_rpkm1),function(x){substr(x,nchar(x),nchar(x)+1)})
coverage=mm_reads[,2]

variance_prop_mat=matrix(ncol=4,nrow=0)
for(i in 1:nrow(Mm_mat_cufflink_coding_rpkm1))
{
    x=as.numeric(Mm_mat_cufflink_coding_rpkm1[i,])
    each=anova(lm(x~stage+rep+coverage))
    total=var(x)*(length(x)-1)
    prop=each[[2]]/total
    variance_prop_mat=rbind(variance_prop_mat,prop)
}

mean_mm=apply(variance_prop_mat,2,mean)
var_mm=apply(variance_prop_mat,2,var)



##Gg
stage=sapply(colnames(Gg_mat_cufflink_coding_rpkm1),function(x){substr(x,1,nchar(x)-2)})
rep=sapply(colnames(Gg_mat_cufflink_coding_rpkm1),function(x){substr(x,nchar(x),nchar(x)+1)})
coverage=gg_reads[,2]

variance_prop_mat=matrix(ncol=4,nrow=0)
for(i in 1:nrow(Gg_mat_cufflink_coding_rpkm1))
{
    x=as.numeric(Gg_mat_cufflink_coding_rpkm1[i,])
    each=anova(lm(x~stage+rep+coverage))
    total=var(x)*(length(x)-1)
    prop=each[[2]]/total
    variance_prop_mat=rbind(variance_prop_mat,prop)
}

mean_gg=apply(variance_prop_mat,2,mean)
var_gg=apply(variance_prop_mat,2,var)



###
anova_test=function(Mat,reads_mat)
{
stage=sapply(colnames(Mat),function(x){substr(x,1,nchar(x)-2)})
rep=sapply(colnames(Mat),function(x){substr(x,nchar(x),nchar(x)+1)})
coverage=reads_mat[,2]

variance_prop_mat=matrix(ncol=4,nrow=0)
for(i in 1:nrow(Mat))
{
    x=as.numeric(Mat[i,])
    each=anova(lm(x~stage+rep+coverage))
    total=var(x)*(length(x)-1)
    prop=each[[2]]/total
    variance_prop_mat=rbind(variance_prop_mat,prop)
}
    return(variance_prop_mat)
}

mm_aov=anova_test(Mm_mat_cufflink_coding_rpkm1,mm_reads)
gg_aov=anova_test(Gg_mat_cufflink_coding_rpkm1,gg_reads)

ps_aov=anova_test(Ps_mat_cufflink_coding_rpkm1,ps_reads)

xl_aov=anova_test(Xl_mat_cufflink_coding_rpkm1,xl_reads)
Xt_aov=anova_test(Xt_mat_cufflink_coding_rpkm1,xt_reads)
dr_aov=anova_test(Dr_mat_cufflink_coding_rpkm1,dr_reads[1:32,])
bf_aov=anova_test(Bf_mat_cufflink_coding_rpkm1,bf_reads[1:23,])
ci_aov=anova_test(Ci_mat_cufflink_coding_rpkm1,ci_reads[1:48,])
cg_aov=anova_test(Cg_mat_cufflink_coding_rpkm1,cg_reads)


mean_mm=apply(mm_aov,2,mean)
sd_mm=apply(mm_aov,2,sd)

mean_gg=apply(gg_aov,2,mean)
sd_gg=apply(gg_aov,2,sd)

mean_ps=apply(ps_aov,2,mean)
sd_ps=apply(ps_aov,2,sd)


mean_xl=apply(xl_aov,2,mean)
sd_xl=apply(xl_aov,2,sd)

mean_xt=apply(Xt_aov,2,mean)
sd_xt=apply(Xt_aov,2,sd)

mean_dr=apply(dr_aov,2,mean)
sd_dr=apply(dr_aov,2,sd)

mean_bf=apply(bf_aov,2,mean)
sd_bf=apply(bf_aov,2,sd)

mean_ci=apply(ci_aov,2,mean)
sd_ci=apply(ci_aov,2,sd)


Mat=Cg_mat_cufflink_coding_rpkm1
stage=c(colnames(Mat)[1:10],rep("Cg_T",5),rep("Cg_ED",2),rep("Cg_D",7),rep("Cg_EU",2),rep("Cg_U",6),rep("Cg_LU",2),rep("Cg_P",2),"Cg_S","Cg_J")
rep=c(rep(1,10),1:5,1:2,1:7,1:2,1:6,1:2,1:2,1,1)
coverage=cg_reads[,2]

variance_prop_mat=matrix(ncol=4,nrow=0)
for(i in 1:nrow(Cg_mat_cufflink_coding_rpkm1))
{
    x=as.numeric(Cg_mat_cufflink_coding_rpkm1[i,])
    each=anova(lm(x~stage+rep+coverage))
    total=var(x)*(length(x)-1)
    prop=each[[2]]/total
    variance_prop_mat=rbind(variance_prop_mat,prop)
}
mean_cg=apply(variance_prop_mat,2,mean)
sd_cg=apply(variance_prop_mat,2,sd)
var_cg=apply(variance_prop_mat,2,var)


mean_all=rbind(mean_mm,mean_gg,mean_ps,mean_xl,mean_xt,mean_dr,mean_bf,mean_ci,mean_cg)
colnames(mean_all)=c("stage","replicate","coverage","residule")
var_all=rbind(var_mm,var_gg,var_ps,var_xl,var_xt,var_dr,var_bf,var_ci,var_cg)
colnames(var_all)=c("stage","replicate","coverage","residule")
sd_all=rbind(sd_mm,sd_gg,sd_ps,sd_xl,sd_xt,sd_dr,sd_bf,sd_ci,sd_cg)
colnames(sd_all)=c("stage","replicate","coverage","residule")

superpose.eb = function (x, y, ebl, ebu = ebl, lh = 0.01, ...) {
    segments(x, y + ebu, x, y - ebl, ...)
    segments(x - lh , y + ebu, x + lh , y + ebu,
    ...)
    segments(x - lh, y - ebl, x + lh , y - ebl, ...)
}

pdf("Figure.variance.pdf")
par(mar = c(5, 5, 4, 2) + 0.1)
bp=boxplot(mean_all,col="royalblue2",border=F,las=2,ylab="Percentage of variance %",xaxt="n")
text(c(1:4),c(0.7,0.3,0.2,0.4),c("stage","replicate","coverage","residule"))
points(rep(1,9),mean_all[,1],pch=19,col="royalblue4")
points(rep(2,9),mean_all[,2],pch=19,col="royalblue4")
points(rep(3,9),mean_all[,3],pch=19,col="royalblue4")
points(rep(4,9),mean_all[,4],pch=19,col="royalblue4")

dev.off()
