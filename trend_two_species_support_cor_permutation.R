args<-commandArgs(TRUE)
sp1<-args[1]
sp2<-args[2]
time<-args[3]


dir<-'/picb/compbio.work/Xuchuan/2015_Dec_linearship/poly/corr_0.7_age_model_both_age_related/'

r_cutoff<-0.7
Mm<-17;Gg<-13;Ps<-14;Xl<-17;Xt<-17;Dr<-15;Bf<-11;Ci<-15;Cg<-19
len1<-get(sp1);len2<-get(sp2)
len1_whole<-2:len1;len2_seq<-seq(1,len2,length.out=20)
data_load<-paste(dir,sp1,'_',sp2,'_age_orth_expression.Rdata',sep='');load(data_load)
func<-function(.x)
 {
 #fit1<-smooth.spline(1:len1,.x[1:len1],df=6)
 y_1<-.x[1:len1];x_1<-1:len1;fit1<-lm(y_1~I(x_1)+I(x_1^2)+I(x_1^3)+I(x_1^4) )
 #fit2<-smooth.spline(1:len2,.x[(len1+1):(len1+len2)],df=ifelse(sp2=='Bf',3,6));sp2_inter_y<-predict(fit2,len2_seq)$y
 y_2<-.x[(len1+1):(len1+len2)];x_2<-1:len2;fit2<-lm(y_2~I(x_2)+I(x_2^2)+I(x_2^3)+I(x_2^4) );sp2_inter_y<-predict(fit2,data.frame(x_2=len2_seq))
 cor_parts<-c()
 for( len1_single in len1_whole)
  {
  len1_seq_partly<-seq(1,len1_single,length.out=20);
  sp1_inter_partly_y<-predict(fit1,data.frame(x_1=len1_seq_partly))
  cor_part<-cor(sp2_inter_y,sp1_inter_partly_y)
  cor_parts<-c(cor_parts,cor_part)
  }
 if( all(cor_parts<r_cutoff) ){flag<-0}else{flag<-which.max(cor_parts)+1}
 flag
 }
rows<-nrow(t)
t_per<-data.frame(t[,1:len1],t[sample(rows),(len1+1):(len1+len2)])
temp<-apply(t_per,1,func)
result<-table(temp)
save(temp,result,file=paste(dir,'permute_data/','cor_permutation_',sp2,'_to_',sp1,'_',time,'.Rdata',sep=''))
