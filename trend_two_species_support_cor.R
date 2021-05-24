sp1<-'Mm'
sps<-c('Gg','Ps','Xl','Xt','Dr','Bf','Ci','Cg')

r_cutoff<-0.7
Mm<-17;Gg<-13;Ps<-14;Xl<-17;Xt<-17;Dr<-15;Bf<-11;Ci<-15;Cg<-19

for(sp2 in sps)
	{
len1<-get(sp1);len2<-get(sp2)
len1_whole<-2:len1;len2_seq<-seq(1,len2,length.out=20)
data_load<-paste(sp1,'_',sp2,'_','age_orth_expression.Rdata',sep='');load(data_load)
func<-function(.x)
 {
 ##fit1<-smooth.spline(1:len1,.x[1:len1],df=6)
  y_1<-.x[1:len1];x_1<-1:len1;fit1<-lm(y_1~I(x_1)+I(x_1^2)+I(x_1^3)+I(x_1^4) )
 ##fit2<-smooth.spline(1:len2,.x[(len1+1):(len1+len2)],df=ifelse(sp2=='Bf',3,6));sp2_inter_y<-predict(fit2,len2_seq)$y
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
 c(flag,max(cor_parts))
 }
katy<-t(apply(t,1,func))
temp<-katy[,1];names(temp)<-rownames(t)
corr<-katy[,2];names(corr)<-rownames(t)
result<-table(temp)
save(corr,temp,result,file=paste('./real_number/',sp2,'_to_',sp1,'.Rdata',sep=''))
  }
