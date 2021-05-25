load("../data/orth_pairwise_gene_list.Rdata")
orth_list=ls()
d=grep("Dm",orth_list)
orth_list=orth_list[-d]
load("../data/1.mean_stage_associated_gene_matrix_rpkm2.Rdata")



for(ni in orth_list)
{
    name1=unlist(strsplit(ni,"_"))[1]
    name2=unlist(strsplit(ni,"_"))[2]
    
    orth=get(ni)
    colnames(orth)=c(name1,name2)
    mat1=get(paste(name1,"_stage_specific_matrix",sep=""))
    mat2=get(paste(name2,"_stage_specific_matrix",sep=""))
    mat1_length=ncol(mat1)
    mat2_length=ncol(mat2)
    
    mat1=mat1[apply(mat1,1,sum)>0,]
    mat2=mat2[apply(mat2,1,sum)>0,]
    mat1_r=data.frame(rownames(mat1),mat1)
    mat2_r=data.frame(rownames(mat2),mat2)
    colnames(mat1_r)[1]=name1
    colnames(mat2_r)[1]=name2
    mat=merge(orth,mat1_r,by.x=name1,by.y=name1)
    mat=merge(mat,mat2_r,by.x=name2,by.y=name2)
    rownames(mat)=mat[,1]
    mat=mat[,c(-1,-2)]

    mat_associate_matrix=matrix(nrow=mat1_length,ncol=mat2_length)
    colnames(mat_associate_matrix)=colnames(mat2)
    rownames(mat_associate_matrix)=colnames(mat1)

    p_mat1=matrix(nrow=mat1_length,ncol=mat2_length)
    colnames(p_mat1)=colnames(mat2)
    rownames(p_mat1)=colnames(mat1)
    for( m in 1:mat1_length)
    {
        for (n in (1+mat1_length):(mat1_length+mat2_length))
        {
        
        
            A=rownames(mat)[mat[,m]==1]
            B=rownames(mat)[mat[,n]==1]
            AB=intersect(A,B)
            all=nrow(mat)
            p= phyper(length(AB),length(B),all-length(B),length(A),lower.tail=F)
            p_BH=p*ncol(mat1)*ncol(mat2)
            p_mat1[m,n-mat1_length]= p
            mat_associate_matrix[m,n-mat1_length]=-log10(p_BH)
            #mat_associate_matrix[m,n-mat1_length]=-log10(p_BH)
        }

    }
    print(ni)
print(max(mat_associate_matrix))
save(mat_associate_matrix,file=paste("../data/",ni,"_mat_associate_matrix_rpkm2.Rdata",sep=""))
save(p_mat1,file=paste("../data/",ni,"_mat_Pvalue_matrix_rpkm2.Rdata",sep=""))



}





