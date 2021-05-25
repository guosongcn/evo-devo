load("orth_pairwise_gene_list.Rdata")
orth_list=ls()
d=grep("Dm",orth_list)
orth_list=orth_list[-d]
#orth_list=orth_list[grep("Mm",orth_list)]

library(RColorBrewer) #to use brewer.pal
library(fields)

trans_mat=function(mat)
{
    mat_new=mat[,rev(1:ncol(mat))]
    mat_new=mat_new[rev(1:nrow(mat_new)),]
    return(mat_new)
    
    
}

polygon_heatmap <- function (x, y, unitcell = 0.5, col = col)
{
    polygon(c(x, x+1, x+1, x,x)-unitcell, c(y,y  ,y +1 ,y +1,y )-unitcell,col = col, border=NA)
}

result=matrix(ncol=5,nrow=0)
#  pdf(paste("../figure/fromstart.pairwise.alignment.plot.rpkm2.d.",d,".pdf",sep=""))
pdf("FigS4.heatmap.path1.rpkm1.pdf",width=16,height=16)
par(mfrow=c(6,6))
par(mar=c(2,2,2,2))

for(ni in  orth_list)
{
    #ni="Mm_Gg"
    load(paste(ni,"_mat_associate_matrix_rpkm1.Rdata",sep=""))
    mat_associate_matrix=trans_mat(mat_associate_matrix)
    
    mat_associate_matrix[mat_associate_matrix>1000]<-500
    sp1=unlist(strsplit(ni,"_"))[1]
    sp2=unlist(strsplit(ni,"_"))[2]
    colnames(mat_associate_matrix)=paste(sp2,colnames(mat_associate_matrix),sep="_")
    rownames(mat_associate_matrix)=paste(sp1,rownames(mat_associate_matrix),sep="_")
    
    m=max(mat_associate_matrix[,ncol(mat_associate_matrix)])
    n=max(mat_associate_matrix[nrow(mat_associate_matrix),])
    if(m>=n)
    {
        S=mat_associate_matrix
        
    }else
    {
        S=t(mat_associate_matrix)
        tmp=sp1
        sp1=sp2
        sp2=tmp
    }
    
    
    #color heatmap plot
    S_scale=t(apply(S,1,scale))
    x <- as.vector(S_scale)
    ColRamp <- designer.colors(n=50, col=brewer.pal(9, "Blues"))
    
    ColorCode <- rep("#FFFFFF", length(x)) #default is all white
    Bins <- seq(min(x, na.rm=T), max(x, na.rm=T), length=length(ColRamp))
    for (ii in 1:length(x))
    {
        if (!is.na(x[ii]))
        {
            ColorCode[ii] <- ColRamp[which.min(abs(Bins-x[ii]))]
        }
    }
    SOM_Rows <- dim(S)[1]
    SOM_Columns <- dim(S)[2]
    
    
    ###color heatmap plot
    
    S_max_total=sum(apply(S,2,max))
    d=(-1)
    F=matrix(nrow=nrow(S)+1,ncol=ncol(S)+1)
    # F[,1]<-seq(0,nrow(S))*d
    #F[1,]<-seq(0,ncol(S))*d
    F[,1]=rep(-1,nrow(S)+1)
    F[1,]=rep(-1,ncol(S)+1)
    
    for(i in 2:nrow(F))
    {
        for(j in 2:ncol(F))
        {
            
            match=F[i-1,j-1]+S[i-1,j-1]
            delete=F[i-1,j]+d
            insert=F[i,j-1]+d
            F[i,j]<-max(c(match,delete,insert))
            
        }
    }
    
    
    
    colnames(F)=c("-",colnames(S))
    
    rownames(F)=c("-",rownames(S))
    
    
    i <-nrow(F)
    j <-ncol(F)
    
    F1=matrix(nrow=nrow(S)+1,ncol=ncol(S)+1)
    F1[,1]<-seq(0,nrow(S))*d
    F1[1,]<-seq(0,ncol(S))*d
    for(i in 2:nrow(F1))
    {
        for(j in 2:ncol(F1))
        {
            
            
            F1[i,j]<-S[i-1,j-1]
            
        }
    }
    
    plot(0,bty='n',pch='',ylab='',xlab='',ylim=c(1,i),xlim=c(1,j),xaxt='n',yaxt='n',main=ni)
  #   sapply(1:i,function(x){text(1:j,rep(x,ncol(F1)),round(F1[i-(x-1),],2),cex=0.5)})
    
 #    axis(side=1,at=c(1:j),label=c("-",colnames(S)),las=3)
 #   axis(side=2,at=c(1:i),label=c(rev(rownames(S)),"-"),las=1)
    
    #polygon(c(1,2,2,1,1)-0.5,c(1,1,2,2,1)-0.5,col="grey",border=T)
    
    
    
    
    for(column in 1:SOM_Columns)
    {
        for(row in 1:SOM_Rows)
        {
            polygon_heatmap(column+1,SOM_Rows-row+1,col=ColorCode[(column-1)*SOM_Rows+row])
        }
    }
    
       i<-as.numeric(which.max(S[,ncol(S)])+1)
    n<-1+(nrow(F)-i)
    AlignmentA <- rownames(F)[i]
    AlignmentB <- colnames(F)[j]
    A=1:nrow(S)
    B=1:ncol(S)
    orderA=i
    orderB=j
    while (i > 1 && j > 1)
    {
        
        Score<-F[i,j]
        ScoreDiag <- F[i - 1, j - 1]
        ScoreUp<-  F[i-1, j ]
        ScoreLeft <-  F[i, j-1]
        if (Score == ScoreDiag + S[i-1,j-1])
        {
            AlignmentA <- paste(rownames(F)[i-1],AlignmentA,sep="_")
            AlignmentB <- paste(colnames(F)[j-1],AlignmentB,sep="_")
            orderA=c(i-1,orderA)
            orderB=c(j-1,orderB)
            if((j-1)!=1 && (n+1) !=1)
            {
                segments(j,n,j-1,n+1)
            }
            n=n+1
            i <-i - 1
            j <- j - 1
            
            
        }
        else if (Score == ScoreLeft+d)
        {
            AlignmentA <- paste(rownames(F)[i], AlignmentA,sep="_")
            # AlignmentB <- paste("-" , AlignmentB,sep="_")
            AlignmentB <- paste(colnames(F)[j-1], AlignmentB,sep="_")
            
            orderA=c(i,orderA)
            orderB=c(j-1,orderB)
            if((j-1)!=1 && n !=1)
            {
                segments(j,n,j-1,n)
            }
            
            j <- j - 1
            
            
            
            
            
        }else if(Score == ScoreUp+d )
        {
            AlignmentA <- paste(rownames(F)[i-1] , AlignmentA,sep="_")
            AlignmentB <- paste(colnames(F)[j], AlignmentB,sep="_")
            orderA=c(i-1,orderA)
            orderB=c(j,orderB)
            if(j!=1 && (n+1) !=1)
            {
                segments(j,n,j,n+1)
            }
            i <- i - 1
            
            #AlignmentA <- paste("-" , AlignmentA,sep="_")
            
            
            n=n+1
        }
        
        
        
    }
    
    A_order=orderA-1
    B_order=orderB-1
    A_order=A_order[-1]
    B_order=B_order[-1]
    A_order=rev(nrow(S)-(A_order-1))
    B_order=rev(ncol(S)-(B_order-1))
    order_matrix=rbind(A_order,B_order)
    rownames(order_matrix)=c(sp1,sp2)
    ss_total=0
    for(n in 1:ncol(order_matrix))
    {
        
        ss_total=ss_total+S[order_matrix[1,n],order_matrix[2,n]]
        
    }
    
    S_diff= S_max_total-ss_total
    dd=c(d,ni,S_max_total,ss_total,S_diff)
    result=rbind(result,dd)
    # save(order_matrix,file=paste("../data/",ni,"_d",d,"_alignment_order.Rdata",sep=""))
}



dev.off()
