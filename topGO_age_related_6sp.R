GO_test=function(test_id,all_gene)
{
    
    onts = c( "MF", "BP", "CC" )
    
    geneIDs = all_gene
    inUniverse = geneIDs %in% geneIDs
    inSelection =  geneIDs %in% test_id
    alg <- factor( as.integer( inSelection[inUniverse] ) )
    names(alg) <- geneIDs[inUniverse]
    
    
    
    ## ----run tests, results='hide'-------------------------------------------
    tab = as.list(onts)
    names(tab) = onts
    ### test all three top level ontologies
    for(i in 1:3){
        
        ## prepare data
        tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=10,
        annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
        
        ## run tests
        resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
        #resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
        allGO = usedGO(tgd)
        ## look at results
        all_results <- GenTable( tgd, Fisher.elim = resultTopGO.elim,topNodes=length(allGO))
        
        tab[[i]]=all_results[p.adjust(as.numeric(all_results$Fisher.elim),method="BH")<0.05,]
        #tab[[i]]=all_results[as.numeric(all_results$Fisher.elim)<0.01,]
        
        #res <- printGenes(tgd, whichTerms = all_results$GO.ID,geneCutOff = 1000)
        
    }
    
    
    # ----process results-----------------------------------------------------
    topGOResults <- rbind.fill(tab)
    return(topGOResults)
    
}

library(org.Mm.eg.db) # ENSEMBL genes
library(topGO)
library(plyr)

load("needleman_fromstat_d1_sp6_cluster6_hclust.Rdata")

load("../data/orth3_9_other_to_Mm_by_prediction_expression_norm_df10.Rdata")

load("../data/age_related_gene_list.Rdata")
sign_mat=orth6_mat_pre_norm[orth6_mat_pre_norm$Mm %in% age_related_gene_list$Mm & orth6_mat_pre_norm$Gg %in% age_related_gene_list$Gg & orth6_mat_pre_norm$Ps %in% age_related_gene_list$Ps & orth6_mat_pre_norm$Xl %in% age_related_gene_list$Xl & orth6_mat_pre_norm$Xt %in% age_related_gene_list$Xt & orth6_mat_pre_norm$Dr %in% age_related_gene_list$Dr,]
rownames(sign_mat)=sign_mat$Mm
sign_mat=sign_mat[,7:ncol(sign_mat)]



sign_mat[is.na(sign_mat)]<-0
all_gene=rownames(sign_mat)



sign_table=matrix(ncol=7,nrow=0)
sign_GO_cluster=list()
for(i in 1:length(sp6_cluster6))
{
    sign_go=GO_test(sp6_cluster6[[i]],all_gene)
    sign_GO_cluster[[i]]=sign_go
    sign_go=cbind(cluster=rep(i,nrow(sign_go)),sign_go)
    sign_table=rbind(sign_table,sign_go)
}

save(sign_GO_cluster,file="6sp_sign_GO_6cluster_age_related_BH.Rdata")
write.table(sign_table,file="6sp_sign_GO_6cluster_BH.csv",quote=F,col.names=T,row.names=F,sep=",")
