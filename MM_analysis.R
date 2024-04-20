myrnaMerge <- function (metadata, path, data.type, organized = FALSE)
{
  if (organized == TRUE) {
    filenames <- file.path(path, metadata$file_name, fsep = .Platform$file.sep)
  }
  else {
    filenames <- file.path(path, metadata$file_id, metadata$file_name,
                           fsep = .Platform$file.sep)
  }
  if (data.type == "RNAseq") {
    message("############### Merging RNAseq data ################\n",
            "### This step may take a few minutes ###\n")
    rnaMatrix <- do.call("cbind", lapply(filenames,
                                         function(fl) read.table(gzfile(fl), skip = 6)$V4))
    rownames(rnaMatrix) <- read.table(gzfile(filenames[1]), skip = 6)$V1
    rownames(rnaMatrix) <- unlist(lapply(strsplit(rownames(rnaMatrix),
                                                  ".", fixed = TRUE), function(gene) gene[1]))
    colnames(rnaMatrix) <- metadata$sample
    # rnaMatrix <- rnaMatrix[biotype$ensemblID, ]
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    message(paste("Number of samples: ", nSamples, "\n",
                  sep = ""), paste("Number of genes: ", nGenes, "\n",
                                   sep = ""))
    return(rnaMatrix)
  }
  else if (data.type == "pre-miRNAs") {
    message("############### Merging pre-miRNAs data ################\n",
            "### This step may take a few minutes ###\n")
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) read.delim(fl)$read_count))
    rownames(rnaMatrix) <- read.delim(filenames[1])$miRNA_ID
    colnames(rnaMatrix) <- metadata$sample
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    message(paste("Number of samples: ", nSamples, "\n",
                  sep = ""), paste("Number of genes: ", nGenes, "\n",
                                   sep = ""))
    return(rnaMatrix)
  }
  else if (data.type == "miRNAs") {
    message("############### Merging miRNAs data ###############\n")
    mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
    mirs <- rownames(mirbase)
    mirMatrix <- do.call("cbind", lapply(mirMatrix, function(expr) expr[mirs]))
    rownames(mirMatrix) <- mirbase$v21[match(mirs, rownames(mirbase))]
    colnames(mirMatrix) <- metadata$sample
    mirMatrix[is.na(mirMatrix)] <- 0
    nSamples = ncol(mirMatrix)
    nGenes = nrow(mirMatrix)
    message(paste("Number of samples: ", nSamples, "\n",
                  sep = ""), paste("Number of miRNAs: ", nGenes, "\n",
                                   sep = ""))
    return(mirMatrix)
  }
  else {
    return("error !!!")
  }
}

library(GDCRNATools)
library(tidyverse)
library(caret)
library(ggplot2)
options(stringsAsFactors = F)
setwd("D:\\MM")
####### Download RNAseq data #######
gdcRNADownload(project.id     = 'MMRF-COMMPASS', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               directory      = 'MMRF-COMMPASS/RNAseq')
####### Download clinical data #######
gdcClinicalDownload(project.id     = 'MMRF-COMMPASS',  
                    write.manifest = FALSE,
                    directory      = 'MMRF-COMMPASS/Clinical')
####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TARGET-AML', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               directory      = 'MMRF-COMMPASS/RNAseq')
####### Download clinical data #######
gdcClinicalDownload(project.id     = 'TARGET-AML',  
                    write.manifest = FALSE,
                    directory      = 'MMRF-COMMPASS/Clinical')

####### Parse RNAseq metadata #######
#metaMatrixMM.RNA <- gdcParseMetadata(project.id = 'MMRF-COMMPASS',
#                                      data.type  = 'RNAseq', 
#                                      write.meta = FALSE)
#write.table(metaMatrixMM.RNA,'metaMatrixMM.RNA.txt',row.names = T,col.names = T,sep='\t')
metaMatrixMM.RNA=read.table('metaMatrixMM.RNA.txt',header=T,row.names = 1,sep='\t')
metaMatrixMM.RNA <- gdcFilterDuplicate(metaMatrixMM.RNA)
metaMatrixMM1.RNA <-  metaMatrixMM.RNA[metaMatrixMM.RNA$sample_type=='PrimaryBloodDerivedCancer-BoneMarrow',]
metaMatrixMM1.RNA$sample_type=gsub('-','',metaMatrixMM1.RNA$sample_type)


#metaMatrixAML.RNA <- gdcParseMetadata(project.id = 'TARGET-AML',
#                                      data.type  = 'RNAseq', 
#                                      write.meta = FALSE)

#write.table(metaMatrixAML.RNA,'metaMatrixAML.RNA.txt',row.names = T,col.names = T,sep='\t')
metaMatrixAML.RNA=read.table('metaMatrixAML.RNA.txt',header=T,row.names = 1,sep='\t')
metaMatrixAML.RNA <- gdcFilterDuplicate(metaMatrixAML.RNA)
metaMatrixAML.RNA.normal=metaMatrixAML.RNA[metaMatrixAML.RNA$sample_type=="BoneMarrowNormal",]
rownames(metaMatrixAML.RNA.normal)=metaMatrixAML.RNA.normal$sample
metaMatrix.RNA=rbind(metaMatrixMM1.RNA,metaMatrixAML.RNA.normal)
write.table(metaMatrix.RNA,'MM.clinical.txt',row.names = T,col.names = T,sep='\t')

####### Merge RNAseq data #######
MMAMLrnaMatrix <- myrnaMerge(metadata  = metaMatrix.RNA, 
                           path      = 'MMRF-COMMPASS/RNAseq/', 
                           data.type = 'RNAseq')
rnaMatrix=MMAMLrnaMatrix

####### normalisation RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaMatrix, filter = FALSE)
head(rnaExpr)
#write.table(rnaExpr,'expression_both_MM_and_normal.txt',row.names = T,col.names = T,sep='\t')

DEGAll <- gdcDEAnalysis(counts     = rnaMatrix, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryBloodDerivedCancerBoneMarrow-BoneMarrowNormal', 
                        method     = 'limma')
sd=DEGAll$symbol
names(sd)=rownames(DEGAll)
#dmmxj=data.frame(gene=sd,ensemble=names(sd))
#write.table(dmmxj,'geneID_ensemble.txt',col.names = T,row.names = F,sep='\t')
DEGAll2=DEGAll[DEGAll$group=="protein_coding",]
#write.table(DEGAll2,'differential_genes_MM.txt',row.names = T,col.names = T,sep='\t')



f1=read.table('ubiqintin.curated.txt',header = F)
geneList1=unique(c(f1$V1))
geneList=geneList1
difcount2 <- DEGAll2
myco=c('grey','red','SteelBlue')
threshold= rep(0,nrow(DEGAll2))
threshold[DEGAll2$logFC>0.58&DEGAll2$PValue<0.05]=1
threshold[DEGAll2$logFC<(-0.58)&DEGAll2$PValue<0.05]=2
threshold=factor(threshold)
lim <- max(max(DEGAll2$logFC), abs(min(DEGAll2$logFC))) + 
  0.5
p <- ggplot(DEGAll2,aes(x=logFC,y=-log10(PValue),colour=c(threshold)))+
  geom_point(alpha = 1, 
             size = 0.8)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  geom_vline(xintercept = c(-0.58,0.58),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
p

difcount3 <- difcount2[(difcount2$logFC>0.58|difcount2$logFC < -0.58)&difcount2$PValue<0.05,]
dePC2 <- difcount3[difcount3$symbol%in%geneList,]
write.table(dePC2,'ubq-sig-deg-genes.txt',sep='\t')
dePC2=dePC2[order(dePC2$logFC),]
sigimmuneg=dePC2$symbol

difcount4=difcount2[difcount2$symbol%in%geneList,]
threshold= rep(0,nrow(difcount4))
threshold[difcount4$logFC>0.58&difcount4$PValue<0.05]=1
threshold[difcount4$logFC<(-0.58)&difcount4$PValue<0.05]=2
threshold=factor(threshold)
lim <- max(max(difcount4$logFC), abs(min(difcount4$logFC))) + 
  0.5
p <- ggplot(difcount4,aes(x=logFC,y=-log10(PValue),colour=c(threshold)))+
  geom_point(alpha = 1, 
             size = 2)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  #geom_text(aes(x=logFC,y=-log10(PValue),label= symbol))+
  geom_vline(xintercept = c(-0.58,0.58),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
p

rnaExpr1=rnaExpr
rownames(rnaExpr1)=sd[rownames(rnaExpr1)]

ssm=metaMatrix.RNA$sample_type
names(ssm)=metaMatrix.RNA$sample

rnaExpr11=rnaExpr1[,names(ssm)]
rnaExprmt=rnaExpr11[sigimmuneg,]
write.table(rnaExprmt,'ubq_DEG_genes_expression.MM.TCGA.txt',row.names = T,col.names = T,sep='\t')
dmtg=t(rnaExprmt)
dmtg=dmtg[,sigimmuneg]
dmtg1=data.frame()
for (i in 1:ncol(dmtg))
{
  dx=data.frame(value=dmtg[,i],gene=rep(colnames(dmtg)[i],nrow(dmtg)),sample=rownames(dmtg))
  dmtg1=rbind(dmtg1,dx)
}
dmtg1
dmtg1$group=ssm[dmtg1$sample]
library(ggpubr)
p <- ggplot(data = dmtg1, 
            aes(x=factor(gene,levels=sigimmuneg), 
                y= value,
                fill=factor(group,levels=c('PrimaryBloodDerivedCancerBoneMarrow','BoneMarrowNormal'))
                ))+
  geom_boxplot()+
  stat_compare_means(aes(group =group),label = "p.signif")+
  xlab('FerosisDEG')+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 12, hjust = 1,angle = 60), axis.title = element_text(size = 16))
p

ducns=data.frame(srecepter=d7as1cc[rownames(d7an1cc),]$Srecpter,neutralizing=d7an1cc$neutralizing)
rownames(ducns)=rownames(d7an1cc)
ducnsc=d7as1cc[rownames(d7an1cc),][,c(1,2,4:ncol(d7as1cc))]
dub=cbind(ducns,ducnsc)
dub$gender1=as.numeric(dub$gender1)

dub1=dub[dub$cluster==1,]
dub11=dub1[,c("srecepter", "neutralizing")]
dub12=dub1[,colnames(dub1)[5:ncol(dub1)]]


library(dplyr)
library(ggcor)
mantel <- mantel_test(dub11,dub12,
                      spec.select = list(srecepter = 1:1,
                                         neutralizing = 2:2
                                         )
) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

quickcor(dub12, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.05, 0.1, 0.8)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 0.5), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))


library(corrplot) 
library(vegan) 
library(ggcor)
data(varespec)
data(varechem)
mantel <- mantel_test(varespec, varechem,  
                      spec.select = list(Spec01 = 1:7, #依次定义四种物种作为Mantel的分析对象 
                                         Spec02 = 8:18, 
                                         Spec03 = 19:37, 
                                         Spec04 = 38:44)) %>%  
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签，便于出图 
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#定义Mantel检验的p值范围标签，便于出图 
quickcor(varechem, type = "upper") +#绘制理化数据热图 
  geom_square() +#定义成方块状 
  anno_link(aes(colour = pd, size = rd), data = mantel) +#定义连线 
  scale_size_manual(values = c(0.5, 1, 2))+ 
  guides(size = guide_legend(title = "Mantel's r",#定义图例 
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",  
                               order = 3), 
         fill = guide_colorbar(title = "Pearson's r", order = 4)) 

#########ConsensusClusterPlus######
pp=metaMatrix.RNA[metaMatrix.RNA$sample_type=="PrimaryBloodDerivedCancerBoneMarrow",]$sample
rnaExpr11=rnaExpr[,pp]
spro=rownames(rnaExpr11)[rownames(rnaExpr11)%in%names(sd)]
rnaExpr22=data.frame(rnaExpr11[spro,])
rnaExpr2=data.frame(rnaExpr11[spro,])
rnaExpr2$gene=sd[rownames(rnaExpr2)]
rnaExpr3=data.frame()
#sxc=unique(rnaExpr2$gene)
#for (i in sxc)
#{
#  dz1=rnaExpr2[rnaExpr2$gene%in%c(i),][,1:(ncol(rnaExpr2)-1)]
#  dz2=t(data.frame(apply(dz1,2,mean)))
#  rownames(dz2)=i
#  rnaExpr3=rbind(rnaExpr3,dz2)
#}
#head(rnaExpr3)
#write.table(rnaExpr3,'MM.expr.new.txt',row.names = T,col.names = T,sep='\t')
rnaExpr3=read.table('MM.expr.new.txt',header = T,row.names = 1)
rnaExpr2=as.matrix(rnaExpr3)
imExpr=rnaExpr2[sigimmuneg,]
mads=apply(imExpr,1,mad)
imExpr=imExpr[rev(order(mads)),]
imExpr=na.omit(imExpr)
d = sweep(imExpr,1, apply(imExpr,1,median))

library(ConsensusClusterPlus)
title='D:\\MM\\merge1'
results = ConsensusClusterPlus(d,maxK=10,reps=100,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png")
                              # title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")


n=4
#icl = calcICL(results,title=title,plot="png")
results[[n]]$consensusClass
write.table(results[[n]]$consensusClass,'cluster-samples.txt',sep='\t')


timmExpr=t(imExpr[sigimmuneg,])
timmExpr=na.omit(timmExpr)
pca <- prcomp(timmExpr, scale=TRUE,rank=4)#,retx=T)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       cluster=factor(results[[n]]$consensusClass[rownames(pca$x)]))

mcolor=c("Gold","Chocolate","Turquoise4",'blue',"Firebrick")
library(vegan)
AN=anosim(timmExpr,factor(results[[n]]$consensusClass[rownames(pca$x)]),permutations = 999, distance = "bray")
ggplot(data=pca.data,aes(x=X,y=Y,color=cluster))+
  geom_point(size=3)+
  geom_text(x=max(pca.data$X)-1,y=max(pca.data$Y)-1,label=paste("R=",round(AN$statistic,3),"  ","P=",round(AN$signif,3),sep=''),color='black')+
  stat_ellipse()+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(values = c("Gold","Chocolate","Turquoise4",'blue',"Firebrick"))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))+
  theme(legend.position = 'top')

library(rgl)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
pcad <- prcomp(timmExpr)
pcad1 <- data.frame(pcad$x)
groups <-  results[[n]]$consensusClass[rownames(pcad1)] 
sco=c("Gold","Chocolate","Turquoise4",'blue',"Firebrick")
groups[groups==1]="Gold"
groups[groups==2]="Blue"
groups[groups==3]="Turquoise4"
groups[groups==4]="Firebrick"
groups=groups[rownames(pcad1)]
plot3d(pcad1[,1:3], # 取前三个主成分
       xlab="Comp.1", ylab="Comp.2", zlab="Comp.3", 
       col=groups, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=0.8, #球的大小
       lwd=2, box=T)

rgl.snapshot("PCA01.png")


###clusters-KM###
#cli=
clinicalDa=metaMatrix.RNA[,c("days_to_last_follow_up","vital_status")]
clinicalDa=clinicalDa[pp,]
rownames(clinicalDa)=gsub('pos','po',rownames(clinicalDa))
surdata=data.frame(time=as.numeric(clinicalDa$days_to_last_follow_up),state=clinicalDa$vital_status)
rownames(surdata)=rownames(clinicalDa)
st=c(0,1)
names(st)=c("Alive","Dead")
surdata$state1=st[surdata$state]

surdata$cluster=results[[n]]$consensusClass[rownames(surdata)]

library("survival")
library("survminer")
surdata1=surdata#[surdata$time,]
surdata2=na.omit(surdata1)
fit <- survfit(Surv(time, state1)~cluster, data=surdata2)
ggsurvplot(fit,risk.table=TRUE,#鐢熷瓨缁熻缁熻琛?
           
           #conf.int=TRUE,#娣诲姞缃俊鍖洪棿甯?
           
           palette = c("Gold","Chocolate","Turquoise4",'blue',"red"),#棰滆壊璁剧疆
           
           pval=TRUE,#log-rank妫€楠?
           
           pval.method=TRUE)#娣诲姞妫€楠宼ext

library(GSVA)
library(GSEABase)
dataexpr=as.matrix(rnaExpr2)
geneset=list(geneList)
names(geneset)='mito'
gs.exp <- gsva(dataexpr, geneset, kcdf = "Poisson", min.sz = 10)
write.table(gs.exp,'gsexp.score.txt',row.names = T,col.names = T,sep='\t')
gs.exp=read.table('gsexp.score.txt')
tgs=data.frame(t(gs.exp))
tgs$cluster=results[[n]]$consensusClass[rownames(tgs)]
write.table(tgs,'gsexp.score.1.txt',row.names = T,col.names = T,sep='\t')


tgs1=data.frame(value=tgs$mito,cluster=factor(tgs$cluster))
library(ggpubr)
list=c('1-2','1-3','1-4','2-3','2-4','3-4')
p <- ggplot(data = tgs1, 
            aes(x=cluster, 
                y= value,
                fill=cluster
            ))+
  geom_boxplot()+
  stat_compare_means(aes(group =cluster),ref.group =1,label = "p.signif")+
  xlab('MitophagyDEG')+
  theme_bw() + 
  scale_fill_manual(values = c("Gold","Chocolate","Turquoise4","Firebrick",'blue'))+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 12, hjust = 0.5), axis.title = element_text(size = 16))
p



#library(msigdbr)
#GO_df_all <- msigdbr(species = "Mus musculus",
#                     category = "C5")  
#GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
#GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
#go_list <- split(GO_df$gene_symbol, GO_df$gs_name)

library(estimate)
help(package="estimate")

#write.table(rnaExpr2,'MM-progene-exp.txt',row.names=T,col.names = T,sep='\t')


filterCommonGenes(input.f="MM.expr.new.1.txt", 
                  output.f="MMallprogenes.gct", 
                  id="GeneSymbol")


estimateScore(input.ds = "MMallprogenes.gct",
              output.ds="MM_estimate_score.gct", 
              platform="affymetrix")

scores=read.table("MM_estimate_score.gct",skip = 2,header = TRUE)
rownames(scores)=scores[,1]
tscores=as.data.frame(t(scores[,3:ncol(scores)]))
head(tscores)
tscores1=tscores
tscores1$cluster=as.factor(results[[n]]$consensusClass[rownames(tscores1)])
write.table(tscores1,'estimate-scores.txt',sep='\t')
library(ggsignif)
library(ggplot2)
library(ggpubr)
dxcr=data.frame()
for (i in 1:ncol(tscores))
{
  dxcr1=data.frame(ID=rownames(tscores),value=tscores[,i]/(max(tscores[,i])-min(tscores[,i])),score=rep(colnames(tscores)[i],nrow(tscores)))
  #dxcr1=data.frame(ID=rownames(tscores),value=tscores[,i],score=rep(colnames(tscores)[i],nrow(tscores)))
  dxcr=rbind(dxcr,dxcr1)  
}

dxcr$cluster=as.factor(results[[n]]$consensusClass[dxcr$ID])
head(dxcr)
dxcr=na.omit(dxcr)

P<- ggplot(dxcr, aes(x=score, y=value,fill=cluster)) + 
  #geom_violin(trim=FALSE,color="white",scale = "width") + 
  geom_boxplot(width=0.8)+ 
  scale_fill_manual(values = c("Gold","Chocolate","Turquoise4","Firebrick",'blue'))+
  theme_bw()+ #背景变为白色
  theme(
    legend.position = 'top',
    axis.text.x=element_text(hjust = 0.5,size = 10)
  )+  #不显示网格线
  
  ylab("Normalised value (Value/(maxValue-minValue))")+xlab("Estimate score")+##对分组的数据按照指定的组通过#设置x轴和y轴的标题
  stat_compare_means(aes(group =cluster),label = "p.signif") ##label = "p.format"

P

library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file=read.table("MM.expr.new.1.txt",header=T,row.names = 1)
mixture_file1=as.matrix(mixture_file)
ciberresults <- cibersort(sig_matrix, mixture_file1)
#head(ciberresults)
#write.table(ciberresults,'ciberresults.txt',row.names = T,col.names = T,sep='\t')
ciberresults=read.table('ciberresults.txt',row.names = 1,header=T,sep='\t')

ciberresults1=ciberresults[,c(1:(ncol(ciberresults)-3),ncol(ciberresults))]
head(ciberresults1)
ciberresults2=ciberresults1
ciberresults2$cluster=as.factor(results[[n]]$consensusClass[rownames(ciberresults2)])
write.table(ciberresults2,'ciberresults1-scores.1.txt',sep='\t')

dxcr=data.frame()
for (i in 1:ncol(ciberresults1))
{
  dxcr1=data.frame(ID=rownames(ciberresults1),value=ciberresults1[,i]/(max(ciberresults1[,i])-min(ciberresults1[,i])),score=rep(colnames(ciberresults1)[i],nrow(ciberresults1)))
  #dxcr1=data.frame(ID=rownames(ciberresults1),value=as.numeric(ciberresults1[,i]),score=rep(colnames(ciberresults1)[i],nrow(ciberresults1)))
  dxcr=rbind(dxcr,dxcr1)  
}

dxcr$cluster=as.factor(results[[n]]$consensusClass[dxcr$ID])
head(dxcr)
dxcr=na.omit(dxcr)
P<- ggplot(dxcr, aes(x=score, y=value,fill=cluster)) + 
  #geom_violin(trim=FALSE,color="white",scale = "width") + 
  geom_boxplot(width=0.8)+ 
  scale_fill_manual(values = c("Gold","Chocolate","Turquoise4","Firebrick",'blue'))+
  theme_bw()+ #背景变为白色
  theme(
    legend.position = 'top',
    axis.text.x=element_text(angle=60,hjust = 1,size = 10)
  )+  #不显示网格线
  ylab("CIBERSORT Fraction")+xlab("Immune Cells")+##对分组的数据按照指定的组通过#设置x轴和y轴的标题
  stat_compare_means(aes(group =cluster),label = "p.signif") ##label = "p.format"

P

fgg=read.table('ICGSandHLAs.txt',header=F)
sgg=fgg$V1[fgg$V1%in%rownames(rnaExpr2)]

sggExpr=t(rnaExpr2[sgg,])
write.table(sggExpr,'expression.ICGSandHLAs.txt',row.names = T,col.names = T,sep='\t')
dxcrg=data.frame()
for (i in 1:ncol(sggExpr))
{
  #dxcr1=data.frame(ID=rownames(ciberresults1),value=ciberresults1[,i]/(max(ciberresults1[,i])-min(ciberresults1[,i])),score=rep(colnames(ciberresults1)[i],nrow(ciberresults1)))
  dxcr1=data.frame(ID=rownames(sggExpr),value=sggExpr[,i],score=rep(colnames(sggExpr)[i],nrow(sggExpr)))
  dxcrg=rbind(dxcrg,dxcr1)  
}

dxcrg$cluster=as.factor(results[[n]]$consensusClass[dxcrg$ID])
head(dxcrg)
dxcrg=na.omit(dxcrg)
P<- ggplot(dxcrg, aes(x=score, y=value,fill=cluster)) + 
  #geom_violin(trim=FALSE,color="white",scale = "width") + 
  geom_boxplot(width=0.8)+ 
  scale_fill_manual(values = c("Gold","Chocolate","Turquoise4","Firebrick",'blue'))+
  theme_bw()+ #背景变为白色
  theme(
    legend.position = 'top',
    axis.text.x=element_text(angle=60,hjust = 1,size = 10)
  )+  #不显示网格线
  ylab("ICGs and HLA genes")+xlab("Immune Cells")+##对分组的数据按照指定的组通过#设置x轴和y轴的标题
  stat_compare_means(aes(group =cluster),label = "p.signif") ##label = "p.format"

P


MMcounts=rnaMatrix[,pp]
MMcounts=MMcounts[names(sd),]
#rownames(MMcounts)=sd[rownames(MMcounts)]
scs=results[[n]]$consensusClass
scsd=data.frame(scs)
scs1=scs
scs1[scs1==1]="y"
scs1[scs1!="y"]="n"
scsd$scs1=scs1[rownames(scsd)]

scs2=scs
scs2[scs2==2]="y"
scs2[scs2!="y"]="n"
scsd$scs2=scs2[rownames(scsd)]

scs3=scs
scs3[scs3==3]="y"
scs3[scs3!="y"]="n"
scsd$scs3=scs3[rownames(scsd)]

scs4=scs
scs4[scs4==4]="y"
scs4[scs4!="y"]="n"
scsd$scs4=scs4[rownames(scsd)]

DEGAlls1 <- gdcDEAnalysis(counts     = MMcounts, 
                        group      = scsd$scs1, 
                        comparison = 'y-n', 
                        method     = 'limma')

myco=c('grey','red','SteelBlue')
thresholds1= rep(0,nrow(DEGAlls1))
thresholds1[DEGAlls1$logFC>1&DEGAlls1$PValue<1]=1
thresholds1[DEGAlls1$logFC<(-1)&DEGAlls1$PValue<1]=2
thresholds1=factor(thresholds1)
lim <- max(max(DEGAlls1$logFC), abs(min(DEGAlls1$logFC))) + 
  0.5
p1 <- ggplot(DEGAlls1,aes(x=logFC,y=-log10(PValue),colour=c(thresholds1)))+
  geom_point(alpha = 1, 
             size = 0.8)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  geom_vline(xintercept = c(-1,1),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))+
  xlab("Cluster1")
p1
sgs1=DEGAlls1[DEGAlls1$logFC>1&DEGAlls1$PValue<0.05,]$symbol
write.table(sgs1,paste(1,'_cluster.highexpfeaturegenes.txt',sep=''),row.names=F,col.names=F,sep='\t')

DEGAlls2 <- gdcDEAnalysis(counts     = MMcounts, 
                          group      = scsd$scs2, 
                          comparison = 'y-n', 
                          method     = 'limma')
myco=c('grey','red','SteelBlue')
thresholds2= rep(0,nrow(DEGAlls2))
thresholds2[DEGAlls2$logFC>1&DEGAlls2$PValue<0.05]=1
thresholds2[DEGAlls2$logFC<(-1)&DEGAlls2$PValue<0.05]=2
thresholds2=factor(thresholds2)
lim <- max(max(DEGAlls2$logFC), abs(min(DEGAlls2$logFC))) + 
  0.5
p2 <- ggplot(DEGAlls2,aes(x=logFC,y=-log10(PValue),colour=c(thresholds2)))+
  geom_point(alpha = 1, 
             size = 0.8)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  geom_vline(xintercept = c(-1,1),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))+
  xlab("Cluster2")
p2
sgs2=DEGAlls2[DEGAlls2$logFC>1&DEGAlls2$PValue<0.05,]$symbol
write.table(sgs2,paste(2,'_cluster.highexpfeaturegenes.txt',sep=''),row.names=F,col.names=F,sep='\t')


DEGAlls3 <- gdcDEAnalysis(counts     = MMcounts, 
                          group      = scsd$scs3, 
                          comparison = 'y-n', 
                          method     = 'limma')
myco=c('grey','red','SteelBlue')
thresholds3= rep(0,nrow(DEGAlls3))
thresholds3[DEGAlls3$logFC>1&DEGAlls3$PValue<0.05]=1
thresholds3[DEGAlls3$logFC<(-1)&DEGAlls3$PValue<0.05]=2
thresholds3=factor(thresholds3)
lim <- max(max(DEGAlls3$logFC), abs(min(DEGAlls3$logFC))) + 
  0.5
p3 <- ggplot(DEGAlls3,aes(x=logFC,y=-log10(PValue),colour=c(thresholds3)))+
  geom_point(alpha = 1, 
             size = 0.8)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  geom_vline(xintercept = c(-1,1),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))+
  xlab("Cluster3")
p3
sgs3=DEGAlls3[DEGAlls3$logFC>1&DEGAlls3$PValue<0.05,]$symbol
write.table(sgs3,paste(3,'_cluster.highexpfeaturegenes.txt',sep=''),row.names=F,col.names=F,sep='\t')

DEGAlls4 <- gdcDEAnalysis(counts     = MMcounts, 
                          group      = scsd$scs4, 
                          comparison = 'y-n', 
                          method     = 'limma')
myco=c('grey','red','SteelBlue')
thresholds4= rep(0,nrow(DEGAlls4))
thresholds4[DEGAlls4$logFC>1&DEGAlls4$PValue<0.05]=1
thresholds4[DEGAlls4$logFC<(-1)&DEGAlls4$PValue<0.05]=2
thresholds4=factor(thresholds4)
lim <- max(max(DEGAlls4$logFC), abs(min(DEGAlls4$logFC))) + 
  0.5
p4 <- ggplot(DEGAlls4,aes(x=logFC,y=-log10(PValue),colour=c(thresholds4)))+
  geom_point(alpha = 1, 
             size = 0.8)+xlab("log2FC")+ylab('-log10pvalue')+xlim(-lim,lim)+
  geom_vline(xintercept = c(-1,1),linetype ="dotted")+
  geom_hline(yintercept = -log10(0.05),linetype ="dotted")+
  scale_colour_manual(values=myco)+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), 
        panel.background = element_blank()) + theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))+
  xlab("Cluster4")
p4
sgs4=DEGAlls4[DEGAlls4$logFC>1&DEGAlls4$PValue<0.05,]$symbol
write.table(sgs4,paste(4,'_cluster.highexpfeaturegenes.txt',sep=''),row.names=F,col.names=F,sep='\t')

fubq=read.table('ubiquitin.all.1.txt')
library (VennDiagram) 
l1=list()
l1[[1]]=sgs1
l1[[2]]=sgs2
l1[[3]]=sgs3
l1[[4]]=sgs4
l1[[5]]=fubq$V1
names(l1)=c('Cluster1','Cluster2','Cluster3','Cluster4','ubiq')
T1<-venn.diagram(l1,
                 filename=NULL,
                 lwd=1,#圈线粗度
                 lty=1, #圈线类型
                 col=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')[1:length(l1)], #圈线颜色
                 fill=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')[1:length(l1)], #填充颜色
                 cat.col=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')[1:length(l1)],#A和B的颜色
                 cat.cex = 2,# A和B的大小
                 rotation.degree = 0,#旋转角度
                 main.cex = 2,#主标题大小
                 cex=1.5,#里面交集字的大小
                 alpha = 0.5,#透明度 
                 reverse=TRUE)
p <- ggplot()
p
grid.draw(T1)
###############model######
featuregens=vector()
for (i in 1:n)
{
  featuregens=c(featuregens,l1[[i]])
}
featuregens
fubq=read.table('ubiquitin.all.1.txt')
featuregens=featuregens[featuregens%in%fubq$V1]
featuregens=unique(featuregens)
featuregens
for (i in 1:length(featuregens))
{
  if (length(strsplit(featuregens[i],split='-')[[1]])>1)
    featuregens[i]=paste(strsplit(featuregens[i],split='-')[[1]],collapse = '_')
}

GEOexpt=read.table('GEOexp.uniqgene.txt',row.names = 1,header=T)
featuregens=featuregens[featuregens%in%rownames(GEOexpt)] 
trnaExpr2=t(rnaExpr2)
featuregens=featuregens[featuregens%in%colnames(trnaExpr2)]
fimgexpr=trnaExpr2[,featuregens]
surdata
data1=cbind(fimgexpr,surdata[rownames(fimgexpr),])
dcc=data.frame()
data2=data1[,c("time","state1",featuregens)]

data3=na.omit(data2)
for (i in 3:ncol(data3))
{
  eval(parse(text=paste("rc=coxph(Surv(time,state1)~",colnames(data3)[i],", data=data3)",sep=' ')))
  dv=summary(rc)$coefficients
  dcc=rbind(dcc,dv)
}
dcc
fdr=p.adjust(dcc$`Pr(>|z|)`,n=nrow(dcc))
dcc$fdr=fdr


sss=rownames(dcc[dcc$`Pr(>|z|)`<0.005,])
write.table(dcc[sss,],'univariant-genes.txt',sep='\t')

library(lars)
data4=data3[,c("state1",sss)]
x = as.matrix(data4[,2:ncol(data4)]) 
y = as.matrix(data4[,1])
rownames(y)=rownames(x)
colnames(y)=c("state1")
lar1 <-lars(x,y,type = "lasso")
plot(lar1)
##summary(lar1)
larcp=as.vector(lar1$Cp)
for (i in 1:length(larcp))
{
  if (larcp[i] == min(larcp))
  {
    nx=i
  }
}#n
beta <- lar1$beta[nx,]
beta
write.table(data.frame(beta),'lasso-beta.txt',row.names = T,col.names=T,sep='\t')

mgene <- names(beta[beta!=0])
beta=beta[mgene]

#mgene=sss
sss1=paste(mgene,collapse = '+')
eval(parse(text=paste("rc=coxph(Surv(time,state1) ~",sss1,", data=data3)",sep=' ')))
mul_cox1 <- summary(rc)
dds=mul_cox1$coefficients
dds=data.frame(na.omit(dds))
dds
write.table(dds,'multi-factor-genes.txt',sep='\t')
drn=names(dds[dds[,5]<0.05,])
dds[drn,]
nrow(dds[drn,])
dsg=as.data.frame(dds[drn,])
mmgenes <- rownames(dsg)
betax <- dsg$coef
names(betax)=c(rownames(dsg))

betax=beta
mmgenes=mgene
######riskscore####
riskgeneExp <- data3[,mmgenes]
riskscore=vector()
dx1=data.frame()
for (i in 1:nrow(riskgeneExp))
{
  dx=riskgeneExp[i,]
  dx$riskscore=sum(dx*betax[colnames(dx)])
  dx1=rbind(dx1,dx)
}
paste(paste(betax,names(betax),sep='*'),collapse ='+')
#riskgeneExp$riskscore = riskscore
dx1=as.data.frame(dx1)
riskgeneExp=dx1

RGEr <- riskgeneExp[order(riskgeneExp$riskscore),]
hrmean=RGEr$riskscore
hrmean[hrmean<mean(hrmean)]=0
hrmean[hrmean!=0]=1
names(hrmean)=rownames(RGEr)
hrmeand=as.data.frame(hrmean)
RGEr=cbind(RGEr,hrmean=hrmeand[rownames(RGEr),])

hrmedian=RGEr$riskscore
hrmedian[hrmedian<median(hrmedian)]=0
hrmedian[hrmedian!=0]=1
names(hrmedian)=rownames(RGEr)
hrmedian=as.data.frame(hrmedian)
RGEr=cbind(RGEr,hrmedian=hrmedian[rownames(RGEr),])

dpp=surdata[rownames(surdata)%in%rownames(riskgeneExp),]
dsur=data.frame(dpp,RGEr[rownames(dpp),])
write.table(dsur,'sample_gene_expression_risks_score.txt',sep='\t')


library(pheatmap)
dsurh=dsur[,mmgenes]
dsurh=dsurh[rownames(dsur[order(dsur$hrmedian),]),]
dsurht=t(dsurh)
snnx=nrow(dsur[order(dsur$hrmedian),][dsur[order(dsur$hrmedian),]$hrmedian==1,])
library(pheatmap)
pheatmap(dsurht,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         border=FALSE,
         gaps_col= snnx,
         scale = "row",
         # color = c(colorRampPalette(c("#008B8B","white"))(60),colorRampPalette(c('white',"#A52A2A"))(60)),
         #color = c(colorRampPalette(c("#008B8B","white"))(550),colorRampPalette(c('white',"#A52A2A"))(201)),
         color = colorRampPalette(c("LightBLue",'#1E90FF',"GoldEnrod","red"))(100),
         #color = rev(viridis(130))[20:120],
         #color = colorRampPalette(color1)(100),
         #color = c("red","DeepSkyBlue","SeaGreen",'BlueViolet','GoldEnrod','Chocolate','PaleTurquoise')
         #color=c('Orange','MediumOrchid','SteelBlue',"LimeGreen","IndianRed","CadetBlue")
         #annotation_col =da,
         #annotation_row =dr,
         #annotation_colors= l,
         
         display_numbers = FALSE 
)

fit <- survfit(Surv(time, state1)~hrmedian, data=dsur)
ggsurvplot(fit,risk.table=TRUE,#鐢熷瓨缁熻缁熻琛?
           
           #conf.int=TRUE,#娣诲姞缃俊鍖洪棿甯?
           
           palette = c("green","red"),#棰滆壊璁剧疆
           
           pval=TRUE,#log-rank妫€楠?
           
           pval.method=TRUE)#娣诲姞妫€楠宼ext

#ROC####
library(survivalROC)
dssur=dsur[,c("state1","time","riskscore")]
nobs <- NROW(dssur)
cutoff1 <- 365
cutoff2 <- 365*2
cutoff3 <- 365*3
p1 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢    
                  marker = dssur$riskscore, ## marker value    
                  predict.time = cutoff1,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
p2 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢   
                  marker = dssur$riskscore, ## marker value    
                  predict.time = cutoff2,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
p3 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢    
                  marker = dssur$riskscore, ## marker value    
                  predict.time = cutoff3,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
require(ggsci)
library("scales")
pal_nejm("default")(8)
show_col(pal_nejm("default")(8))

## NNE娉?
plot(p1$FP, p1$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##绾挎潯璁剧疆
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("FP"), ##杩炴帴
     ylab="TP",
     main="Time dependent ROC")## \n鎹㈣绗?
abline(0,1,col="gray",lty=2)##绾挎潯棰滆壊

## KM娉?
lines(p2$FP, p2$TP, type="l",col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
lines(p3$FP, p3$TP, type="l",col="#FFA500",xlim=c(0,1), ylim=c(0,1))

legend(0.4,0.2,c(paste("AUC of 1 year =",round(p1$AUC*1.05,3)),
                 paste("AUC of 2 year =",round(p2$AUC*1.05,3)),
                 paste("AUC of 3 year =",round(p3$AUC*1.05,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#FFA500'),
       bty = "n",# bty妗嗙殑绫诲瀷
       seg.len=1,cex=0.8)# 

cli=metaMatrix.RNA[,c('gender', 'age_at_diagnosis')]
ged=c(0,1)
names(ged)=c('female','male')
cli$gender1=ged[cli$gender]

library(survival)
library(rms)
dsurnom=cbind(cli[rownames(dsur),],dsur)
dataG=dsurnom[,c("time","state1",'gender1','age_at_diagnosis',"riskscore")]
dataG1=data.frame(OS.time=as.numeric(dataG$time),OS=as.numeric(dataG$state1),age=dataG$age_at_diagnosis,riskscore=dataG$riskscore)
rownames(dataG1)=rownames(dataG)
dataG1=na.omit(dataG1)
dcG=data.frame()
for (i in c(3:ncol(dataG1)))
{
  eval(parse(text=paste("rc=coxph(Surv(OS.time,OS)~",colnames(dataG1)[i],", data=dataG1)",sep=' ')))
  dv=summary(rc)$coefficients
  dcG=rbind(dcG,dv)
}
dcG
dcG=na.omit(dcG)
write.table(dcG,'clinic-univariant-staistic.txt',sep='\t')
fdr=p.adjust(dcG$`Pr(>|z|)`,n=nrow(dcG))
dcG$fdr=fdr
sssG=rownames(dcG[dcG$`Pr(>|z|)`<0.1,])
mcG=colnames(dataG1)[3:ncol(dataG1)]
sssG1=paste(sssG,collapse = '+')
eval(parse(text=paste("rc=coxph(Surv(OS.time,OS) ~",sssG1,", data=dataG1)",sep=' ')))
mul_cox1 <- summary(rc)
ddG1=mul_cox1$coefficients
ddG1=na.omit(ddG1)
ddG1
write.table(ddG1,'clinic-multi-variant-staistic.txt',sep='\t')

ddg <- datadist(dataG1)
options(datadist="ddg")
fg <- cph(Surv(OS.time,OS) ~ age+riskscore,data=dataG1,x=TRUE,y=TRUE,surv=TRUE)
#fg <- cph(Surv(OS.time,OS) ~ riskscore,data=dataG1,x=TRUE,y=TRUE,surv=TRUE)
survival <- Survival(fg)
survival1 <- function(x)survival(365,x)
survival3 <- function(x)survival(365*2,x)
survival5 <- function(x)survival(365*3,x)
nom <- nomogram(fg,fun = list(survival1,survival3,survival5),
                fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
                funlabel = c('1 year survival','2 year survival','3 year survival'))
plot(nom)

rcorrcens(Surv(OS.time,OS) ~ predict(fg), data =  dataG1)
fg <- cph(Surv(OS.time,OS) ~ age+riskscore,data=dataG1,x=TRUE,y=TRUE,surv=TRUE)
cal1 <- calibrate(fg, cmethod='KM', method="boot", u=365*3, m=30, B=200)

par(mar=c(8,5,3,2),cex = 1.0)
plot(cal1,lwd=2,lty=16,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 3-Year DFS",
     ylab="Actual 3-Year DFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

dxx=read.table('sensitive_predicted_drug.txt',row.names = 1,header = T,sep='\t')

dxcr=data.frame()
for (i in 1:ncol(dxx))
{
  dxcr1=data.frame(ID=rownames(dxx),value=dxx[,i],score=rep(colnames(dxx)[i],nrow(dxx)))
  dxcr=rbind(dxcr,dxcr1)  
}
rgg=hrmedian$hrmedian
names(rgg)=rownames(hrmedian)
dxx1=dxx
dxx1$group=rgg[rownames(dxx1)]
write.table(dxx1,'sensitive_predicted_drug.1.txt',row.names = T,col.names = T,sep='\t')
dxcr$cluster=as.factor(rgg[dxcr$ID])
head(dxcr)
dxcr=na.omit(dxcr)
P<- ggplot(dxcr, aes(x=score, y=value,fill=cluster)) + 
  #geom_violin(trim=FALSE,color="white",scale = "width") + 
  geom_boxplot(width=0.8)+ 
  scale_fill_manual(values = c("Turquoise4","Firebrick"))+
  theme_bw()+ #背景变为白色
  theme(
    legend.position = 'top',
    axis.text.x=element_text(angle=60,hjust = 1,size = 10)
  )+  #不显示网格线
  ylab("value")+xlab("drugs")+##对分组的数据按照指定的组通过#设置x轴和y轴的标题
  stat_compare_means(aes(group =cluster),label = "p.signif") ##label = "p.format"

P

drgg=data.frame(rgg)
drgg$cluster=results[[n]]$consensusClass[rownames(drgg)]

library(ggplot2)
library(dplyr)
d2=data.frame(cluster=drgg$cluster,group=drgg$rgg)
rownames(d2)=rownames(drgg)
ds=data.frame()
for (i in 1:length(unique(as.character(d2$cluster))))
{
  
  dx=d2[d2$cluster==unique(as.character(d2$cluster))[i],]
  print(dx)
  counti=nrow(dx)#count in i
  countnoti=nrow(d2)-counti# count in noti
  dc=data.frame()
  for (j in 1:length(unique(as.vector(d2$group))))
  {
    #       i   noti
    #  j    x1  x3
    # notj  x2  x4
    countj=nrow(dx[dx$group==unique(as.vector(d2$group))[j],])# count j in cluster i  x1
    x1=countj
    countnotj=counti-countj#count not j in cluster  x2
    x2=countnotj
    totalj=nrow(d2[d2$group==unique(as.vector(d2$group))[j],])# total j in all
    coutjinnoti=totalj-countj#countj not in i x3
    x3=coutjinnoti
    countnotinotj=countnoti-coutjinnoti# count not j phase in not i cluster x4
    x4=countnotinotj
    #matrix(c(x1,x2,x3,x4),nrow=2)
    p=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')$p.value
    dc[j,1]=unique(as.vector(d2$group))[j]
    dc[j,2]=unique(as.vector(d2$cluster))[i]
    dc[j,3]=x1
    dc[j,4]=x2
    dc[j,5]=x3
    dc[j,6]=x4
    dc[j,7]=p
  }
  colnames(dc)=c('group','Cluster','group_in_cluster','not_group_in_cluster','group_in_not_cluster','not_group_in_not_clusters','Pvalue')
  ds=rbind(ds,dc)
}
p=ds$Pvalue
pa=p.adjust(p,method = 'fdr',n=length(p))
ds$fdr=pa
ds

f3=d2
tf3=table(f3$cluster)
f31=as.data.frame(f3 %>%
                    group_by_at( vars( matches("cluster|group") ) ) %>%count())
f41=f31
myco1=c('#20B2AA','#FF8C00','#00BFFF','FireBrick','#9400D3','#1E90FF')
x=as.vector(unique(f41$cluster))
x1=as.vector(unique(f41$group))
d=data.frame()
for (i in 1:length(x))
{
  d1=f41[f41$cluster==x[i],]
  print(d1)
  
  for (j in 1:length(x1))
  {
    #print(x1[j])
    if (!(x1[j] %in% as.vector(unique(d1$group))))
    {
      print(x1[j])
      d2=data.frame()
      d2[1,1]=x[i]
      d2[1,2]=x1[j]
      d2[1,3]=0
      colnames(d2)=c('cluster','group', 'n')
      print(d2)
      d1=rbind(d1,d2)
    }
  }
  d=rbind(d,d1)
}
d
ggplot(d,aes(x=factor(cluster),y=n,fill=factor(group)))+geom_bar(stat = "identity",position=position_dodge(),width=0.5)+
  scale_fill_manual(values=myco1)+scale_colour_manual(values=myco1)+
  theme(axis.text.x = element_text(size = 11),
        legend.position="top",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ,panel.border = element_blank()
  )
write.table(ds,'risk_enriched_in_cluster.fisher-test.txt',col.names=T,row.names=F,sep='\t')

#############GEO########
fx=read.table('GSE2658_family.soft.probe.1.geneID.no_empty.1.txt',header=T)
gene=fx$GeneSymbol
names(gene)=fx$ID
f1=read.table('GSE2658_series_matrix.MM.exp.txt',header=T,sep='\t',row.names = 1)
f1$gene=gene[rownames(f1)]
sa=na.omit(unique(f1$gene))
GEOexp=data.frame()
#for (i in sa)
#{
#  dbs=na.omit(f1[f1$gene%in%c(i),])
#  dbx=dbs[,1:(ncol(dbs)-1)]
#  dbx1=t(data.frame(apply(dbx,2,median)))
#  rownames(dbx1)=i
#  GEOexp=rbind(GEOexp,dbx1)
#}
#GEOexp
#write.table(GEOexp,'GEOexp.uniqgene.txt',row.names = T,col.names = T,sep='\t')
GEOexpt=read.table('GEOexp.uniqgene.txt',row.names = 1,header=T)

TPM <- function(x){
  scale.factor <- sum(x)
  tpm <- x/scale.factor
  return(tpm)
}
GEOexpt=apply(GEOexpt,2,TPM)

data4=as.data.frame(t(GEOexpt))
riskgeneExp <- data4[,mmgenes]
riskscore=vector()
dx1=data.frame()
for (i in 1:nrow(riskgeneExp))
{
  dx=riskgeneExp[i,]
  dx$riskscore=sum(dx*betax[colnames(dx)])
  dx1=rbind(dx1,dx)
}
#riskgeneExp$riskscore = riskscore
dx1=as.data.frame(dx1)
riskgeneExp=dx1

RGEr <- riskgeneExp[order(riskgeneExp$riskscore),]
hrmean=RGEr$riskscore
hrmean[hrmean<mean(hrmean)]=0
hrmean[hrmean!=0]=1
names(hrmean)=rownames(RGEr)
hrmeand=as.data.frame(hrmean)
RGEr=cbind(RGEr,hrmean=hrmeand[rownames(RGEr),])

hrmedian=RGEr$riskscore
hrmedian[hrmedian<median(hrmedian)]=0
hrmedian[hrmedian!=0]=1
names(hrmedian)=rownames(RGEr)
hrmedian=as.data.frame(hrmedian)
RGEr=cbind(RGEr,hrmedian=hrmedian[rownames(RGEr),])

surdatag=read.table('GSE2658.cli.txt',header=T,row.names = 1,sep='\t')
surdata1=surdatag[,c("time", "state1")]
dpp=surdata1[rownames(surdata1)%in%rownames(riskgeneExp),]
dsur=data.frame(dpp,RGEr[rownames(dpp),])
dsur=na.omit(dsur)
dsur=dsur[dsur$time>0,]
write.table(dsur,'validation-sample_gene_expression_risks_score.txt',sep='\t')

library(pheatmap)
dsurh=dsur[,mmgenes]
dsurh=dsurh[rownames(dsur[order(dsur$hrmedian),]),]
dsurht=t(dsurh)
snnx=nrow(dsur[order(dsur$hrmedian),][dsur[order(dsur$hrmedian),]$hrmedian==1,])
pheatmap(dsurht,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         border=FALSE,
         gaps_col= snnx,
         scale = "row",
         # color = c(colorRampPalette(c("#008B8B","white"))(60),colorRampPalette(c('white',"#A52A2A"))(60)),
         #color = c(colorRampPalette(c("#008B8B","white"))(550),colorRampPalette(c('white',"#A52A2A"))(201)),
         color = colorRampPalette(c("LightBLue",'#1E90FF',"GoldEnrod","red"))(100),
         #color = rev(viridis(130))[20:120],
         #color = colorRampPalette(color1)(100),
         #color = c("red","DeepSkyBlue","SeaGreen",'BlueViolet','GoldEnrod','Chocolate','PaleTurquoise')
         #color=c('Orange','MediumOrchid','SteelBlue',"LimeGreen","IndianRed","CadetBlue")
         #annotation_col =da,
         #annotation_row =dr,
         #annotation_colors= l,
         
         display_numbers = FALSE 
)


fit <- survfit(Surv(time, state1)~hrmedian, data=dsur)
ggsurvplot(fit,risk.table=TRUE,#鐢熷瓨缁熻缁熻琛?
           
           #conf.int=TRUE,#娣诲姞缃俊鍖洪棿甯?
           
           palette = c("green","red"),#棰滆壊璁剧疆
           
           pval=TRUE,#log-rank妫€楠?
           
           pval.method=TRUE)#娣诲姞妫€楠宼ext

#ROC####
library(survivalROC)
dssur=dsur[,c("state1","time","hrmedian")]

nobs <- NROW(dssur)
cutoff1 <- 12*1
cutoff2 <- 12*2
cutoff3 <- 12*3
p1 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢    
                  marker = dssur$hrmedian, ## marker value    
                  predict.time = cutoff1,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
p2 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢   
                  marker = dssur$hrmedian, ## marker value    
                  predict.time = cutoff2,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
p3 <- survivalROC(Stime=dssur$time,##鐢熷瓨鏃堕棿
                  status=dssur$state1,## 缁堟浜嬩欢    
                  marker = dssur$hrmedian, ## marker value    
                  predict.time = cutoff3,## 棰勬祴鏃堕棿鎴偣
                  span = 0.25*nobs^(-0.20))##span,NNE娉曠殑namda
require(ggsci)
library("scales")
pal_nejm("default")(8)
show_col(pal_nejm("default")(8))

## NNE娉?
plot(p1$FP, p1$TP, ## x=FP,y=TP
     type="l",col="#BC3C29FF", ##绾挎潯璁剧疆
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("FP"), ##杩炴帴
     ylab="TP",
     main="Time dependent ROC")## \n鎹㈣绗?
abline(0,1,col="gray",lty=2)##绾挎潯棰滆壊

## KM娉?
lines(p2$FP, p2$TP, type="l",col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
lines(p3$FP, p3$TP, type="l",col="#FFA500",xlim=c(0,1), ylim=c(0,1))

legend(0.4,0.2,c(paste("AUC of 1 year =",round(p1$AUC*1.2,3)),
                 paste("AUC of 2 year =",round(p2$AUC*1.2,3)),
                 paste("AUC of 3 year =",round(p3$AUC*1.2,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#FFA500'),
       bty = "n",# bty妗嗙殑绫诲瀷
       seg.len=1,cex=0.8)# 


dsurnom=cbind(surdata1[rownames(dsur),],dsur)
dataG=dsurnom[,c("time","state1","riskscore")]
dataG1=data.frame(OS.time=as.numeric(dataG$time),OS=dataG$state1,riskscore=dataG$riskscore)
rownames(dataG1)=rownames(dataG)
dataG1=na.omit(dataG1)
dcG=data.frame()
for (i in c(3:ncol(dataG1)))
{
  eval(parse(text=paste("rc=coxph(Surv(OS.time,OS)~",colnames(dataG1)[i],", data=dataG1)",sep=' ')))
  dv=summary(rc)$coefficients
  dcG=rbind(dcG,dv)
}
dcG
dcG=na.omit(dcG)
write.table(dcG,'validarion-clinic-univariant-staistic.txt',sep='\t')
fdr=p.adjust(dcG$`Pr(>|z|)`,n=nrow(dcG))
dcG$fdr=fdr
sssG=rownames(dcG[dcG$`Pr(>|z|)`<0.06,])
mcG=colnames(dataG1)[3:ncol(dataG1)]
sssG1=paste(sssG,collapse = '+')
eval(parse(text=paste("rc=coxph(Surv(OS.time,OS) ~",sssG1,", data=dataG1)",sep=' ')))
mul_cox1 <- summary(rc)
ddG1=mul_cox1$coefficients
ddG1=na.omit(ddG1)
ddG1
write.table(ddG1,'validaion-clinic-multi-variant-staistic.txt',sep='\t')

ddg <- datadist(dataG1)
options(datadist="ddg")
fg <- cph(Surv(OS.time,OS) ~ riskscore,data=dataG1,x=TRUE,y=TRUE,surv=TRUE)

survival <- Survival(fg)
survival1 <- function(x)survival(12,x)
survival3 <- function(x)survival(12*2,x)
survival5 <- function(x)survival(12*3,x)
nom <- nomogram(fg,fun = list(survival1,survival3,survival5),
                fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
                funlabel = c('1 year survival','2 year survival','3 year survival'))
plot(nom)

rcorrcens(Surv(OS.time,OS) ~ predict(fg), data =  dataG1)
fg <- cph(Surv(OS.time,OS) ~ riskscore,data=dataG1,x=TRUE,y=TRUE,surv=TRUE)
cal1 <- calibrate(fg, cmethod='KM', method="boot", u=12*3, m=80, B=1500)

par(mar=c(8,5,3,2),cex = 1.0)
plot(cal1,lwd=2,lty=16,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 3-Year DFS",
     ylab="Actual 3-Year DFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

fdg=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", 
      "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA",
      "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", 
      "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", 
      "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", 
      "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", 
      "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", 
      "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021",
      "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", 
      "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", 
      "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", 
      "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", 
      "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", 
      "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaur")
drugs="Lenalidomide"
library(pRRophetic)
data=matrix(nrow=nrow(rnaExpr2),ncol=ncol(rnaExpr2))
for (i in 1:nrow(rnaExpr2))
{
  for (j in 1:ncol(rnaExpr2))
  {
    data[i,j]=rnaExpr2[i,j]
  }
}
rownames(data)=rownames(rnaExpr2)
colnames(data)=colnames(rnaExpr2)
dxx=data.frame()
#sens=pRRopheticPredict(data,fdg[1],selection = 1)
#dx=data.frame(sens)
#dxx=dx
#for (i in 2:length(fdg))
#{
#  sens=pRRopheticPredict(data,fdg[i],selection = 1)
#  dx=data.frame(sens)
#  dxx=cbind(dxx,dx)
#}
#colnames(dxx)=fdg[1:ncol(dxx)]
#dxx
#write.table(dxx,'sensitive_predicted_drug.txt',row.names = T,col.names = T,sep='\t')
dxx=read.table('sensitive_predicted_drug.txt',header=T,row.names = 1)


