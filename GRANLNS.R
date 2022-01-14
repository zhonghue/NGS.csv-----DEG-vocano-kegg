#从csv格式的二代测序counts中读取到Rstudio中
exp=GSE129985_deSeq2_counts
exp <- as.data.frame(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:5]

#group是GSM（GEO sample）中，测序样本的不同处理。
#比如GRASLND shRNA（用GRASLND shRNA处理）。
#是factor格式。
group=c('scrambled shRNA',
        'scrambled shRNA',
        'GRASLND shRNA',
        'GRASLND shRNA')
group <- factor(group)
#index1=c(rowSums(exp>0))
#class(index1)
#exp2=exp[index1,]

#sample是每个样本的名字，比如21天用空shRNA处理的细胞。
#建立一个表格f，使得：
#（1）行名为‘样本名字’；
#（2）列名为‘condition’；
#（3）内容为‘group’中的内容。
sample=c(colnames(exp2))
f=as.data.frame(group)
rownames(f)=c(sample)
colnames(f)='condition'

#使用DESeq2处理表格。
library('DESeq2')
dds=DESeqDataSetFromMatrix(countData = exp, #countData的内容是csv数据中的counts
                           colData = f,#colData数据由f表格提供
                           design = ~ condition)#design以f表格中的condition提供。
dep=DESeq(dds)
res = results(dep, contrast = c("condition",#contrast是对比的意思，即分子/分母。
                                'GRASLND shRNA',#实验组是分子。
                                'scrambled shRNA'))#对照组是分母。

DEG1 <- as.data.frame(res)#DEG是“差异表达基因”。内有p值和Foldchange。

write.csv(DEG1, file="results.csv")

#将差异表达基因矩阵以p从低到高（从差异到无明显差异排序）
resOrdered=res[order(res$pvalue),]

sum(res$padj<0.1, na.rm = T)#反馈一下p.adj<0.1的数目。
sum(res$pvalue<0.05, na.rm = T)#反馈一下p.value<0.05的数目。
sum(!is.na(res$pvalue))

resSig=subset(resOrdered, padj<0.1)#挑选出p.adj<0.1的DEG子集，取名DEG2。
DEG2=as.data.frame(resSig)
write.csv(resSig, file='results1.csv')

#绘制热图
library(pheatmap)
choose_gene=head(rownames(DEG2),50)#取top 50的基因名字。
choose_matrix=exp[choose_gene,]#我这已经没有exp2了
choose_matrix=t(scale(t(choose_matrix)))#转置、归一化处理、转置回原始表格。
pheatmap(choose_matrix, filename = "DEG_top50_heatmap.png")

#绘制火山图
DEG1=na.omit(DEG1)
logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   
DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))


library(ggplot2)
g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   #这里将pvalue取负对数
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     #绘制点图
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +    #轴标签
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   #设定颜色
ggsave(g,filename='volcano.png')
# 如下图，一般会重点关注两边的，偏上方的基因点

library(clusterProfiler)
mydata <- read.table('results1.csv',header = T,
                     sep=',', stringsAsFactors = F)
genelist = mydata$X

BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
# 安装加载人类go注释包
go.all <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',
                   pAdjustMethod = 'BH', pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2, keyType = "SYMBOL") #keyType分为SYMBOL、ENTREZID、ENSEMBL
dim(go.all)
head(go.all)
go.all.df <- as.data.frame(go.all)
barplot(go.all, showCategory = 50, drop=T)
dotplot(go.all, showCategory = 20)

BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
go.BP <- enrichGO(genelist, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,
                  keyType = 'SYMBOL')
plotGOgraph(go.BP)

go.MF <- enrichGO(genelist, 
                  OrgDb = org.Hs.eg.db, 
                  ont='MF',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,
                  keyType = 'SYMBOL')
plotGOgraph(go.MF)

go.cc <- enrichGO(genelist, 
                  OrgDb = org.Hs.eg.db, 
                  ont='CC',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,
                  keyType = 'SYMBOL')
plotGOgraph(go.cc)

