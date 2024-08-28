rm(list = ls()) #清空
setwd('D:\\Project\\AU\\a') #设置工作文件夹

# Step1 数据下载-------------------------------------------------------------------
library(GEOquery) 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
options(timeout = 300) #设置超时时间
gse = "GSE37265" # 设置要下载的数据集ID
eSet <- getGEO(gse,
               destdir = '.',
               getGPL = F)  #下载指定GSE ID的数据集，并存储在当前目录，不下载GPL文件


exp <- exprs(eSet[[1]]) #从下载的数据集中提取表达矩阵

fivenum(exp) #查看表达数据的五数概括
summary(exp) #查看表达数据的概要统计信息
boxplot(exp) #绘制表达数据的箱线图
exp[1:4,1:4] #查看表达矩阵的前4行4列数据

# exp <- log2(exp+1) #对表达数据进行log2转换以减少数据的偏态
# fivenum(exp) #查看转换后的五数概括
# summary(exp) #查看转换后的概要统计信息
# boxplot(exp) #绘制转换后的箱线图

pd <- pData(eSet[[1]]) #提取样本信息数据
gpl <- eSet[[1]]@annotation #提取平台信息
save(gse,exp,pd,gpl,file = "step2_output.Rdata")


# Step2 数据预处理-------------------------------------------------------------------
rm(list = ls()) #清空当前的工作环境
load("step2_output.Rdata") #加载保存的数据
library(stringr) #install.packages("stringr")
library(dplyr) #install.packages("dplyr")
summary(exp) 
boxplot(exp) 
View(pd)  #查看临床信息
options(timeout = 3000)  #设置超时选项

# 条件语句，如果条件为真，则执行内部代码块
if(T){
  a = getGEO(gpl,destdir = ".") #从 GEO 数据库获取平台数据
  b = a@dataTable@table         #获取平台数据表格
  colnames(b)                   #查看表格的列名
  ids2 = b[,c("ID","Gene Symbol")]  #提取 ID 和 Gene Symbol 列，并重命名列名
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]  #过滤掉 symbol 列为空或包含 "///" 的行
}

#将表达矩阵转换为数据框，并添加探针 ID 作为一列
exp1 <- data.frame(exp)
exp1 <- mutate(exp1, probe_id = rownames(exp))

#将表达数据框与注释数据框按 probe_id 合并
exp1 <- merge(x = exp1, y = ids2, by = "probe_id")

#移除 symbol 列中重复的行，确保每个基因只保留一个表达值
exp1 <- exp1[!duplicated(exp1$symbol),]

#将行名设置为 symbol，并移除 probe_id 和 symbol 列
row.names(exp1) <- exp1$symbol
exp1 <- exp1[,-1]          #移除 probe_id 列
exp1 <- exp1[,-ncol(exp1)] #移除 symbol 列
save(exp1,pd, file = "step3_output.Rdata")


# Step3 数据分类-------------------------------------------------------------------
rm(list = ls())
load("step3_output.Rdata")
rm(exp, ids2, group_list)
pd$group <- NA
pd$group <- ifelse(pd$`condition:ch1` == "Ulcer", "Ulcer","Normal")
table(pd$group)
pd$group <- factor(pd$group, levels = c("Normal", "Ulcer"))

pd <- pd[order(pd$group), ]
group_list <- pd$group

match_positions <- match(pd$geo_accession, colnames(exp1))
exp2 <- exp1[, match_positions]

save(group_list,exp2,pd,pd, file = "step4_output.Rdata")



# Step4 PCA分析-------------------------------------------------------------------
rm(list = ls()) #清空
load("step4_output.Rdata")
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggsci)
cors <- pal_lancet()(5)

dat=as.data.frame(t(exp2))
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", #使用点来表示样本
                         col.ind = group_list, #根据分组进行着色
                         palette = cors,###色彩颜色根据分组个数决定
                         addEllipses = TRUE, #添加椭圆，以表示各组样本的分布区域
                         legend.title = "Groups") #图例标题
print(pca_plot)

save(pca_plot,file = 'pca.Rdata')


# Step5 差异分析-------------------------------------------------------------------
rm(list = ls())
load("step4_output.Rdata")
library(limma)
library(tidyverse)
exp2 <- normalizeBetweenArrays(as.matrix(exp2), method = "quantile")
design=model.matrix(~group_list)
fit=lmFit(data.frame(exp2), design)
fit=eBayes(fit)
deg=topTable(fit, coef=2, number = Inf)
deg <- mutate(deg,probe_id=rownames(deg))
deg <- deg[order(deg$logFC, decreasing = TRUE), ]

logFC_t= 1.5
change=ifelse(deg$P.Value>0.05,'Stable',
              ifelse(abs(deg$logFC) < logFC_t,'Stable',
                     ifelse(deg$logFC >= logFC_t,'Up','Down') ))
table(change)
deg <- mutate(deg, change)
table(deg$change)

# write.csv(data.frame(deg),'degs.csv')

library(ggsci)
cors <- pal_lancet()(3)
fivenum(deg$logFC)
fivenum(-log10(deg$P.Value))

library(ggplot2)
volcano_plot <- ggplot(deg,aes(logFC,
               -log10(P.Value)))+
  geom_point(size = 2.5, 
             alpha = 0.8, 
             aes(color = change),
             show.legend = T)+
  scale_color_manual(values = c('#00468BFF','gray','#ED0000FF'))+
  ylim(0, 15)+
  xlim(-10, 10)+
  labs(x = 'LogFC',y = '-Log10(P.Value)')+
  geom_hline(yintercept = -log10(0.05),
             linetype = 2,
             color = 'black',lwd = 0.8)+
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = 2, 
             color = 'black', lwd = 0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
volcano_plot

ggsave(filename = "volcano_plot.tiff", plot = volcano_plot, 
       device = "tiff", dpi = 300,
       width = 9, height = 8)


deg_gene <- deg$probe_id[deg$change %in% c('Up','Down')]
# write.csv(data.frame(deg_gene),'deg_gene.csv')

#Heatmap Plot
library(pheatmap)#install.packages("pheatmap")
library(ggplot2)
deg_genes <- deg %>%
  filter(change %in% c('Up', 'Down'))

exp_deg <- exp2[rownames(exp2) %in% deg_genes$probe_id,]
table(deg_genes$change)

library(ComplexHeatmap)

exp_hot <- apply(exp_deg, 1, scale)
rownames(exp_hot) <- colnames(exp_deg)
exp_hot <- t(exp_hot)

library(circlize)
col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("#00468BFF", "white", "#ED0000FF")
)
table(group_list)
heatmap_plot <- Heatmap(exp_hot, name = "Expression", cluster_rows = T, row_km = 2, 
                        col = col_fun ,show_column_names = F , show_row_names = F, 
                        row_gap = unit(2, "mm"),  # 关闭行之间的间隙
                        column_gap = unit(2, "mm"),  # 关闭列之间的间隙
                        clustering_distance_columns = 'pearson', row_title  = "Genes", 
                        row_title_side = "left" , column_km = 2, column_title = NULL ,
                        top_annotation = HeatmapAnnotation(
                          foo = anno_block(# 设置填充色
                            gp = gpar(fill = c("#cD0000FF","#00468BFF")),
                            # 设置文本标签
                            labels = c("Ulcer", "Normal"), 
                            # 文本标签样式
                            labels_gp = gpar(col = "white", fontsize = 10)
                          )))
print(heatmap_plot)

tiff(filename = "heatmap_plot.tiff", width = 9, height = 8, units = "in", res = 300)  # Set resolution as needed
draw(heatmap_plot)
dev.off()
ggsave(filename = "heatmap_plot.tiff", plot = heatmap_plot, device = "tiff", width = 9, height = 8)
save(deg,heatmap_plot,file = 'volcano_heatmap.Rdata')
save(deg,deg_gene,file = 'enrich_analysis.Rdata')
save(exp2,pd,group_list,deg, file = "step5_output.Rdata")


# Step6 富集分析-------------------------------------------------------------------
rm(list = ls())
load('enrich_analysis.Rdata')
library(clusterProfiler)
library(org.Hs.eg.db)#将Symbol转换为基因符号Entrez ID
library(dplyr)
library(DOSE)
library(ggplot2)
library(tidyr)
gene <- deg_gene

diff_entrez<-bitr(
  gene,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Hs.eg.db'
)

head(diff_entrez)

###GO
go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",#p值矫正方法
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)

go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')

go_geo<- simplify(go_enrich, cutoff=0.7, by="p.adjust",
                  select_fun=min)

go_result<-go_geo@result
head(go_result)

#dotplot
go_dot <- dotplot(go_geo,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory=10,
        split='ONTOLOGY',
        label_format = Inf)+#不换行
  #分面
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
            )
go_dot
ggsave("go_dot.tiff", plot = go_dot, device = "tiff", 
       width = 16,  height = 9, 
       units = "in",dpi = 300)

#barplot
go_bar <- barplot(go_geo,
        x = "Count",
        color = "p.adjust",
        showCategory=10,
        split='ONTOLOGY',
        label_format = Inf)+
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )
go_bar
ggsave("go_bar.tiff", plot = go_bar, device = "tiff", 
       width = 16,  height = 9, 
       units = "in",dpi = 300)

###KEGG富集
KEGG_enrich <- clusterProfiler::enrichKEGG(gene = diff_entrez$ENTREZID,
                                           organism = "hsa", #物种Homo sapiens 
                                           pvalueCutoff = 0.05,#pvalue阈值
                                           qvalueCutoff = 0.05,#qvalue阈值
                                           pAdjustMethod = "BH",#p值矫正方法，one of "holm", 
                                           #"hochberg", "hommel", 
                                           #"bonferroni", "BH", 
                                           #"BY", "fdr", "none"
                                           minGSSize = 10,#富集分析中考虑的最小基因集合大小
                                           maxGSSize = 500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbol
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')
#提取KEGG富集结果表格                 
KEGG_result<-KEGG_enrich@result
head(KEGG_result)

#dotplot
kegg_dot <- dotplot(KEGG_enrich,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10)#展示top10通路
kegg_dot
ggsave("kegg_dot.tiff", plot = kegg_dot, device = "tiff", 
       width = 8,  height = 6, 
       units = "in",dpi = 300)

#barplot
kegg_bar <- barplot(KEGG_enrich,
        x = "Count",
        color = "p.adjust",
        showCategory = 10)
kegg_bar

ggsave("kegg_bar.tiff", plot = kegg_bar, device = "tiff", 
       width = 8,  height = 6, 
       units = "in",dpi = 300)

###GSEA

symbol <- rownames(deg)
entrez <- bitr(symbol, 
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = 'org.Hs.eg.db'
)
head(symbol)


gene_list <- deg$logFC
names(gene_list) <- rownames(deg)

gene_list <- gene_list[names(gene_list) %in% entrez[,1]]
names(gene_list) <- entrez[match(names(gene_list), entrez[, 1]), 2]
length(gene_list)
head(gene_list)

gene_list <- sort(gene_list, decreasing = T)
head(gene_list)


#KEGG-GSEA
KEGG_gse <- gseKEGG(geneList = gene_list, 
                    organism = "hsa", 
                    minGSSize = 10, 
                    maxGSSize = 500, 
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH", 
                    verbose = FALSE, 
                    eps = 0)
KEGG_gse <- setReadable(KEGG_gse, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID")
KEGG_gse_result <- KEGG_gse@result
# write.csv(KEGG_gse_result, file = c('gsea(KEGG).csv'))
#可视化
library(ggridges)
library(ggplot2)
library(enrichplot)

###山峦图
kegg1 <- ridgeplot(KEGG_gse, 
                   showCategory = 10, 
                   fill = "p.adjust", 
                   decreasing = T)
kegg1

###气泡图
kegg2 <- dotplot(KEGG_gse,
                 x = "GeneRatio",
                 color = "p.adjust",
                 showCategory = 10)
kegg2

#ES图绘制(gseaplot2函数)
#单个基因集可视化：
gsea_kegg <- gseaplot2(KEGG_gse,
                       geneSetID = 1,
                       color = "red",
                       rel_heights = c(1.5, 0.5, 1), #子图高度
                       subplots = 1:3, #显示哪些子图
                       pvalue_table = F, #是否显示pvalue表
                       title = KEGG_gse$Description[1],
                       ES_geom = "line") #"dot"将线转换为点
gsea_kegg

ggsave("gsea_kegg.tiff", plot = gsea_kegg, device = "tiff", 
       width = 10,  height = 6, 
       units = "in",dpi = 300)

library(ggplot2)
cors <- pal_lancet("lanonc",alpha = 0.6)(9)
cors
#同时可视化多个基因集：
p2 <- gseaplot2(KEGG_gse,
                geneSetID = 1:3, #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = cors[1:3],
                pvalue_table = TRUE,
                ES_geom = "line")
p2


# Stet7 WGCNA分析-------------------------------------------------------------------
rm(list = ls())
load("step5_output.Rdata")

library(tidyverse)
exp_wgcna <- exp2[,colnames(exp2) %in% pd$geo_accession]
samples <- pd[,c("geo_accession","group")]
table(samples$group)
samples$group <- ifelse(samples$group== "Ulcer","1","0")
rownames(samples) <- seq_len(nrow(samples))

library(WGCNA)
#筛选基因
dataExpr <- exp_wgcna
m.mad <- apply(dataExpr, 1, mad)
dataExprVar <- dataExpr[which(m.mad >
                                quantile(m.mad, probs=0.25)),]
dataExpr <- as.data.frame(t(dataExprVar))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

if(!gsg$allOK){
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)

sampleTree = hclust(dist(dataExpr), method = "average")
par(cex = 0.9);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(dataExpr)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 110”和“cutHeight = 110”中换成你的cutoff
  abline(h = 125, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 125, minSize = 10)
keepSamples = (clust==1)
datExpr = dataExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr=as.data.frame(datExpr)
sample_to_delete <- c("GSM1337307","GSM1337308")
samples <- samples[!samples$geo_accession %in% sample_to_delete,]

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers,
                        verbose = 5 )
tiff("1Threshold.tiff",width = 9*300, height = 5*300, res = 300)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft$powerEstimate
#构建网络，找出gene module
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 200,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE, 
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)
tiff("2module.tiff",width = 10*300, height = 5*300, res = 300)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

#把gene module输出到文件
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),
             datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

#表型与模块的相关性
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)

group_list <- ifelse(samples$group == '0', 'Normal', 'Ulcer')
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=rownames(datExpr)
modTraitCor = cor(MEsWW,design, use = "p")

colnames(MEsWW)
modlues=MEsWW
nSamples <- ncol(datExpr)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")

dim(textMatrix) = dim(modTraitCor)

tiff("3Module-trait.tiff",width = 6*300, height = 6*300, res = 300)
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(design), 
               yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), 
               colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

ACC = as.data.frame(samples$group)
names(ACC) = "ACC"
modNames = substring(names(MEs), 3)
modNames

# 计算基因表达数据 (datExpr) 和模块特征基因 (MEs) 之间的相关性
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# 计算基因模块成员相关性的 p 值
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), ncol(datExpr)))

# 给 geneModuleMembership 和 MMPvalue 添加列名，表示模块成员和相应的 p 值
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# 计算基因表达数据 (datExpr) 和表型数据 (ACC) 之间的相关性
geneTraitSignificance = as.data.frame(cor(datExpr, ACC, use = "p"))

# 计算基因表型相关性的 p 值
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ncol(datExpr)))

# 给 geneTraitSignificance 和 GSPvalue 添加列名，表示基因显著性和相应的 p 值
names(geneTraitSignificance) = paste("GS.", names(ACC), sep="")
names(GSPvalue) = paste("p.GS.", names(ACC), sep="")

# table(mergedColors)
colnames(modlues)
colnames(geneModuleMembership)

# 选择 turquoise 模块中的基因
moduleGenes = moduleColors=="yellow"
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
tiff("4Module membership.tiff",width = 6*300, height = 6*300, res = 300)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, "MM5"]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in MMyellow module"),
                   ylab = "Gene significance for Ulcer status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "yellow")
#abline(h=0.4,v=0.8,col="red",lwd=1.5)#筛选线，用于筛选模块中具有高度相关性和显著性的基因
dev.off()



# Step8 机器学习 ------------------------------------------------------------------
rm(list = ls())
load("step4_output.Rdata")
gsea <- read.csv("gsea(KEGG).csv")
genes <- gsea$core_enrichment[1]
gene_list <- strsplit(genes, "/")[[1]]

# 提取交集基因的表达数据
exp_ml <- exp2[rownames(exp2) %in% gene_list,]
ml <- data.frame(t(exp_ml))

samples <- pd[,c("geo_accession","group")]
table(samples$group)
samples$group <- ifelse(samples$group== "Ulcer","1","0")
rownames(samples) <- seq_len(nrow(samples))              

# 构建包含样本分组信息的数据集 x
x <- data.frame(cbind(group = samples$group, ml))
# 将分组信息转换为因变量 y，并将 x 转换为矩阵形式以用于模型训练
y=data.matrix(x$group)
y <- as.factor(y)
x <- as.matrix(x[,-1])# 去除掉第一列的分组信息，因为它只是标签而非特征

###Lasso
library(glmnet)
# 设置随机种子
set.seed(110000)
fit=glmnet(x,y,family = "binomial",maxit = 100000, nfold = 10)
tiff("Lasso1.tiff",width = 6*300, height = 6*300, res = 300)
plot(fit,xvar="lambda",label = TRUE)
dev.off()

cvfit = cv.glmnet(x,y,family="binomia",maxit = 100000, nfold = 10)
tiff("Lasso2.tiff",width = 6*300, height = 6*300, res = 300)
plot(cvfit)
dev.off()


coef=coef(fit,s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)#查看模型的相关系数geneCoef
geneCoef

lassoGene <- lassoGene[-1]
actCoef<- actCoef[-1]
lassoGene
###randomForest
library(randomForest)
library(caret)
set.seed(110)
ctrl <- trainControl(method = "cv", number = 10)
rf <- randomForest(y~.,  data = x , ntree = 200 , trainControl = ctrl)

# rf <- randomForest(y~.,  data = x , ntree = 200)
tiff("RF1.tiff",width = 6*300, height = 6*300, res = 300)
plot(rf, main = 'Random Forest', lwd = 2, ylim = c(0,1))
dev.off()

optionTrees = which.min(rf$err.rate[, 1])
#rf2 = randomForest(y~., data = x, ntree = optionTrees, importance = T)
rf2 = randomForest(y~., data = x, ntree = optionTrees)
importance = importance(x = rf2)

tiff("RF2.tiff",width = 6*300, height = 6*300, res = 300)
varImpPlot(rf2, main = 'Feature Importance')
dev.off()
rfGenes = importance[order(importance[, 'MeanDecreaseGini'], decreasing = T), ]

## 使用中位数作为阈值选择重要特征
median_gini = median(importance[, 'MeanDecreaseGini'])
rfGenes = names(rfGenes[rfGenes > median_gini])
rfGenes

both <- intersect(lassoGene, rfGenes)
both

save(both, exp_ml, samples,group_list, file = 'machine_learn.Rdata')

#####Boruta##
# library(Boruta) ##install.packages("Boruta")
# set.seed(1124)
# x1$group <- as.factor(x1$group)
# 
# Var.Selec<-Boruta(
#   group~.,
#   x1,
#   pValue = 0.01, #confidence level. 可信水平
#   mcAdj = TRUE, #是否使用Bonferroni调整
#   #if set to TRUE, a multiple comparisons adjustment using the Bonferroni method will be applied.
#   maxRuns = 500, #迭代最多次数
#   doTrace =0,#可以选0-3，运行结果显示的详细程度，0不显示轨迹
#   holdHistory = TRUE, #如果设置为TRUE，则存储完整的重要性历史记录，并将其作为结果的ImpHistory元素返回。
#   getImp = getImpRfZ #用于获取属性重要性的函数。默认值是 getImpRfZ，它从 ranger 包运行随机森林并收集平均降低精度测量的 Z 分数。
# )
# Var.Selec
# 
# #（绿色是重要的变量，红色是不重要的变量，蓝色是影子变量，黄色是Tentative变量）。
# tiff("Boruta1.tiff",width = 6*300, height = 6*300, res = 300)
# plotImpHistory(Var.Selec,
#                ylab="Z-Scores",
#                las=2)
# dev.off()
# 
# tiff("Boruta2.tiff",width = 6*300, height = 6*300, res = 300)
# plot(Var.Selec,
#      whichShadow=c(F,F,F),
#      xlab="",
#      ylab="Z-Scores",
#      las=2)
# dev.off()
# 
# getConfirmedFormula(Var.Selec)
# getNonRejectedFormula(Var.Selec)
# boruta_genes <- getSelectedAttributes(Var.Selec,withTentative=FALSE)
# boruta_genes


#Veen
rm(list = ls())
load('machine_learn.Rdata')
genes <- both
roc_exp <- as.matrix(exp_ml[genes, ])
table(group_list)
roc_exp<-exp_ml[match(genes,rownames (exp_ml)),]
samples$group <- ifelse(samples$group == "0", 'Normal', 'Ulcer')
samples$group <- factor(samples$group)
library(ggsci)
cors <- pal_lancet()(length(genes))

library(pROC)
roc_exp1 <- as.data.frame(t(roc_exp))

roc_exp1$response = samples$group
auc_values <- sapply(roc_exp1[, -ncol(roc_exp1)], function(x) auc(roc(roc_exp1$response, x)))
sorted_indices <- order(auc_values, decreasing = TRUE)

tiff("ROC.tiff",width = 6*300, height = 6*300, res = 300)
plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), 
     xlab="1 - Specificity (FPR)", ylab="Sensitivity (TPR)", 
     main="ROC Curves", 
     cex.lab=1.2, cex.axis=1.1, cex.main=1.3)
grid(lty = "dotted", col = "lightgray")
abline(a=0, b=1, lty=2, col="darkgray")
aucText <- sapply(sorted_indices, function(i) {
  roc1 <- roc(roc_exp1$response, as.vector(roc_exp1[,i]))
  lines(roc1, col=cors[i], lwd=2)
  paste0(colnames(roc_exp1)[i], ", AUC=", sprintf("%0.3f", auc(roc1)))
})
legend("bottomleft", aucText, lwd=2, bty="n", col=cors[sorted_indices], 
       cex=0.8, text.font=3)
dev.off()
