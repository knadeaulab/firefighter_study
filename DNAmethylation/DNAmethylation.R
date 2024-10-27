setwd("R_codes/DNAmethylation")
library(EnhancedVolcano)
library(limma)
library(circlize)
library(annotatr)
library(GenomicRanges)
library(regioneR)
library(reshape)
library(randomcoloR)
library(RColorBrewer)
library(dplyr)
options(bedtools.path = "/opt/homebrew/bin/")
library(bedtoolsr)
library(bios2mds)
library(magrittr)
library(ggpubr)
library(gridBase)
library(grid)
library(ComplexHeatmap)
library("viridis")
library(minfi)
library(dplyr)
library(PCAtools)
library(ENmix)
load("DNAmethylation_V3.RData") ## Loading pre-processed dataset ##
#### Section-1: PCA analysis ####
colnames(Mdata) == rownames(metadata) ## Mdata variable contrins M-values for each CpG site
x_vars = Mdata
y_var = metadata[colnames(x_vars),c("BMI","age","Race","Group")]
rownames(y_var) == colnames(x_vars) ## should be all TRUE
p <- PCAtools::pca(x_vars, metadata = as.data.frame(y_var), removeVar = 0.1)
out = eigencorplot(p, components = getComponents(p, seq_len(5)),
                   metavars = c('BMI','age', 'Race', 'Group'), ##
                   main  = "PCA analysis", scale = F, corMultipleTestCorrection = "BH")

#### Section-2: Differential methylation analysis ####
phenoData <- new("AnnotatedDataFrame", data=metadata)
eSet <- Biobase::ExpressionSet(assayData=Mdata, phenoData=phenoData)
eSet$Group <- factor(eSet$Group)
tmp = ifelse(eSet$Race=="Caucasian", "Caucasian", "other")
tmp[is.na(tmp)]<- "other"
eSet$Race2 <- factor(tmp)
design <- model.matrix(~ eSet$Group + eSet$age)
fit <- lmFit(eSet, design)
fit <- eBayes(fit)
res = topTable(fit, n = Inf)
res$logFC = res$eSet.Group1
## volcano plot 
vp = EnhancedVolcano::EnhancedVolcano(toptable = res, x = 'logFC', y = 'adj.P.Val', lab = rownames(res), 
                                      labSize=0,  FCcutoff = 0.5, pCutoff = 0.2, xlim = c(-2, 2.3), ylim = c(0,15),
                                      pointSize = 4)

res_sig = dplyr::filter(res, abs(logFC) > 0.5 & adj.P.Val < 0.20 )

#### Section-3: Manhattan plot ####
mm = match(rownames(res), locinfo$probes)
res$chr = locinfo$chr[mm]
res$start = locinfo$start[mm]
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
         'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
         'chr18','chr19','chr20','chr21','chr22')

colrs = randomcoloR::randomColor(22, luminosity = "dark")
res$chr = factor(res$chr, levels = chrs )
manhat = ggplot(res, aes(x = chr, y = -log10(adj.P.Val), color = as.factor(chr))) +
  geom_jitter(alpha = 0.75, width = 0.4, size = 1.2) +
  scale_color_manual(values = colrs) +
  geom_hline(yintercept = -log10(0.01), color = "red", linetype = "dashed") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.35))

#### Section-4: MDS plot ####
res_sig = dplyr::filter(res, abs(logFC) > 0.5 & adj.P.Val < 0.10 )
##mds  = t(FF.M[rownames(res_sig),]) %>% dist() %>% cmdscale()
mds  = t(Mdata[rownames(res_sig),]) %>% dist() %>% cmdscale()
colnames(mds) <- c("Dim.1", "Dim.2")
mds = as.data.frame(mds)
mm = match(rownames(mds), metadata$ID)
mds$Group = metadata$Group[mm]
mds$female = metadata$female[mm]

mdsplot = ggplot(mds, aes(x = Dim.1, y = Dim.2, color = as.factor(Group))) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c('#4363d8','#f58231')) +
  theme_classic(base_size = 16)


#### Section-5: GREAT tool for CpG site-> gene-ID association analysis ####
## Note: GREAT (http://great.stanford.edu/public/html/index.php) 
## can be executed using web-browser. The following section of code 
## will generate two bed files that can be used as input for GREAT 
## web-browser. The output will include distance from TSS also.
## If you want, you can skip this section, please directly move to next section

set.seed(123)
dd <- toGRanges(locinfo)
annots = c('hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries')
annotations = build_annotations(genome = 'hg19', annotations = annots)
dm_annotated = annotate_regions(
                  regions = dd,
                  annotations = annotations,
                  ignore.strand = TRUE,
                  quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
dm_annsum = summarize_annotations(annotated_regions = dm_annotated,quiet = TRUE)
print(dm_annsum)
bedout = data.frame()
for (feature in unique(df_dm_annotated$annot.type)){
  message(feature)
  tmp = dplyr::filter(df_dm_annotated, annot.type == feature)
  tmp = tmp[,c("seqnames","start", "end")] %>% dplyr::distinct()
  ## merge regions that fall under 1000 bp window and make it a single entry
  tet= bt.cluster(i = tmp, d = 1000)
  a = tet[,c("V1","V2","V4")] %>% group_by(V1,V4) %>% dplyr::summarise_all(min)
  b = tet[,c("V1","V2","V4")] %>% group_by(V1,V4) %>% dplyr::summarise_all(max)
  tmp = cbind(a[,c("V1","V2")],a[,"V2"] + 1)
  colnames(tmp) = c("chr","start","end")
  tmp$group = gsub("hg19_", "", feature)
  bedout = rbind(bedout, tmp)
}
mm = match(rownames(res_sig), df_dm_annotated$probes)
res_sig$chr = df_dm_annotated$seqnames[mm]
res_sig$start = df_dm_annotated$start[mm]
res_sig$end = df_dm_annotated$end[mm]
a = res_sig[,c("chr", "start", "end")] ## foreground for pathway analysis 
b = locinfo[,c("chr", "start", "end")] ## background 
colnames(a) = NULL
colnames(b) = NULL
write.table(a, file = "significantSite.bed", sep="\t", quote = F, row.names = F)
write.table(b, file = "background.bed", sep="\t", quote = F, row.names = F)
## the above two files were used as input for http://great.stanford.edu/public/html/index.php
## from the result page of GREAT, I extracted the region-geneID association excel file- GREAT_output_association_table.xlsx
res_sig = dplyr::filter(res, adj.P.Val < 0.05 ) # using only q < 0.01 as threshold #
mm = match(rownames(res_sig), df_dm_annotated$probes)
res_sig$chr = df_dm_annotated$seqnames[mm]
res_sig$start = df_dm_annotated$start[mm]
res_sig$end = df_dm_annotated$end[mm]
a = res_sig[,c("chr", "start", "end")] ## foreground for new pathway analysis 
colnames(a) = NULL
write.table(a, file = "significantSite_2.bed", sep="\t", quote = F, row.names = F)
## the above file was also used as input for http://great.stanford.edu/public/html/index.php
## from the result page of GREAT, I extracted the region-geneID association excel file- GREAT_output_association_table.xlsx
GL = read.table("../Files/DataF2C.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "Gene_Description.tsv", sep="\t", quote = F, row.names = F)

GL1 = read.table("../Files/DataF2D.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "Gene_Description2.tsv", sep="\t", quote = F, row.names = F)


#### Section-6: upload gene list from association file to metascape (https://metascape.org/gp/index.html)  ####
CpGassociatedGenes = read.table("../Files/DataF2C.txt", sep="\t")
CpGassociatedGenes2 = read.table("../Files/DataF2D.txt", sep="\t")
write.table(file="CpGassociatedGenes4Metascape.txt", CpGassociatedGenes$V1, quote = F, row.names = F)
write.table(file="CpGassociatedGenes4Metascape_withQ_less_than_0.05.txt", CpGassociatedGenes2$V1, quote = F, row.names = F)
## upload these gene list to metascape web-browser to replicate figures and supplementary info 


#### Section-7: circos plot ####
FFid = rownames(metadata)[metadata$Group==1]
ctrlid = rownames(metadata)[metadata$Group!=1]
a = rowMedians(Mdata[,FFid])
b = rowMedians(Mdata[,ctrlid])
a1 = data.frame(Methsite = rownames(Mdata), median = a)
b1 = data.frame(Methsite = rownames(Mdata), median = b)
mm = match(a1$Methsite, locinfo$probes)
a1 = cbind(locinfo[mm,c("chr","start","end")], a1)
mm = match(b1$Methsite, locinfo$probes)
b1 = cbind(locinfo[mm,c("chr","start","end")], b1)
a1$Methsite = NULL
colnames(a1)  = c("chr","start","end","value1")
b1$Methsite = NULL
colnames(b1)  = c("chr","start","end","value1")
target = c("genes_promoters","genes_intergenic","genes_introns","genes_exons","genes_3UTRs","genes_5UTRs" )
colrs=randomcoloR::randomColor(count = length(target), luminosity="dark" )
colrs = c('#e6194b','#3cb44b','#f032e6','#4363d8','#f58231','#911eb4')
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.clear()
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
#circos.genomicHeatmap(a1, col = col_fun, side = "outside",connection_height = NULL, heatmap_height = 0.07) ## FF 
circos.genomicHeatmap(b1, col = col_fun, side = "outside",connection_height = NULL, heatmap_height = 0.07) ## controls
for ( i in 1:length(colrs)){
  gene_type = dplyr::filter(bedout,group %in% c(target[i]))
  cnt = length(unique(gene_type$group))
  gene_type$group = as.numeric(factor(gene_type$group, labels = cnt))
  circos.genomicDensity(gene_type, col = colrs[i], track.height = 0.05)
}
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
upViewport()
h = dev.size()[2]
lgd_meth = Legend(title = "Methylation", col_fun = col_fun)
lgd_list = packLegend(lgd_meth, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")
circos.clear()

#### Section-8: PFAS analysis ####
centrix1 = masterlist$Sentrix_Position.x
centrix0 = masterlist$Sentrix_ID
colms = paste(centrix0,"_",centrix1, sep="")
masterlist$ID = colms
common  = intersect(colnames(FF.M),colms) ## all IDs matched 
metadata00 = masterlist[masterlist$ID %in% common,] ## complete metadata of all samples available in methylation data I have. 
metadata00$Group = ifelse(metadata00$protocol.name %in% c("Firefighter"),"FireFighter", "Control")
rownames(metadata00) = metadata00$ID

pfas3 = pfas2[pfas2$X %in% metadata00$PPID,]
mm = match(pfas3$X, metadata00$PPID)
pfas3$Status = metadata00$Group[mm]

X = FF.M[,metadata00$ID]
colnames(X) = metadata00$PPID

com = intersect(colnames(X) , pfas3$X)
mm = match(com,pfas3$X)
pfas3 = pfas3[mm,]
mm = match(com,colnames(X)) ## if all TRUE we are good to go 
X2 = X[mm]
colnames(X2)==pfas3$X

TPFOA = pfas3$Total.PFOA
TPFOS = pfas3$Total.PFOS


getcorrPFOA <- function(x2){
  
  corr = cor.test(TPFOA,as.numeric(x2), method = "s")
  rho = corr$estimate
  pval = corr$p.value
  result <- correlationBF(TPFOA,as.numeric(x2))
  res = describe_posterior(result)
  out = 0
  if(res$ROPE_Percentage < 0.10 & res$BF > 1.0){
    if(pval < 0.05){
      out = res$Median
      out = rho
    }
  }
  return(out)
}
getcorrPFOS <- function(x2){
  
  corr = cor.test(TPFOS,as.numeric(x2), method = "s")
  rho = corr$estimate
  pval = corr$p.value
  out = 0
  result <- correlationBF(TPFOS,as.numeric(x2))
  res = describe_posterior(result)
  if(res$ROPE_Percentage < 0.10 & res$BF > 1.0){
    if(pval < 0.05){
      out = res$Median
      out = rho
    }
  }
  return(out)
}

## following bayesian analysis is extremely time-consuming. however Outcome variables (corrmat1,corrmat2) are already loaded 
#corrmat1 = apply(X2, MARGIN = 1, FUN=function(x2) getcorrPFOA(x2))
#corrmat2 = apply(X2, MARGIN = 1, FUN=function(x2) getcorrPFOS(x2))

PFOA_epigen = corrmat1[abs(corrmat1) >= 0.60]
PFOS_epigen = corrmat2[abs(corrmat2) >= 0.60] 
edgelist1 = as.data.frame(PFOA_epigen)
edgelist1$PFAS = "PFOA"
edgelist1$from = rownames(edgelist1)
colnames(edgelist1) = c("weight","from","to")

edgelist2 = as.data.frame(PFOS_epigen)
edgelist2$PFAS = "PFOS"
edgelist2$from = rownames(edgelist2)
colnames(edgelist2) = c("weight","from","to")
edgelist = rbind(edgelist1,edgelist2)
rownames(edgelist) =NULL
edgelist$weight = NULL


## Circos plot pfos in SE or firefighters ##
colr = randomcoloR::randomColor(count = 2, luminosity = c("random"))
names(colr) = c("PFOS","PFOS")

colr = c(PFOS="red",PFOA="darkblue")

circos.clear()
par(cex = 0.9, mar = c(0, 0, 0, 0))
#circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
chordDiagram(edgelist, transparency = 0, grid.col = colr, 
             annotationTrack = c("grid"), preAllocateTracks = list(
               track.height = 0.040,
               col = colr,
               track.margin = c(0.00, 0.00)))


## bed file for GREAT analysis ##
## for those cpG sites that are significantly 
## associated pfoa pfos 
mm = match(edgelist1$to, locinfo$probes)
bed1 = locinfo[mm,c("chr", "start", "end")]
colnames(bed1) = NULL
mm = match(edgelist2$to, locinfo$probes)
bed2 = locinfo[mm,c("chr", "start", "end")]
colnames(bed2) = NULL
write.table(bed1, file = "significantSite_PFOA.bed", sep="\t", quote = F, row.names = F)
write.table(bed2, file = "significantSite_PFOS.bed", sep="\t", quote = F, row.names = F)
## use these bed files represent CpG sites significantly associated with pfas (total pfoa & total pfos)
## along with previously created background.bed file to upload in great webserver ###
## also get thier distance from TSS and gene-ID association ##
## upload the gene list from gene-ID association file to metascape webserver 
## to get affected biological processes/pathways 
## Lets create a scatter plot between most  correlated sites and Total PFOA ##
BX2= ENmix::M2B(X2) ## CONVERT M values to B values for plotting 
ggdf = as.data.frame(PFOA_epigen)
ggdf$MethSite = rownames(ggdf)
ggdf$cor = abs(ggdf$PFOA_epigen)
ord = order(ggdf$cor, decreasing = T)
ggdf = ggdf[ord,]
methlation = as.data.frame(t(BX2[ggdf$MethSite[1], ]))
pfas3$X == rownames(methlation) ## if all TRUE we are good to go 
tmpsource = data.frame()
plots = list()
for (i in 1:length(ggdf$MethSite)){
  message(i)
  #si = "cg04432046_BC11"
  si = ggdf$MethSite[i]
  methlation = as.data.frame(t(BX2[si, ]))
  mydf = data.frame(X = methlation[,1], Y =  pfas3$Total.PFOA, G = "Firefighter")
  mod = cor.test(x = mydf$X, y = mydf$Y, method = "p")
  pval = mod$p.value
  mydf$CpGSite = si
  tmpsource = rbind(tmpsource, mydf)
  
  g = ggplot(mydf, aes(x = X, y = Y, fill = G)) + 
    geom_point(size = 2.0, aes(color= G)) + 
    geom_smooth(method = "lm", aes(fill=G), color = "#4D4D4D") + 
    xlab(label = colnames(methlation)[1] ) + ylab(label = "Total PFOA") +
    scale_fill_manual(values = c('#f58231')) +
    scale_color_manual(values = c('#f58231')) +
    theme(line= element_blank()) + theme_bw(base_size = 11)
  plots[[i]] = g
}
tmpsource$PFAS = "PFOA"
pdf("ScatterPlot_PFOA_Methylation.pdf", width = 2.8, height = 2.0)
for(pl in plots){
  print(pl)
}
dev.off()

### Lets create a scatter plot between most correlated sites and Total PFOS ###
ggdf = as.data.frame(PFOS_epigen)
ggdf$MethSite = rownames(ggdf)
ggdf$cor = abs(ggdf$PFOS_epigen)
ord = order(ggdf$cor, decreasing = T)
ggdf = ggdf[ord,]

methlation = as.data.frame(t(BX2[ggdf$MethSite[1], ]))
pfas3$X == rownames(methlation) ## if all TRUE we are good to go 

plots = list()

for (i in 1:length(ggdf$MethSite)){
  message(i)
  #si = "cg04432046_BC11"
  si = ggdf$MethSite[i]
  methlation = as.data.frame(t(BX2[si, ]))
  mydf = data.frame(X = methlation[,1], Y =  pfas3$Total.PFOS, G = "Firefighter")
  mod = cor.test(x = mydf$X, y = mydf$Y, method = "p")
  pval = mod$p.value
  mydf$CpGSite = si
  mydf$PFAS = "PFOS"
  tmpsource = rbind(tmpsource, mydf)
  
  g = ggplot(mydf, aes(x = X, y = Y, fill = G)) + 
    geom_point(size = 2.0, aes(color= G)) + 
    geom_smooth(method = "lm", aes(fill=G), color = "#4D4D4D") + 
    xlab(label = colnames(methlation)[1] ) + ylab(label = "Total PFOS") +
    scale_fill_manual(values = c('#f58231')) +
    scale_color_manual(values = c('#f58231')) +
    theme(line= element_blank()) + theme_bw(base_size = 11)
  plots[[i]] = g
}

pdf("ScatterPlot_PFOS_Methylation.pdf", width = 2.8, height = 2.0)
for(pl in plots){
  print(pl)
}
dev.off()


#### Section-9: FF experience vs DNAmethylation analysis ####
metadataFF = dplyr::filter(metadata, PPID %in% Experience$PID)
centrix1 = metadataFF$Sentrix_Position.x
centrix0 = metadataFF$Sentrix_ID
metadataFF$ID = paste(centrix0,"_",centrix1, sep="")
mm = match(metadataFF$PPID, Experience$PID)

metadataFF$Experience = Experience$Years.of.exposure[mm]
FFset = Mdata[,metadataFF$ID]

# Check if IDs match
if (!all(colnames(FFset) == metadataFF$ID)) {
  stop("Column names of FFset do not match IDs in metadataFF")
}

results <- data.frame(CpG_site = rownames(FFset), 
                      Estimate = numeric(nrow(FFset)),
                      Std_Error = numeric(nrow(FFset)),
                      t_value = numeric(nrow(FFset)),
                      p_value = numeric(nrow(FFset)))



# Run GLM for each CpG site
for (i in 1:nrow(FFset)) {
  #message(i)
  cpg_data <- data.frame(M_value = as.numeric(FFset[i, ]), Experience = metadataFF$Experience)
  
  # Fit GLM
  glm_model <- glm(M_value ~ Experience, data = cpg_data, family = gaussian())
  
  # Extract results
  model_summary <- summary(glm_model)
  coef_summary <- coef(model_summary)
  
  results$Estimate[i] <- coef_summary['Experience', 'Estimate']
  results$Std_Error[i] <- coef_summary['Experience', 'Std. Error']
  results$t_value[i] <- coef_summary['Experience', 't value']
  results$p_value[i] <- coef_summary['Experience', 'Pr(>|t|)']
}
results$fdr = p.adjust(results$p_value, method = "fdr")
results2 = dplyr::filter(results, p_value < 0.05 )
#write.csv(results2, "cpg_glm_results.csv", row.names = FALSE)

## check common sites ##
res_sig = dplyr::filter(res, adj.P.Val < 0.05)
com = intersect(rownames(res_sig), results2$CpG_site)
per = length(com)/nrow(res_sig) #24.42% (n = 107) of the sites assocaites with experience are also differentially methylated
a = res_sig[com,]

dim(a)
sum(a$eSet.Group1 < 0) ## 39 higher methylation (hyper-methylation) in firefighter
sum(a$eSet.Group1 > 0) ## 68 lower methylation (hypo-methylation) in firefighter

# xx=as.data.frame(Mdata["cg01901100_BC21",])
# mm = match(rownames(xx),metadata$ID)
# xx$group = metadata$Group[mm]
# colnames(xx)[1] = "m"
# mean(xx$m[xx$group == 0])
# mean(xx$m[xx$group == 1])

mm = match(com, df_dm_annotated$probes)
a$chr = df_dm_annotated$seqnames[mm]
a$start = df_dm_annotated$start[mm]
a$end = df_dm_annotated$end[mm]
a = a[,c("chr","start","end")]
colnames(a) = NULL
write.table(a, file = "significantSite_experience.bed", sep="\t", quote = F, row.names = F)
## use this bed file to represent CpG sites significantly Differentially methylated sites 
## associated with years of exposure, along with previously created background.bed file 
## to upload in great webserver to get their distance from TSS and gene-ID association ##
## upload the gene list from gene-ID association file to metascape webserver 
## to get affected biological processes/pathways 

GL1 = read.table("../Files/DataF2D.txt", sep="\t")
# Here is the output of GREAT, lets see names of the genes
gene_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = GL$V1, 
                                  columns = c("GENENAME"), 
                                  keytype = "SYMBOL")
write.table(gene_ids, file = "Gene_Description2.tsv", sep="\t", quote = F, row.names = F)



