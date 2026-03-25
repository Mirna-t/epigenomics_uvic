library(minfi)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)#    install.packages("BiocManager")

library(heatmap.2x) #devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")
library(gplots)
library(DMRcate)
library(dbplyr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)
library(missMethyl)

#devtools::install_version("dbplyr", version = "2.3.4")


setwd("/Users/agomez/Desktop/UVIC2018/GSE115382_RAW")
idat.folder <- "/Users/agomez/Desktop/UVIC2018/GSE115382_RAW" ###change for wour path dir

targets <- read.metharray.sheet(base=idat.folder)

###loading data

rgset <- read.metharray.exp(targets = targets,verbose = T)


####quality control first day

#The class of RGSet is a RGChannelSet object. 
#This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. 
#Note that this object contains the intensities of the internal control probes as well. 


phenoData <- pData(rgset)

##A MethylSet objects contains only the methylated and unmethylated signals

MSet <- preprocessRaw(rgset) 


#The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
qc <- getQC(preprocessRaw(rgset) )


plotQC(qc)

#To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color the densities by group:
densityPlot(preprocessRaw(rgset) , sampGroups = pData(MSet)$Status)

#Sex prediction

#By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X)med(X), 
#and the median total intensity of the Y-chromosome-mapped probes, denoted med(Y)med(Y),
#one can observe two different clusters of points corresponding to which gender the samples belong to.
#To predict the gender, minfi separates the points by using a cutoff on log2med(Y)



####detection pvals

detP<-detectionP(rgset)
#Plotting the mean detection p-value for each sample will allow us to gauge whether any samples have many failed probes 
#- this will be indicated by a large mean detection p-value. Samples with mean detection p-values exceeding a cutoff 
#such as 0.05 can be excluded from further analysis.

barplot(colMeans(detP),col=factor(pData(rgset)$Status),
        las=2,cex.names=0.8,main="Mean detection p-values")
###in case you want to remove samples

keep<-colMeans(detP) < 0.01
rgset<-rgset[,keep]
targets <- targets[keep,]
detP <- detP[,keep]


#####analysis


keep <- rowSums(detP < 0.01) == ncol(rgset)
table(keep)
#not.keep<-names(keep[keep=="FALSE"])


rgset <- rgset[keep,]


###Normalization

gRatioSet.illumina <- preprocessIllumina(rgset)##Illumina
gRatioSet.quantile <- preprocessQuantile(rgset,verbose=T)##Illumina

#####

#

ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann850k<-as.data.frame(ann850k)
#keep <- !(featureNames(gRatioSet.quantile) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])

#table(keep)

#gRatioSet.quantile <- gRatioSet.quantile[keep,]

## remove probes with SNPs at CpG or SBE site
#####get snp info
#snp.info<-getSnpInfo(gRatioSet.quantile)
####mapping

#gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)

####get values to proceed with analysis
betas<-getBeta(gRatioSet.quantile)

###M values for analysis
Mval<-getM(gRatioSet.quantile)

Mval <- rmSNPandCH(Mval, dist=2, 
                   mafcut=0.05,
                   rmcrosshyb = T,
                   rmXY = T) ### get out SNPs of my samples

##DMP finding
####LIMMA
library(limma)


design <- model.matrix(~pData(gRatioSet.quantile)$Status)##
colnames(design)[2]<-"Comp"
####generate adjusted beta matrix
fit <- lmFit(Mval, design)

fit2 <- eBayes(fit)

results <- topTable(fit2,
                    coef="Comp",
                    num=dim(fit2)[1],
                    sort.by="P",
                    adjust.method = "BH") ####Pval list 1

res<-subset(results,results$adj.P.Val<.05)


##dmpfinder option
group<-as.factor(pData(gRatioSet.quantile)$Status)
group<-relevel(group,ref="WT")


dmp <- dmpFinder(Mval, pheno = group, type = "categorical")

dmp<-subset(dmp,dmp$qval<.05)

##base R
dmp <- merge(dmp,
             ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")],
             by="row.names")
##dplyr

res<- res %>%
      mutate(Name=rownames(.)) %>%
      left_join(ann850k %>%
                  as.data.frame() %>% 
                  dplyr::select(Name,UCSC_RefGene_Name,Relation_to_Island,UCSC_RefGene_Group)
                  )

####distribution of logFC

hist(res$logFC)

####distribution of adj pval

hist(-log10(res$adj.P.Val))

###distribution of logFC

hist(res$logFC)


## Outputs
proms<-res[grep("TSS1500|TSS200|5'UTR|1stExon",res$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]


proms.island.filt <- proms.island %>%
                     dplyr::filter(logFC > 2.5 | logFC < -2.5)

write.csv(proms.island.filt,file="Results_proms_islands_logFCHigh.csv",row.names=T)

####Heatmap

colnames(betas)<-targets$Status

##subset data
my.data<- betas %>%
        as.data.frame() %>%
        dplyr::filter(rownames(.) %in% proms.island.filt$Name)

samples<-ifelse(pData(gRatioSet.quantile)$Status=="KO","red","blue1")

#####scale data

data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


spcol<-samples
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2x(as.matrix(data), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", 
           ColSideColors=spcol,
           trace="none", 
           dendrogram="column", 
           cexRow=1, cexCol=.7,
           main="WT vs KO Proms + Islands",
           # labCol=NA,
           labRow=NA, 
           density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
)


##GO enrichment
#####GO enrichment analysis


UP<- proms.island.filt  %>%
     dplyr::filter(logFC>0)
 
DOWN<- proms.island.filt  %>%
  dplyr::filter(logFC<0)    

genesid<-DOWN$UCSC_RefGene_Name

genesid<- strsplit(as.character(genesid),';')
genesid<-unique(unlist(genesid))

genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
genesid <- genesid[!is.na(genesid)]
#####GO enrichment analysis



eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
)


####KEGG

ego2 <-enrichKEGG(
  eg$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)



dotplot(ego2, title=" Prom islands Cpgs ",showCategory=20)

##results

ego2@result %>%
            dplyr::filter(grepl("carbon",Description))
########missMethyl

enrichment_GO <- gometh(proms.island.filt$Name,
                        all.cpg = rownames(Mval),
                        collection = "GO", 
                        array.type = "EPIC",
                        plot.bias = T,
                        prior.prob = T,
                        equiv.cpg = T,
                        anno = ann850k) 

enrichment_GO %>%
              dplyr::filter(ONTOLOGY=="BP" & FDR<.05)


enrichment_GO %>%
  dplyr::filter(ONTOLOGY=="BP" & grepl("carbon",TERM))

enrichment_GO<-enrichment_GO[enrichment_GO$FDR<.05,]



#DMRs


#####bumphunter
#library(bumphunter)
#pheno <- pData(gRatioSet.quantile)$Status
#designMatrix <- model.matrix(~ pheno)
#Run the algorithm with a large number of permutations, say B=1000:

# dmrs <- bumphunter(gRatioSet.quantile, design = designMatrix, cutoff = 0.01, B=1000, type="Beta")

# library(doParallel)
#registerDoParallel(cores = 3)
#The results of bumphunter are stored in a data frame with the rows being the different differentially methylated regions (DMRs):
#  

pheno <- as.factor(pData(gRatioSet.quantile)$Status)
pheno<-relevel(pheno,ref="WT")
designMatrix <- model.matrix(~ pheno)

myannotation <- cpg.annotate("array", 
                             Mval, 
                             what="M", 
                             arraytype = "EPIC",
                             analysis.type="differential", 
                             design=designMatrix, 
                             coef=2)
###DMR finding

dmrcoutput <- dmrcate(myannotation, 
                      lambda=1000, 
                      C=2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")


dmrs.res<-as.data.frame(results.ranges)
###filter those with less than 5 CpGs inside DMR
my.dmrs<-subset(dmrs.res,dmrs.res$no.cpgs >= 5)
my.dmrs<-subset(my.dmrs,my.dmrs$Fisher <.05)


write.csv(my.dmrs,"dmrs.proms.csv")

###select specifical

dmrs.res[grep("IDH3",dmrs.res$overlapping.genes),]

#GFRA1 and GSTM2,SEPT9
#GSTM2

groups <- c(KO="magenta", WT="forestgreen")
type<-pheno
cols <- groups[as.character(type)]

DMR.plot(ranges=results.ranges, dmr=14269, CpGs=ilogit2(myMs), phen.col=cols, genome="hg19")

