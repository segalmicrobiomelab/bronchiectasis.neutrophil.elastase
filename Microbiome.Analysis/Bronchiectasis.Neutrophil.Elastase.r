#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(vegan)
library(dplyr)
library(MetaboSignal)
library(qiime2R)
library(phyloseq)

#Set Theme
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
theme2<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.01

#Load Data
physeq<-qza_to_phyloseq("no-miss-table-dada2.qza","rooted-tree_quality.qza","taxonomy.qza", "BX_Sputum_convert.txt")

# Remove taxa with 0 abundance
OTU.Table = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
OTU.Rel.Table = transformSampleCounts(OTU.Table, normalizeSample)

# # Create phyllum and order tables (do it after normalization and out of the relative table)
Phylum.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Phylum")
Class.rel.table  = tax_glom(OTU.Rel.Table, taxrank = "Class")
Order.rel.table  = tax_glom(OTU.Rel.Table, taxrank = "Order")
Family.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Family")
Genus.rel.table  = tax_glom(OTU.Rel.Table, taxrank = "Genus")

#Create Genus Table (Raw Data)
Genus.table  = tax_glom(OTU.Table, taxrank = "Genus")
# Remove taxa with 0 abundance
#Genus.Table.2 = subset_taxa(physeq, rowSums(otu_table(Genus.table)) != 0)

#Relative Genus Table
#Genus.Rel.Table = transformSampleCounts(Genus.Table.2, normalizeSample)

#Elastase
#Genus.Rel.Table.Elastase = subset_samples(Genus.rel.table, Sample_Project=="Elastase")
Genus.Table.Elastase = subset_samples(Genus.table, Neutrophil_Elastase_Cat %in% c("Low_aNE","High_aNE"))
OTU.Table.Elastase = subset_samples(OTU.Table, Neutrophil_Elastase_Cat %in% c("Low_aNE","High_aNE"))

sample_data(Genus.Table.Elastase)$aNE <- ifelse(sample_data(Genus.Table.Elastase)$NEUTROPHIL_ELASTASE_SPUTUMugml<20, "Low_aNE", "High_aNE")
Genus.Table.Elastase.Pseudo = subset_samples(Genus.Table.Elastase, Chronic_infection_bug %in% c("P._aeruginosa"))
Genus.Table.Elastase.Pseudo = subset_samples(Genus.Table.Elastase, PSA_PCR %in% c("1"))
Genus.Table.Elastase.Hemo = subset_samples(Genus.Table.Elastase, Chronic_infection_bug %in% c("H._influenzae"))
Genus.Table.Elastase.Hemo = subset_samples(Genus.Table.Elastase, HFLU_PCR %in% c("1"))

Genus.Rel.Table.Elastase = subset_samples(Genus.rel.table, Neutrophil_Elastase_Cat %in% c("Low_aNE","High_aNE"))
sample_data(Genus.Rel.Table.Elastase)$aNE <- ifelse(sample_data(Genus.Rel.Table.Elastase)$NEUTROPHIL_ELASTASE_SPUTUMugml<20, "Low_aNE", "High_aNE")

Genus.Rel.Table.Elastase.Pseudo = subset_samples(Genus.Rel.Table.Elastase, Chronic_infection_bug %in% c("P._aeruginosa"))
Genus.Rel.Table.Elastase.Pseudo = subset_samples(Genus.Rel.Table.Elastase, PSA_PCR %in% c("1"))
Genus.Rel.Table.Elastase.Hemo = subset_samples(Genus.Rel.Table.Elastase, Chronic_infection_bug %in% c("H._influenzafdssde"))
Genus.Rel.Table.Elastase.Hemo = subset_samples(Genus.Rel.Table.Elastase, HFLU_PCR %in% c("1"))


OTU.Rel.Table.Elastase = subset_samples(OTU.Rel.Table, Neutrophil_Elastase_Cat %in% c("Low_aNE","High_aNE"))


#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================
#PCOA of Sputum Data Comparing High and Low NeutroPhil Elastase
#Create Distance Matrix
vegdist   = vegdist(t(otu_table(Genus.Rel.Table.Elastase)), method="bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(Genus.Rel.Table.Elastase), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~aNE,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="aNE",suffixes=c("",".centroid"))

pdf("16S_Neutrophil_Elastase_Cat_BRAY.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=aNE)) +
    geom_point(size=5,alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#EA3323","#296218")) +
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=aNE), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=aNE))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "NTM", "Other.Culture.Pos")), size=10) +
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=aNE), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme2
    #theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#Calculate PERMANOVA 
adonis(vegdist ~ sample_data(Genus.Rel.Table.Elastase)$aNE)
Terms added sequentially (first to last)

                                           Df SumsOfSqs MeanSqs F.Model      R2
sample_data(Genus.Rel.Table.Elastase)$aNE   1     2.994 2.99419  13.513 0.06876
Residuals                                 183    40.548 0.22158         0.93124
Total                                     184    43.542                 1.00000
                                          Pr(>F)
sample_data(Genus.Rel.Table.Elastase)$aNE  0.001 ***
Residuals
Total
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
adonis(formula = vegdist ~ sample_data(Genus.Rel.Table.Elastase)$aNE)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                                                               Df SumsOfSqs
sample_data(Genus.Rel.Table.Elastase)$Neutrophil_Elastase_Cat   1     2.994
Residuals                                                     183    40.548
Total                                                         184    43.542
                                                              MeanSqs F.Model
sample_data(Genus.Rel.Table.Elastase)$Neutrophil_Elastase_Cat 2.99419  13.513
Residuals                                                     0.22158
Total
                                                                   R2 Pr(>F)
sample_data(Genus.Rel.Table.Elastase)$Neutrophil_Elastase_Cat 0.06876  0.001
Residuals                                                     0.93124
Total                                                         1.00000

sample_data(Genus.Rel.Table.Elastase)$Neutrophil_Elastase_Cat ***
Residuals
Total
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



pdf("16S_Elastase_Neutrophil_Elastase_Cat_with_Value_BRAY.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color=Neutrophil_Elastase_Cat)) +
    #geom_point(size=5) +
    geom_point(aes(size=as.numeric(as.character(NEUTROPHIL_ELASTASE_SPUTUMugml))),alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    scale_color_manual(values=c("#EA3323","#296218")) +
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Neutrophil_Elastase_Cat), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Neutrophil_Elastase_Cat))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "NTM", "Other.Culture.Pos")), size=10) +
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Neutrophil_Elastase_Cat), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme2
    #theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#////////////////
#Alpha-Diversity
#////////////////
#Calcultes Shannon Diversity
sample_data(Genus.Rel.Table.Elastase)$Shannon = diversity(otu_table(Genus.Rel.Table.Elastase), index = "shannon", MARGIN = 2, base = exp(1))
#Convert to data frame for ggplot
shannon = as.data.frame(sample_data(Genus.Rel.Table.Elastase))
#Remove any zero values
shannon[shannon==0] <- NA

#Set Order Of Figure
shannon$or <-ifelse(shannon$Neutrophil_Elastase_Cat=="High_aNE", 1,NA)
shannon$or <-ifelse(shannon$Neutrophil_Elastase_Cat=="Low_aNE",2 ,shannon$or)

pdf("16S_Neutrophil_Elastase_Cat_SHANNON.pdf", height = 7, width = 5)
    ggplot(shannon, aes(x=reorder(Neutrophil_Elastase_Cat,+or), y=Shannon, fill=Neutrophil_Elastase_Cat)) + 
    stat_boxplot(geom ='errorbar', width=0.1)+
    geom_boxplot(outlier.shape = NA, width=0.5)+
    #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
    geom_jitter(shape=1, position=position_jitter(0.2))+
    scale_fill_manual(values=c("#EA3323","#296218")) +
    #scale_x_discrete(labels = c("NTM", "Negative"))+ 
    #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("Shannon Diversity") + 
    xlab("Neutrophil_Elastase_Cat")+
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme
dev.off()

wilcox.test(shannon$Shannon ~ shannon$Neutrophil_Elastase_Cat)
	Wilcoxon rank sum test with continuity correction

data:  shannon$Shannon by shannon$Neutrophil_Elastase_Cat
W = 2720, p-value = 8.579e-05
alternative hypothesis: true location shift is not equal to 0


#---------------
#---------------
#Run The DESEQ2
#---------------
#---------------

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Genus.Table.Elastase, ~Neutrophil_Elastase_Cat)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset
diagdds$Neutrophil_Elastase_Cat <- droplevels(diagdds$Neutrophil_Elastase_Cat)

#Relevel Data
diagdds$Neutrophil_Elastase_Cat <- relevel(diagdds$Neutrophil_Elastase_Cat, ref ="Low_aNE")

#Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
diagdds<- DESeq(diagdds)
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#=========================================================
//////////////////////TABLES/////////////////////////////
#=========================================================
#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Genus.Table.Elastase)[rownames(res), ], "matrix"))
#res.bal = cbind(as(res.bal, "data.frame"), as(tax_table(Species.bal.table)[rownames(res.bal), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$OTU)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))

#Make the full trail the First Column
res$names <- res$Gene.symbol
res$Gene.symbol <- res$row2

######get abundance data - mean relative - use otu.relative.table - THIS CODE ADDS REL ABUNDANCE AS THE DOT SIZE 
{
#decide what otu to save 
otu.to.save <-as.character(res$names)

#from relative table we should get the mean across the row of the otu table
#OW.OTU.rel.table.df <- data.frame(otu_table(Genus.Table.Elastase))
OW.OTU.rel.table.df <- data.frame(otu_table(Genus.Rel.Table.Elastase))
OW.OTU.rel.table.df.meanRA <- rowMeans(OW.OTU.rel.table.df)

#need to subset AND reorder just the otus that we have 
OW.OTU.rel.table.df.meanRA.save <- OW.OTU.rel.table.df.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- OW.OTU.rel.table.df.meanRA.save
}


#Keep only the variables you need for pathway analysis
res.1 <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.2 <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
write.table(res.1,file=  "16S_Neutrophil_Elastase_Cat.txt", sep="\t", col.names = NA, row.names = TRUE)
write.table(res.2,file="16S_Neutrophil_Elastase_Cat_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)


#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================
# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res$sig <- -log10(res$adj.P.Val)
sum(is.infinite(res$sig))

#If there infinite value set a maximum
res[is.infinite(res$sig),"sig"] <- 350

##Set the colors for your volcano plat
cols <- densCols(res$logFC, res$sig)
cols[res$pvalue ==0] <- "purple"
cols[res$logFC > 0 & res$adj.P.Val < alpha ] <- "red"
cols[res$logFC < 0 & res$adj.P.Val < alpha ] <- "green"
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

#Final Figure
pdf(file="16S_Neutrophil_Elastase_Cat_Volcano_Plot_rel_Abundance_Small_FDR_0.2.pdf", width=5, height=5)
    ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
    geom_point(color=cols, size = ifelse(res$adj.P.Val < alpha, 200 * res$abundance, 2), alpha=0.5) + #Chose Colors and size for dots
    geom_text_repel(aes(label=ifelse(res$adj.P.Val < alpha & res$Genus!="g__", as.character(res$Genus),'')),size=3,force=25, segment.colour="grey",segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
    xlab("Effect size: log2(fold-change)") + #label X Axis
    ylab("-log10(adjusted p-value)") + #label Y Axis
    #ylim(0,20)+
    theme #Set Theme
dev.off() 

res$abun = res$abundance * 200    
abundplot <-    
    ggplot(res, aes(x = logFC, y = sig,label=Gene.symbol)) +
    geom_point(aes(size=abundance*200),fill="white",pch=21)+
    #geom_point(color=cols, aes(size = ifelse(adj.P.Val < alpha, 10^-2 * abundance, 2)), alpha=0.5) +
    #guides(size="legend") #Chose Colors and size for dots
    scale_size(name   = "Rel Abundance",
    breaks = fivenum(res$abundance)*200,
    labels = fivenum(res$abundance))


legend <- get_legend(abundplot)

pdf("Volcano_Plot_Size_Legend.pdf", height = 5, width = 5)
grid.draw(legend)
dev.off()

#=========================================================
///////////////BARPLOT OF TAXA////////////////
#=========================================================

#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])
keepOTUs = rownames(res[res$padj < 0.5, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.Table.Elastase))), ntaxa(Genus.Rel.Table.Elastase)), Genus.Rel.Table.Elastase)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# agglomerate taxa
#glom <- tax_glom(phy, taxrank = 'Genus')

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comaring Status to a character vector from a factor
dat$Location <- as.character(dat$Neutrophil_Elastase_Cat)


# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
#dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")

#If you want the trail also do this
#dat$o[is.na(dat$o)] <- paste(as.character(dat$Order[is.na(dat$o)]),"(u.f.)", sep=" ")
#dat$o <- ifelse(is.na(dat$Order), dat$Order, dat$o)
#dat$o[is.na(dat$o)] <- paste(as.character(dat$Class[is.na(dat$o)]),"(u.c.)", sep=" ")
#dat$o <- ifelse(is.na(dat$Class), dat$Class, dat$o)
#dat$o <- gsub("o__", "", dat$o)

dat$Genus <- dat$o
dat <- dat[!is.na(dat$Genus),]


#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]
data <- setDT(dat)[,list(Abundance=as.numeric(mean(Abundance, na.rm=TRUE)), iqr=as.numeric(sd(Abundance, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "High_aNE", "Low_aNE"))


#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]

data$Location <-ifelse(data$Location=="Negative","pNegative",data$Location)
dat$Location <-ifelse(dat$Location=="Negative","pNegative",dat$Location)

#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#EA3323","#296218")) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c("#EA3323","#296218")) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("Genus") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESeq2_Neutrophil_Elastase.pdf", height = 10, width = 10)
print(p)
dev.off()





#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.Table.Elastase))), ntaxa(Genus.Rel.Table.Elastase)), Genus.Rel.Table.Elastase)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)


#Prune Data For BiPlot
BAL.OTU.Rel.wh1 = genefilter_sample(x10, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(x10))
BAL.OTU.Rel.table1B = prune_taxa(BAL.OTU.Rel.wh1, x10)

#Ordinate NMDS and BRAY
BAL.ords <- ordinate(BAL.OTU.Rel.table1B , "NMDS", "bray")
#Gather Data For Biplot from Ordination
pord <- plot_ordination(BAL.OTU.Rel.table1B, BAL.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape = "Neutrophil_Elastase_Cat")
#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$Strain)
#from relative table we should get the mean across the row of the otu table
BAL.OTU.rel.table.df <- data.frame(otu_table(x10))
BAL.OTU.rel.table.df.meanRA <- rowMeans(BAL.OTU.rel.table.df)
#need to subset AND reorder just the otus that we have 
BAL.OTU.rel.table.df.meanRA.save <- BAL.OTU.rel.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- BAL.OTU.rel.table.df.meanRA.save



#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Neutrophil_Elastase_Cat,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="Neutrophil_Elastase_Cat",suffixes=c("",".centroid"))
#Create Taxa name with Genus and OTU
ord$taxa <- paste(ord$Genus,ord$Strain,sep="_")
ord$taxa2 <- ifelse(ord$taxa=="Samples_NA",ord$taxa,as.character(ord$Genus))

pdf("16S_Neutrophil_Elastase_BiPlot.pdf", height = 15, width = 20)
    ggplot(ord, aes(PC1, PC2, color=Neutrophil_Elastase_Cat)) +
    geom_point(size= ifelse(ord$Neutrophil_Elastase_Cat=="Taxa", 200 * ord$abundance, 2),alpha=0.7) +    
    geom_point(data=subset(ord,Neutrophil_Elastase_Cat!="Taxa"),size=5,alpha=0.7) +
    #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
    xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
    ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=subset(centroids,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1, y=PC2, color=Neutrophil_Elastase_Cat), size=0) +
    geom_segment(data=subset(ord,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Neutrophil_Elastase_Cat))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
    geom_label_repel(data=subset(centroids,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1, y=PC2, label=c("High_aNE", "Low_aNE")), size=10) +
    geom_text_repel(data=subset(ord,Neutrophil_Elastase_Cat=="Taxa"),aes(label=taxa2), size=4)+
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme2
    #theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#Bray
BAL.OTU.Rel.dist = distance(BAL.OTU.Rel.table1B, "bray")

######NMDS 
BAL.OTU.Rel.mds <- metaMDS(BAL.OTU.Rel.dist, trace = FALSE)
BAL.OTU.Rel.mds

#Vectors
BAL.table <- sample_data(BAL.OTU.Rel.table1B)
#Create categorical columns for each variable
BAL.table$bug <- ifelse(BAL.table$Chronic_infection_bug==0,"Negative",
                        ifelse(BAL.table$Chronic_infection_bug=="A._xylosoxidans","A._xylosoxidans",
                        ifelse(BAL.table$Chronic_infection_bug %in% c("E._coli","E.coli"),"E._coli",
                        ifelse(BAL.table$Chronic_infection_bug=="E._coli_ESBL","E._coli_ESBL",
                        ifelse(BAL.table$Chronic_infection_bug=="H._influenzae","H._influenzae",
                        ifelse(BAL.table$Chronic_infection_bug=="K._pneumoniae","K._pneumoniae",
                        ifelse(BAL.table$Chronic_infection_bug=="K._pneumoniae_ESBL-","K._pneumoniae_ESBL-",
                        ifelse(BAL.table$Chronic_infection_bug=="M._abscessus","A._xylosoxidans",
                        ifelse(BAL.table$Chronic_infection_bug=="A._xylosoxidans","M._abscessus",
                        ifelse(BAL.table$Chronic_infection_bug=="MRSA","MRSA",
                        ifelse(BAL.table$Chronic_infection_bug=="MSSA","MSSA",
                        ifelse(BAL.table$Chronic_infection_bug=="P._aeruginosa","P._aeruginosa",
                        ifelse(BAL.table$Chronic_infection_bug=="P._mirabilis","P._mirabilis",
                        ifelse(BAL.table$Chronic_infection_bug=="S._maltophilia","S._maltophilia",
                        ifelse(BAL.table$Chronic_infection_bug=="S._pneumoniae","S._pneumoniae",
                        "Polymicrobial")))))))))))))))

BAL.table$ntm <-ifelse(BAL.table$Chronic_infection_bug=="M._abscessus",1,0)
BAL.table$Klebsiella <-ifelse(BAL.table$Chronic_infection_bug=="K._pneumoniae",1,0)
BAL.table$Klebsiella_ESBL <-ifelse(BAL.table$Chronic_infection_bug=="K._pneumoniae_ESBL-",1,0)
BAL.table$EColi <-ifelse(BAL.table$Chronic_infection_bug=="E._coli",1,0)
BAL.table$EColi_ESBL <-ifelse(BAL.table$Chronic_infection_bug=="E._coli_ESBL",1,0)
BAL.table$Achromobacter <-ifelse(BAL.table$Chronic_infection_bug=="A._xylosoxidans",1,0)
BAL.table$Pseudomonas <-ifelse(BAL.table$Chronic_infection_bug=="P._aeruginosa",1,0)
BAL.table$Haemophilus <-ifelse(BAL.table$Chronic_infection_bug=="H._influenzae",1,0)
BAL.table$Negative <-ifelse(BAL.table$Chronic_infection_bug=="Negative",1,0)
BAL.table$MSSA <-ifelse(BAL.table$Chronic_infection_bug=="MSSA",1,0)
BAL.table$MRSA <-ifelse(BAL.table$Chronic_infection_bug=="MRSA",1,0)


#Change the different categories for 0 and 1 for each of the different columns 

BAL.table <- subset(BAL.table, select=c("ntm","Klebsiella","Klebsiella_ESBL","EColi","EColi_ESBL","Achromobacter","Pseudomonas","Haemophilus","Negative","MSSA","MRSA"))
BAL.OTU.Rel.ef <- envfit(BAL.OTU.Rel.mds, BAL.table, permu = 999)
BAL.OTU.Rel.ef 

pdf(file = "BAL.Vectorgram.pruned.pdf", width = 20, height = 15)
plot(BAL.OTU.Rel.mds, type = "t", display = "sites")
plot(BAL.OTU.Rel.ef , p.max=0.999)
dev.off()



#Use this file with taxa in different colors so you can edit the figure in illustrator accurately
pdf("16S_BAL_NTM_v_OtherCulture_BiPlot.pdf", width = 15, height = 20)
    plot_ordination(BAL.noNeg.OTU.Rel.table1B, BAL.noNeg.ords, type = "biplot", color = "FINAL_BRONCHIECTASIS_CULTURE_Simple", title = "Lower Airway Biplot", shape = "FINAL_BRONCHIECTASIS_CULTURE_Simple", label= "Genus") + 
    geom_point(size = 2) +
    #geom_text(mapping = aes(label = Genus), size = 4, vjust = 1.5)+
    scale_shape_manual(values = c(19, 1, 6))+
    scale_color_manual(values = c("#EA3323","#296218", "#1805F0")) + 
    #stat_ellipse(type = "norm", linetype = 2) +
    #stat_ellipse(type = "t") +
    theme
dev.off()

vector <- as.data.frame(scores(BAL.OTU.Rel.ef, display="vectors"))
vector$name <-rownames(vector)
vector <-vector[vector$name!="Pseudomonas.1",]
vector$name <-ifelse(vector$name=="ntm","NTM",vector$name)
vector$name <-ifelse(vector$name=="Stenotroph","Stenotrophomonas",vector$name)
vector$name <-ifelse(vector$name=="Strep","Streptococcus",vector$name)

pdf(file = "BAL.Vectorgram.pruned.pdf", width = 20, height = 15)
    ggplot(ord, aes(PC1, PC2, color=Neutrophil_Elastase_Cat)) +
    geom_point(size= ifelse(ord$Neutrophil_Elastase_Cat=="Taxa", 200 * ord$abundance, 2),alpha=0.7) +    
    geom_point(data=subset(ord,Neutrophil_Elastase_Cat!="Taxa"),size=5,alpha=0.7) +
    #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
    xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
    ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
    #plot ellipse
    #stat_ellipse(type = "t") + 
    #plot point and lines from centroid
    geom_point(data=subset(centroids,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1, y=PC2, color=Neutrophil_Elastase_Cat), size=0) +
    geom_segment(data=subset(ord,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Neutrophil_Elastase_Cat))+ 
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids 
    #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
    geom_label_repel(data=subset(centroids,Neutrophil_Elastase_Cat!="Taxa"), aes(x=PC1, y=PC2, label=c("High_aNE", "Low_aNE")), size=10) +
    geom_text_repel(data=subset(ord,Neutrophil_Elastase_Cat=="Taxa"),aes(label=taxa2), size=4)+
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    geom_segment(data=vector,aes(x=0, y=0, xend=4*NMDS1, yend=4*NMDS2), color="blue",arrow = arrow(length = unit(0.1,"cm")))+
    geom_label(data=vector, aes(x=4.5*NMDS1, y=4.5*NMDS2,  label = name), color="blue",label.size=NA,fill="white")+
    theme2
    #theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()
