setwd("~/Dropbox/tsa_2014/rnaseq/deseq2_host/") # change to your working directory

library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("pheatmap") # for pretty heatmaps
library("tidyverse") # for wrangling, plotting
library("VennDiagram") # for making Venn diagram
library("vegan") # for distance matrix calculations
library("ape") # for PCoA

# Load data ---------------------------------------------------------------
countdata <- read.table("../allcounts_Mcav_CladeC.txt",header=TRUE,row.names = 1) 
head(countdata)
tail(countdata)
nrow(countdata)
# 36106 genes mapped

# simplify sample names
names(countdata)
names(countdata) <- gsub("MC","",names(countdata)) # replace the "MC" with nothing
names(countdata) <- gsub(".fq.sam.counts","",names(countdata)) # replace the ".fq.sam.counts" with nothing
names(countdata)

# make conditions table in excel and load it here
# or make it here if it's a simple design
conds <- read.csv("../conds.csv")
head(conds)
str(conds)

# VERY IMPORTANT!
# The sample names in the "conds" table must match the sample names in the counts matrix (column names)
# The order must be the same. Be sure to check it!

# First, check that the names in the countdata file exist as sample names in conds
table(names(countdata) %in% conds$sam)
# Yes! All of the names in count data are also in the conds table

# Next, check that the order matches
table(names(countdata) == conds$sam)
# No! The order is totally off. Let's fix it.

# Reorder the samples conds table to match the order in the counts matrix
conds <- conds[match(names(countdata), conds$sam),]
head(conds)

# check that it worked
table(names(countdata) == conds$sam)
# yes!

# now that everything matches, lets make our lives easier by giving our samples more informative names
head(conds)
conds$name <- paste(conds$sam,conds$bank,conds$pheno,sep="_")
head(conds)

names(countdata) <- conds$name
head(countdata)

# SEPARATE HOST AND SYMBIONT ------------------------
length(grep("Mcavernosa",row.names(countdata)))
# 21748 host genes
length(grep("isogroup",row.names(countdata)))
# 14358 symbiont genes

hostCounts <- countdata[grep("Mcavernosa",row.names(countdata)),]
head(hostCounts)
tail(hostCounts)
nrow(hostCounts)

symCounts <- countdata[grep("isogroup",row.names(countdata)),]
head(symCounts)
tail(symCounts)
nrow(symCounts)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS -------------------------------------------
# Change here depending on if you are analyzing HOST or SYMBIONT
countdata <- hostCounts
# countdata <- symCounts

#----------Total counts?
totalCounts <- colSums(countdata)
min(totalCounts) # 72909
mean(totalCounts) # 327146.5
max(totalCounts)  # 541846

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = conds,
  design = ~ pheno + bank)

save(conds, countdata, dds, file="dds_Start.Rdata")
load("dds_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(countdata,1,mean)
table(means>3)
# FALSE  TRUE 
# 13594  8154   

means3 <- names(means[means>3])
head(means3)
length(means3)
#8154

countFilt <- countdata[row.names(countdata) %in% means3,]
head(countFilt)

totalCountsFilt <- colSums(countFilt)
totalCountsFilt

min(totalCountsFilt) #71638
max(totalCountsFilt) #528713
mean(totalCountsFilt) #319547.4

# check sample order one more time... just in case!
table(names(countFilt) == as.vector(conds$name))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = countFilt,
  colData = conds,
  design = ~ pheno + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("pheno", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out2 <- c("11A_e_s", "17A_w_s", "3A_w_r", "7A_e_s")
barelyout2 <- c("18B_w_s")
out1 <- c("1B_w_s")

# remove out3 and out2
dim(countFilt)
countdata.out <- countFilt %>% select(-one_of(c(out2)))
dim(countdata.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!name %in% c(out2))
dim(conds.out)

# START HERE TO ADD GROWTH DATA TO CONDS ------------
# update "conds.out" with calcification rate and tissue growth normalized to surface area
# add those two new columns to conds.out (merge function)
# save(countdata.out, conds.out, file="counts_conds_filtered.RData")
load("counts_conds_filtered.RData")
load("../../tissueGrowth/growth.RData")

head(conds.out)
head(polyps)

# correct genets in polyps 
conds.out[!conds.out$oldgeno==conds.out$newgeno,]
# 7 and 29 not in data set
polyps$genet <- gsub("^21", "2117", polyps$genet)
polyps$genet <- gsub("^17", "2117", polyps$genet)
polyps$genet

# summary stats
mean(polyps$calcification) #0.000369189 g / cm2 / day

# annual rates in g
mean(polyps$calcification*365) # 0.134754 g / cm / year
sd(polyps$calcification*365) #  0.0609493 g / cm / year
hist(polyps$calcification*365) #  0.0609493 g / cm / year

# daily rates
mean(polyps$calcification*1000) # 0.369189 g / cm / day

# polyp over the year (total number, not normalized to surface area)
mean(polyps$tissuegrowth) #10.94828
sd(polyps$tissuegrowth) #6.697531
range(polyps$tissuegrowth) # 0 - 30

mean(polyps$normtissue) #2.384791
sd(polyps$normtissue) #1.341659
range(polyps$normtissue) # 0 - 6.683003

# first, make a per-genet average of calcification rate (g/cm2/day) 
# and normtissue (polyp/cm2/year)
avgGrowth <- polyps %>% group_by(genet) %>%
  rename(newgeno = genet) %>%
  summarize(meanCalc = round(mean(calcification*365*1000),0),
            meanPolyp = round(mean(normtissue),0))
hist(avgGrowth$meanCalc)
hist(avgGrowth$meanPolyp)

test<- merge(conds.out, avgGrowth, by = "newgeno", all.x=T)
head(test)

conds.out <- test
head(conds.out)

# subset for no NA
conds.growth <- conds.out %>% filter(!is.na(meanCalc))
head(conds.growth)

head(countdata.out)
counts.growth <- countdata.out[,names(countdata.out) %in% conds.growth$name]
head(counts.growth)
dim(counts.growth)

table(names(counts.growth) %in% conds.growth$name)
table(names(counts.growth) == conds.growth$name)

test <- conds.growth[match(names(counts.growth), conds.growth$name),]
head(test)
table(names(counts.growth) == test$name)
conds.growth <- test

# Reconstruct data object (filtered and outliers removed) ------------------------------
# START HERE TO MAKE NEW MODELS THAT INCLUDE DIFFERENT CONDITIONS TO TEST -------
# change the design

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = countdata.out,
  colData = conds.out,
  design =  ~ pheno + bank)

ddsCalcification <- DESeqDataSetFromMatrix(
  countData = counts.growth,
  colData = conds.growth,
  design =  ~ meanCalc)

ddsPolyp <- DESeqDataSetFromMatrix(
  countData = counts.growth,
  colData = conds.growth,
  design =  ~ meanPolyp)


# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds.pheno <- DESeq(ddsFiltOut)

deds.calc <- DESeq(ddsCalcification)

deds.polyp <- DESeq(ddsPolyp)

#---Results
resPheno <- results(deds.pheno, independentFiltering = F, contrast=c("pheno","r","s"))
resPheno
# log2 fold change (MLE): pheno r vs s 
# Wald test p-value: pheno r vs s 
# here "resistant" is the numerator and "susceptible" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed 
# in resistant compared to susceptible genets
table(resPheno$padj<0.1)

resBank <- results(deds.pheno, independentFiltering = F, contrast=c("bank", "e", "w"))
resBank
# log2 fold change (MLE): bank e vs w 
# Wald test p-value: bank e vs w 
# here "east" is the numerator and "west" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed 
# in east bank compared to west bank genets
table(resBank$padj<0.1)

resPolyp <- results(deds.polyp, independentFiltering = F)
resPolyp
# log2 fold change (MLE): meanPolyp 
# Wald test p-value: meanPolyp
table(resPolyp$padj<0.1)


resCalc <- results(deds.calc, independentFiltering = F)
resCalc
# log2 fold change (MLE): meanCalc 
# Wald test p-value: meanCalc 
table(resCalc$padj<0.05)

# make a new VSD file (this is a useful format for making heatmaps later)
vsd.pheno.bank <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)
head(assay(vsd.pheno.bank))

vsd.calc <- varianceStabilizingTransformation(ddsCalcification, blind=TRUE)
head(assay(vsd.growth))

vsd.polyp <- varianceStabilizingTransformation(ddsPolyp, blind=TRUE)
head(assay(vsd.polyp))

# Write results for making heatmaps and other downstream analyses ---------------------------------------

###--------------Get pvals
head(resPheno)
valsPheno <- cbind(resPheno$pvalue, resPheno$padj)
head(valsPheno)
colnames(valsPheno)=c("pval.pheno", "padj.pheno")

head(resBank)
valsBank <- cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.bank", "padj.bank")

head(resCalc)
valsCalc <- cbind(resCalc$pvalue, resCalc$padj)
head(valsCalc)
colnames(valsCalc)=c("pval.calc", "padj.calc")

head(resPolyp)
valsPolyp <- cbind(resPolyp$pvalue, resPolyp$padj)
head(valsPolyp)
colnames(valsPolyp)=c("pval.polyp", "padj.polyp")

#Make vsd data and pvals table
vsdpvals.pheno.bank <- as.data.frame(cbind(assay(vsd.pheno.bank),valsPheno, valsBank))
head(vsdpvals.pheno.bank)
dim(vsdpvals.pheno.bank)
# 8154   54

vsdpvals.calc <- as.data.frame(cbind(assay(vsd.calc),valsCalc))
head(vsdpvals.calc)
dim(vsdpvals.calc)
# 8154    43

vsdpvals.polyp <- as.data.frame(cbind(assay(vsd.polyp),valsPolyp))
head(vsdpvals.polyp)
dim(vsdpvals.polyp)
# 8154   43

# SAVE and LOAD (start here for downstream analyses) -------
save(avgGrowth,conds.out, countdata.out, ddsFiltOut,conds.growth, counts.growth,
     deds.calc, deds.polyp, deds.pheno,
     resPheno, resPolyp, resCalc, resBank,
     vsdpvals.calc, vsdpvals.pheno.bank, vsdpvals.polyp,
     file="ddsCoralFiltOut.Rdata")

load("ddsCoralFiltOut.Rdata")

# Write results for BLAST
sigCalc2BLAST <- row.names(resCalc[resCalc$padj<0.05,])

# Plot some interesting genes

head(vsdpvals.calc)

vsd2plot <- vsdpvals.calc %>%
  rownames_to_column("gene") %>%
  filter(gene == "Mcavernosa00597") %>%
  select(-c(gene,pval.calc,padj.calc)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("name")
vsd2plot
pheno2plot <- conds.growth %>% select(name,meanCalc)
pheno2plot

vsdpheno2plot <- merge(vsd2plot,pheno2plot,by="name")
head(vsdpheno2plot)

vsdpheno2plot %>% ggplot(aes(x = V1, y = meanCalc)) +
  geom_point() +
  geom_smooth(method=lm, col="black")+
  xlab("IHH expression")+
  ylab("Calcification (mg/cm2/year)")+
  theme_classic()

# Write results for supplementary table
sigPheno4Table <- cbind(resPheno[resPheno$padj <0.05,],"comparison"=rep("resistant vs. susceptible",1))
sigBank4Table <- resBank[resBank$padj <0.05,] # no DEGs
sigCalc4Table <- cbind(resCalc[resCalc$padj <0.05,],"comparison"=rep("calcification",70))
sigPolyp4Table <- cbind(resPolyp[resPolyp$padj <0.05,],"comparison"=rep("polyp",1))

sigAll4Table <- rbind(sigPheno4Table,sigBank4Table,sigCalc4Table,sigPolyp4Table)
write.table(sigAll4Table, file = "res_all_sig.txt", quote=F, sep="\t")

# Look at the results -------

# how many DE genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resPheno$padj < 0.05)
resPheno[resPheno$padj <0.05,]
# Mcavernosa25949 52.8801998470291 3.20094674281573 0.51907468857887 6.16664001009051 6.97563173166543e-10 5.68793011399999e-06

resPheno[resPheno$padj <0.15,]


table(resBank$padj < 0.05)
# none

table(resCalc$padj<0.05)
# FALSE  TRUE 
# 8084    70 
resCalc[resCalc$padj <0.01,]

table(resPolyp$padj<0.05)
# 1
resPolyp[resPolyp$padj <0.05,]

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds.pheno, main="Dispersion Plot")
plotDispEsts(deds.calc, main="Dispersion Plot")
plotDispEsts(deds.polyp, main="Dispersion Plot")

#MA plot
plotMA(resPheno, ylim = c(-2, 2), main="MA Plot Pheno") 
plotMA(resBank, ylim = c(-2, 2), main="MA Plot Bank")
plotMA(resCalc, ylim = c(-0.01, 0.01), main="MA Plot Growth") # the LFCs are tiny because they are PER UNIT (mg per cm2 per year)
plotMA(resPolyp, ylim = c(-2, 2), main="MA Plot Bank") # per unit (polyp per cm)

# Venn diagram ---------
sig.pheno <-  row.names(vsdpvals.pheno.bank[vsdpvals.pheno.bank$padj.pheno<0.05 & !is.na(vsdpvals.pheno.bank$padj.pheno),])
sig.bank <- row.names(vsdpvals.pheno.bank[vsdpvals.pheno.bank$padj.bank<0.05 & !is.na(vsdpvals.pheno.bank$padj.bank),])
sig.calc <- row.names(vsdpvals.calc[vsdpvals.calc$padj.calc<0.05 & !is.na(vsdpvals.calc$padj.calc),])
sig.polyp <- row.names(vsdpvals.polyp[vsdpvals.polyp$padj.polyp<0.05 & !is.na(vsdpvals.polyp$padj.polyp),])

candidates <- list("Phenotype"=sig.pheno,"Calcification" = sig.calc, "Polyp" = sig.polyp)

prettyvenn <- venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "forestgreen", "magenta"),
  alpha = 0.5,
  label.col = c("black"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("black"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

# Make a gene expression heatmaps ------------
load("ddsCoralFiltOut.Rdata")

# subset for just the expression data
expg.calc <- vsdpvals.calc[,1:41]
head(expg.calc)

# Make p-value cut-offs
sig.calc <- row.names(vsdpvals.calc[vsdpvals.calc$padj.calc<0.01 & !is.na(vsdpvals.calc$padj.calc),])

# subset expression data for significant DEGs
exp.calc <- expg.calc[row.names(expg.calc) %in% sig.calc,]
nrow(exp.calc) #70 at 0.05; 9 at 0.01

# Load annotations
gg <- read.delim("/Volumes/wrightlab/genomes/Mcav_genome/Mcav_gene2annot.tab", header=F)
head(gg)

# phenotype heatmap ---------
# how many are annotated?
table(row.names(exp.calc) %in% gg$V1)

# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.calc)){
  s=subset(gg,V1==i)
  gnames=append(gnames,paste(s$V2[1],i,sep="."))
  expg=rbind(expg,exp.calc[i,])
} 
row.names(expg)=gnames
expl=expg
means <- apply(expl,1,mean) # means of rows
explc <- expl-means # subtracting them
head(explc)
row.names(explc)

# shorten the really long names
row.names(explc)[3] <- "Intercellular signal essential for a variety of patterning events during development.Mcavernosa00597"
row.names(explc)[35] <- "Catalyzes ATP-dependent carboxylation of biotin and carboxyl transfer to pyruvate.Mcavernosa20920"  
row.names(explc)[50] <- "Histone lysine demethylase and a ribosomal histidine hydroxylase.Mcavernosa26817"
row.names(explc)[59] <- "DNA-dependent RNA polymerase.Mcavernosa30823"

# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.25)(100)

# cluster plot
pheatmap(as.matrix(explc), color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         clustering_distance_rows="correlation", cluster_cols=T)

# sort by calcification rate
sort.by.calc <- conds.growth %>% arrange(meanCalc) %>% select(name, meanCalc)
explc.sorted <- explc[c(match(sort.by.calc$name,names(explc)))]

# plot sorted
quartz()
pheatmap(as.matrix(explc.sorted), color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         clustering_distance_rows="correlation", cluster_cols=F)
quartz()
barplot(sort.by.calc$meanCalc,ylab = "mg per cm2 per year", col="black")

# BLAST sig genes to check annots
resCalc[resCalc$padj<0.01,]
#                 baseMean       log2FoldChange                lfcSE              stat               pvalue                 padj
# Mcavernosa04698 26.2318035356699   0.0047656661692057 0.000910343287160094  5.23502093817011 1.64966161805633e-07 0.000336283520840784
# Mcavernosa07002 41.5559521676303 -0.00516341847466237 0.000897607837849515 -5.75242133249735  8.7974209355306e-09 7.17341703083165e-05
# Mcavernosa07810 50.8939431450941  0.00411750463481326 0.000831433707970492  4.95229456701242 7.33434678581644e-07 0.000996737728192454
# Mcavernosa23221 1003.95986912731 -0.00929382859471799  0.00193653669250793 -4.79920087787333 1.59299958474164e-06  0.00185561694485476
# Mcavernosa26244 24.8278571870337 -0.00529674378092343  0.00115776657868846 -4.57496690474834 4.76295272029861e-06  0.00431523516459054
# Mcavernosa28697 32.8905579744309 -0.00711560651720192  0.00140125706629949 -5.07801650984224 3.81395560651743e-07 0.000621979880310863
# Mcavernosa29491 70.0411101485753 -0.00536067876522525  0.00101898085228895  -5.2608238449069  1.4341137796543e-07 0.000336283520840784
# Mcavernosa31935 9.95455790722664  -0.0165214137038463  0.00303523177962117 -5.44321320525589 5.23279526855823e-08 0.000213341063099119
# Mcavernosa34054 22.9499130278549  -0.0206834338270276  0.00448896785193616 -4.60761460301093 4.07314717616008e-06  0.00415155525930116

# find gene of interest
resCalc[row.names(resCalc) == "Mcavernosa12226",]

# not enough DEGs to make pheatmap for others, but what are the genes?
resPheno[resPheno$padj<0.05,]
# log2 fold change (MLE): pheno r vs s 
#                   baseMean          log2FoldChange    lfcSE               stat               pvalue                 padj
# Mcavernosa25949   52.8801998470291  3.20094674281573  0.51907468857887    6.16664001009051 6.97563173166543e-10 5.68793011399999e-06
# NCBI BLAST result
# PREDICTED: Orbicella faveolata lymphocyte antigen 6H-like (LOC110039806), transcript variant X2, mRNA
# Score	          Expect	Identities	  Gaps	      Strand
# 514 bits(278)	  2e-141	390/446(87%)	0/446(0%)	  Plus/Plus

resBank[resBank$padj<0.05,]
# no DEGs

resPolyp[resPolyp$padj<0.05,]
# Mcavernosa00313 22.2027188216228 -0.368021211962868 0.0762307520036241 -4.82772637406715 1.38100664898065e-06 0.0112607282157883
# NCBI BLAST result
# PREDICTED: Orbicella faveolata uncharacterized LOC110062166 (LOC110062166), partial mRNA
#Score	          Expect	  Identities	Gaps	Strand
#1306 bits(707)	  0.0	      939/1053(89%)	8/1053(0%)	Plus/Plus
# Uniprot BLAST
# ANK_REP_REGION domain-containing protein Pdam
# E 1.8e-155 score 1335 ident 52.9%

# PCoA (for bank + pheno)---------
load("ddsCoralFiltOut.Rdata")

# get variance stabilized expression data
head(vsdpvals.pheno.bank)
names(vsdpvals.pheno.bank)
# subset for just the expression data
exp <- vsdpvals.pheno.bank[c(1:50)]
head(exp)

# make sure condition data match expression data
table(conds.out$name == names(exp))
# NO

# But it should be in the exp data right? 
table(conds.out$name %in% names(exp))
# Yes!

# Reorder the conds. out to match the expression data
conds.out <- conds.out[match(names(exp), conds.out$name),]
head(conds.out)

# check that it worked
table(conds.out$name == names(exp))
# yes!

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis(t(exp)~pheno+bank,
                    data=conds.out,method="manhattan")
adonisRes
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# pheno      1  10564329 10564329  1.3147 0.02635  0.040 * 
# bank       1  12751723 12751723  1.5870 0.03180  0.004 **

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA (bank+pheno)------
margin <- .25

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA for mid by site type
pdf(file = "fig_pcoa_pheno_bank_host.pdf", width = 5, height = 5)
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = "Coral Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5)
# plot "spiders" connecting the samples from each bank
ordiellipse(scores,conds.out$pheno,label=F,cex=1)

points(scores[conds.out$pheno=="r" & conds.out$bank=="w",xaxis],
       scores[conds.out$pheno=="r" & conds.out$bank=="w",yaxis],
       col="black", pch=15, cex=1.5) +
points(scores[conds.out$pheno=="r" & conds.out$bank=="e",xaxis],
         scores[conds.out$pheno=="r" & conds.out$bank=="e",yaxis],
         col="black", pch=16, cex=1.5) +
points(scores[conds.out$pheno=="s" & conds.out$bank=="w",xaxis],
         scores[conds.out$pheno=="s" & conds.out$bank=="w",yaxis],
         col="grey60", pch=15, cex=1.5) +
points(scores[conds.out$pheno=="s" & conds.out$bank=="e",xaxis],
         scores[conds.out$pheno=="s" & conds.out$bank=="e",yaxis],
         col="grey60", pch=16, cex=1.5) 

# legend of sites 
legend("topright", 
       c("Resistant", "Susceptible"),
       pch=c(8,8), 
       col=c("black","grey60"), cex=1.5, bty = "n")
legend("bottomright", 
       c("West Bank", "East Bank"),
       pch=c(0,1), 
       col=c("black","black"), cex=1.5, bty = "n")

#insert p value 
legend("topleft",
       paste("Phenotype p = ",adonisRes$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  

legend("bottomleft", 
       paste("Bank p = ",adonisRes$aov.tab$`Pr(>F)`[2], sep=" "), 
       cex=1.5, bty='n')
dev.off()

# PCoA (growth)-----------------
# get variance stabilized expression data for growth
head(vsdpvals.calc)
names(vsdpvals.calc)
# subset for just the expression data
exp.g <- vsdpvals.calc[c(1:41)]
head(exp.g)
names(exp.g)
# make sure condition data match expression data
table(conds.growth$name == names(exp.g))

# compute dissimilarity indices
dd.veg.g <- vegdist(t(exp.g), "manhattan")
div.dd.veg.g <- dd.veg.g/1000
head(div.dd.veg.g)

# perform PERMANOVA  
set.seed(1)
adonisRes.g <- adonis(t(exp.g)~meanCalc,
                    data=conds.growth,method="manhattan")
adonisRes.g
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# meanCalc   1  11100334 11100334  1.5366 0.03791  0.011 *
# Residuals 39 281740503  7224115         0.96209         
# Total     40 292840837                  1.00000          

# compute principal coordinate decomposition
dd.pcoa.g <- pcoa(div.dd.veg.g)
head(dd.pcoa.g)
scores.g <- dd.pcoa.g$vectors

# plotting PCoA (growth)----------------
pdf(file = "fig_pcoa_calc_host.pdf", width = 5, height = 5)
plot(scores.g[,xaxis], scores.g[,yaxis],type="n", 
     main = "Coral Gene Expression",
     xlim=c(min(scores.g[,xaxis])-margin,max(scores.g[,xaxis])+margin),
     ylim=c(min(scores.g[,2])-margin,max(scores.g[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa.g$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa.g$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5)

# plot the points (size == meanCalc)
points(scores.g[,xaxis],
       scores.g[,yaxis],
       col="black", pch=15, cex=c(conds.growth$meanCalc/100))

#insert p value 
legend("topright",
       paste("Calc p = ",adonisRes.g$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  
dev.off()

# Write results for GO/KOG analysis by negative log pvalue-----------------------------------
load("ddsCoralFiltOut.Rdata")

# write GO results by phenotype, log-transformed signed pvalue
head(resPheno)
logs <- data.frame(cbind("gene"=row.names(resPheno),
                         "logP"=round(-log(resPheno$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resPheno$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#3913 4241  
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOPheno_logP.csv",sep=",")

# write GO results by phenotype, log-fold change
head(resPheno)
lfc <- data.frame(cbind("gene"=row.names(resPheno),
                         "LFC"=round(resPheno$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOPheno_LFC.csv",sep=",")

# write GO results by bank, log-transformed signed pvalue
head(resBank)
logs <- data.frame(cbind("gene"=row.names(resBank),
                         "logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#4189 3965  
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOBank_logP.csv",sep=",")

# write GO results by bank, log-fold change
head(resBank)
lfc <- data.frame(cbind("gene"=row.names(resBank),
                        "LFC"=round(resBank$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOBank_LFC.csv",sep=",")

# write GO results by calc rate, log-transformed signed pvalue
head(resCalc)
logs <- data.frame(cbind("gene"=row.names(resCalc),
                         "logP"=round(-log(resCalc$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resCalc$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#4082 4072
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOCalc_logP.csv",sep=",")

# write GO results by calc rate, log-fold change
head(resCalc)
lfc <- data.frame(cbind("gene"=row.names(resCalc),
                        "LFC"=round(resCalc$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOCalc_logP.csv_LFC.csv",sep=",")

# write GO results by polyp growth, log-transformed signed pvalue
head(resPolyp)
logs <- data.frame(cbind("gene"=row.names(resPolyp),
                         "logP"=round(-log(resPolyp$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resPolyp$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#3873 4281 
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOPolyp_logP.csv",sep=",")

# write GO results by polyp growth, log-fold change
head(resPolyp)
lfc <- data.frame(cbind("gene"=row.names(resPolyp),
                        "LFC"=round(resPolyp$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOPolyp_LFC.csv",sep=",")