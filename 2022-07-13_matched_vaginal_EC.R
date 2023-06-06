source("LDM_fun.R") # LDM, permanovaFL, adjust.data.by.covariates
library(vegan)      # vegdist
library(GUniFrac)   # Rarefy, GUniFrac
library(prodlim)   
library("devtools")
library("RJSONIO")
library(survival)
library(parallel)
library(matrixStats)

#------------------------
# OTU table
#------------------------

otu.tab <- read.table("./DATA/Fiji vaginal and endocervical abundance  data_V4_VC.txt", sep="\t", header=TRUE, row.names=1)
dim(otu.tab) # 182 270
otu.tab[1:3,1:3]
# Acidaminococcus_intestini Acinetobacter_baumannii Actinobacillus_pleuropneumoniae
# 20V                          59                       4                               0
# 20C                           5                       8                               0
# 20F3V                         0                       4                               0
summary(rowSums(otu.tab))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1265   443627  1789983  3636499  4132926 40099996 

otu.tab.ID = rownames(otu.tab)
otu.tab.sam = sapply(strsplit(otu.tab.ID, split = "[VC]"), function(x) x[[1]])
table(otu.tab.sam)
table(table(otu.tab.sam))
# 2 
# 91 
otu.tab.VC = sapply(otu.tab.ID, function(x) substr(x, nchar(x), nchar(x)))
table(otu.tab.VC)
# C  V 
# 91 91 

#------------------------
# meta data
#------------------------

metadata = read.table("./DATA/Fiji vaginal and endocervical metadata_V5_VC.txt", sep="\t", header=TRUE, as.is=TRUE, na.strings=c("N/A", "NA", "."))
dim(metadata) # 182 20
metadata[1:3,]
# Sl.. Sample_ID Age Age.category Clinic                Ethnicity Ct.Status TimePoint Trich..Lab. Trich..Metaphlan. BV..Lab. Candida..Lab.
# 1    1       20V  21        18-24    oxf iTaukei Pacific Islander  Negative  Baseline    Negative          Negative Positive      Negative
# 2    2       20C  21        18-24    oxf iTaukei Pacific Islander  Negative  Baseline    Negative          Negative Positive      Negative
# 3    3     20F3V  21        18-24    oxf iTaukei Pacific Islander  Positive  Followup    Negative          Negative Negative      Negative
# NG..Lab. HPV..HPViewer. HPV.Types MG..VIRGO. subCST     score New.CST  X
# 1 Positive       Negative  Negative   Negative   IV-B 0.7209466   IV-D2 NA
# 2 Positive       Negative  Negative   Negative   IV-B 0.6866705    IV-B NA
# 3 Negative       Negative  Negative   Negative   IV-B 0.6313116    IV-B NA

all(metadata$Sample_ID==otu.tab.ID) # TRUE


table(metadata$Age.category)
# >24 18-24 
# 129    53 

table(metadata$New.CST)
# I-A   I-B    II III-A III-B  IV-A  IV-B IV-C1 IV-C3 IV-D0 IV-D1 IV-D2  IV-E     V 
# 3     7     1    66    24     5    56     3     3    37    16    13    18     6 

metadata$Ethnicity = gsub(" ", "", metadata$Ethnicity)
table(metadata$Ethnicity, useNA="ifany")
# Indo-Fijian iTaukei Pacific Islander                    Other   Other Pacific Islander                     <NA> 
#     65                      148                       17                       26                        2 

metadata$BV..Lab. = gsub(" ", "", metadata$BV..Lab.)
table(metadata$BV..Lab., useNA="ifany")
# Negative Positive     <NA> 
#     174       39       45 

metadata$NG..Lab. = gsub(" ", "", metadata$NG..Lab.)
table(metadata$NG..Lab., useNA="ifany")
# Negative Positive     <NA> 
#     242       10        6 

table(metadata$MG..VIRGO., useNA="ifany")
# Negative Positive 
# 249        9 

# HPV

metadata$Candida..Lab. = gsub(" ", "", metadata$Candida..Lab.)
table(metadata$Candida..Lab., useNA="ifany")
# Negative Positive     <NA> 
#     194       13       51 

#------------------------
# Variable of interest
# Confounders: can be multiple confounders
#------------------------


dim(otu.tab)
dim(metadata)
w <- which(metadata$HPV.Types %in% c("Highrisk"))  # "Highrisk", "Lowrisk", "Negative"    #  Ct.Status: c("Negative", "Positive"))
otu.tab <- otu.tab[w,]
metadata <- metadata[w,]
otu.tab.sam <- otu.tab.sam[w]
otu.tab.ID <- otu.tab.ID[w]
otu.tab.VC <- otu.tab.VC[w]
dim(otu.tab)
dim(metadata)
length(otu.tab.VC)



var = otu.tab.VC ###### change to any "var" that is of interest: Age.category
var.name = "V_or_C"   ############################### change to the name of "var"
var.type = 1   # 1: categorical, 0: continuous， 2: survival
if (var.type==0) var = as.numeric(var)
if (var.type==1) var = as.factor(var)

dim(metadata) # 182  20
dim(otu.tab)  # 182 270

if (var.type==1) {
    var.level = sort(unique(var))
    var.nlevel = length(var.level)
    table(var)
}


#------------------------
#------------------------
# Alpha diversity
#------------------------
#------------------------

w_V = seq(1, length(otu.tab.sam), by=2)
w_C = seq(2, length(otu.tab.sam), by=2)
all(otu.tab.sam[w_V]==otu.tab.sam[w_C]) # TRUE
all(otu.tab.VC[w_V]=="V") # TRUE
all(otu.tab.VC[w_C]=="C") # TRUE

log.chao1=log(estimateR(otu.tab)[2,])
shannon=diversity(otu.tab)
evenness=shannon/log(specnumber(otu.tab))        
# evenness is a scaled version of shannon, scaled by the number of species present, 
# so that all variation relates to distribution within taxa not number of taxa

# boxplots and p-values, adjusting for confounders

pdf(paste("plot_alpha_", var.name, "_HRHPV.pdf", sep=""), height=9, width=9)
par(mfcol=c(2,2), pty="s")

pvalue = signif(wilcox.test(log.chao1[w_V], log.chao1[w_C], paired=TRUE)$p.value, 2) # wilcox.test for two-level var; krustal.test for multi-level var
boxplot(log.chao1 ~ var, main="log Chao1", ylab="log Chao1", xlab=var.name)
legend(x="bottomright", legend=paste("p =", pvalue), bty="n",cex = 1.0)
boxplot(log.chao1[w_V]-log.chao1[w_C], main="log Chao1", ylab="Diff of vaginal and endocervical samples", xlab=var.name)
legend(x="bottomright", legend=paste("p =", pvalue), bty="n",cex = 1.0)
abline(h=0)

pvalue = signif(wilcox.test(shannon[w_V], shannon[w_C], paired=TRUE)$p.value, 2)
boxplot(shannon ~ var, main="Shannon", ylab="Shannon", xlab=var.name) 
legend(x="bottomright", legend=paste("p =", pvalue), bty="n",cex = 1.0)
boxplot(shannon[w_V]-shannon[w_C], main="Shannon", ylab="Diff of vaginal and endocervical samples", xlab=var.name)
legend(x="bottomright", legend=paste("p =", pvalue), bty="n",cex = 1.0)
abline(h=0)

dev.off()



#------------------------
#------------------------
# Beta diversity
#------------------------
#------------------------


res.jaccard <- jaccard.mean( otu.tab ) 
distJ.mean.sq <- res.jaccard$jac.mean.sq.o2

(p.permanova.BC <- signif(permanovaFL(otu.tab  ~ var, data = metadata, dist.method="bray", seed=67817,
                                      cluster.id=otu.tab.sam, perm.between.type="none", perm.within.type="free")$p.permanova, 2)) # Age+Sex+Collect_Site
(p.permanova.J <- signif(permanovaFL(distJ.mean.sq  ~ var, data=metadata, seed=67817, square.dist=FALSE,
                                     cluster.id=otu.tab.sam, perm.between.type="none", perm.within.type="free")$p.permanova, 2))  # parameters for requesting PERMANOVA-D2-A2
# for Age, 0.1 and 0.032

# Adjusting for confounders)
######################################### use "~ 1" in formula if there is no confounder
distBC <- adjust.data.by.covariates(formula= ~ 1, data=metadata, otu.table=otu.tab, dist.method="bray")$adj.dist  # Age+Sex+Collect_Site
distJ <- adjust.data.by.covariates(formula= ~ 1, data=metadata, dist=distJ.mean.sq, square.dist=FALSE)$adj.dist
##########################################

distBC.eigen <- eigen(distBC, symmetric=TRUE)
distJ.eigen <- eigen(distJ, symmetric=TRUE)
PC.percent.BC = signif(distBC.eigen$values[1:2]/sum(distBC.eigen$values),3)*100
PC.percent.J = signif(distJ.eigen$values[1:2]/sum(distJ.eigen$values),3)*100

# PCoA plot

col.palette = c("red","blue","black","gray","green4","orange","purple")
pch.palette = rep(19, 7)
col = as.character(factor(var, labels=col.palette[1:var.nlevel])) # assign colors to different levels of var
pch = as.numeric(as.character(factor(var, labels=pch.palette[1:var.nlevel]))) # assign symbols to different levels of var
level.index = c(1:var.nlevel)
w = which(var %in% var.level[level.index])

pdf(paste("plot_beta_", var.name, ".pdf", sep=""), height=9, width=18)

par(mfrow=c(1,2),pty="s")
plot(distBC.eigen$vectors[,1],distBC.eigen$vectors[,2], main="Bray-Curtis", 
     xlab=paste("PC1 (", PC.percent.BC[1], "%)", sep=""), ylab=paste("PC2 (", PC.percent.BC[2], "%)", sep=""), col=col[w], pch=pch[w])
text(distBC.eigen$vectors[,1],distBC.eigen$vectors[,2]+0.01, labels=otu.tab.sam, cex= 0.7)
ordiellipse(distBC.eigen$vectors, var, conf=0.9, col=col.palette[level.index])
legend(x="topleft", legend=var.level[level.index], col=col.palette[level.index], pch=pch.palette[level.index])
####################################3 remove "| confounder" in permanovaFL(.) if there is no confounder
legend(x="bottomright", legend=paste("p =", p.permanova.BC), bty="n",cex=1.0)

plot(distJ.eigen$vectors[,1],distJ.eigen$vectors[,2], main="Jaccard", 
     xlab=paste("PC1 (", PC.percent.J[1], "%)", sep=""), ylab=paste("PC2 (", PC.percent.J[2], "%)", sep=""), col=col[w], pch=pch[w])
text(distJ.eigen$vectors[,1],distJ.eigen$vectors[,2]+0.005, labels=otu.tab.sam, cex= 0.7)
ordiellipse(distJ.eigen$vectors, var, conf=0.9, col=col.palette[level.index])
##################################### remove "| confounder" in permanovaFL(.) if there is no confounder
legend(x="bottomright", legend=paste("p =", p.permanova.J), bty="n",cex=1.0)

dev.off()



#------------------------
#------------------------
# LDM (relative abundance)
#------------------------
#------------------------

w = which(colSums(otu.tab>0)<5)
otu.tab = otu.tab[,-w]
dim(otu.tab) # 182 244

fdr.nominal = 0.1

res.ldm <- ldm(otu.tab ~ var, data = metadata, dist.method="bray", n.rej.stop=100, seed=123,
               cluster.id=otu.tab.sam, perm.between.type="none", perm.within.type="free")
res.ldm$p.global.omni 

res.ldm$detected.otu.omni 

### order the detected OTUs by their p-values
w1 = match(res.ldm$detected.otu.omni[[1]], names(res.ldm$q.otu.omni))
o = w1[order(res.ldm$p.otu.omni[w1])]
# o = order(res.ldm$p.otu.omni)[1:10]
summary.tab = data.frame(raw.pvalue=signif(res.ldm$p.otu.omni[o],3),
                         adj.pvalue=signif(res.ldm$q.otu.omni[o],3),
                         mean.freq=signif(res.ldm$mean.freq[o],3),
                         #direction=t(ifelse(res.ldm$beta>0, "+", "-"))[o],
                         otu.name=names(res.ldm$q.otu.omni)[o], # can use taxonomy assignment
                         row.names=NULL)
#colnames(summary.tab)[4] = paste("direction.", rownames(res.ldm$beta[1,]), sep="")
summary.tab

adj.data.1 <- adjust.data.by.covariates(formula= ~ 1, # "~ 1" if no cov # Age.category+Ethnicity
                                        data=metadata, otu.table=otu.tab,
                                        center.otu.table=FALSE) # suppress centering
# to retain the original range
n.detected.otu.omni <- length(o)
pdf(paste("boxplot_", var.name, "pos_new.pdf", sep=""), height=15, width=15) 
par(mfrow=c(6,6), pty="s", mar=c(3,2.7,4,2))
for (i in 1:n.detected.otu.omni) {
    boxplot(adj.data.1$y.freq[,o[i]]~var,
            main=names(res.ldm$q.otu.omni)[o[i]], ylab="Relative Abundance", xlab="")
    legend(x="topright", legend=paste("p =", signif(res.ldm$p.otu.omni[o[i]],3)), bty="n", cex=1.3)
}
dev.off()


#-----------------------------------------------
# LDM-A (presence-absence)
#-----------------------------------------------

res.ldmA <- ldm(formula=otu.tab ~ var, 
                data=metadata, seed=123, n.rej.stop=1000, 
                n.rarefy="all",
                cluster.id=otu.tab.sam, perm.between.type="none", perm.within.type="free") # parameter for requesting LDM-A
res.ldmA$p.global.pa

res.ldmA$detected.otu.pa  


w1 = match(res.ldmA$detected.otu.pa[[1]], names(res.ldmA$q.otu.pa))
o = w1[order(res.ldmA$p.otu.pa[w1])]
# o = order(res.ldmA$p.otu.pa)[1:10]
summary.ldmA.tab = data.frame(raw.pvalue=signif(res.ldmA$p.otu.pa[o],3),
                              adj.pvalue=signif(res.ldmA$q.otu.pa[o],3),
                              prob.presence=signif(colMeans(res.ldmA$phi)[o],3),
                              #direction=t(ifelse(res.ldmA$beta>0, "+", "-"))[o],
                              otu.name=names(res.ldmA$q.otu.pa)[o],
                              row.names=NULL)
#colnames(summary.ldmA.tab)[4] = paste("direction.", rownames(res.ldmA$beta), sep="")
summary.ldmA.tab

adj.phi = t(t(lm(res.ldmA$phi~1, data=metadata)$resid) #"~ 0" if no cov # Age.category+Ethnicity
            + colMeans(res.ldmA$phi))
n.detected.otu.pa <- length(o)
pdf(paste("boxplot_pa_", var.name, "pos_new.pdf", sep=""), height=15, width=15) 
par(mfrow=c(6,6), pty="s", mar=c(3,2.7,4,2))
for (i in c(1:n.detected.otu.pa)) {
    boxplot(adj.phi[,o[i]] ~ var,
            main=names(res.ldmA$q.otu.pa)[o[i]], ylab="Chance of Presence", xlab=var.name)
    legend(x="topright", legend=paste("p =", signif(res.ldmA$p.otu.pa[o[i]],3)), bty="n",cex=1.3)
}
dev.off()




