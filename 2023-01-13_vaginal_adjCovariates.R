source("LDM_fun.R") # LDM, permanovaFL, adjust.data.by.covariates, v5.0
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

otu.tab <- read.table("./DATA/Fiji vaginal and endocervical abundance  data_V3_V.txt", sep="\t", header=TRUE, row.names=1)
dim(otu.tab) # 229 270
otu.tab[1:3,1:3]
# Acidaminococcus_intestini Acinetobacter_baumannii Actinobacillus_pleuropneumoniae
# 16V                         175                       0                               0
# 16F3V                         0                       2                               0
# 20V                          59                       4                               0 
summary(rowSums(otu.tab))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1265   973495  2608518  4884610  5453379 47820195 

# pdf("plot_libsize.pdf", height=4.5, width=4.5)
# par(mfcol=c(1,1), pty="s")
# plot(sort(rowSums(otu.tab)))
# dev.off()

otu.tab.ID = rownames(otu.tab)
otu.tab.sam = sapply(strsplit(otu.tab.ID, split = "[A-Z]+"), function(x) x[[1]])
length(table(otu.tab.sam)) # 111
table(table(otu.tab.sam))
# 1  2  3 
# 10 84 17 



#------------------------
# meta data
#------------------------

metadata = read.table("./DATA/Fiji vaginal and endocervical metadata_V5_V.txt", sep="\t", header=TRUE, as.is=TRUE, na.strings=c("N/A", "NA", "."))
dim(metadata) # 229  19
metadata[1:3,]
# Sl.. Sample_ID Age Age.category Clinic                Ethnicity Ct.Status TimePoint Trich..Lab. Trich..Metaphlan. BV..Lab. Candida..Lab. NG..Lab.
# 1    1       16V  29          >24    oxf iTaukei Pacific Islander  Positive  Baseline    Negative          Negative Positive      Negative Negative
# 2    2     16F3V  29          >24    oxf iTaukei Pacific Islander  Negative  Followup    Negative          Negative Negative      Positive Negative
# 3    3       20V  21        18-24    oxf iTaukei Pacific Islander  Negative  Baseline    Negative          Negative Positive      Negative Positive
# HPV..HPViewer. HPV.types MG..VIRGO. subCST     score New.CST
# 1       Positive   Lowrisk   Negative  III-B 0.6172243   III-B
# 2       Positive  Highrisk   Negative  III-A 0.9907228   III-A
# 3       Negative  Negative   Negative   IV-B 0.7209466   IV-D2

all(metadata$Sample_ID==otu.tab.ID) # TRUE

metadata$otu.tab.sam = otu.tab.sam

table(metadata$Age.category, useNA="ifany")
# >24 18-24 
# 158    71 

metadata$Ethnicity = gsub(" ", "", metadata$Ethnicity)
table(metadata$Ethnicity, useNA="ifany")
# Indo-Fijian iTaukei Pacific Islander                    Other   Other Pacific Islander                     <NA> 
#     64                    126                     16                     21                      2 

table(metadata$New.CST, useNA="ifany")
# I-A   I-B III-A III-B  IV-A  IV-B IV-C1 IV-C3 IV-D0 IV-D1 IV-D2  IV-E     V 
# 3     7    59    21     5    47     2     3    33    15    13    15     6 
metadata$New.CST.pooled = metadata$New.CST
metadata$New.CST.pooled[metadata$New.CST.pooled %in% c("I-A", "I-B")] = "I"
metadata$New.CST.pooled[metadata$New.CST.pooled %in% c("IV-C1", "IV-C3")] = "IV-C"
table(metadata$New.CST.pooled, useNA="ifany")
# I III-A III-B  IV-A  IV-B  IV-C IV-D0 IV-D1 IV-D2  IV-E     V 
# 10    59    21     5    47     5    33    15    13    15     6 

metadata$BV..Lab. = gsub(" ", "", metadata$BV..Lab.)
table(metadata$BV..Lab., useNA="ifany")
# Negative Positive     <NA> 
#    174       39       16 

metadata$NG..Lab. = gsub(" ", "", metadata$NG..Lab.)
table(metadata$NG..Lab., useNA="ifany")
# Negative Positive     <NA> 
#     217       10        2 

table(metadata$MG..VIRGO., useNA="ifany")
# Negative Positive 
# 222        7 

table(metadata$HPV.types, useNA="ifany")
# Highrisk  Lowrisk Negative 
# 86       55       88 

metadata$Candida..Lab. = gsub(" ", "", metadata$Candida..Lab.)
table(metadata$Candida..Lab., useNA="ifany")
# Negative Positive     <NA> 
#     194       13       22


#------------------------
# Variable of interest
# Confounders: can be multiple confounders
#------------------------


### Sankhya: change the variable in the following two lines to any variable after "###"
metadata$var = metadata$Ethnicity ### Age.category, Ethnicity, New.CST.pooled, BV..Lab., NG..Lab., MG..VIRGO., HPV.types, Candida..Lab., Ct.Status
var.name = "Ethnicity"   
if (var.name == "Ethnicity") {
    metadata$var[metadata$var=="iTaukeiPacificIslander"] = "iTaukei"
    metadata$var[metadata$var=="OtherPacificIslander"] = "Other Pacific Islanders"
}

var.type = 1   # 1: categorical, 0: continuous， 2: survival
if (var.type==0) metadata$var = as.numeric(metadata$var)
if (var.type==1) metadata$var = as.factor(metadata$var)

# Missing values
w = which(is.na(metadata$var))
if (length(w)>0) {
    metadata = metadata[-w,]
    otu.tab = otu.tab[-w,]
}

if (var.type==2) {
    metadata$var = as.numeric(metadata$var)
    metadata$var.status = metadata$Overal_Survival
}
dim(metadata) # 229 20
dim(otu.tab)  # 229 270

if (var.type==1) {
    var.level = sort(unique(metadata$var))
    var.nlevel = length(var.level)
    table(metadata$var)
}



### create the strata variable

unique.ind = unique(metadata$otu.tab.sam)
metadata$strata = rep(NA, length(metadata$otu.tab.sam))
for (i in 1:length(unique.ind)) {
    w = which(metadata$otu.tab.sam==unique.ind[i])
    metadata$strata[w] = length(w)
}
table(metadata$strata)
# strata
# 1   2   3 
# 10 168  51 


table(metadata$Ct.Status, useNA="ifany")
# Negative Positive 
# 146       83
table(metadata$Ct.Status, metadata$var)
metadata$Ct.Status1 = matrix(1*(metadata$Ct.Status=="Positive"), ncol=1)
ldm(Ct.Status1 ~ var,   
    data=metadata, seed=67817, scale.otu.table=FALSE, freq.scale.only=TRUE,
    cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free", n.perm.max=1000)$p.global.freq


#------------------------
#------------------------
# Alpha diversity
#------------------------
#------------------------

log.chao1=log(estimateR(otu.tab)[2,])
shannon=diversity(otu.tab)
evenness=shannon/log(specnumber(otu.tab))        
# evenness is a scaled version of shannon, scaled by the number of species present, 
# so that all variation relates to distribution within taxa not number of taxa

# boxplots and p-values, adjusting for confounders

if (var.type %in% c(0, 1)) { 
    
    ######################################## use "~ 1" (intercept only) in lm(.) if there is no confounder
    log.chao1.r = lm(log.chao1~Age.category+Ct.Status, data=metadata)$residuals+mean(log.chao1)  # Age.category+Ct.Status
    shannon.r = lm(shannon~Age.category+Ct.Status, data=metadata)$residuals+mean(shannon)
    #################################################################
    
    pdf(paste("plot_alpha_", var.name, "_adjAgeCt.pdf", sep=""), height=4.5, width=9)
    par(mfcol=c(1,2), pty="s")
    
    if (var.type==1) {
        boxplot(log.chao1.r ~ metadata$var, main="log Chao1", ylab="log Chao1", xlab=var.name, xaxt="n") # xaxt="n"
        text(x=1:var.nlevel, par("usr")[3], cex=0.7, labels = var.level, srt = 10, pos = 1, xpd = TRUE)
        # pvalue = signif(kruskal.test(log.chao1.r ~ metadata$var)$p.value, 2) # wilcox.test for two-level var; krustal.test for multi-level var
        log.chao1 = matrix(log.chao1, ncol=1)
        # pvalue1 <- ldm(log.chao1 | Age.category+Ct.Status ~ var,   # Age.category+Ethnicity
        #                      data=metadata, seed=67817,
        #                      scale.otu.table=FALSE, freq.scale.only=TRUE,
        #                      cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free")$p.global.freq
        legend(x="topright", legend=paste("p =", signif(pvalue1,3)), bty="n",cex = 1.0)
        
        boxplot(shannon.r ~ metadata$var, main="Shannon", ylab="Shannon", xlab=var.name, xaxt="n") 
        text(x=1:var.nlevel, par("usr")[3], cex=0.7, labels = var.level, srt = 10, pos = 1, xpd = TRUE)
        # pvalue = signif(kruskal.test(shannon.r ~ metadata$var)$p.value, 2)
        shannon = matrix(shannon, ncol=1)
        # pvalue2 <- ldm(shannon | Age.category+Ct.Status ~ var,   # Age.category+Ethnicity
        #               data=metadata, seed=67817,
        #               scale.otu.table=FALSE, freq.scale.only=TRUE,
        #               cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free")$p.global.freq
        legend(x="topright", legend=paste("p =", signif(pvalue2,3)), bty="n",cex = 1.0)
    } else if (var.type==0) {
        plot(y=log.chao1.r, x=metadata$var, main="log Chao1", ylab="log Chao1", xlab=var.name)                # scatter plot, to replace boxplot
        res <- summary(lm(log.chao1.r~metadata$var))       # fit linear regression
        abline(res$coefficients[,1])                               # add the fitted line to the scatter plot
        pvalue = signif(res$coefficient[2,4], 2)
        legend (x="bottomleft", legend=paste("p =", pvalue), bty="n",cex = 1.0) # add p-value for testing the slope=0
        
        plot(y=shannon.r, x=metadata$var, main="Shannon", ylab="Shannon", xlab=var.name)                        
        res <- summary(lm(shannon.r~metadata$var))       
        abline(res$coefficients[,1])                                                      
        pvalue = signif(res$coefficient[2,4], 2)
        legend (x="bottomleft", legend=paste("p =", pvalue), bty="n",cex = 1.0) 
    } 
    
    dev.off()

} else if (var.type==2) {
    (cox.fit.chao1 = coxph(Surv(var, var.status) ~ Age+Sex+Collect_Site+log.chao1, data=metadata, ties = "efron", control = coxph.control(timefix = FALSE)))
    (cox.fit.shannon = coxph(Surv(var, var.status) ~ Age+Sex+Collect_Site+shannon, data=metadata, ties = "efron", control = coxph.control(timefix = FALSE)))
}

#------------------------
#------------------------
# Beta diversity
#------------------------
#------------------------

if (var.type==2) {
    cox.fit = coxph(Surv(var, var.status) ~ Age+Sex+Collect_Site, data=metadata, ties = "efron", control = coxph.control(timefix = FALSE))
    metadata$var = residuals(cox.fit, type='martingale') # deviance
}

res.jaccard <- jaccard.mean( otu.tab ) 
distJ.mean.sq <- res.jaccard$jac.mean.sq.o2

(p.permanova.BC <- signif(permanovaFL(otu.tab | Age.category+Ct.Status  ~ var, 
                                      data = metadata, dist.method="bray", seed=67817,
                                      cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free")$p.permanova, 2)) # Age.category+Ethnicity
(p.permanova.J <- signif(permanovaFL(distJ.mean.sq | Age.category+Ct.Status ~ var, 
                                     data=metadata, seed=67817, square.dist=FALSE,
                                     cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free")$p.permanova, 2))  # parameters for requesting PERMANOVA-D2-A2
# for Age, 0.1 and 0.032

if (var.type==1) {
    
    w = which(colSums(otu.tab>0)==0)
    if (length(w)>0) {
        otu.tab = otu.tab[,-w]
    }
    dim(otu.tab)
    
    # Adjusting for confounders)
    ######################################### use "~ 1" in formula if there is no confounder
    distBC <- adjust.data.by.covariates(formula= ~ Age.category+Ct.Status, data=metadata, otu.table=otu.tab, dist.method="bray")$adj.dist  # Age.category+Ethnicity
    distJ <- adjust.data.by.covariates(formula= ~ Age.category+Ct.Status, data=metadata, dist=distJ.mean.sq, square.dist=FALSE)$adj.dist
    ##########################################
    
    distBC.eigen <- eigen(distBC, symmetric=TRUE)
    distJ.eigen <- eigen(distJ, symmetric=TRUE)
    PC.percent.BC = signif(distBC.eigen$values[1:2]/sum(distBC.eigen$values),3)*100
    PC.percent.J = signif(distJ.eigen$values[1:2]/sum(distJ.eigen$values),3)*100
    
    # PCoA plot
    
    col.palette = c("red","blue","green4","orange","black","gray","purple","palegreen", "cyan3", "plum", "peachpuff")
    pch.palette = rep(19, 11)
    col = as.character(factor(metadata$var, labels=col.palette[1:var.nlevel])) # assign colors to different levels of var
    pch = as.numeric(as.character(factor(metadata$var, labels=pch.palette[1:var.nlevel]))) # assign symbols to different levels of var
    level.index = c(1:var.nlevel)
    w = which(metadata$var %in% var.level[level.index])
    
    pdf(paste("plot_beta_", var.name, "_adjAgeCt.pdf", sep=""), height=4.5, width=11.5)
    
    par(mfrow=c(1,2), mar=c(5, 2.5, 4, 7.5), pty="s", xpd=TRUE)
    plot(distBC.eigen$vectors[,1],distBC.eigen$vectors[,2], main="Bray-Curtis", 
         xlab=paste("PC1 (", PC.percent.BC[1], "%)", sep=""), ylab=paste("PC2 (", PC.percent.BC[2], "%)", sep=""), 
         col=col[w], pch=pch[w]
         ,xlim=c(-0.16, 0.2), #xlim=c(min(distBC.eigen$vectors[,1])*, max(distBC.eigen$vectors[,1])*1.5), 
         ylim=c(-0.22, 0.3), # ylim=c(min(distBC.eigen$vectors[,2])*1.5, max(distBC.eigen$vectors[,2])*1.5)
         )
    ordiellipse(distBC.eigen$vectors, metadata$var, conf=0.9, col=col.palette[level.index])
    legend(x="bottomright", legend=paste("p =", p.permanova.BC), bty="n",cex=1.0)
    
    plot(distJ.eigen$vectors[,1],distJ.eigen$vectors[,2], main="Jaccard", 
         xlab=paste("PC1 (", PC.percent.J[1], "%)", sep=""), ylab=paste("PC2 (", PC.percent.J[2], "%)", sep=""), 
         col=col[w], pch=pch[w]
         ,xlim=c(min(distJ.eigen$vectors[,1])*1.5, max(distJ.eigen$vectors[,1])*1.5), 
         ylim=c(min(distJ.eigen$vectors[,2])*1.5, max(distJ.eigen$vectors[,2])*1.5)
         )
    ordiellipse(distJ.eigen$vectors, metadata$var, conf=0.9, col=col.palette[level.index])
    legend(x="bottomright", legend=paste("p =", p.permanova.J), bty="n",cex=1.0)
    
    legend(x="topright", inset=c(-0.75,0), bty="n", legend=var.level[level.index], col=col.palette[level.index], pch=pch.palette[level.index])
    
    dev.off()
}


#------------------------
#------------------------
# LDM (relative abundance)
#------------------------
#------------------------

w = which(colSums(otu.tab>0)<5)
otu.tab = otu.tab[,-w]
dim(otu.tab) # 229 248

fdr.nominal = 0.1

res.ldm <- ldm(otu.tab ~ var, data = metadata, dist.method="bray", n.rej.stop=100, seed=123, 
               cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free") # Age.category+Ethnicity
res.ldm$p.global.omni # 0.734

res.ldm$detected.otu.omni 

### order the detected OTUs by their p-values
w1 = match(res.ldm$detected.otu.omni[[1]], names(res.ldm$q.otu.omni))
o = w1[order(res.ldm$p.otu.omni[w1])]
#o = order(res.ldm$p.otu.omni)[1:10] ### Sankhya: uncomment this line if there is species detected and you want information of the species with the top 10 smallest p-values
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
pdf(paste("boxplot_", var.name, "_new.pdf", sep=""), height=6, width=15) 
par(mfrow=c(2,5), pty="s")
for (i in 1:n.detected.otu.omni) {
    boxplot(adj.data.1$y.freq[,o[i]]~metadata$var,
            main=names(res.ldm$q.otu.omni)[o[i]], ylab="Relative Abundance", xlab="", xaxt="n")
    text(x=1:var.nlevel, par("usr")[3], cex=0.7, labels = var.level, srt = 45, pos = 1, xpd = TRUE, offset=2.4)
    legend(x="topright", legend=paste("p =", signif(res.ldm$p.otu.omni[o[i]],3)), bty="n", cex=1.1)
}
dev.off()


#-----------------------------------------------
# LDM-A (presence-absence)
#-----------------------------------------------

res.ldmA <- ldm(formula=otu.tab ~ var, # Age.category+Ethnicity
                data=metadata, seed=123, n.rej.stop=100, 
                n.rarefy="all",
                cluster.id=otu.tab.sam, strata=metadata$strata, perm.within.type="free", perm.between.type="free") # parameter for requesting LDM-A
res.ldmA$p.global.pa # 0.834

res.ldmA$detected.otu.pa  


w1 = match(res.ldmA$detected.otu.pa[[1]], names(res.ldmA$q.otu.pa))
o = w1[order(res.ldmA$p.otu.pa[w1])]
# o = order(res.ldmA$p.otu.pa)[1:10] ### Sankhya: uncomment this line if there is species detected and you want information of the species with the top 10 smallest p-values
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
pdf(paste("boxplot_pa_", var.name, "_new.pdf", sep=""), height=6, width=15) 
par(mfrow=c(2,5), pty="s")
for (i in c(1:n.detected.otu.pa)) {
    boxplot(adj.phi[,o[i]] ~ metadata$var,
    main=names(res.ldmA$q.otu.pa)[o[i]], ylab="Chance of Presence", xlab="", xaxt="n")
    text(x=1:var.nlevel, par("usr")[3], cex=0.7, labels = var.level, srt = 45, pos = 1, xpd = TRUE, offset=2.4)
    legend(x="topright", legend=paste("p =", signif(res.ldmA$p.otu.pa[o[i]],3)), bty="n",cex=1.1)

}
dev.off()



