# connected barplot of testdaat split on drug corr into sensitive and resistant populations

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(ggpubr)
library(ggalluvial)
library(rstatix)
library(cowplot)

obj <- obj.integrated
                                                  
unique(obj$Replicate)



################################################################################
# Does normalized enrichment for four neftel states show correlation with same 
# normalized enrichment in vivo?

inVivoSummaryData <- read.csv(file = "inVivoNeftelShiftSummary2.csv")
inVivoSummaryData



################################################################################


# sensResMarkers <- FindMarkers(obj, ident.1 = "Resistant", ident.2 = "Sensitive", test.use = "MAST", only.pos = F, slot = "scale.data")
# topResMarkers <- sensResMarkers %>% top_n(50, wt = avg_diff)
# topResMarkers

# # gsea via clusterProfiler...
# organism <- "org.Hs.eg.db"
# library(organism, character.only = T)
# 
# # Prepare gene list input...
# Res_gene_list <- sensResMarkers$avg_diff
# class(Res_gene_list)
# names(Res_gene_list) <- rownames(sensResMarkers)
# gene_list <- na.omit(Res_gene_list)
# gene_list <- sort(gene_list, decreasing = T)
# 
# gseResistance <- gseGO(geneList = gene_list,
#                        ont = "ALL",
#                        keyType = "SYMBOL",
#                        pvalueCutoff = 0.05,
#                        verbose = T,
#                        OrgDb = org.Hs.eg.db,
#                        pAdjustMethod = "BH",
#                        eps = 0)
# 
# require(DOSE)
# 
# library(EnhancedVolcano)
# EnhancedVolcano(toptable = sensResMarkers, 
#                 lab = rownames(sensResMarkers), 
#                 x = 'avg_diff', y = "p_val_adj",
#                 # ylim = c(0, max(-log10(sensResMarkers[["p_val_adj"]]), na.rm=TRUE) + 0.2)
#                 )
# 
# 
# 
# clusterProfiler::dotplot(gseResistance, showCategory = 9, split = ".sign") + facet_grid(.~.sign)
# 

#########################################################################
#########################################################################
obj@meta.data
groupingCounts <- table(obj@meta.data[,"singscorebin"],
                        obj@meta.data[,"CellCycleIdents"])

groupingCountsSplit <- as.data.frame(table(obj@meta.data[,"singscorebin"],
                                           obj@meta.data[,"CellCycleIdents"],
                                           obj@meta.data[,"Replicate"]))

groupingCountsSplit.2 <- as.data.frame(table(obj@meta.data[,"singscorebin"],
                                             obj@meta.data[,"CellCycleIdents"],
                                             obj@meta.data[,"Replicate"]))

# fix below...
colnames(groupingCountsSplit.2) <- c("State", "Idents", "patient_ID", "Freq") 

# relabel here... i.e. combine NPCs and MESs

for (i in 1:length(rownames(groupingCountsSplit.2))){
  print(groupingCountsSplit.2$State[i])
  if (groupingCountsSplit.2$State[i] %in% c("NPC1", "NPC2")){
    groupingCountsSplit.2$State2[i] <- "NPC"
  } else if (groupingCountsSplit.2$State[i] %in% c("MES1", "MES2")){
    groupingCountsSplit.2$State2[i] <- "MES"
  } else {
    groupingCountsSplit.2$State2[i] <- as.character(groupingCountsSplit.2$State[i])
  }
}


groupingCountsSplit.2$ptID_sens_state2 <- paste0(groupingCountsSplit.2$patient_ID, "_", groupingCountsSplit.2$Idents, "_", groupingCountsSplit.2$State2)
groupingCountsSplit.2.2 <- groupingCountsSplit.2[c("Freq", "ptID_sens_state2")]

groupingCountsSplit.2.2 <- as.data.frame(
  groupingCountsSplit.2.2 %>% 
    group_by(ptID_sens_state2) %>% 
    mutate(across(.fns = sum)) %>% 
    distinct()
)

# now re-expand column to ptID, sensitivity, state2.... 

for (i in 1:length(rownames(groupingCountsSplit.2.2))){
  groupingCountsSplit.2.2$patient_ID[i] <- unlist(strsplit(groupingCountsSplit.2.2$ptID_sens_state2[i], split = "_", fixed = T))[1]
  groupingCountsSplit.2.2$Idents[i] <- unlist(strsplit(groupingCountsSplit.2.2$ptID_sens_state2[i], split = "_", fixed = T))[2]
  groupingCountsSplit.2.2$State2[i] <- unlist(strsplit(groupingCountsSplit.2.2$ptID_sens_state2[i], split = "_", fixed = T))[3]
}

groupingCountsSplit.2.2$ptID_sens_state2 <- NULL

groupingCountsSplit.2 <- groupingCountsSplit.2.2

rm(groupingCountsSplit.2.2)

# groupingCountsSplit.2$Freq <- groupingCountsSplit.2$Freq + 1

#########################################################################
#########################################################################

sensitiveSum <- sum(groupingCounts[,"alisertib"])
resistantSum <- sum(groupingCounts[,"vehicle"])

#########################################################################
# Need sums for individual patients to calculate proportions as %

colnames(groupingCountsSplit) <- c("groupBy", "sensitivity", "ptID", "Freq")
groupingCountsSplit$ptID_sensitivity <- paste0(groupingCountsSplit$ptID, "_", groupingCountsSplit$sensitivity)
splitSums <- as.data.frame(groupingCountsSplit %>% 
                             dplyr::group_by(ptID_sensitivity) %>% 
                             # dplyr::group_by(sensitivity) %>% # added this...
                             dplyr::summarise(Freq = sum(Freq)))
rownames(splitSums) <- splitSums$ptID_sensitivity
for (i in 1:length(rownames(splitSums))){
  splitSums$ptID[i] <- unlist(strsplit(splitSums$ptID_sensitivity[i], split = "_", fixed = T))[1]  
}

splitSums

#########################################################################
#########################################################################
# Convert frequencies to % of whole
groupingCounts2 <- as.data.frame(groupingCounts)

colnames(groupingCounts2) <- c("groupBy", "sensitivity", "Freq","Replicate") 

for (i in 1:length(groupingCounts2$Freq)){
  if (groupingCounts2$sensitivity[i] == "alisertib"){
    groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/sensitiveSum*100
  }
  if (groupingCounts2$sensitivity[i] == "vehicle"){
    groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/resistantSum*100
  }
}

# for four states...

groupingCounts2.2 <- as.data.frame(groupingCounts)
colnames(groupingCounts2.2) <- c("State", "Sensitivity", "Freq") 

# relabel here... i.e. combine NPCs and MESs

for (i in 1:length(rownames(groupingCounts2.2))){
  print(groupingCounts2.2$State[i])
  if (groupingCounts2.2$State[i] %in% c("NPC1", "NPC2")){
    groupingCounts2.2$State2[i] <- "NPC"
  } else if (groupingCounts2.2$State[i] %in% c("MES1", "MES2")){
    groupingCounts2.2$State2[i] <- "MES"
  } else {
    groupingCounts2.2$State2[i] <- as.character(groupingCounts2.2$State[i])
  }
}

##########################################################################
#

splitSums
splitSums$ptID_sensitivity <- as.character(splitSums$ptID_sensitivity)
for (i in 1:length(rownames(groupingCountsSplit))){
  ptID <- as.character(groupingCountsSplit$ptID[i])
  print(ptID)
  ptID_senssum <- splitSums[paste0(ptID, "_alisertib"),]$Freq
  ptID_ressum <- splitSums[paste0(ptID, "_vehicle"),]$Freq
  print(ptID_senssum)
  print(ptID_ressum)
  if (groupingCountsSplit$sensitivity[i] == "alisertib"){
    groupingCountsSplit$Percent[i] <- groupingCountsSplit$Freq[i]/ptID_senssum*100
  }
  if (groupingCountsSplit$sensitivity[i] == "vehicle"){
    groupingCountsSplit$Percent[i] <- groupingCountsSplit$Freq[i]/ptID_ressum*100
  }
}
splitSums1 <- splitSums[splitSums$Freq != 0, ]
groupingCountsSplit.2n <- groupingCountsSplit.2[groupingCountsSplit.2$Freq !=0, ]
for (i in 1:length(rownames(groupingCountsSplit.2n))){
  ptID <- as.character(groupingCountsSplit.2n$patient_ID[i])
  print(ptID)
  ptID_senssum <- splitSums1[paste0(ptID, "_alisertib"),]$Freq
  ptID_ressum <- splitSums1[paste0(ptID, "_vehicle"),]$Freq
  print(ptID_senssum)
  print(ptID_ressum)
  if (groupingCountsSplit.2n$Idents[i] == "alisertib"){
    groupingCountsSplit.2n$Idents[i] <- groupingCountsSplit.2n$Freq[i]/ptID_senssum*100
  }
  if (groupingCountsSplit.2n$Sensitivity[i] == "vehicle"){
    groupingCountsSplit.2n$Percent[i] <- groupingCountsSplit.2n$Freq[i]/ptID_ressum*100
  }
}

groupingCountsSplit.2 <- groupingCountsSplit.2[groupingCountsSplit.2$State2 %in% c("NPC", "OPC", "AC", "MES"),]
groupingCountsSplit.2
##########################################################################
##########################################################################
# alluvial plots for individual patients... 
# then make a box plot of delta proportions using each patient as separate datapoint for stats...
# does the cut-off need to be reset for each individual patient? I.e. to apply 50/50 sens vs res?
# does individual patient data need to be scaled independently...? hope not. 
groupingCountsSplit$ptID <- as.character(groupingCountsSplit$ptID)
alluvialList <- list()
for (i in 1:length(unique(groupingCountsSplit$ptID))){
  print(i)
  ptID_search <- as.character(unique(groupingCountsSplit$ptID)[i])
  plotdf <- subset(groupingCountsSplit, ptID == as.character(ptID_search))
  plotdf$groupBy <- as.factor(plotdf$groupBy)
  plotdf$sensitivity <- factor(plotdf$sensitivity, levels = c("alisertib", "vehicle"))
  alluvsplit <- ggplot(plotdf, 
                       aes(x = sensitivity, stratum = groupBy, alluvium = groupBy,
                           y = Percent, fill = groupBy, label = groupBy)) + 
    scale_fill_brewer(type = "qual", palette = "Spectral") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgrey") +
    geom_stratum() + 
    theme(legend.position = "bottom") + theme_minimal() + 
    ggtitle(paste(unique(groupingCountsSplit$ptID)[i], "Transcriptional State Shift"))
  alluvialList[[i]] <- alluvsplit
}

plot_grid(plotlist = alluvialList, labels = "AUTO")

alluvialList.2 <- list()
for (i in 1:length(unique(groupingCountsSplit.2$patient_ID))){
  print(i)
  ptID_search <- as.character(unique(groupingCountsSplit.2$patient_ID)[i])
  plotdf <- subset(groupingCountsSplit.2, patient_ID == as.character(ptID_search))
  plotdf$State2 <- as.factor(plotdf$State2)
  plotdf$Sensitivity <- factor(plotdf$Sensitivity, levels = c("alisertib", "vehicle"))
  
  totalsum_alluv <- sum(plotdf$Freq)
  sensSum_alluv <- sum(plotdf[plotdf$Idents == "alisertib",]$Freq)
  resSum_alluv <- sum(plotdf[plotdf$Idents == "vehicle",]$Freq)
  respct_alluv <- resSum_alluv/totalsum_alluv
  senspct_alluv <- sensSum_alluv/totalsum_alluv
  
  alluvsplit <- ggplot(plotdf, 
                       aes(x = Sensitivity, stratum = State2, alluvium = State2,
                           y = Percent, fill = State2, label = State2)) + 
    scale_fill_brewer(type = "qual", palette = "Spectral") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgrey", width = 1/3) +
    geom_stratum(width = as.numeric(c(senspct_alluv, senspct_alluv, 
                                      senspct_alluv, senspct_alluv, 
                                      respct_alluv, respct_alluv, 
                                      respct_alluv, respct_alluv))) + 
    theme(legend.position = "bottom") + theme_minimal() + NoLegend()
  ggtitle(paste(unique(groupingCountsSplit.2$patient_ID)[i], "Transcriptional State Shift"))
  alluvialList.2[[i]] <- alluvsplit
}

pdf(file = "merge_individualPatient_FourState_alluvial.pdf", width = 10)
plot_grid(plotlist = alluvialList.2, labels = "AUTO")
dev.off()
##########################################################################
##########################################################################

# How do we plot change in proportion for each state, with datapoints for each individual patient...?
# Do we plot % resistant - % senstiive?
# or % gain or loss of proportion?
# maybe here just make a general plot, showing avg proportions of each starting...?

groupingCountsSplit
df2 <- data.frame()
for (i in 1:length(unique(groupingCountsSplit$ptID))){
  print(paste("i = ", i))
  ptID_search <- as.character(unique(groupingCountsSplit$ptID)[i])
  plotdf <- subset(groupingCountsSplit, ptID == as.character(ptID_search))
  print(plotdf$groupBy)
  for (j in 1:length(unique(plotdf$groupBy))){
    print(j)
    groupBy_search <- unique(plotdf$groupBy)[j]
    plotdf_state <- subset(plotdf, groupBy == groupBy_search)
    print(plotdf_state)
    rownames(plotdf_state) <- plotdf_state$sensitivity
    plotdf_state
    if (0 %in% plotdf_state$Percent){
      print("zero here...")
    } else {
      # extract resistant percent
      respct <- plotdf_state["vehicle", "Percent"]
      resFreq <- plotdf_state["vehicle", "Freq"]
      resFreq
      senspct <- plotdf_state["alisertib", "Percent"]
      sensFreq <- plotdf_state["alisertib", "Freq"]
      df1 <- data.frame(groupBy = plotdf_state$groupBy[1], 
                        ptID = plotdf_state$ptID[1], 
                        deltapct = respct - senspct)
      # df1.1 <- data.frame(groupBy = plotdf_state$State2[1],
      #                     ptID = plotdf_state$ptID[1],
      #                     detapct = respct - senspct)
      df2 <- rbind(df2, df1)
    }
  }
  # resistance - sensitive
}

df2 <- subset(df2, groupBy != "Non-Neoplastic")
df2

percShiftBP <- ggbarplot(df2, x = "groupBy", y = "deltapct", title = "ISOSCELES shift", 
                         ylab = "Change in % Proportion", xlab = "Transcriptional State", color = "groupBy", fill = "groupBy",
                         add = c("mean_se", "jitter"), position = position_dodge(.8)) + 
  scale_color_manual(values = c(OPC = "black", NPC1 = "black", NPC2 = "black", AC = "black", MES1 = "black", MES2 = "black"))

pdf(file = "merge_stateShiftByPatient.pdf")
percShiftBP
dev.off()

################################################################################

# Four states...

groupingCountsSplit.2$Percent <- groupingCountsSplit.2$Percent + 1
df2.2 <- data.frame()
for (i in 1:length(unique(groupingCountsSplit.2$patient_ID))){
  print(paste("i = ", i))
  ptID_search <- as.character(unique(groupingCountsSplit.2$patient_ID)[i])
  plotdf <- subset(groupingCountsSplit.2, patient_ID == as.character(ptID_search))
  plotdf$State <- NULL
  print(plotdf)
  print(plotdf$State2)
  for (j in 1:length(unique(plotdf$State2))){
    print(j)
    groupBy_search <- unique(plotdf$State2)[j]
    plotdf_state <- subset(plotdf, State2 == groupBy_search)
    print(plotdf_state)
    rownames(plotdf_state) <- plotdf_state$Sensitivity
    plotdf_state
    if (0 %in% plotdf_state$Percent){
      print("zero here...")
    } else {
      # extract resistant percent
      respct <- plotdf_state["vehicle", "Percent"]
      resFreq <- plotdf_state["vehicle", "Freq"]
      resFreq
      senspct <- plotdf_state["alisertib", "Percent"]
      sensFreq <- plotdf_state["alisertib", "Freq"]
      df1.2 <- data.frame(State2 = plotdf_state$State2[1], 
                          ptID = plotdf_state$patient_ID[1], 
                          deltapct = respct - senspct)
      # df1.1 <- data.frame(groupBy = plotdf_state$State2[1],
      #                     ptID = plotdf_state$ptID[1],
      #                     detapct = respct - senspct)
      df2.2 <- rbind(df2.2, df1.2)
    }
  }
  # resistance - sensitive
}

df2.2 <- subset(df2.2, State2 != "Non-Neoplastic")

colnames(df2.2)

df2.2
################################################################################
# ANOVA
################################################################################

res.aov <- aov(formula = deltapct~State2, data = df2.2)
res.aov

tukey <- TukeyHSD(res.aov)
tukey

gh_res <- games_howell_test(data = df2.2, formula = deltapct~State2, conf.level = 0.95, detailed = F)
gh_res$comp <- paste0(gh_res$group1, "_", gh_res$group2)
gh_res$inv.p.adj <- -log(gh_res$p.adj)

pdf(file = "merge_inv.p.adj.barplot.pdf")
ggbarplot(gh_res[order(gh_res$inv.p.adj, decreasing = T),], y = "inv.p.adj", x = "comp", ggtheme = theme_minimal(), fill = "comp", xlab = "State Comparison", ylab = "Probability (-log(p.adj))") + 
  NoLegend()
dev.off()

################################################################################
# T tests...
################################################################################
# 1 MES_NPC

MES_NPC <- subset(df2.2, df2.2$State2 %in% c("MES", "NPC"))


MES_NPC_ttest_res <- t.test(data = MES_NPC, deltapct ~ State2 , paired = F)
MES_NPC_ttest_res$p.value

# paired
MES_NPC_2 <- reshape(MES_NPC, direction = "wide", idvar = "ptID", timevar = "State2" )


MES_NPC_ttest_paired_res <- t.test(Pair(deltapct.MES, deltapct.NPC) ~ 1, data = MES_NPC_2)
MES_NPC_ttest_paired_res$p.value

################################################################################
# 2 MES_OPC

MES_OPC <- subset(df2.2, df2.2$State2 %in% c("MES", "OPC"))


MES_OPC_ttest_res <- t.test(data = MES_OPC, deltapct ~ State2 , paired = F)
MES_OPC_ttest_res$p.value

# paired
MES_OPC_2 <- reshape(MES_OPC, direction = "wide", idvar = "ptID", timevar = "State2" )


MES_OPC_ttest_paired_res <- t.test(Pair(deltapct.MES, deltapct.OPC) ~ 1, data = MES_OPC_2)
MES_OPC_ttest_paired_res$p.value

################################################################################
# 3 MES_AC

MES_AC <- subset(df2.2, df2.2$State2 %in% c("MES", "AC"))


MES_AC_ttest_res <- t.test(data = MES_AC, deltapct ~ State2 , paired = F)
MES_AC_ttest_res$p.value

# paired
MES_AC_2 <- reshape(MES_AC, direction = "wide", idvar = "ptID", timevar = "State2" )


MES_AC_ttest_paired_res <- t.test(Pair(deltapct.MES, deltapct.AC) ~ 1, data = MES_AC_2)
MES_AC_ttest_paired_res$p.value

################################################################################
# 4 NPC_OPC

NPC_OPC <- subset(df2.2, df2.2$State2 %in% c("NPC", "OPC"))


NPC_OPC_ttest_res <- t.test(data = NPC_OPC, deltapct ~ State2 , paired = F)
NPC_OPC_ttest_res$p.value

# paired
NPC_OPC_2 <- reshape(NPC_OPC, direction = "wide", idvar = "ptID", timevar = "State2" )


NPC_OPC_ttest_paired_res <- t.test(Pair(deltapct.NPC, deltapct.OPC) ~ 1, data = NPC_OPC_2)
NPC_OPC_ttest_paired_res$p.value

################################################################################
# 5 AC_OPC

AC_OPC <- subset(df2.2, df2.2$State2 %in% c("AC", "OPC"))


AC_OPC_ttest_res <- t.test(data = AC_OPC, deltapct ~ State2 , paired = F)
AC_OPC_ttest_res$p.value

# paired
AC_OPC_2 <- reshape(AC_OPC, direction = "wide", idvar = "ptID", timevar = "State2" )


AC_OPC_ttest_paired_res <- t.test(Pair(deltapct.AC, deltapct.OPC) ~ 1, data = AC_OPC_2)
AC_OPC_ttest_paired_res$p.value

################################################################################
# 6 NPC_AC

NPC_AC <- subset(df2.2, df2.2$State2 %in% c("NPC", "AC"))


NPC_AC_ttest_res <- t.test(data = NPC_AC, deltapct ~ State2 , paired = F)
NPC_AC_ttest_res$p.value

# paired
NPC_AC_2 <- reshape(NPC_AC, direction = "wide", idvar = "ptID", timevar = "State2" )


NPC_AC_ttest_paired_res <- t.test(Pair(deltapct.NPC, deltapct.AC) ~ 1, data = NPC_AC_2)
NPC_AC_ttest_paired_res$p.value

################################################################################
# List and turn into a dataframe to then put into illustrator...

MES_NPC_ttest_res$p.value
MES_NPC_ttest_paired_res$p.value
MES_OPC_ttest_res$p.value
MES_OPC_ttest_paired_res$p.value
MES_AC_ttest_res$p.value
MES_AC_ttest_paired_res$p.value
NPC_AC_ttest_res$p.value
NPC_AC_ttest_paired_res$p.value
NPC_OPC_ttest_res$p.value
NPC_OPC_ttest_paired_res$p.value
AC_OPC_ttest_res$p.value
AC_OPC_ttest_paired_res$p.value

################################################################################

percShiftBP.2 <- ggbarplot(df2.2, x = "State2", y = "deltapct", title = "ISOSCELES                                                             shift", 
                           ylab = "Change in % Proportion", xlab = "Transcriptional State", color = "State2", fill = "State2", order = c("AC", "MES", "NPC", "OPC"),
                           add = c("mean_se", "jitter"), position = position_dodge(.8)) + 
  scale_color_manual(values = c(OPC = "black", NPC = "black", AC = "black", MES = "black"))

pdf(file = "merge_stateShiftByPatient_FourStates.pdf")
percShiftBP.2
dev.off()

percShiftBP.3 <- ggbarplot(df2.2, x = "State2", y = "deltapct", title = "ISOSCELES shift", size = 1,
                           ylab = "Change in % Proportion", xlab = "Transcriptional State", color = "State2", fill = "State2", order = c("AC", "MES", "NPC", "OPC"),
                           add = c("mean_se"), position = position_dodge(.8)) + 
  scale_color_manual(values = c(OPC = "black", NPC = "black", AC = "black", MES = "black")) + 
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", size  = 1, color = "black") + 
  theme(axis.line=element_line(size = 1))

pdf(file = "merge_stateShiftByPatient_FourStates_v2.pdf")
percShiftBP.3
dev.off()
##########################################################################
##########################################################################

## Plot FC instead...

groupingCountsSplit
df3 <- data.frame()
for (i in 1:length(unique(groupingCountsSplit$ptID))){
  print(paste("i = ", i))
  ptID_search <- as.character(unique(groupingCountsSplit$ptID)[i])
  plotdf <- subset(groupingCountsSplit, ptID == as.character(ptID_search))
  print(plotdf$groupBy)
  for (j in 1:length(unique(plotdf$groupBy))){
    print(j)
    groupBy_search <- unique(plotdf$groupBy)[j]
    plotdf_state <- subset(plotdf, groupBy == groupBy_search)
    print(plotdf_state)
    rownames(plotdf_state) <- plotdf_state$sensitivity
    plotdf_state
    if (0 %in% plotdf_state$Percent){
      print("zero here...")
    } else {
      # extract resistant percent
      respct <- plotdf_state["vehicle", "Percent"]
      resFreq <- plotdf_state["vehicle", "Freq"]
      resFreq
      senspct <- plotdf_state["alisertib", "Percent"]
      print(plotdf_state)
      df1 <- data.frame(groupBy = plotdf_state$groupBy[1], 
                        ptID = plotdf_state$patient_ID[1], 
                        fc = (resFreq - sensFreq) / sensFreq)
      df3 <- rbind(df3, df1)
    }
  }
  # resistance - sensitive
}

df3 <- subset(df3, groupBy != "Non-Neoplastic")

ggbarplot(df3, x = "groupBy", y = "fc", title = "ISOSCELES shift", 
          ylab = "Change in % Proportion", xlab = "Transcriptional State", color = "groupBy", fill = "groupBy",
          add = c("mean_se", "jitter"), position = position_dodge(.8)) + 
  scale_color_manual(values = c(OPC = "black", NPC1 = "black", NPC2 = "black", AC = "black", MES1 = "black", MES2 = "black"))

##########################################################################
##########################################################################

groupingCounts2$groupBy <- as.factor(groupingCounts2$groupBy)
groupingCounts2$sensitivity <- factor(groupingCounts2$sensitivity, levels = c("alisertib", "vehicle"))

library(ggalluvial)
is_alluvia_form(groupingCounts2)
alluv <- ggplot(groupingCounts2,
                aes(x = sensitivity, stratum = groupBy, alluvium = groupBy,
                    y = Percent,
                    fill = groupBy, label = groupBy)) +
  scale_fill_brewer(type = "qual", palette = "Spectral") +
  # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") + theme_minimal() +
  ggtitle("Transcriptional State Shift")

pdf("alisertib_cutOff_0_alluvial.pdf")
alluv
dev.off()

################################################################################
# compare johnson et al isosceles output to in vivo experiment...
################################################################################

inVivoGroupingCounts <- read.csv("inVivoAli_freqPerc_shiftAlluvialData_2.csv", row.names = 1)

groupingCounts3 <- groupingCounts2
colnames(groupingCounts3) <- c("State", "Sensitivity", "Freq", "Percent")

for (i in 1:length(rownames(groupingCounts3))){
  print(groupingCounts3$State[i])
  if (groupingCounts3$State[i] %in% c("NPC1", "NPC2")){
    groupingCounts3$State2[i] <- "NPC"
  } else if (groupingCounts3$State[i] %in% c("MES1", "MES2")){
    groupingCounts3$State2[i] <- "MES"
  } else {
    groupingCounts3$State2[i] <- as.character(groupingCounts3$State[i])
  }
}

alluv2 <- ggplot(groupingCounts3,
                 aes(x = Sensitivity, stratum = State2, alluvium = State,
                     y = Percent,
                     fill = State2, label = State2)) +
  scale_fill_brewer(type = "qual", palette = "Spectral") +
  # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") + theme_minimal() +
  ggtitle("Transcriptional State Shift ")

pdf("FourState_merge_alluvial.pdf")
alluv2
dev.off()

################################################################################
# Scatter plot for correlation?

# First, merge data for comparison...

inVivoGroupingCounts$State <- NULL 


inVivoGroupingCounts$arm_state2 <- paste0(inVivoGroupingCounts$Arm, "_", inVivoGroupingCounts$State2)


inVivoGroupingCounts <- inVivoGroupingCounts[c("Freq", "Percent", "arm_state2")]

inVivoGroupingCounts <- as.data.frame(
  inVivoGroupingCounts %>% 
    group_by(arm_state2) %>% 
    mutate(across(.fns = sum)) %>% 
    distinct()
)

# now re-expand column to ptID, sensitivity, state2.... 

for (i in 1:length(rownames(inVivoGroupingCounts))){
  inVivoGroupingCounts$Arm[i] <- unlist(strsplit(inVivoGroupingCounts$arm_state2[i], split = "_", fixed = T))[1]
  inVivoGroupingCounts$State2[i] <- unlist(strsplit(inVivoGroupingCounts$arm_state2[i], split = "_", fixed = T))[2]
}

inVivoGroupingCounts$arm_state2 <- NULL

# experimental results...
# Calculate actual shift here... i.e.
df2.3 <- data.frame()
for (j in 1:length(unique(inVivoGroupingCounts$State2))){
  print(j)
  groupBy_search <- unique(inVivoGroupingCounts$State2)[j]
  print(groupBy_search)
  plotdf_state <- subset(inVivoGroupingCounts, State2 == groupBy_search)
  print(plotdf_state)
  rownames(plotdf_state) <- plotdf_state$Arm
  print(plotdf_state)
  if (0 %in% plotdf_state$Percent){
    print("zero here...")
  } else {
    # extract resistant percent
    respct <- plotdf_state["alisertib", "Percent"]
    resFreq <- plotdf_state["alisertib", "Freq"]
    resFreq
    senspct <- plotdf_state["vehicle", "Percent"]
    sensFreq <- plotdf_state["vehicle", "Freq"]
    df1.3 <- data.frame(groupBy = plotdf_state$State2[1], 
                        # ptID = plotdf_state$ptID[1], 
                        deltapct = respct - senspct)
    # df1.1 <- data.frame(groupBy = plotdf_state$State2[1],
    #                     ptID = plotdf_state$ptID[1],
    #                     detapct = respct - senspct)
    df2.3 <- rbind(df2.3, df1.3)
  }
}

df2.2mean <- df2.2
df2.2mean$ptID <- NULL
df2.2mean <- as.data.frame(
  df2.2mean %>% 
    group_by(State2) %>% 
    mutate(across(.fns = sum)) %>% 
    distinct()
)

compdf <- merge(df2.2mean, df2.3, by.x = "State2", by.y = "groupBy")
colnames(compdf) <- c("State", "Johnson", "InVivo")

compdf2 <- merge(df2.2, df2.3, by.x = "State2", by.y = "groupBy")
colnames(compdf2) <- c("State", "ptID", "Johnson", "InVivo")

InVivoShift <- ggscatter(df2.2, x = "state", y = "ptID", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                                  conf.int = T) + stat_cor(method = "pearson")

johnsonVsInVivoShift2 <- ggscatter(compdf2, x = "Johnson", y = "InVivo", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                                   conf.int = T) + stat_cor(method = "pearson")

pdf("mergeVsInVivoShiftScatter.pdf")
InVivoShift
dev.off()

pdf("mergeVsInVivoShiftScatter_IndividualPoints.pdf")
johnsonVsInVivoShift2
dev.off()
# switch groupingCounts3 to percent shift... i.e. resist - sens?


# Initialize df2.3
df2.3 <- data.frame()

# Loop through unique State2 values
for (j in 1:length(unique(groupingCounts3$State2))) {
  groupBy_search <- unique(groupingCounts3$STATE2)[j]
  plotdf_state <- subset(groupingCounts3, State2 == groupBy_search)
  
  if (0 %in% plotdf_state$Percent) {
    # Skip calculations when percentage is zero
    next
  } else {
    # extract resistant percent
    respct <- plotdf_state["alisertib", "Percent"]
    senspct <- plotdf_state["vehicle", "Percent"]
    df1.3 <- data.frame(groupBy = plotdf_state$State2[1], deltapct = respct - senspct)
    df2.3 <- rbind(df2.3, df1.3)
  }
}


# Print the resulting data frame

df2.3mean <- df2.3
df2.3mean$ptID <- NULL
df2.3mean <- as.data.frame(
  df2.3mean %>% 
    group_by(groupBy) %>% 
    mutate(across(.fns = sum)) %>% 
    distinct()
)

compdf <- merge(df2.3mean, df2.3, by.x = "State2", by.y = "groupBy")
colnames(compdf) <- c("State", "Johnson", "InVivo")

compdf2 <- merge(df2.2, df2.3, by.x = "State2", by.y = "groupBy")
colnames(compdf2) <- c("State", "ptID", "Johnson", "InVivo")

InVivoShift <- ggscatter(gc, x = "State2", y = "patient_ID", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                         conf.int = T) + stat_cor(method = "pearson")

johnsonVsInVivoShift2 <- ggscatter(compdf2, x = "Johnson", y = "InVivo", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                                   conf.int = T) + stat_cor(method = "pearson")

##


#gc <- groupingCountsSplit.2[groupingCountsSplit.2$Freq !=0, ]
gc <- groupingCountsSplit[groupingCountsSplit$Freq !=0, ]

gc <- groupingCounts3
combine_groups <- function(group) {
  if (group %in% c("MES1", "MES2")) {
    return("MES")
  } else if (group %in% c("NPC1", "NPC2")) {
    return("NPC")
  } else if (group == "AC") {
    return("AC")
  }
  else if (group == "OPC") {
    return("OPC")
  }
}


# Print the updated data frame
print(groupingCounts2.2)



gc$State2 <- sapply(gc$groupBy, combine_groups)

combined_df <- gc %>%
  group_by(ptID, State2) %>%
  summarise(
    Freq = mean(Freq),
    Percent = mean(Percent),
    .groups = 'drop'
  )

combined_df <- combined_df %>%
  mutate(State2 = ifelse(State2 %in% c('MES1', 'MES2'), 'MES', State2))

# Similarly, combine the 'NPC1' and 'NPC2' rows into one row with 'State2' as 'NPC'
combined_df <- combined_df %>%
  mutate(State2 = ifelse(State2 %in% c('NPC1', 'NPC2'), 'NPC', State2))

# Print the resulting data frame
print(combined_df)



gc1 <- as.data.frame(combined_df)

# Print the updated dataframe
print(gc)
dev.off()

# Load the necessary libraries
library(ggplot2)
library(ggpubr)

gg <- ggplot(data = gc, aes(x = State2, y = Percent, color = ptID)) +
  geom_jitter(width = 0.2) +
  labs(x = "State2", y = "Percent") +
  scale_color_discrete(name = "patient_ID")

gg1 <- ggplot(data = gc, aes(x = State2, y = Percent, fill = ptID)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "State2", y = "Percent") +
  scale_fill_discrete(name = "patient_ID") +
  theme_minimal()

gg1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))



metadata <- obj.integrated@meta.data
cell_counts <- table(metadata$orig.ident, metadata$singscorebin)
print(cell_counts)

table(obj.integrated$Replicate)
