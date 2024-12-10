#read in your data. 
#The sample names in the "featureinfo" file must match sample names in the node table

featTab<-read.csv("featureinfo.csv", header = TRUE)
node<- read.csv("node.csv", header=TRUE)
hitDat<-read.csv("BioactivityData.csv", header=TRUE)

# formatting --------------------------------------------------------------
colnames(hitDat)<-c("sampleid", "Bioactivity")
tab<-merge(hitDat, featTab, by="sampleid")

# Bioactive table codes. -------------------------------------------------
if(any(is.na(tab[,2]))) tab <- tab[!is.na(tab[,2]),] 
tab2 <- tab

num_feat<-ncol(tab2)-2
correlation_data <- data.frame(features = rep(NA, num_feat), correlation_Coeff = rep(NA, num_feat), pvalue = rep(NA, num_feat))
correlation_data$features<-colnames(tab2)[-c(1,2)]

for (i in 1:nrow(correlation_data)) {  
  feat <- correlation_data$features[i]
  chemical_data <- data.frame(Bioactivity = tab2$Bioactivity, Feat_ab = tab2[[feat]])
  chemical_data <- subset(chemical_data, Feat_ab != 0)
  
  if (nrow(chemical_data) <= 3) {
    #if statements for if there are not enough data points for correlation
    correlation_data$pvalue[i] <- NA
    correlation_data$correlation_Coeff[i] <- NA
  } else {
    chemical_data$Bioactivity <- as.numeric(chemical_data$Bioactivity)
    chemical_data$Feat_ab <- as.numeric(chemical_data$Feat_ab)
    
    correlation_test <- cor.test(chemical_data$Bioactivity, chemical_data$Feat_ab, method = "spearman")
    
    correlation_data$pvalue[i] <- correlation_test$p.value
    correlation_data$correlation_Coeff[i] <- correlation_test$estimate
  }
}

#make pval-adjust column
correlation_data$pvalues_adjusted <- p.adjust(correlation_data$pvalue, method = "fdr", n = length(na.omit(correlation_data$pvalue)))

#make shared name column
correlation_data$`shared id` <- sub(".*_(\\d+)$", "\\1", correlation_data$features)

# rational library building and comparisons--------------------------------------------------------------------
Ratnod<-as.data.frame(cbind(node$componentindex, node$UniqueFileSources))
colnames(Ratnod)<-c("GNPSFam", "samples")

Ratnod <- subset(Ratnod, GNPSFam != "-1") #remove singlets
Ratnod$samples<-gsub(".mzML", "", Ratnod$samples)

library(dplyr)
node_concat <- Ratnod %>%
  group_by(GNPSFam) %>%
  summarize(samples = paste(samples, collapse = "|")) %>%
  ungroup()
node_concat_cleaned  <- node_concat %>%
  mutate(samples = sapply(strsplit(samples, "\\|"), function(x) paste(unique(x), collapse = "|")))
df<-node_concat_cleaned
chosensamps <- c()
percent_div_reached<-c()
while(nrow(df) > 0) {
  most_frequent_sample <- names(sort(table(unlist(strsplit(df$samples, "\\|"))), decreasing = TRUE))[1]
  chosensamps <- c(chosensamps, most_frequent_sample)
  df <- df %>%
    filter(!grepl(most_frequent_sample, samples))
  percent_div_reached<-c(percent_div_reached,  100-((nrow(df)/nrow(node_concat_cleaned))*100)  )
}
chosensamps_dat<-as.data.frame(cbind(chosensamps, percent_div_reached))
chosensamps_dat$chosensamps<-gsub(".mzML", "", chosensamps_dat$chosensamps)


#how many of the chosen are taken
samps80<- chosensamps[1:sum(chosensamps_dat$percent_div_reached < 80)]
samps95<-chosensamps[1:sum(chosensamps_dat$percent_div_reached < 95)]
samps100<-chosensamps

# sigfeats ----------------------------------------------------------------
sig_subset<-subset(correlation_data, correlation_data$pvalues_adjusted <.05)#pval cuttoff
sig_subset<-subset(sig_subset, sig_subset$correlation_Coeff >0.5) #coefficient cuttoff

signif_feats<-sig_subset$features
# formatting --------------------------------------------------------------

#rearrange tab2
ret_table<- tab2
row.names(ret_table)<-ret_table[,1]
ret_table<-ret_table[ -c(1,2)]

row.names(ret_table)<-sub(" ", "", row.names(ret_table))

#subsetbucket
ratlib80<-subset(ret_table, row.names(ret_table)  %in% samps80)
ratlib95<- subset(ret_table, row.names(ret_table) %in% samps95)
ratlib100<- subset(ret_table, row.names(ret_table) %in% samps100)

#remove feats not in rational libraries
sums <- colSums(ratlib80)
zero_sum_columns <- names(sums[sums == 0])
ratlib80 <- ratlib80[, !(colnames(ratlib80) %in% zero_sum_columns)]

sums <- colSums(ratlib95)
zero_sum_columns <- names(sums[sums == 0])
ratlib95 <- ratlib95[, !(colnames(ratlib95) %in% zero_sum_columns)]

sums <- colSums(ratlib100)
zero_sum_columns <- names(sums[sums == 0])
ratlib100 <- ratlib100[, !(colnames(ratlib100) %in% zero_sum_columns)]

rat80<- colnames(ratlib80)
rat95<- colnames(ratlib95)
rat100<- colnames(ratlib100)

# feature retention analysis -----------------------------------------------------------
#signif_feats has the p-adjusted sig features from full library
intersect_r80_sigfeats <- intersect(rat80, signif_feats)
intersect_r95_sigfeats <- intersect(rat95, signif_feats)
intersect_r100_sigfeats <- intersect(rat100, signif_feats)

# exporting data ----------------------------------------------------------
write.csv(correlation_data, "correlationData.csv")

#printing significant feature data
cat("Total significant features:", length(signif_feats), "\n")
cat("retained in 80% scaffold library:", length(intersect_r80_sigfeats), "\n")
cat("retained in 95% scaffold library:", length(intersect_r95_sigfeats), "\n")
cat("retained in 100% scaffold library:", length(intersect_r100_sigfeats), "\n")

