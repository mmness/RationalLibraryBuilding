#input node and bucket table from the classical molecular networking job
node<- read.csv("node.csv", header=TRUE)
bucket<-read.csv("bucket.csv", header=TRUE)

#bioactivity data
hitDat<-read.csv("BioactivityDataExample.csv", header=TRUE)

#this code uses a pearson correlation with a boneferri correction method. It filters out features negatively correlated together you can adjust to what fits your data

# formatting --------------------------------------------------------------
Ratnod<-node
node$ID<-paste("X",node$precursor.mass, round(node$RTMean,1), node$cluster.index, sep = "_")
node<-as.data.frame(cbind(node$cluster.index, node$ID))
colnames(node)<-c("cluster.index", "ID")

colnames(bucket)<-sub("X.OTU.ID", "cluster.index", colnames(bucket))
bucket<-merge(node, bucket, by="cluster.index", all=FALSE)  ##change to all=TRUE once this issue gets worked out

bucket <- bucket[, -which(names(bucket) == "cluster.index")]
t_buck <- as.data.frame(t(as.matrix(bucket)))
colnames(t_buck)<-bucket$ID
t_buck <- t_buck[!(rownames(t_buck) == "ID"), ]
t_buck[] <- lapply(t_buck, as.numeric)

#add in bioactivity data
colnames(hitDat)<-c("sampleid", "inf")
t_buck<-cbind(sampleid = rownames(t_buck), t_buck)
t_buck<-merge(hitDat, t_buck, by="sampleid", all=TRUE)
colnames(t_buck)<-sub("sampleid", "Sample_name", colnames(t_buck))
colnames(t_buck)<-sub("inf", "Bioactiv", colnames(t_buck))

# Bioactive table codes. -------------------------------------------------
tab<-t_buck
if(any(is.na(tab[,2]))) tab <- tab[!is.na(tab[,2]),] # Take out blank rows in the table
tab2 <- tab
tab2[,-c(1:2)] <- t(apply(tab2[,-c(1:2)], 1, function(x) (x+1)/sum((x+1)))) #+1 to all and normalize
tab2[,-1] <- lapply(tab2[,-1], as.numeric) 

ct <- t(sapply(3:ncol(tab2), function(x) unlist(cor.test(scale(tab2[, 2])[, 1], scale(tab2[, x])[, 1], method = "pearson")[c("estimate", "p.value")])))

ct <- rbind(c("cor"," p_value"), c(0,0), ct)
tab3 <- rbind(t(ct),  as.matrix(tab2))
rownames(tab3) <- NULL

new = t(tab3)
colnames(new) = new[1,]
new = new[-1,]
new = cbind(0:(nrow(new)-1), rownames(new), new)
rownames(new) <- NULL
colnames(new)[1:2] <- c("shared name", "IDs")
new[1,1] <- ""

# sigfeats ----------------------------------------------------------------
# Show the features ID with correlation coefficient
nm <- colnames(tab)
new <- cbind(new[,1:4], c(0, p.adjust(as.numeric(ct[-c(1:2),2]), method = "bonferroni")), new[,-c(1:4)]) #can change p value correction method 
colnames(new)[5] <- "p_value_corrected"

#MN-->subset with significant features
new_subset<-as.data.frame(new[-1,])
new_subset<-new_subset[new_subset$p_value_corrected < 0.05, ]
new_subset<-new_subset[new_subset$cor > 0, ] #remove negative correlations #can adjust if %growth was measured instead

signif_feats<-new_subset$IDs

# rational library building and comparisons--------------------------------------------------------------------
Ratnod<-as.data.frame(cbind(Ratnod$componentindex, Ratnod$UniqueFileSources))
colnames(Ratnod)<-c("GNPSFam", "samples")

Ratnod <- subset(Ratnod, GNPSFam != "-1") #remove singlets
Ratnod$samples <- gsub(".mzML", "", Ratnod$samples) #removes the .mzML from the sample names (not needed)

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

#how many of the chosen do we take?
samps80<- chosensamps[1:sum(chosensamps_dat$percent_div_reached < 80)]
samps95<-chosensamps[1:sum(chosensamps_dat$percent_div_reached < 95)]
samps100<-chosensamps


# formatting --------------------------------------------------------------
#rearrange t_buck 
table<-t_buck
row.names(table)<-table[,1]
table<-table[,-c(1,2)]
#subsetbucket
bucket80<-subset(table, row.names(table) %in% samps80)
bucket95<- subset(table, row.names(table) %in% samps95)
bucket100<- subset(table, row.names(table) %in% samps100)

#remove feats not in rational libraries
sums <- colSums(bucket80)
zero_sum_columns <- names(sums[sums == 0])
bucket80 <- bucket80[, !(colnames(bucket80) %in% zero_sum_columns)]

sums <- colSums(bucket95)
zero_sum_columns <- names(sums[sums == 0])
bucket95 <- bucket95[, !(colnames(bucket95) %in% zero_sum_columns)]

sums <- colSums(bucket100)
zero_sum_columns <- names(sums[sums == 0])
bucket100 <- bucket100[, !(colnames(bucket100) %in% zero_sum_columns)]

rat80<- colnames(bucket80)
rat95<- colnames(bucket95)
rat100<- colnames(bucket100)

# feature retention analysis -----------------------------------------------------------
#signif_feats has the p-adjusted sig features from full library
intersect_r80_sigfeats <- intersect(rat80, signif_feats)
intersect_r95_sigfeats <- intersect(rat95, signif_feats)
intersect_r100_sigfeats <- intersect(rat100, signif_feats)

# exporting data ----------------------------------------------------------

#writing data
write.csv(new, "matrix_with_correlation_data.csv", row.names=FALSE)
write.csv(signif_feats, "significant_features.csv")

#printing significant feature data
cat("Total significant features:", length(signif_feats), "\n")
cat("retained in 80% scaffold library:", length(intersect_r80_sigfeats), "\n")
cat("retained in 95% scaffold library:", length(intersect_r95_sigfeats), "\n")
cat("retained in 100% scaffold library:", length(intersect_r100_sigfeats), "\n")


