#input your node table and bucket table CSV files from your Classical molecular Networking job
node<-read.csv("node.csv", header=T)
bucketTab<-read.csv("bucket.csv", header=T)

#install.packages("dplyr") #if needed
#install.packages("vegan") #if needed
library("vegan")
library("dplyr")

#formatting node table
node_forbuck<-node
node<-as.data.frame(cbind(node$componentindex, node$UniqueFileSources))
colnames(node)<-c("GNPSFam", "samples")
node <- subset(node, GNPSFam != "-1") #remove singlets
node$samples <- gsub(".mzML|.mzXML", "", node$samples) #removes the .mzML from the sample names

#aggregating by scaffold family
node_concat <- node %>%
  group_by(GNPSFam) %>%
  summarize(samples = paste(samples, collapse = "|")) %>%
  ungroup()
node_concat_cleaned  <- node_concat %>%
  mutate(samples = sapply(strsplit(samples, "\\|"), function(x) paste(unique(x), collapse = "|")))

#Rational Library Build
#loop continues until 100% scaffold diversity is reached
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

#formatting bucket table, merging by families
node_forbuck<-as.data.frame(cbind(node_forbuck$cluster.index, node_forbuck$componentindex))
colnames(node_forbuck)<-c("ID", "GNPS.Fam")
node_forbuck$GNPS.Fam <- paste("Fam", node_forbuck$GNPS.Fam, sep = "")
colnames(bucketTab)<-sub("X.OTU.ID", "ID", colnames(bucketTab))
bucketTab<-merge(node_forbuck, bucketTab, by="ID", all=FALSE)
bucketTab <- bucketTab[, -which(names(bucketTab) == "ID")]
bucketTab <- subset(bucketTab, bucketTab$GNPS.Fam != "Fam-1") #removes singets
famMergedBucket <- aggregate(. ~ GNPS.Fam, data = bucketTab, FUN = sum)

#formating
row.names(famMergedBucket)<-famMergedBucket$GNPS.Fam
famMergedBucket <- famMergedBucket[, -which(names(famMergedBucket) == "GNPS.Fam")]
t_mergedBuck<-as.data.frame(t(famMergedBucket))
t_mergedBuck<- t_mergedBuck %>%mutate_all(~ ifelse(. > 0, 1, .))

# Making FAC --------------------------------------------------------------------
png("scaffold_accumulation_curve.png", width = 800, height = 600)  

totalfeatures <- ncol(t_mergedBuck)
par(mar = c(5, 5, 4, 5) + 0.1)  
accum_curve_rands <- specaccum(t_mergedBuck, method = "random", permutations = 50)
plot(accum_curve_rands,
     xlab='Number of Extracts',
     ylab='Percent of Total Scaffolds',
     col='orange2', lwd=2, 
     main="Scaffold Accumulation Curve",
     yaxt = "n", 
     xaxt = "n", 
     cex.lab=1.15, 
     cex.main=1.3)
abline(h=(totalfeatures*0.8), col = "black", lwd=1.5, lty = 2)
abline(h=(totalfeatures*0.95), col = "black", lwd=1.5, lty = 2)
abline(h=(totalfeatures*1), col = "black", lwd=1.5, lty = 2)

custom_y_labels <- c(0, 10, 20, 30, 40, 50, 60,70,80,90,100)
axis(2, label = custom_y_labels, at = ((custom_y_labels/100)*totalfeatures), cex.axis = 0.8)  

custom_x_ticks <- seq(0, 1400, by = 100)  # tick marks every 100
axis(1, at = custom_x_ticks, cex.axis = 0.8, labels = FALSE)  

custom_x_labels <- seq(0, 1400, by = 200)  #labels every 200
axis(1, label = custom_x_labels , at = custom_x_labels, cex.axis = 0.8)  

#Rational Library Curve
subset_data_ness <- t_mergedBuck[chosensamps, ]
accum_curve_opt_ness <- specaccum(subset_data_ness, method = "collector")
plot(accum_curve_opt_ness, col = "blue", add = TRUE, lwd = 5)

#legend
legend("bottomright", 
       legend = c("Random Selections", "Rational Library Method"), 
       col = c("orange2", "blue"), 
       lwd = 3, 
       cex=1.3)

dev.off()