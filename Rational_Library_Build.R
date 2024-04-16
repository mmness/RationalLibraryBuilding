#input your node table from the Classical Molecular Networking Job
node<-read.csv("node.csv", header=T) 

#install.packages("dplyr") #if needed
library(dplyr)

#formatting
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

#formatting
max_length<-length(chosensamps) 
RL80<- chosensamps[1:sum(chosensamps_dat$percent_div_reached < 80)]
RL95<-chosensamps[1:sum(chosensamps_dat$percent_div_reached < 95)]
RatLibrary_100scaffold_diversity<-chosensamps
RatLibrary_80scaffold_diversity <- c(RL80, rep(NA, max_length - length(RL80)))
RatLibrary_95scaffold_diversity <- c(RL95, rep(NA, max_length - length(RL95)))
summary<- cbind(RatLibrary_100scaffold_diversity,RatLibrary_95scaffold_diversity,RatLibrary_80scaffold_diversity)

#export
write.csv(summary, "RationalLibraries.csv", na="" )





