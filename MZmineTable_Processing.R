# OVERSPLIT REMOVAL -------------------------------------------------------
feat<-read.csv("quant.csv") 
mz_tol<-0.0005
rt_tol<-0.2
oversplit_num<-5

feat$row.m.z<-round(feat$row.m.z, 4) 
feat$row.retention.time<-round(feat$row.retention.time, 2) 
mz_table<-as.data.frame(  table(feat$row.m.z)  )
potential_oversplits<-mz_table[mz_table$Freq >= oversplit_num,]
potential_mz_forremov<-as.vector(potential_oversplits$Var1)

CIs_to_drop<-numeric()
for (i in 1:length(potential_mz_forremov)) {
  mz<-subset(feat, feat$row.m.z<=as.numeric(potential_mz_forremov[i])+mz_tol & feat$row.m.z>=as.numeric(potential_mz_forremov[i])-mz_tol)
  rt_sorted<-mz[order(mz$row.retention.time),]
  for (j in 2:nrow(rt_sorted)){
    dif<-rt_sorted$row.retention.time[j]-rt_sorted$row.retention.time[j-1]
    if (dif < rt_tol){
      CIs_to_drop<-c(CIs_to_drop, rt_sorted$row.ID[j])}}}
CIs_to_drop<-unique(CIs_to_drop)
cleaned_data<- subset(feat, !(row.ID %in% CIs_to_drop))
cleaned_data <- cleaned_data[, names(cleaned_data) != "X"]

length(CIs_to_drop) 
write.csv(cleaned_data, "quant_OSremoved.csv", row.names = FALSE) 

# BLANK REMOVAL, by fold changes -----------------------------------------------------------
tab<-read.csv("quant_OSremoved.csv", header=TRUE) 
l<-5 
blanknames<-"blank|Blank|H9|H10|H11|H12" #H9-12 contained the cheerio media extract

tab <- tab[, names(tab) != "X"]
fcount <- c()
for(i in 1:nrow(tab)) {
  bmax <- max(tab[i, grep(blanknames, colnames(tab))])
  if(bmax==0) {
    fcount <- c(fcount, 1 )} else {
    fcount <- c(fcount, sum(tab[i, -grep(blanknames, colnames(tab))][,-1|-2|-3] > l*bmax) ) }}
num_removed<-nrow(tab)-nrow(tab[-which(fcount==0),])

cat("number of features removed:", num_removed, "\n") 
write.csv(tab[-which(fcount==0),], "quant_OSremov_5blankremoved.csv", row.names=FALSE)

# TIC NORM ----------------------------------------------------------------
table<-read.csv("quant_OSremov_5blankremoved.csv")

norm_table<- apply(table[, -(1:3)], 2, function(x) x / sum(x) )
norm_table<-cbind(table[,1:3], norm_table)
norm_table <- norm_table[, names(norm_table) != "X"]                   

write.csv(norm_table, "quant_OSremov_5blankremoved_TICnorm.csv", row.names = FALSE) #name what you want but keep .csv
# FORMATTING --------------------------------------------------------
table<-read.csv("quant_OSremov_5blankremoved_TICnorm.csv")

table$row.m.z<-round(table$row.m.z, 4) 
table$row.retention.time<-round(table$row.retention.time, 2) 
sampleid <-numeric()
for(i in 1:nrow(table)){
  id<-c(table$row.m.z[i], table$row.retention.time[i], table$row.ID[i])
  nums<-c("X", paste(id, collapse = "_"))
  sampleid[i]<-paste(nums, collapse = "")
}
tab.w.id<-cbind(sampleid, table[, -(1:3)])
names(tab.w.id) <- sub(".mzML.Peak.area|.mzXML.Peak.area", "", names(tab.w.id))
trans<-t(tab.w.id)
rownames(trans)[1] <- "sampleid"
trans<-cbind(row.names(trans), trans)
colnames(trans)<-trans[1,]
trans<-trans[-1,]
trans[is.na(trans)] <- 0    
trans <- trans[rownames(trans) != "X", ]

write.csv(trans, "quant_OSremov_5blankremoved_TICnorm_formatted.csv", row.names = FALSE, na="NA")

