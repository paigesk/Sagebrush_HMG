# Paige Skinner 2024
# Random Forest -- For each gene find which motifs differentiate Sagebrush and Arabidopsis using only the paralog genes

##Load packages
install.packages("vioplot")
library(varSelRF)
library(ranger)
library(vioplot)

setwd("/media/Artemisia/HMG_Protein/07_HMG_Promoter")
path <- "/media/Artemisia/HMG_Protein/07_HMG_Promoter"
setwd("/path/to/your/project")

#read in the frequency table showing the frequency of each motif within the promoters of genes you want to analyze
data <- read.csv("Frequency_table_motifs.csv",header = T)
head(data)

#add genes to row names
rownames(data) <- data[,1]

#remove column with genes 
data2 <- data[,-1]

#make object of all the gene names
genez <- rownames(data)

#create objects containing the gene names of the paralogs only 
Artemisiaz <- genez[grepl("Art",genez)]
Arabidopsisz <- c(gsub("Art","",Artemisiaz), "HMGB13", "HMGB1")

#Subset to only include paralogs from the two species
data_simp <- subset(data2, row.names(data2) %in% Artemisiaz |row.names(data2) %in% Arabidopsisz)
head(data_simp)

#Add column to assign the taxon
data_simp$Taxon <- NA
data_simp$Taxon <- row.names(data_simp)

#assign random numbers to indicate sagebrush or arath taxon which is in the last column
data_simp[!grepl("Art",Arabidopsisz),ncol(data_simp)] <- 1239458389457986304578
data_simp[grepl("Art",genez),ncol(data_simp)] <- 2131234123412342



#determine the number of features
n_features <- length(colnames(data_simp))-1

#run the random forest on only the motifs, not including column with taxon
message("Selecting variables by successive feature elimination and OOB error")
message("...using 50,000 trees at each iteration")
varselrf_comp <- varSelRF::varSelRF(data_simp[,-ncol(data_simp)],Class = as.factor(data_simp$Taxon),recompute.var.imp=F,whole.range=T,ntree=50000,ntreeIterat=50000, keep.forest = TRUE,verbose=TRUE,vars.drop.num=1,vars.drop.frac=NULL)
SelectedVars <- as.factor(varselrf_comp$selected.vars)
  sh <- varselrf_comp$selec.history
  message(sprintf("%s features were found as important",length(SelectedVars)))
  print(SelectedVars)
  #paralogs only: "ngatt"  "tatcca"
  #all: aatccaa ctctt gataa ggata ngatt tatcca
  data_simp$Taxon <- as.factor(data_simp$Taxon)
  

  NewRF <- ranger::ranger(as.formula(sprintf("Taxon ~ %s",paste(SelectedVars, collapse = " + "))), data = data_simp, num.trees = 50000 ,importance = "permutation")
  pz <- ranger::importance_pvalues(NewRF,method="altmann",num.permutations = 100,formula=as.formula(sprintf("Taxon ~ %s",paste(SelectedVars, collapse = " + "))),data = data_simp)
  Forest <- as.data.frame(pz)
  #p value says NA

  Forest[,c(2,3)] <- Forest[,c(1,2)]
  Forest[,1] <- row.names(Forest)
  row.names(Forest) <- NULL
  colnames(Forest) <- c("Variable","Importance","p-Value")

  #add column for wilcox test p value 
  Forest$Wilcox_pvalue <- NA
  
#run wilcox test and plot the selectedVars that were found

for(i in 1:length(SelectedVars)){
OneVar <- data_simp[which(data_simp$Taxon == 2131234123412342), which(colnames(data_simp) == SelectedVars[i])]
TwoVar <- data_simp[which(data_simp$Taxon == 1239458389457986304578), which(colnames(data_simp) == SelectedVars[i])]
 Wilcox1 <- wilcox.test(as.numeric(OneVar),as.numeric(TwoVar))
row <- which(Forest$Variable == SelectedVars[i])
   Forest$Wilcox_pvalue[row] <- Wilcox1$p.value
   

 #pvalue is 0.0009894 (paralog)

 pdf(paste0(path,"/Random_Forest/sagebrushvarath_vioplot" ,i,".pdf"))
  vioplot::vioplot(OneVar, side = "left", plotCentre = "line", ylim=c(.8*min(OneVar,TwoVar),max(OneVar,TwoVar)+.2*max(OneVar,TwoVar)),col = "red", xaxt= "n",main=SelectedVars[i])# = TRUE)
 vioplot::vioplot(TwoVar, side = "right", plotCentre = "line", col = "blue", add = TRUE )
 legend("topright", legend = c("Sagebrush", "A. thaliana"), fill = c("red","blue"))
dev.off()
}
  
  head(Forest)
 
#write variable importance and p-value table
  write.csv(Forest,"Random_Forest/Forest_sagebrushvarath_paraloggenes_res.csv",row.names = F)
  
  
  
  
  
  
  
  
  OneVar <- data_simp[which(data_simp$Taxon == 2131234123412342), which(colnames(data_simp) == SelectedVars[i])]
TwoVar <- data_simp[which(data_simp$Taxon == 1239458389457986304578), which(colnames(data_simp) == SelectedVars[i])]
 Wilcox2 <- wilcox.test(as.numeric(OneVar),as.numeric(TwoVar))
 #p value is 0.002116 (paralog)
 
 #All genes p values: 
 #aatccaa(0.001888) ctctt(0.0003809) gataa(0.02081) ggata(0.008091) ngatt(0.0002597) tatcca(0.0002597)
 
 

 pdf(paste0("Random_Forest/Cluster/sagebrushvarath_allgenes_vioplot2.pdf"))
 vioplot::vioplot(OneVar, side = "left", plotCentre = "line", ylim=c(.8*min(OneVar,TwoVar),max(OneVar,TwoVar)+.2*max(OneVar,TwoVar)),col = "red", xaxt= "n",ylab=SelectedVars[2])# = TRUE)
 vioplot::vioplot(TwoVar, side = "right", plotCentre = "line", col = "blue", add = TRUE )
dev.off()















