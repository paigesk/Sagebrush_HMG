# Paige Skinner 2024
### Extract the Promoter and Gene sequences for each HMG gene ####

#load the necessary packages 
library(seqinr)
library(Biostrings)


## Run the code for Arabidopsis and sagebrush -- the genomes must be saved in separate folders ##

#set path to project directory 
path <- ("path/to/project/directory/")

######## Arabidopsis ###########

#read in Gene Dataframe from ncbi summary that contains the accession number for the chromosome, start, stop and strand infromation
#with column added for the gene names (eg. HMGB1)
TargetGPs <- read.csv(paste0(path, "Genome/Arabidopsis_Genome/ncbi_Arath_HMG.csv"), header = T)
head(TargetGPs)

#Add column for promoter and gene sequence
TargetGPs$PromoterSeq <- NA
TargetGPs$GeneSeq <- NA

#Separate multifasta file into individual oneliners for each chromosome 
fasta <- readDNAStringSet(paste0(path,"Genome/Arabidopsis_Genome/Arath_Genome.fasta"),format="fasta") 

for(i in 1:5){
  write.fasta(as.character(fasta[i]),fasta@ranges@NAMES[i], paste0(path,"Genome/Arabidopsis_Genome/Scaffolds/",substring(fasta@ranges@NAMES[i],1,10),".fasta"), open = "w", nbchar = 10000000000, as.string = FALSE)
}




#####Loop to Mine Gene and Promoter Sequences####

# Create a vector w/ scaffold IDs

ScaffID <- unique(TargetGPs$Accession)

#First loop acts on scaffolds
for(i in 1:length(ScaffID)){
  #ID in TargetGPs rows with target scaffoldID
  row2edit <- which(TargetGPs$Accession == ScaffID[i])
  #Open FASTA file of target scaffold
  print("Opens scaffold sequence, be patient!")
  ScaffSeq <- seqinr::read.fasta(file = paste0(path, "Genome/Arabidopsis_Genome/Scaffolds/", ScaffID[i], ".fasta"), as.string = TRUE)
  #Second loop gets seq data (gene and promoter) for all rows assigned to scaffID in TargetGPs 
  pb <- txtProgressBar(min = 0, max = length(row2edit), style = 3)
  for(j in 1:length(row2edit)){
    #First check strand
    if(TargetGPs$Orientation[row2edit[j]] == "plus"){
      #Promoter DNA sequence
      TargetGPs$PromoterSeq[row2edit[j]] <- substr(ScaffSeq, (TargetGPs$Begin[row2edit[j]]-1500), (TargetGPs$Begin[row2edit[j]]-1)) 
      TargetGPs$GeneSeq[row2edit[j]] <- substr(ScaffSeq, (TargetGPs$Begin[row2edit[j]]), (TargetGPs$End[row2edit[j]]))
      }else{
      #code for - strand
      #Get Promoter DNA sequence
      tmp2 <- substr(ScaffSeq, (TargetGPs$End[row2edit[j]]+1), (TargetGPs$End[row2edit[j]]+1500)) 
      #get reverse compliment
      tmp2_comp <- rev(comp(s2c(tmp2)))
      #paste a string of the reverse compliment seq into the targetGPs file
      TargetGPs$PromoterSeq[row2edit[j]] <-  paste(tmp2_comp, collapse = "")
      
      #Get gene DNA sequence
      tmp1 <- substr(ScaffSeq, (TargetGPs$Begin[row2edit[j]]), (TargetGPs$End[row2edit[j]])) 
      #get reverse compliment
      tmp1_comp <- rev(comp(s2c(tmp1)))
      
      #paste a string of the reverse compliment seq into the targetGPs file
      TargetGPs$GeneSeq[row2edit[j]] <-  paste(tmp1_comp, collapse = "")
    }
    #update pb
    setTxtProgressBar(pb, j)
  }
  close(pb)
}

head(TargetGPs)
write.csv(TargetGPs,paste0(path, "Sequences/Arath_HMG_promoter_geneseq.csv"), row.names = F)

#write promoter sequences into a fasta

#make a dataframe of promoter sequences
Prom_seq <- as.data.frame(cbind(paste0(TargetGPs$Locus.tag ,"_",TargetGPs$X),TargetGPs$PromoterSeq))

#export to a fasta
write.fasta(as.list(Prom_seq$V2),Prom_seq$V1,paste0(path,"Sequences/Prom_seq_Arath.fasta"))

#make a dataframe of gene seq
gene_seq <- as.data.frame(cbind(paste0(TargetGPs$Locus.tag ,"_",TargetGPs$X),TargetGPs$GeneSeq))

#export to a fasta
write.fasta(as.list(gene_seq$V2),gene_seq$V1,paste0(path,"Sequences/Gene_seq_Arath_extracted.fasta"))



####### Sagebrush #######

#read in Gene Dataframe extracted from gff file that contains the gene and scaffold names for sagebrush in the IDT3 genome
TargetGPs <- read.csv(paste0(path, "Genome/IDT3_Genome/Sagebrush_HMG_domains.csv"), header = T)
head(TargetGPs)

##scaffolds are already separated out within the /Scaffolds folder

#Add column for promoter sequence
TargetGPs$PromoterSeq <- NA
TargetGPs$GeneSeq <- NA

#####Loop to Mine Gene and Promoter Sequences####

# Create a vector w/ scaffold IDs
# - It will serve to subset data and get sequences
ScaffID <- unique(TargetGPs$ScaffID)

#First loop acts on scaffolds

for(i in 1:length(ScaffID)){
  #ID in TargetGPs rows with target scaffoldID
  row2edit <- which(TargetGPs$ScaffID == ScaffID[i])
  #Open FASTA file of target scaffold
  print("Opens scaffold sequence, be patient!")
  ScaffSeq <- seqinr::read.fasta(file = paste0(path, "Genome/IDT3_Genome/Scaffolds/", ScaffID[i], ".fa"), as.string = TRUE)
  #Second loop gets seq data (gene and promoter) for all rows assigned to scaffID in TargetGPs 
  pb <- txtProgressBar(min = 0, max = length(row2edit), style = 3)
  for(j in 1:length(row2edit)){
    #First check strand
    if(TargetGPs$Strand[row2edit[j]] == "+"){
      #Promoter DNA sequence
      TargetGPs$PromoterSeq[row2edit[j]] <- substr(ScaffSeq, (TargetGPs$Start[row2edit[j]]-1500), (TargetGPs$Start[row2edit[j]]-1)) 
      TargetGPs$GeneSeq[row2edit[j]] <- substr(ScaffSeq, (TargetGPs$Start[row2edit[j]]), (TargetGPs$End[row2edit[j]]))
    }else{
      #code for - strand
      #Get Promoter DNA sequence
      tmp2 <- substr(ScaffSeq, (TargetGPs$End[row2edit[j]]+1), (TargetGPs$End[row2edit[j]]+1500)) 
      #get reverse compliment
      tmp2_comp <- rev(comp(s2c(tmp2)))
      #paste a string of the reverse compliment seq into the targetGPs file
      TargetGPs$PromoterSeq[row2edit[j]] <-  paste(tmp2_comp, collapse = "")
      
      #Get gene DNA sequence
      tmp1 <- substr(ScaffSeq, (TargetGPs$Start[row2edit[j]]), (TargetGPs$End[row2edit[j]])) 
      #get reverse compliment
      tmp1_comp <- rev(comp(s2c(tmp1)))
      
      #paste a string of the reverse compliment seq into the targetGPs file
      TargetGPs$GeneSeq[row2edit[j]] <-  paste(tmp1_comp, collapse = "")
    }
    #update pb
    setTxtProgressBar(pb, j)
  }
  close(pb)
}

write.csv(TargetGPs,paste0(path, "Sequences/Sagebrush_HMG_domains_promoter_geneseq.csv"), row.names = F)


#remove duplicate rows
Prom <- TargetGPs[!duplicated(TargetGPs$PromoterSeq),]

#export fasta files 
Prom_seq <- as.data.frame(cbind(Prom$ProtID,Prom$PromoterSeq))

write.fasta(as.list(Prom_seq$V2),Prom_seq$V1,paste0(path,"Sequences/Sagebrush_Prom_seq.fasta"))

gene_seq <- as.data.frame(cbind(Prom$ProtID,Prom$GeneSeq))

write.fasta(as.list(gene_seq$V2),gene_seq$V1,paste0(path,"Sequences/Sagebrush_Gene_seq.fasta"))

#### Now you can enter your fasta files of the promoter sequences into NewPLACE! ##















