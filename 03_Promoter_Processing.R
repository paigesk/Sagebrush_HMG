#Paige Skinner and Sven Buerki 2024

#### Process the NewPLACE output files for each promoter -> Find the unique motifs, import the keywords and replace with curated lexicon keyword file #####
#NEWPLACE link -> https://www.dna.affrc.go.jp/PLACE/?action=newplace

#working directory - Set to project file containing Keyword_Lexicon, NewPLACE_Output, and Lexicon_Output files
## Keyword_Lexicon contains lexicon and place.seq files, NewPLACE_Output contains output files acquired from NEWPLACE database, Lexicon_Output is empty 
setwd("/path/to/your/project")


######################################## 1. Get keywords for each motif ########################################

#File Inputs
# Open document with SITE pubs information (from NEWPLACE)
place_seq <- readLines("Keyword_Lexicon/place_seq.txt")


# Open NEWPLACE lexicon file (From Supporting Table 4)
NEWPLACE_lexicon <- read.csv('Keyword_Lexicon/Supporting_Table_4.csv', header = T)


#List all the promoter sequence files
# - OUTPUT of NEWPLACE web search converted to csv format
prom_seq_files_ext <- list.files( path ="NEWPLACE_Output/", pattern = ".csv", full.names = F)

#remove file extension
prom_seq_files <- tools::file_path_sans_ext(prom_seq_files_ext) 

length(prom_seq_files)


#first loop opens the file and reformats the data 

for(i in 1:length(prom_seq_files)){
  #import file
  seq_mat <- read.csv(paste0( "NEWPLACE_Output/",prom_seq_files[i], ".csv"))
  #Add the length and location of the stop for each motif
  seq_mat$Length <- nchar(seq_mat$Signal_Sequence)
  seq_mat$Loc_Stop <- seq_mat$Loc_Start + seq_mat$Length
  #create a string of the start, stop and strand for each motif 
  seq_mat$Unique_Site <- paste(seq_mat$Loc_Start, seq_mat$Loc_Stop, seq_mat$Strand, sep = "|")
  
  # Create a list of unique sites
  unique_sites <- unique(seq_mat$Unique_Site)
  
  # Create reduced matrix with unique sites
  seq_red <- as.data.frame(matrix(ncol = 6, nrow = length(unique_sites)))
  colnames(seq_red) <- c("Unique_Site_ID", "NewPLACE_Signal_Sequence","SITE_Pubs","NEWPLACE_Words","ACCEPTED_word","NewPLACE_Motif_Element")
  # Assign unique sites ID
  seq_red$Unique_Site_ID <- unique_sites
  
#second loop gets information from imported data to reduced matrix with unique sites
  for(j in 1:nrow(seq_red)){
    #Subset seq_mat based on Unique_Site_ID
    tmp <- subset(seq_mat, seq_mat$Unique_Site == seq_red$Unique_Site_ID[j])
    # Get newPLACE sequence
    seq_red$NewPLACE_Signal_Sequence[j] <- tmp$Signal_Sequence[1]
    # Get Pubs
    seq_red$SITE_Pubs[j] <- paste(tmp$SITE_., collapse = "|")

  
  #get the keywords associated with each publication from the newplace database
  
    #Use the pub IDs to find the keywords
    tmp2 <- seq_red$SITE_Pubs[j]
    tmp2 <- strsplit(tmp2, split = "[|]")[[1]]
    
    #Finds where the pub ID is located, corresponds to AC in place_seq file
    toFind <- match(paste("AC   ", tmp2, sep = ""), place_seq)
    #Find, and extract unique words associated with publications
    LinesToExtract <- as.vector(sapply(toFind, function(x) grep("KW", place_seq[x:as.numeric(x+30)])))
    if(length(unlist(LinesToExtract)) > length(toFind)){
      if(class(LinesToExtract) == "list"){
        #print("list")
        #Update toFind to match LinesToExtract
        toFind <- rep(toFind, lengths(LinesToExtract))
        seq_red$NEWPLACE_Words[j] <- paste(unique(strsplit(paste(gsub(";", "|", gsub(" ", "", gsub("KW   ", "|", place_seq[toFind+unlist(LinesToExtract)-1]))), collapse = ""), split = "[|]")[[1]]), collapse = "|")
      }
      if(class(LinesToExtract) == "integer"){
        #Update toFind to match LinesToExtract
        toFind <- rep(toFind, length(LinesToExtract))
        seq_red$NEWPLACE_Words[j] <- paste(unique(strsplit(paste(gsub(";", "|", gsub(" ", "", gsub("KW   ", "|", place_seq[toFind+unlist(LinesToExtract)-1]))), collapse = ""), split = "[|]")[[1]]), collapse = "|")
      }
    }else{
      seq_red$NEWPLACE_Words[j] <- paste(unique(strsplit(paste(gsub(";", "|", gsub(" ", "", gsub("KW   ", "|", place_seq[toFind+LinesToExtract-1]))), collapse = ""), split = "[|]")[[1]]), collapse = "|")
    }
    #}
    
    
    #Replace NewPLACE Words with Accepted words from curated lexicon file
  
      
      #split individual Newplace keywords
      tmp3 <- seq_red$NEWPLACE_Words[j]
      tmp4 <- strsplit(tmp3, split = "[|]")[[1]]
      
      #match keywords to lexicon file
      rows <- match(tmp4, NEWPLACE_lexicon$NEWPLACE_word)
      
      accepted_word <- NEWPLACE_lexicon$ACCEPTED_word[rows]
      
      #find unique keywords
      unique <- unique(accepted_word)
      
      #put into the reducex matrix
      seq_red$ACCEPTED_word[j] <- paste(unique, collapse = "|")

    
   
    #Convert seq to lower case
   
      seq_red$NewPLACE_Motif_Element[j] <- tolower(seq_red$NewPLACE_Signal_Sequence[j])
  }
    
    #output the file
    write.csv(seq_red, paste0("Lexicon_Output/",prom_seq_files[i],'_KW.csv'), row.names = F)}
    
    
    
    ##### 2. Find the frequency of each word from lexicon and each motif within each gene and output merged file #####
 
#Make a newplace lexicon that is reduced to only the accepted words and categories
 
  Red_NEWPLACE_Lexicon <- NEWPLACE_lexicon[!duplicated(NEWPLACE_lexicon$ACCEPTED_word), -2]


  
  #Create a merged file with all the words and categories from the files we just made 
lex_files_ext <- list.files( path ="Lexicon_Output/", pattern = ".csv", full.names = F)

#remove file extension
lex_seq_files <- tools::file_path_sans_ext(lex_files_ext) 

length(lex_seq_files)

#Merge files to make frequency table

#Merge files where each line is one keyword

#create an empty objects to enter data into within the loop
all <- NULL
freq <- NULL

for(l in 1:length(lex_seq_files)){
  
  lex <- read.csv(paste0( "Lexicon_Output/",lex_seq_files[l], ".csv"))
  
  #Add column showing the protein ID and copy from the file name   
  
  string <- paste(lex_seq_files[l])
  string1 <- substring(string, 1,9)
  
  lex$ProtID <- string1
  
  #Merge files to make frequency table
  freq <- rbind(lex,freq)
  
  #Merge files where each line is one keyword
  for(m in 1:nrow(lex)){
    
    tmp5 <- unlist(strsplit(lex$ACCEPTED_word[m], split = "[|]"))
    #get rid of first position that is empty
    df <- as.data.frame(matrix(tmp5[-1]))
    #add the motif name associated with the keywords
    df$Motif <- lex$NewPLACE_Motif_Element[m]
    df$ProtID <- lex$ProtID[m]
    
    
    for(h in 1:nrow(df)){
      
      #find the row that matches the lexicon word
      rows<- match(df$V1[h],Red_NEWPLACE_Lexicon$ACCEPTED_word)
      
      
      #add the category and type from lexicon
      df$Category[h] <- Red_NEWPLACE_Lexicon$Category[rows]
      df$Type[h] <- Red_NEWPLACE_Lexicon$Types[rows]
    na.omit(df)
    
    }
    #merge
    all <- rbind(na.omit(df),all)
    }
}

head(all)
nrow(all)
head(freq)
#Manually changed the names from protein ID in sagebrush and locus ID in A. thaliana to gene name
all$Gene <- NA

colnames(all)[1] <- "Word"

for(j in 1:nrow(all)){
  if(all$ProtID[j] == "AT1G04880"){ all$Gene[j] <- "HMGB15"}
  else if(all$ProtID[j] == "AT1G20693"){ all$Gene[j] <- "HMGB2"}
  else if(all$ProtID[j] == "AT1G20696"){ all$Gene[j] <- "HMGB3"}
  else if(all$ProtID[j] == "AT1G55650"){ all$Gene[j] <- "HMGB11"}
  else if(all$ProtID[j] == "AT1G76110"){ all$Gene[j] <- "HMGB9"}
  else if(all$ProtID[j] == "AT2G17560"){ all$Gene[j] <- "HMGB4"}
  else if(all$ProtID[j] == "AT2G34450"){ all$Gene[j] <- "HMGB14"}
  else if(all$ProtID[j] == "AT3G13350"){ all$Gene[j] <- "HMGB10"}
  else if(all$ProtID[j] == "AT3G28730"){ all$Gene[j] <- "SSRP1"}
  else if(all$ProtID[j] == "AT3G51880"){ all$Gene[j] <- "HMGB1"}
  else if(all$ProtID[j] == "AT4G11080"){ all$Gene[j] <- "HMGB13"}
  else if(all$ProtID[j] == "AT4G23800"){ all$Gene[j] <- "HMGB6"}
  else if(all$ProtID[j] == "AT4G35570"){ all$Gene[j] <- "HMGB5"}
  else if(all$ProtID[j] == "AT5G23405"){ all$Gene[j] <- "HMGB12"}
  else if(all$ProtID[j] == "AT5G23420"){ all$Gene[j] <- "HMGB7"}
  else if(all$ProtID[j] == "ANN37363_"){ all$Gene[j] <- "ArtHMGB1a"}
  else if(all$ProtID[j] == "ANN37364_"){ all$Gene[j] <- "ArtHMGB1b"}
  else if(all$ProtID[j] == "ANN36225_"){ all$Gene[j] <- "ArtHMGB7"}
  else if(all$ProtID[j] == "ANN32443_"){ all$Gene[j] <- "ArtSSRP1"}
  else if(all$ProtID[j] == "ANN34248_"){ all$Gene[j] <- "ArtHMGB11"}
  else if(all$ProtID[j] == "ANN10184_"){ all$Gene[j] <- "ArtHMGB15"}
  else if(all$ProtID[j] == "ANN36958_"){ all$Gene[j] <- "ArtHMGB13a"}
  else if(all$ProtID[j] == "ANN36988_"){ all$Gene[j] <- "ArtHMGB13b"}
  else if(all$ProtID[j] == "ANN00312_"){ all$Gene[j] <- "ArtHMGB10"}
}

head(all)
tail(all)
#output merged file that contains all information 

write.csv(all, "Merged_motifs_words_cats.csv", row.names = F)

#type2 was manually edited to separate abiotic_stimuli into specific categories 

#type3 was manually edited to separate stress into specific categories

##See Supporting Table 6 for Type2 and Type3 categories
####Copy in information from supporting table 6
adj <- read.csv("Supporting_File_6.csv", header =T)

all$Type2 <- NA
all$Type3 <- NA

for(b in 1:nrow(all)){
  
  row <- match(all$Word[b], adj$Word)
  
  all$Type2[b] <- adj$Type2[row]
  all$Type3[b] <- adj$Type3[row]
  
}

write.csv(all,"Merged_motifs_words_cats_adj.csv", row.names = F)


#manually change gene names to export a frequency table that will be used in random forest analysis
freq$Gene <- NA

for(j in 1:nrow(freq)){
  if(freq$ProtID[j] == "AT1G04880"){ freq$Gene[j] <- "HMGB15"}
  else if(freq$ProtID[j] == "AT1G20693"){ freq$Gene[j] <- "HMGB2"}
  else if(freq$ProtID[j] == "AT1G20696"){ freq$Gene[j] <- "HMGB3"}
  else if(freq$ProtID[j] == "AT1G55650"){ freq$Gene[j] <- "HMGB11"}
  else if(freq$ProtID[j] == "AT1G76110"){ freq$Gene[j] <- "HMGB9"}
  else if(freq$ProtID[j] == "AT2G17560"){ freq$Gene[j] <- "HMGB4"}
  else if(freq$ProtID[j] == "AT2G34450"){ freq$Gene[j] <- "HMGB14"}
  else if(freq$ProtID[j] == "AT3G13350"){ freq$Gene[j] <- "HMGB10"}
  else if(freq$ProtID[j] == "AT3G28730"){ freq$Gene[j] <- "SSRP1"}
  else if(freq$ProtID[j] == "AT3G51880"){ freq$Gene[j] <- "HMGB1"}
  else if(freq$ProtID[j] == "AT4G11080"){ freq$Gene[j] <- "HMGB13"}
  else if(freq$ProtID[j] == "AT4G23800"){ freq$Gene[j] <- "HMGB6"}
  else if(freq$ProtID[j] == "AT4G35570"){ freq$Gene[j] <- "HMGB5"}
  else if(freq$ProtID[j] == "AT5G23405"){ freq$Gene[j] <- "HMGB12"}
  else if(freq$ProtID[j] == "AT5G23420"){ freq$Gene[j] <- "HMGB7"}
  else if(freq$ProtID[j] == "ANN37363_"){ freq$Gene[j] <- "ArtHMGB1a"}
  else if(freq$ProtID[j] == "ANN37364_"){ freq$Gene[j] <- "ArtHMGB1b"}
  else if(freq$ProtID[j] == "ANN36225_"){ freq$Gene[j] <- "ArtHMGB7"}
  else if(freq$ProtID[j] == "ANN32443_"){ freq$Gene[j] <- "ArtSSRP1"}
  else if(freq$ProtID[j] == "ANN34248_"){ freq$Gene[j] <- "ArtHMGB11"}
  else if(freq$ProtID[j] == "ANN10184_"){ freq$Gene[j] <- "ArtHMGB15"}
  else if(freq$ProtID[j] == "ANN36958_"){ freq$Gene[j] <- "ArtHMGB13a"}
  else if(freq$ProtID[j] == "ANN36988_"){ freq$Gene[j] <- "ArtHMGB13b"}
  else if(freq$ProtID[j] == "ANN00312_"){ freq$Gene[j] <- "ArtHMGB10"}
}
head(freq)
#Create a frequency table
##turn into a frequency table showing the number of occurrences of each motif in each gene promoter
table <- table(freq$Gene, freq$NewPLACE_Motif_Element)
head(table)

# Convert the contingency table to a data frame
df_freq <- as.data.frame.matrix(table)
head(df_freq)

#save
write.csv(df_freq, "Frequency_table_motifs.csv", row.names = T)
