## Sven Buerki 2023
## Edited by Paige Skinner 2023

###~~~
# CREATE NEWPLACE WORD LEXICON
###~~~
#set path to project directory 
path <- ("path/to/project/directory/")
##~~
# Open document with SITE pubs information (from NEWPLACE)
#NEWPLACE link -> https://www.dna.affrc.go.jp/PLACE/?action=newplace
##~~
place_seq <- readLines("Keyword_Lexicon/place_seq.txt")

##~~
# Find and extract keywords
##~~
# Keywords are identified by "KW" @ line start
KeyWords <- strsplit(paste(gsub(";", "|", gsub(" ", "", gsub("KW   ", "", place_seq[grep("KW   ", place_seq)]))), collapse = "|"), split = "[|]")

##~~
# Pivot table
##~~
# Produce a pivot table with word and frequency in NEWPLACE
pivKeyWords <- sort(table(KeyWords), decreasing = T)
#Remove data ""
pivKeyWords <- pivKeyWords[-which(names(pivKeyWords) == "")]

#wordCount <- data.frame(word = names(pivKeyWords), freq = as.vector(pivKeyWords))
#wordcloud2::wordcloud2(data=wordCount, size=1, color='random-dark')

##~~
# Produce lexicon
##~~

NEWPLACELexicon <- data.frame(NEWPLACE_word = names(pivKeyWords), 
                              NEWPLACE_word_freq = as.vector(pivKeyWords), 
                              ACCEPTED_word = names(pivKeyWords),
                              Category = rep("NA", length(pivKeyWords)),
                              Types = rep("NA", length(pivKeyWords)))
write.csv(NEWPLACELexicon, file = "Keyword_Lexicon/NEWPLACE_lexicon.csv", row.names = F)


##### Refining Lexicon #####

NEWPLACE_lexicon <- read.csv(paste0(path, 'NEWPLACE_lexicon.csv'), header = T)

#Create file with unique accepted words and export for editing

Red_NEWPLACE_Lexicon <- NEWPLACE_lexicon[!duplicated(NEWPLACE_lexicon$ACCEPTED_word), -2]

write.csv(Red_NEWPLACE_Lexicon,paste0(path, 'Red_NEWPLACE_lexicon.csv'), row.names = F)


#import updated lexicon file 

Red_NEWPLACE_lexicon_Refined <- read.csv(paste0(path, 'Red_NEWPLACE_lexicon_Update_6_1.csv'), header = T)

###Update lexicon file with refined words
NEWPLACE_lexicon_V2 <- NEWPLACE_lexicon[,-8:-11]
NEWPLACE_lexicon_V2$ACCEPTED_word_refined <- NA
NEWPLACE_lexicon_V2$Category_refined <- NA
NEWPLACE_lexicon_V2$Types_refined <- NA

for(i in 1:nrow(NEWPLACE_lexicon_V2)){
  
 row <- match(NEWPLACE_lexicon_V2$ACCEPTED_word[i], Red_NEWPLACE_lexicon_Refined$ACCEPTED_word)
 
   
   NEWPLACE_lexicon_V2$ACCEPTED_word_refined[i] <- Red_NEWPLACE_lexicon_Refined$ACCEPTED_word_refined[row]
   
   NEWPLACE_lexicon_V2$Category_refined[i] <- Red_NEWPLACE_lexicon_Refined$Category[row]
    
  NEWPLACE_lexicon_V2$Types_refined[i] <- Red_NEWPLACE_lexicon_Refined$Types[row]
 }
#remove old names
NEWPLACE_lexicon_V2_final <- NEWPLACE_lexicon_V2[,-3:-5]  
#rename columns 
colnames(NEWPLACE_lexicon_V2_final) <- c("NEWPLACE_word","NEWPLACE_word_freq","Information","Source", "ACCEPTED_word","Category","Types")

#correct names so there are not misspelled repeats

unique(NEWPLACE_lexicon_V2_final$Category)
unique(NEWPLACE_lexicon_V2_final$Types)

for(k in 1:nrow(NEWPLACE_lexicon_V2_final)){
  
  if(NEWPLACE_lexicon_V2_final$Category[k] == "Where "){
    
    NEWPLACE_lexicon_V2_final$Category[k] <- paste("Where")}
  
  else if(NEWPLACE_lexicon_V2_final$Category[k] == "Other"){
    NEWPLACE_lexicon_V2_final$Category[k] <- paste("Others")}
  
  else if(NEWPLACE_lexicon_V2_final$Types[k] == "Hormone"){
    NEWPLACE_lexicon_V2_final$Types[k] <- paste("Hormones")}
  
  else if(NEWPLACE_lexicon_V2_final$ACCEPTED_word[k] =="#NAME?"){
    NEWPLACE_lexicon_V2_final$ACCEPTED_word[k] <- paste("-S")
    NEWPLACE_lexicon_V2_final$NEWPLACE_word[k] <- paste("-S")
    NEWPLACE_lexicon_V2_final$Category[k] <- paste("Others")
    NEWPLACE_lexicon_V2_final$Types[k] <- paste("Others")
  }
  }
    
 

write.csv(NEWPLACE_lexicon_V2_final,paste0(path, 'NEWPLACE_lexicon_V2.csv'), row.names = F)
