#Paige Skinner 2024

############# Create heat plots showing the number of occurrences of each keyword within the promoters ##########


#working directory - Set to project file
setwd("/path/to/your/project")

#load necessary packages
library(pheatmap)

#import data that has been adjusted with type2 and type3 added
Merged <- read.csv("Merged_motifs_words_cats_adj.csv", header = T)


#remove Trash, binidng site, others and TF from dataframe since these do not tell you information on the function
all_red <- subset(Merged, Merged$Type3 != "Trash" & Merged$Type3 != "Binding_site" & Merged$Type3 != "Others" & Merged$Type3 != "TF" & Merged$Type3 != "Flower"
                  & Merged$Type3!= "Leaf" & Merged$Type3 != "Shoot" & Merged$Type3 != "Root")


#Heatplot with only the specific stress categories
specific_stress <- subset(Merged, Merged$Type3 == "Drought" | Merged$Type3 == "Cold" | Merged$Type3 == "Oxygen" | Merged$Type3 == "Low-CO2" | Merged$Type3 == "Starvation"
                  | Merged$Type3 == "Heat" | Merged$Type3 == "Stress")

#create a frequency table to show the number of occurrences of each keyword
tab_s <- table(specific_stress$Type3,specific_stress$Gene)

# sort the genes and categories
order_ind <- c( "Drought", "Cold","Heat","Oxygen" ,"Low-CO2","Starvation", "Stress")
order_genes <- c("ArtHMGB1a", "ArtHMGB1b", "HMGB1", "ArtHMGB7", "HMGB7", "ArtHMGB10","HMGB10", "ArtHMGB11","HMGB11","ArtHMGB13a",
                 "ArtHMGB13b","HMGB13", "ArtHMGB15", "HMGB15","ArtSSRP1", "SSRP1", "HMGB2", "HMGB3", "HMGB4", "HMGB5", "HMGB6", "HMGB9", "HMGB12", "HMGB14" )

#apply sorting to table
tab_s2 <- tab_s[order_ind,order_genes]

#create heatmap 
pdf("Heatmaps/Specific_stress.pdf")
pheatmap(
  tab_s2, xlabs = "Keyword", ylabs = "Gene", main = "", fontsize_row=12,fontsize_col = 12, legend = T,
  cluster_rows = F,  # Disable row clustering
  cluster_cols = FALSE,  # Disable column clustering
  show_rownames = T,  # Do not show row names
  show_colnames = T,  # Do not show column names
  border_color = NA, # Remove border color
  cellheight=20, cellwidth = 15, 
  
)

dev.off()

#heatplot with general categories

general <- subset(Merged, Merged$Type3 == "Biotic_stimuli" | Merged$Type3 == "Hormones" | Merged$Type3 == "Growth" | Merged$Type3 == "Ontology" 
                          | Merged$Type3 == "Abiotic_stimuli" | Merged$Type3 == "Light")



#create a frequency table to show the number of occurences of each keyword
tab_g <- table( general$Type3,general$Gene)
# sort the genes and categories
order_ind <- c("Abiotic_stimuli","Light","Biotic_stimuli","Hormones","Growth","Ontology")
order_genes <- c("ArtHMGB1a", "ArtHMGB1b", "HMGB1", "ArtHMGB7", "HMGB7", "ArtHMGB10","HMGB10", "ArtHMGB11","HMGB11","ArtHMGB13a",
                 "ArtHMGB13b","HMGB13", "ArtHMGB15", "HMGB15","ArtSSRP1", "SSRP1", "HMGB2", "HMGB3", "HMGB4", "HMGB5", "HMGB6", "HMGB9", "HMGB12", "HMGB14" )

#apply sorting to table
tab_g2 <- tab_g[order_ind,order_genes]

#create heatmap 
pdf("Heatmaps/General_stress.pdf")
pheatmap(
  tab_g2, xlabs = "Keyword", ylabs = "Gene", main = "", fontsize_row=12,fontsize_col = 12, legend = T,
  cluster_rows = F,  # Disable row clustering
  cluster_cols = FALSE,  # Disable column clustering
  show_rownames = T,  # Do not show row names
  show_colnames = T,  # Do not show column names
  border_color = NA, # Remove border color
  cellheight=20, cellwidth = 15, 

  
)
dev.off()











