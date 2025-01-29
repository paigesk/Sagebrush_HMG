This repository contains the R scripts used for the analyses in the paper:

"An in silico Approach to Identify Candidate Genes in a Stack of Genes: Do High Mobility Group Proteins Contribute to Phenotypic Plasticity in a Keystone Shrub?"

Authors: Paige Skinner, John M. A. Wojahn, Adam Renfrow, Sven Buerki

Journal: 

DOI: 

The analyses include extraction of promoter sequences 1500 bps upstream of the start codon in each gene of interest, analyses of promoter function through categorization to compare the model species (Arabidopsis thaliania) and the non-model species of interest (Artemisia tridentata). This repository is designed to ensure the reproducibility of all results presented in the paper.


## Prerequisites

This code was developed using:
- R version:  R version 4.3.2 (2023-10-31)
- R libraries:
  - ggplot2: Plotting.
  - pheatmap: Plotting.
  - seqinr: Writing fasta files.
  - Biostrings: Reading fasta files.


Install the required libraries using:
R
install.packages(c( "ggplot2”,”pheatmap”,”seqinr”,”Biostrings”))

### Data

- Data from publication is stored in zenudo: 10.5281/zenodo.14218603


- Directory structure required for code (files listed must be added prior to running code):
	- Project/
		- Genome/ : Contains genomes of interest.
			- Arabipdosis_Genome/: fasta file of full Arabidopsis thaliana genome.
				- ncbi_Arath_HMG.csv : gff file that has been subset to include information for only the genes of interest. 
				- Scaffolds/: fasta file is separated into individual chromosomes and output here
				- IDT3_Genome/
				- Sagebrush_HMG_domains.csv : gff file that has been subset to include information for only the genes of interest. 
				- Scaffolds/: fasta file of each individual scaffold in the sagebrush (Artemisia tridentata subsp. tridentata IDT3 genome)
		- Heatmaps/: Output of heatmap plots are placed here.
		- Keyword_Lexicon/: 
			- place_seq.txt: downloaded from the NewPLACE database
			- supporting_table_4.csv: Curated lexicon file found in doi: 
		- NEWPLACE_Output/: Contains the output files from the NewPLACE database after entering the promoter sequences into the database and converting into csv files with the following column names:  	
							Factor_or_Site_Name, Loc_Start, Strand, Signal_Sequence, SITE_#
		- Sequences/: Promoter and gene sequences extracted from the scaffold fasta files will be output here. 



## License

This repository is licensed under the GNU GPL License. See LICENSE.txt for details.

## Citation

If you use this code or data, please cite the following paper:

[Author(s)]. ([Year]). *[Title of the Paper]*. [Journal Name]. DOI: [Link]

## Contact

For questions or feedback, please contact:
- Paige Skinner (mailto:paigeskinner@u.boisestate.edu)

## Example repository: bolded are outputs

           °--07_HMG_Promoter                                    
               ¦--Genome                                         
               ¦   ¦--Arabidopsis_Genome                         
               ¦   ¦   ¦--Arath_Genome.fasta                     
               ¦   ¦   ¦--ncbi_Arath_HMG.csv                     
               ¦   ¦   °--Scaffolds                              
               ¦   ¦       ¦--CP002684.1.fasta                   
               ¦   ¦       ¦--CP002685.1.fasta                   
               ¦   ¦       ¦--CP002686.1.fasta                   
               ¦   ¦       ¦--CP002687.1.fasta                   
               ¦   ¦       °--CP002688.1.fasta                   
               ¦   °--IDT3_Genome                                
               ¦       ¦--Sagebrush_HMG_domains.csv              
               ¦       °--Scaffolds                              
               ¦           ¦--Scaffold_1.fa                      
               ¦           ¦--Scaffold_2.fa                      
               ¦           ¦--Scaffold_3.fa                      
               ¦           ¦--Scaffold_4.fa                      
               ¦           ¦--Scaffold_5.fa                      
               ¦           ¦--Scaffold_6.fa                      
               ¦           ¦--Scaffold_7.fa                      
               ¦           ¦--Scaffold_8.fa                      
               ¦           °--Scaffold_9.fa                      
               ¦--Heatmaps                                       
               ¦   ¦--General_stress.pdf                         
               ¦   °--Specific_stress.pdf                        
               ¦--Keyword_Lexicon                                
               ¦   ¦--NEWPLACE_lexicon_V3.csv                    
               ¦   °--place_seq.txt                              
               ¦--Lexicon_Output                                 
               ¦   ¦--ANN00312_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN10184_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN32443_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN34248_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN36225_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN36958_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN36988_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN37363_IDT3_PromOutput_KW.csv            
               ¦   ¦--ANN37364_IDT3_PromOutput_KW.csv            
               ¦   ¦--AT1G04880_Arath_PromOutput_KW.csv          
               ¦   ¦--AT1G20693_Arath_PromOutput_KW.csv          
               ¦   ¦--AT1G20696_Arath_PromOutput_KW.csv          
               ¦   ¦--AT1G55650_Arath_PromOutput_KW.csv          
               ¦   ¦--AT1G76110_Arath_PromOutput_KW.csv          
               ¦   ¦--AT2G17560_Arath_PromOutput_KW.csv          
               ¦   ¦--AT2G34450_Arath_PromOutput_KW.csv          
               ¦   ¦--AT3G13350_Arath_PromOutput_KW.csv          
               ¦   ¦--AT3G28730_Arath_PromOutput_KW.csv          
               ¦   ¦--AT3G51880_Arath_PromOutput_KW.csv          
               ¦   ¦--AT4G11080_Arath_PromOutput_KW.csv          
               ¦   ¦--AT4G23800_Arath_PromOutput_KW.csv          
               ¦   ¦--AT4G35570_Arath_PromOutput_KW.csv          
               ¦   ¦--AT5G23405_Arath_PromOutput_KW.csv          
               ¦   °--AT5G23420_Arath_PromOutput_KW.csv          
               ¦--Merged_motifs_words_cats_adj.csv               
               ¦--Merged_motifs_words_cats.csv
               ¦--Frequency_table_motifs.csv                   
               ¦--NEWPLACE_Output                                
               ¦   ¦--ANN00312_IDT3_PromOutput.csv               
               ¦   ¦--ANN10184_IDT3_PromOutput.csv                
               ¦   ¦--ANN32443_IDT3_PromOutput.csv               
               ¦   ¦--ANN34248_IDT3_PromOutput.csv               
               ¦   ¦--ANN36225_IDT3_PromOutput.csv               
               ¦   ¦--ANN36958_IDT3_PromOutput.csv               
               ¦   ¦--ANN36988_IDT3_PromOutput.csv               
               ¦   ¦--ANN37363_IDT3_PromOutput.csv               
               ¦   ¦--ANN37364_IDT3_PromOutput.csv               
               ¦   ¦--AT1G04880_Arath_PromOutput.csv             
               ¦   ¦--AT1G20693_Arath_PromOutput.csv             
               ¦   ¦--AT1G20696_Arath_PromOutput.csv             
               ¦   ¦--AT1G55650_Arath_PromOutput.csv             
               ¦   ¦--AT1G76110_Arath_PromOutput.csv             
               ¦   ¦--AT2G17560_Arath_PromOutput.csv             
               ¦   ¦--AT2G34450_Arath_PromOutput.csv             
               ¦   ¦--AT3G13350_Arath_PromOutput.csv             
               ¦   ¦--AT3G28730_Arath_PromOutput.csv             
               ¦   ¦--AT3G51880_Arath_PromOutput.csv             
               ¦   ¦--AT4G11080_Arath_PromOutput.csv             
               ¦   ¦--AT4G23800_Arath_PromOutput.csv             
               ¦   ¦--AT4G35570_Arath_PromOutput.csv             
               ¦   ¦--AT5G23405_Arath_PromOutput.csv             
               ¦   °--AT5G23420_Arath_PromOutput.csv                           
               °--Sequences                                      
                   ¦--Arath_HMG_promoter_geneseq.csv             
                   ¦--Gene_seq_Arath_extracted.fasta             
                   ¦--Prom_seq_Arath.fasta                       
                   ¦--Sagebrush_Gene_seq.fasta                   
                   ¦--Sagebrush_HMG_domains_promoter_geneseq.csv 
                   °--Sagebrush_Prom_seq.fasta  
