#Load the following packages
library(ggplot2)
library(readxl)
library(vegan)
library(adespatial)
library(ade4)
library(ggpubr)
library(readr)
library(tidyverse)
library(dplyr)
library(stats)
library(scales)

#import datasets
OTU_V_final <- read_csv("OTU_V_final.csv")


#import color scheme for each bacteria
Species_colors<- read_xlsx("Species_colors.xlsx")
Species_col<- as.data.frame(Species_colors$Species, stringsAsFactors = TRUE)
Species_color_scheme <- Species_colors$Color_hex

#assign levels list as names of colors list
names(Species_color_scheme) <- levels(Species_col$`Species_colors$Species`) 


#Create the final figure
dataset <- OTU_V_final %>% group_by(species, function_name, Ct_status, category) %>% 
                            summarize(average_gene = mean(count_of_species)) %>% 
                            filter(!category == "Poorly Characterized", !function_name == "Nuclear structure")

dataset %>% filter(average_gene < 10) %>% mutate(species = "Other")

dataset <- rbind( (dataset %>% filter(average_gene >=10)) , (dataset %>% filter(average_gene < 10) %>% mutate(species = "Other")) )

ggplot(dataset, aes(y=function_name, x= average_gene, fill=species)) + 
  geom_bar(position="stack", stat="identity", show.legend = TRUE)+
  scale_fill_manual(values = Species_color_scheme, limits = force)+
  xlab("Average Gene Count")+
  ylab("")+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  facet_grid(rows = vars(category, Ct_status), scale = "free_y", space = "free_y")+
  guides(fill=guide_legend(ncol=1))

