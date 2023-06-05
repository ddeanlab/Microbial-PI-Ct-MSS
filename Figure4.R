library(ggplot2)
library(readxl)
library(vegan)
library(adespatial)
library(ade4)
library(ggpubr)
library(readr)
library(tidyverse)
library(dplyr)

#import datasets
OTU_pos_neg <- read_xlsx("Final_OTU_pos_to_neg.xlsx", sheet = "Pairwise")
OTU_pos_pos <- read_xlsx("Final_OTU_pos_pos.xlsx", sheet = "Pairwise")
OTU_neg_neg <- read_xlsx("Final_OTU_neg_neg.xlsx", sheet = "Pairwise")
OTU_neg_neg <- OTU_neg_neg %>% filter(TimePoint != "Second Followup")

#import dataset to run function
dataset <- OTU_pos_neg

#new rows for each species
raw_abun_perc <- dataset %>% gather(species, count_of_sample, 12:281)

#rank the species by read count
raw_abun_perc_25_pos <- raw_abun_perc %>% group_by(species) %>%
  summarize(avg_density = mean(count_of_sample)) %>%
  arrange(desc(avg_density)) %>%
  mutate(rank = row_number())%>% 
  dplyr::select(rank, species, avg_density) %>%
  ungroup() %>% 
  top_n(25) %>% 
  rename(avg_density_pos = avg_density)


testing_df <- OTU_pos_neg %>% gather(species, per_sample_100, 12:281)

#filter OTU table for data on only the top 25 species data
testing_df_top25 <- testing_df %>% 
  filter(species %in% raw_abun_perc_25_pos$species)

#change all species not part of top 25 to be in "Others" category
testing_df_not25 <- testing_df %>% 
  filter(!species %in% raw_abun_perc_25_pos$species) %>% 
  mutate(species = "Others")
testing_df <- rbind(testing_df_top25,testing_df_not25)

#change species name format
testing_df$species <- gsub("_", " ", testing_df$species)

#set sample ID order for plotting
sample_order <- unique(paste0(as.numeric(substr(OTU_pos_neg$Sample, 
                                                1,
                                                nchar(OTU_pos_neg$Sample)-1)), 
                              "V"))

#make sure this is in alphabetical order
Species_colors<- read_xlsx("Species_colors.xlsx")
Species_colors$Species <- gsub("_", " ", Species_colors$Species)
#Species as factors table
Species_col<- as.data.frame(Species_colors$Species, stringsAsFactors = TRUE)

#colors list
Species_color_scheme <- Species_colors$Color_hex

#assign levels list as names of colors list
names(Species_color_scheme) <- levels(Species_col$`Species_colors$Species`) 

#plot microbiome for each patient at baseline and followup
ggplot(testing_df , aes(x = TimePoint, y = per_sample_100, fill = species))+
  geom_col(position = "fill", show.legend = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = Species_color_scheme, limits = force)+
  theme_bw()+
  facet_wrap(~factor(Sample,levels=sample_order), ncol = 15)+
  theme(axis.title= element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())
