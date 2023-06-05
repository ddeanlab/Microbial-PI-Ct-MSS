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
#I used summary_abundance_allvag.csv file

CST_centroids_012920 <- read_csv('CST_centroids_012920.csv')

#import dataset to run function
dataset <- CST_centroids_012920 # centroids file

#new rows for each species
raw_abun_perc <- dataset %>% gather(species, count_of_sample, 2:200) # change the number here

#rank the species by read count
raw_abun_perc_25_pos <- raw_abun_perc %>% group_by(species) %>%
  summarize(avg_density = mean(count_of_sample)) %>%
  arrange(desc(avg_density)) %>%
  mutate(rank = row_number())%>% 
  dplyr::select(rank, species, avg_density) %>%
  ungroup() %>% 
  top_n(25) %>% 
  rename(avg_density_pos = avg_density)


testing_df <- CST_centroids_012920 %>% gather(species, per_sample_100, 2:200)  # centroids file, number

#filter OTU table for data on only the top 25 species data
testing_df_top25 <- testing_df %>% 
  filter(species %in% raw_abun_perc_25_pos$species)

#change all species not part of top 25 to be in "Others" category
testing_df_not25 <- testing_df %>% 
  filter(!species %in% raw_abun_perc_25_pos$species) %>% 
  mutate(species = "Others")
testing_df <- rbind(testing_df_top25,testing_df_not25)

#change species name format
#testing_df$species <- gsub("_", " ", testing_df$species)

#set sample ID order for plotting
sample_order <- unique(paste0(as.numeric(substr(CST_centroids_012920$sub_CST, 
                                                1,
                                                nchar(CST_centroids_012920$sub_CST)-1)), 
                              "V"))  # centroids file
#write.csv(raw_abun_perc_25_pos, 'top25_centroids.csv')

cutoff <- 0.01
cst_name <- "IV-B"
plot_cutoff <- 0.001
Species_colors <- read_csv("top25_centroids_colorcode.csv")
#select data from one specific centroid
`%notin%` <- Negate(`%in%`)
temp <- raw_abun_perc%>% filter(sub_CST == cst_name)
colnames(temp) <- c("sub_CST", "species", "per_sample_100")
#temp <- temp_raw[temp_raw$species != 'Others']
#not_temp <- temp_raw[temp_raw$species %notin% testing_df_top25$species]
#temp <- temp[order(temp$species),]
#less_than_cutoff <- temp %>% filter(per_sample_100 < cutoff)
#temp <- temp %>% filter(per_sample_100 >= cutoff)
not_temp <- temp[temp$species %notin% Species_colors$Species,]
temp <- temp[temp$species %in% Species_colors$Species,]
not_not_in <- temp %>% filter(per_sample_100 < plot_cutoff)
temp <- temp %>% filter(per_sample_100 >= plot_cutoff)
temp[nrow(temp)+1, ] <- list(cst_name,'Others', sum(not_not_in$per_sample_100) + sum(not_temp$per_sample_100))
temp <- temp[order(temp$species),]

#make sure this is in alphabetical order
#import Species_colors.csv file

#Species_colors$Species <- gsub("_", " ", Species_colors$Species)
#Species as factors table
Species_colors <- Species_colors[Species_colors$Species %in% temp$species, ]
Species_colors <- Species_colors[order(Species_colors$Species),]
Species_col<- as.data.frame(Species_colors$Species, stringsAsFactors = TRUE)

#colors list
Species_color_scheme <- Species_colors$Color_hex

#assign levels list as names of colors list
names(Species_color_scheme) <- levels(Species_col$`Species_colors$Species`) 

library(ggplot2)
library(ggrepel)
# Use this to plot the centroids pie chart
temp_label <- temp
temp_label$per_sample_100[temp_label$per_sample_100 < 0.01] <- 0

slices <- temp_label$per_sample_100
lbls <- temp_label$species
pct <- slices
pct <- round(slices*10000)/100 #add percents to labels
# ad % to labels
lbldf <- data_frame(lbls, pct)
pct <- paste0(pct, "%")
pct[pct == "0%"] <-""


df2 <- temp %>% 
  mutate(csum = rev(cumsum(rev(per_sample_100))), 
         pos = per_sample_100/2 + lead(csum, 1),
         pos = if_else(is.na(pos), per_sample_100/2, pos))
df2$per_sample_100[df2$per_sample_100 < cutoff] <- ""

labels <- temp$species
result <- ggplot(temp, aes(x="", y=per_sample_100, fill=species))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  geom_label_repel(data = df2,
                  aes(y = pos, label = pct, segment.size = 0.2),
                  size = 6, nudge_x = 1 ,nudge_y = 0.9, show.legend = FALSE, box.padding = 0.1) +
  scale_fill_manual(values = Species_color_scheme, limits = force) + theme_void()
result<- result + theme(legend.position = "none")
dir <- paste("CST\ centroids\ pie\ chart/",cst_name)
ggsave(paste(dir, ".png"), width = 20, height = 20, units = "cm")
