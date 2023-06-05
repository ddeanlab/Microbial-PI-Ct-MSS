library(ggplot2)
library(readxl)
library(vegan)
library(adespatial)
library(ade4)
library(ggpubr)
library(readr)
library(tidyverse)
library(dplyr)

data <- read_csv("New CST mean abundance.csv")
Species_colors <- read_csv("Species_colors.csv")


cutoff <- 0.01
# Replace here for different new subCST
cst_name <- "IV-E"
plot_cutoff <- 0.001
#select data from one specific centroid
`%notin%` <- Negate(`%in%`)
selected <- select(data, Species, cst_name)
temp <- selected[selected$Species %in% Species_colors$Species,]
not_temp <- selected[selected$Species %notin% Species_colors$Species,]

#temp <- temp_raw[temp_raw$species != 'Others']
#not_temp <- temp_raw[temp_raw$species %notin% testing_df_top25$species]
#temp <- temp[order(temp$species),]
#temp <- temp %>% filter(per_sample_100 >= cutoff)
temp[nrow(temp)+1, ] <- list('Others', sum(not_temp$`IV-E`))
#less_than_cutoff <- temp %>% filter(`IV-E` < plot_cutoff)
#other <-sum(less_than_cutoff$`IV-E`) + temp$`IV-E`[temp$Species == "Others"]
#temp <- temp %>% filter(`IV-E` >= plot_cutoff)
#temp[nrow(temp)+1, ] <- list('Others', other)
#temp$`IV-E`[temp$Species == "Others"]<- other
temp <- temp[order(temp$Species),]

#make sure this is in alphabetical order
#import Species_colors.csv file

#Species_colors$Species <- gsub("_", " ", Species_colors$Species)
#Species as factors table
Species_colors <- Species_colors[Species_colors$Species %in% temp$Species, ]
Species_colors <- Species_colors[order(Species_colors$Species),]
Species_col<- as.data.frame(Species_colors$Species, stringsAsFactors = TRUE)

#colors list
Species_color_scheme <- Species_colors$Color_hex

#assign levels list as names of colors list
names(Species_color_scheme) <- levels(Species_col$`Species_colors$Species`) 

library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra) 
# Use this to plot the centroids pie chart
temp_label <- temp
temp_label$`IV-E`[temp_label$`IV-E`< 0.01] <- 0

slices <- temp_label$`IV-E`
lbls <- temp_label$Species
pct <- slices
pct <- round(slices*10000)/100 #add percents to labels
# ad % to labels
lbldf <- data_frame(lbls, pct)
pct <- paste0(pct, "%")
pct[pct == "0%"] <-""


df2 <- temp %>% 
  mutate(csum = rev(cumsum(rev(`IV-E`))), 
         pos = `IV-E`/2 + lead(csum, 1),
         pos = if_else(is.na(pos), `IV-E`/2, pos))
df2$`IV-E`[df2$`IV-E` < cutoff] <- ""

labels <- temp$Species
result <- ggplot(temp, aes(x="", y=`IV-E`, fill=Species))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  geom_label_repel(data = df2,
                  aes(y = pos, label = pct, segment.size = 0.2),
                  size = 4, nudge_x = 1 ,nudge_y = 0.9, show.legend = FALSE, box.padding = 0.1) +
  scale_fill_manual(values = Species_color_scheme, limits = force) + theme_void()
result <-result + ggtitle(paste("CST ",cst_name)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
result <- result + theme(legend.position = "none")

dir <- paste("CST\ centroids\ pie\ chart/",cst_name)
ggsave(paste(dir, ".png"))
legend <- cowplot::get_legend(result)
grid.newpage()
grid.draw(legend)

