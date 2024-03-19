##%#########################################################################%##
#                                                                             #
#                    Dissertation script - Zoja Manček Páli                   #
#                              Started: 12.3.2024                             #
#                                                                             #
##%#########################################################################%##

#WD
setwd("~/") #erases previously set WDs
setwd("Personal repo - zmancekpali/Dissertation") #sets a new one
getwd() #check that it's worked

#Libraries
library(ggmap)
library(ggspatial)
library(gridExtra)
library(tidyverse)

#Set up the google maps connection
ggmap::register_google(key = "AIzaSyDnersipSvcXuK4tCDbr8NOpa-qsrYf9pc", 
                       write = TRUE) #register your own Google API Key here

#Data
leaves <- read.csv("traits_analysis2.csv")


leaves <- leaves %>% 
  select("type", "code", "latin_name", "long", "lat") %>%  #select the relevant columns
  mutate(type = recode(type, "Alien" = "Alien species",
                       "Invasive" = "Invasive species", 
                       "Naturalised" = "Naturalised species", 
                       "Native" = "Native species")) %>%  #recode the invasion type names
  distinct(long, lat, .keep_all = TRUE) #remove multiple rows (avoids overplotting)

(leaves_counts <- leaves %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))


#Maps
rbge_map <- get_googlemap("Royal Botanic Gardens Edinburgh", zoom = 16, maptype = "satellite")
#Simple RBGE map: ----
(rbge_simple_map <- ggmap(rbge_map) +
    geom_point(data = leaves, aes(x = long, y = lat, color = type, shape = type), 
               size = 3) +
    scale_color_manual(values = c("#5EA8D9", "#CC168F", "green", "#EEC900"),
                       name = "Invasion type") +
    scale_shape_manual(values = c(17, 15, 18, 16), name = "Invasion type") +
    ylab("Latititude") +
    xlab("Longitude") +
    theme_void() +
    theme(legend.position = c(0.85, 0.87),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.size = unit(1.5, "line")) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                   line_col = 'white',
                                                                   fill = 'white')))
ggsave("rbge_map_simple.jpg", rbge_simple_map, path = "Plots", units = "cm", 
       width = 20, height = 20)


#Map with species abbreviations: ----
(abbr_map <- ggmap(rbge_map) +
  geom_point(data = leaves, aes(x = long, y = lat, color = type, shape = type), 
             size = 3) +
  scale_color_manual(values = c("#5EA8D9", "#CC168F", "green", "#EEC900"),
                     name = "Invasion type") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Invasion type") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_void() +
  theme(legend.position = c(0.85, 0.87),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.key.size = unit(1.5, "line")) +
  ggrepel::geom_label_repel(data = leaves, aes(x = long, y = lat, label = code),
                            max.overlaps = 20, box.padding = 0.5, 
                            point.padding = 0.1, segment.color = "white", 
                            size = 3, fontface = "plain") +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                 line_col = 'white',
                                                                 fill = 'white')))
ggsave("rbge_map_with_abbreviations.jpg", abbr_map, path = "Plots", units = "cm", 
       width = 20, height = 20)


#Map with species names: ----
(names_map <- ggmap(rbge_map) +
    geom_point(data = leaves2, aes(x = long, y = lat, color = type, shape = type), 
               size = 3) +
    scale_color_manual(values = c("#5EA8D9", "#CD6090", "#2CB82E", "#EEC900"),
                       name = "Invasion type") +
    scale_shape_manual(values = c(16, 17, 18, 15), name = "Invasion type") +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_void() +
    theme(legend.position = c(0.85, 0.87),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.size = unit(1.5, "line")) +
    ggrepel::geom_label_repel(data = leaves, aes(x = long, y = lat, label = latin_name),
                              max.overlaps = 50, box.padding = 0.5, point.padding = 0.1, 
                              segment.color = "white", size = 2.5, fontface = "italic") +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                   line_col = 'white',
                                                                   fill = 'white')))
ggsave("rbge_map_with_names.jpg", names_map, path = "Plots", units = "cm", 
       width = 20, height = 20)

#Edinburgh map: ----
edinburgh_map <- get_googlemap("Edinburgh", zoom = 12, maptype = "satellite")

(edi_10zoom <- ggmap(edinburgh_map) +
    geom_point(data = NULL, aes(x = -3.209664, y = 55.965140), color = "red", 
               size = 3, shape = 17) +
    theme_void() +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                   line_col = 'white',
                                                                   fill = 'white')))
ggsave("edi_map(10).jpg", edi_10zoom, path = "Plots", units = "cm", 
       width = 20, height = 20)


(edi_11zoom <- ggmap(edinburgh_map) +
    geom_point(data = NULL, aes(x = -3.209664, y = 55.965140), color = "red", 
                                        size = 3, shape = 17) +
    theme_void() +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                   line_col = 'white',
                                                                   fill = 'white')))
ggsave("edi_map(11).jpg", edi_11zoom, path = "Plots", units = "cm", 
       width = 20, height = 20)

(edi_12zoom <- ggmap(edinburgh_map) +
    geom_point(data = NULL, aes(x = -3.209664, y = 55.965140), color = "red", 
               size = 3, shape = 17) +
    theme_void() +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                   line_col = 'white',
                                                                   fill = 'white')))
ggsave("edi_map(12).jpg", edi_12zoom, path = "Plots", units = "cm", 
       width = 20, height = 20)



#Scotland map ----
scotland_map <- get_googlemap("Scotland", zoom = 7, maptype = "satellite")

(scotland <- ggmap(scotland_map) +
  geom_point(data = NULL, aes(x = -3.209664, y = 55.965140), color = "red", 
             size = 3, shape = 17) +
  theme_void() +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering (text_col = 'white',
                                                                 line_col = 'white',
                                                                 fill = 'white')))

ggsave("scotland_map.jpg", scotland, path = "Plots", units = "cm", 
       width = 20, height = 20)

