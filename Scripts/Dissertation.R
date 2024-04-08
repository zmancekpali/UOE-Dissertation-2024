##%#########################################################################%##
#                                                                             #
#                    Dissertation script - Zoja Manček Páli                   #
#                              Started: 30.9.2023                             #
#                                                                             #
##%#########################################################################%##

#WD
setwd("~/") #erases previously set WDs
setwd("UOE Dissertation 2024") #sets a new one
getwd() #check that it's worked


#Libraries
library(ape)
library(cowplot)
library(dunn.test)
library(ecodist)
library(ellipse)
library(e1071)
library(ggfortify)
library(ggpubr)
library(ggrepel)
library(ggsignif)
library(goeveg)
library(grid)
library(gridExtra)
library(lme4)
library(lmtest)
library(MASS)
library(multcomp)
library(nortest)
library(rstatix)
library(sjPlot)
library(stats)
library(tidyverse)
library(vegan)

#Data
trees <- read.csv("Data/traits_analysis2.csv")
trees <- trees %>% 
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  mutate(code_two = recode(code_two,
                           "CB" = "C. bullatus",
                           "R. ponticum" = "Invasive")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus'))) %>%  #rearranges the categories in this order
  mutate(latin_name = recode(latin_name,
                             "Malus pumila" = "Malus x domestica")) %>% #changing the latin name to match Preston et al. (2002)
  mutate(type = recode(type, RR = "Native", PA = "Naturalised")) %>%  #slight changes to the classifications
  rename("LMA" = "lma") %>% 
  rename("LDMC" = "ldcm") %>% 
  rename("Chl" = "chl") %>% 
  rename("Rleaf" = "Dark_resp")

trees_pos <- trees %>% 
  filter(A >= 0) #for physiological measurements

(trees_counts <- trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))
#1 invasive, 12 naturalised, 20 native, 1 alien


nns <- trees %>% 
  filter(type %in% c('Native', 'Naturalised', "Invasive")) %>% #excluding the alien group for initial analysis
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  arrange(type = factor(type, levels = c('Native', 'Naturalised', 'Invasive'))) #rearranges the categories in this order

nns_pos <- nns %>% 
  filter(A >= 0) 

(nns_counts <- nns %>%
  group_by(type) %>%
  summarise(unique_species = n_distinct(code)))
#1 invasive, 12 naturalised, 20 native

traits.palette <- c("Invasive" = "#CD6090", 
                    "Native" = "#698B69",
                    "Naturalised" = "#EEC900", 
                    "C. bullatus" = "#5EA8D9")    #defining 3 colours


cn_trees <- read.csv("Data/cn_analysis.csv")
cn_trees <- cn_trees %>% 
  dplyr::mutate(canopy_pos = recode(canopy_pos, 
                             "M" = "Lower",
                             "U" = "Upper"), #recode canopy positions
         code_two = recode(code_two,
                           "CB" = "C. bullatus", 
                           "R. ponticum" = "Invasive")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus'))) %>% #rearranges the categories in this order
  mutate(c_n = C/N) %>% 
  rename("CN" = "c_n") %>% 
  mutate(type = recode(type, RR = "Native", PA = "Naturalised"))  #slight changes to the classifications

cn_nns <- cn_trees %>% 
  filter(type %in% c('Native', 'Naturalised', "Invasive")) %>% #excluding the alien group for initial analysis
  mutate(canopy_pos = recode(canopy_pos, 
                             "M" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  arrange(type = factor(type, levels = c('Native', 'Naturalised', 'Invasive'))) %>%   #rearranges the categories in this order
  mutate(c_n = C/N)

(cn_counts <- cn_trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))

#Mean trait values for each group - for invasion index
(means_trees <- nns_pos %>% 
  group_by(type) %>% 
  summarise(mean_lma = mean(LMA),
            mean_ldmc = mean(LDMC), 
            mean_chl = mean(Chl),
            mean_A = mean(A),
            mean_E = mean(E),
            mean_g = mean(g),
            mean_dr = mean(DR)))

(means_cn_trees <- cn_trees %>% 
  group_by(type) %>% 
  summarise(mean_cn = mean(C/N)))

#Step 1 - NMDS ----
merged_trees_nns <- merge(nns_pos, cn_nns[, c("code", "canopy_pos", "CN")], 
                          by = c("code", "canopy_pos"))

numeric_cols_nns <- colnames(merged_trees_nns)[sapply(merged_trees_nns, 
                                                      is.numeric)] 
numeric_data_nns <- merged_trees_nns[, numeric_cols_nns]
numeric_data_nns <- numeric_data_nns %>% select(Chl, LMA, LDMC, A, 
                                                E, g, CN, Rleaf)

#finding the lowest stress for up to 6 dimensions:
dimcheckMDS(numeric_data_nns,
            distance = "euclidean",
            k = 6) #goeveg package
#generally accepted that stress < 0.2 is a fair fit for ordination, so will use 2 dimensions (stress = 0.073) 

#2-dimensional NMDS (nns)----
nmds_nns <- metaMDS(numeric_data_nns, distance = "euclidean", k = 2)
nmds_coords_nns <- as.data.frame(scores(nmds_nns, "sites"))
nmds_coords_nns$type <- merged_trees_nns$type

plot(nmds_nns, type = "t") #base r NMDS plot
stressplot(nmds_nns) #stressplot; linear R^2 = 0.995; non-linear R^2 = 0.98
(stress_nns <- nmds_nns$stress) #0.0731529

#ggplot NMDS
hull.data <- data.frame()
for (i in unique(nmds_coords_nns$type)) {
  temp <- nmds_coords_nns[nmds_coords_nns$type == i, ][chull(nmds_coords_nns[nmds_coords_nns$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_nns_plot <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords_nns, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()))
ggsave("nmds_nns_2d.jpg", nmds_nns_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 


diss_matrix_nns <- vegdist(numeric_data_nns, method = "euclidean")
anosim(diss_matrix_nns, merged_trees_nns$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.09491), indicating a slight but significant difference between the groups

#drivers of the variation:
en = envfit(nmds_nns, numeric_data_nns, permutations = 999, na.rm = TRUE)

plot(nmds_nns)
plot(en)
#The arrow(s) point to the direction of most rapid change in the variable (direction of the gradient)
#The length of the arrow(s) is proportional to the correlation between ordination and environmental variable (strength of the gradient)

plot(nmds_nns)
plot(en, p.max = 0.05) #can also only plot the significant ones this way


hull.data <- data.frame()
for (i in unique(nmds_coords_nns$type)) {
  temp <- nmds_coords_nns[nmds_coords_nns$type == i, ][chull(nmds_coords_nns[nmds_coords_nns$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

(drivers_nmds <- ggplot(data = nmds_coords_nns, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.7) + #add polygons for invasive type
    geom_point(data = nmds_coords_nns, aes(colour = type), size = 3) + 
    scale_colour_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) + 
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) + 
    theme_classic() +
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.text = element_text(size = 9, colour = "grey30"),
          legend.position = c(0.93, 0.94), legend.direction = "vertical", legend.title = element_blank()) +
    geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                 size =1, alpha = 0.5, colour = "grey30") +
    geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
              fontface = "bold", label = row.names(en_coord_cont)) +
    labs(colour = "type"))

ggsave("drivers_nmds_nns_2d.jpg", drivers_nmds, path = "Plots", units = "cm", 
       width = 20, height = 20)

#importance of each leaf trait:
en
#lma and g have high r2 values (close to 1), indicating strong correlations with both NMDS1 and NMDS2.
#This suggests that these traits play a significant role in differentiating between invasive and native species in terms of their leaf characteristics.
#Additionally, the significant p-values (***) indicate that these correlations are statistically robust.
#chl and c_n also show significant correlations with the NMDS axes, but to a lesser extent compared to lma and g (still ***)
#ldcm and A have relatively low r2 values and (but are still significany, * and ***, respectively), suggesting weaker correlations with the NMDS axes.
#E and DR have very low r^2 values and are insignificant





#2-dimensional NMDS (with alien - don't know if I will use) ----
merged_trees <- merge(trees_pos, cn_trees[, c("code", "canopy_pos", "CN")], by = c("code", "canopy_pos"))

numeric_cols <- colnames(merged_trees)[sapply(merged_trees, is.numeric)] 
numeric_data <- merged_trees[, numeric_cols]
numeric_data <- numeric_data %>% select(Chl, LMA, LDMC, A, E, g, CN, DR)

#finding the lowest stress for up to 6 dimensions:
dimcheckMDS(numeric_data,
            distance = "euclidean",
            k = 6) #goeveg package

nmds <- metaMDS(numeric_data, distance = "euclidean", k = 2) 
nmds_coords <- as.data.frame(scores(nmds, "sites"))
nmds_coords$type <- merged_trees$type

stressplot(nmds) #non-metric R^2 = 0.995, linear R^2 = 0.979
(stress <- nmds$stress) #0.07331421

hull.data <- data.frame()
for (i in unique(nmds_coords$type)) {
  temp <- nmds_coords[nmds_coords$type == i, ][chull(nmds_coords[nmds_coords$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_plot <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.7) + #add polygons for invasive type
    geom_polygon(data = hull.data[hull.data$type == "Alien", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", 
                                  "Invasive" = "#CD6090", 
                                  "Naturalised" = "#EEC900", 
                                  "Alien" = "#5EA8D9")) +
    scale_fill_manual(values = c("Native" = "#698B69", 
                                 "Invasive" = "#CD6090", 
                                 "Naturalised" = "#EEC900", 
                                 "Alien" = "#5EA8D9")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), 
          legend.direction = "vertical", 
          legend.title = element_blank()))

ggsave("combined_nmds_2d.jpg", nmds_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix <- vegdist(numeric_data, method = "euclidean")
anosim(diss_matrix, merged_trees$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.1405), indicating a slight but significant difference between the groups

#drivers of the variation:
en_total = envfit(nmds, numeric_data, permutations = 999, na.rm = TRUE)

plot(nmds)
plot(en_total)
#The arrow(s) point to the direction of most rapid change in the variable (direction of the gradient)
#The length of the arrow(s) is proportional to the correlation between ordination and environmental variable (strength of the gradient)

plot(nmds)
plot(en_total, p.max = 0.05) #can also only plot the significant ones this way


hull.data1 <- data.frame()
for (i in unique(nmds_coords$type)) {
  temp <- nmds_coords[nmds_coords$type == i, ][chull(nmds_coords[nmds_coords$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data1 <- rbind(hull.data, temp)
}

en_coord_cont_total = as.data.frame(scores(en_total, "vectors")) * ordiArrowMul(en_total)

(drivers_nmds <- ggplot(data = nmds_coords, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data = hull.data1[hull.data1$type != "Invasive", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data1[hull.data1$type == "Invasive", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.7) + #add polygons for invasive type
    geom_polygon(data = hull.data1[hull.data1$type == "Alien", ], 
                 aes(x = NMDS1, y = NMDS2, group = type, fill = type), 
                 alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords, aes(colour = type), size = 3) + 
    scale_color_manual(values = c("Native" = "#698B69", 
                                  "Invasive" = "#CD6090", 
                                  "Naturalised" = "#EEC900", 
                                  "Alien" = "#5EA8D9")) +
    scale_fill_manual(values = c("Native" = "#698B69", 
                                 "Invasive" = "#CD6090", 
                                 "Naturalised" = "#EEC900", 
                                 "Alien" = "#5EA8D9")) +
    theme_classic() +
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.key = element_blank(), 
          legend.text = element_text(size = 9, colour = "grey30"),
          legend.position = c(0.93, 0.93), 
          legend.direction = "vertical", 
          legend.title = element_blank()) +
    geom_segment(data = en_coord_cont, aes(x = 0, y = 0, 
                                           xend = NMDS1, 
                                           yend = NMDS2), 
                 size =1, alpha = 0.5, colour = "grey30") +
    geom_text(data = en_coord_cont_total, aes(x = NMDS1, y = NMDS2), 
              colour = "black", fontface = "bold", 
              label = row.names(en_coord_cont_total)) +
    labs(colour = "type")) #can't get these arrows to work - idk why

ggsave("drivers_nmds_nns_2d.jpg", drivers_nmds, path = "Plots", units = "cm", 
       width = 20, height = 20)

#importance of each leaf trait:
en_total
#lma and g have high r2 values (close to 1), indicating strong correlations with both NMDS1 and NMDS2.
#This suggests that these traits play a significant role in differentiating between invasive and native species in terms of their leaf characteristics.
#Additionally, the significant p-values (***) indicate that these correlations are statistically robust.
#chl and c_n also show significant correlations with the NMDS axes, but to a lesser extent compared to lma and g (still ***)
#A has a relatively low r2 value and (but is still significant ***), suggesting weaker correlations with the NMDS axes.
#E, LDMC, and DR have very low r^2 values and are insignificant




#Step 2: NNS boxplots + post-hoc (if applicable) ----
#LMA ----
lma_mod <- lm(LMA ~ type, data = nns)
autoplot(lma_mod)
shapiro.test(resid(lma_mod)) #residuals not distributed normally
bartlett.test(LMA ~ type, data = nns) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox <- boxcox(LMA ~ 1, data = nns) #the λ is the highest point on the curve
(lma_lambda <- lma_boxcox$x[which.max(lma_boxcox$y)]) #λ = -0.1818182
nns <- nns %>% mutate(transformed_lma = (LMA ^ (lma_lambda - 1)) / lma_lambda) #Box-Cox transformation applied in a new column

lma_mod_trans <- lm(transformed_lma ~ type, data = nns)
autoplot(lma_mod_trans)
shapiro.test(resid(lma_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_lma ~ type, data = nns) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kruskal <- nns %>% kruskal_test(LMA ~ type)) #n = 198; df = 2; p = 0.000542
(lma_effect <- nns %>% kruskal_effsize(LMA ~ type)) #effect size = 0.0669; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0669
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 6.69% of the variance in LMA

#Dunn post-hoc test
(dunn_lma1 <- nns %>% dunn_test(LMA ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #native/invasive and native/naturalised differ significantly

(lma_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                       'Invasive')), #reorders the types 
                           y = LMA, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g cm"^-2*")"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"),
          legend.position = "none"))

ggsave("lma_boxplot1.jpg", lma_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 


#Average chlorophyll ----
chl_mod <- lm(Chl ~ type, data = nns)
autoplot(chl_mod)
shapiro.test(resid(chl_mod)) #residuals not distributed normally
bartlett.test(Chl ~ type, data = nns) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
chl_boxcox <- boxcox(nns$Chl ~ 1)
(chl_lambda <- chl_boxcox$x[which.max(chl_boxcox$y)]) #λ = -0.1010101
nns <- nns %>% mutate(transformed_chl = (Chl ^ (chl_lambda - 1)) / chl_lambda) #Box-Cox transformation applied in a new column

chl_mod_trans <- lm(transformed_chl ~ type, data = nns)
autoplot(chl_mod_trans)
shapiro.test(resid(chl_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_chl ~ type, data = nns) #heteroscedascity           

#Transformation did not work, moving on to non-parametric alternative:
(chl_kruskal <- nns %>% kruskal_test(Chl ~ type)) #n = 198; df = 2; p = 0.00037 
(chl_effect <- nns %>% kruskal_effsize(Chl ~ type)) #effect size = 0.0708; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0708
#this value indicates the % of variance in the dependent variable (chl) explained by the invasion status
#so this explains 7.08% of the variance in chl


#Dunn post-hoc test
(dunn_chl1 <- nns %>% dunn_test(Chl ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #native/invasive and naturalised/invasive differ significantly

(chl_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                       'Invasive')), #reorders the types 
                           y = Chl, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("Chlorophyll content (SPAD)"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("chl_boxplot.jpg", chl_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 

#Assimilation rate ----
a_mod <- lm(A ~ type, data = nns_pos)
autoplot(a_mod)
shapiro.test(resid(a_mod)) #residuals distributed normally
bartlett.test(A ~ type, data = nns_pos) #heteroscedascity
anova(a_mod) #p = 0.0006566; df = 2

tukey_a <- aov(a_mod)
TukeyHSD(tukey_a, conf.level = 0.95) #native/invasive and naturalised/invasive differ significantly

(a_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                     'Invasive')), #reorders the types 
                         y = A, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("A (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none")) 

ggsave("a_boxplot.jpg", a_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 

#LDCM ----
ldcm_mod <- lm(LDMC ~ type, data = nns)
autoplot(ldcm_mod)
shapiro.test(resid(ldcm_mod)) #residuals distributed normally
bartlett.test(LDMC ~ type, data = nns) #homoscedascity
anova(ldcm_mod) #p = 0.001018; df = 2

tukey_ldcm <- aov(ldcm_mod)
TukeyHSD(tukey_ldcm, conf.level = 0.95) #native/invasive and naturalised/invasive differ significantly

(ldmc_boxplot <- ggplot(nns, 
                        aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                        'Invasive')), #reorders the types 
                            y = LDMC, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("LDMC (g" ~ "g"^-1~")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("ldmc_boxplot.jpg", ldmc_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 

#Transpiration rate ----
e_mod <- lm(E ~ type, data = nns_pos)
autoplot(e_mod)
shapiro.test(resid(e_mod)) #residuals distributed normally
bartlett.test(E ~ type, data = nns_pos) #homoscedascity
anova(e_mod) #NS; p-value = 0.5231

#no post-hoc test because no significant differences found overall

(e_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                     'Invasive')), #reorders the types 
                         y = E, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(paste("E (", mu, "mol H"[2]*"O" ~ "m"^-2*"s"^-1, ")"))) +    
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("e_boxplot.jpg", e_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 

#Dark respiration ----
dr_mod <- lm(Rleaf ~ type, data = nns_pos)
autoplot(dr_mod)
shapiro.test(resid(dr_mod)) #residuals not distributed normally
bartlett.test(Rleaf ~ type, data = nns_pos) #heteroscedascity (just barely)

#cannot do a boxcox transformation bc values are negative       

(dr_kruskal <- nns_pos %>% kruskal_test(Rleaf ~ type)) #n = 196; df = 2; p = 0.288  
(dr_effect <- nns_pos %>% kruskal_effsize(Rleaf ~ type)) #effect size = 0.00253; small magnitude
#report as: small effect size is detected, eta2[H] = 0.00253
y <- expression(paste("R"[leaf], " (", mu, "mol CO"[2], "m"^-2*"s"^-1, ")"))

(dr_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                     'Invasive')), #reorders the types 
                         y = Rleaf, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = y) +    
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("dr_boxplot.jpg", dr_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 

#GH20 ----
g_mod <- lm(g ~ type, data = nns_pos)
autoplot(g_mod)
shapiro.test(resid(g_mod)) #residuals not distributed normally
bartlett.test(g ~ type, data = nns_pos) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
g_boxcox <- boxcox(nns_pos$g ~ 1)
(g_lambda <- g_boxcox$x[which.max(g_boxcox$y)]) #λ = -0.6666667
nns_pos <- nns_pos %>% mutate(transformed_g = (g ^ (g_lambda - 1)) / g_lambda) #Box-Cox transformation applied in a new column

g_mod_trans <- lm(transformed_g ~ type, data = nns_pos)
autoplot(g_mod_trans)
shapiro.test(resid(g_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_g ~ type, data = nns_pos) #homoscedascity           

#Transformation did not work, moving on to non-parametric alternative:
(g_kruskal <- nns_pos %>% kruskal_test(g ~ type)) #n = 196; df = 2; p = 0.0256  
(g_effect <- nns_pos %>% kruskal_effsize(g ~ type)) #effect size = 0.0276; small magnitude
#report as: small effect size is detected, eta2[H] = 0.0276 
#this value indicates the % of variance in the dependent variable (g) explained by the invasion status
#so this explains 2.76% of the variance in g

#Dunn post-hoc test
(dunn_g1 <- nns_pos %>% dunn_test(g ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #invasive/native and invsive/naturalised differ significantly

(g_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                     'Invasive')), #reorders the types 
                         y = g, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(paste("gs (", mu, "mol H"[2]*"O" ~ "m"^-2*"s"^-1, ")"))) +    
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("g_boxplot.jpg", g_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 


#C:N ----
cn_mod <- lm(CN ~ type, data = cn_nns)
autoplot(cn_mod)
shapiro.test(resid(cn_mod)) #residuals not distributed normally
bartlett.test(CN ~ type, data = cn_nns) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox <- boxcox(CN ~ 1, data = cn_nns) #the λ is the highest point on the curve
(cn_lambda <- cn_boxcox$x[which.max(cn_boxcox$y)]) #λ = -0.7070707
cn_nns <- cn_nns %>% mutate(transformed_cn = (CN ^ (cn_lambda - 1)) / cn_lambda) #Box-Cox transformation applied in a new column

cn_mod_trans <- lm(transformed_cn ~ type, data = cn_nns)
autoplot(cn_mod_trans)
shapiro.test(resid(cn_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_cn ~ type, data = cn_nns) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kruskal <- cn_nns %>% kruskal_test(CN ~ type)) #n = 96; df = 2; p = 0.00108   
(cn_effect <- cn_nns %>% kruskal_effsize(CN ~ type)) #effect size = 0.125  ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.125  
#this value indicates the % of variance in the dependent variable (C/N) explained by the invasion status
#so this explains 12.5% of the variance in C/N

#Dunn post-hoc test
(dunn_cn1 <- cn_nns %>% dunn_test(CN ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #invasive/native differ significantly

(cn_boxplot <- ggplot(cn_nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 
                                                     'Invasive')), #reorders the types 
                         y = CN, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("C/N ratio"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("cn_boxplot.jpg", cn_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 



#Step 3: compare alien species ----
#LMA ----
lma_mod2 <- lm(LMA ~ type, data = trees)
autoplot(lma_mod2)
shapiro.test(resid(lma_mod2)) #residuals not distributed normally
bartlett.test(LMA ~ type, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox2 <- boxcox(LMA ~ 1, data = trees) #the λ is the highest point on the curve
(lma_lambda2 <- lma_boxcox2$x[which.max(lma_boxcox2$y)]) #λ = -0.1414141
trees <- trees %>% mutate(transformed_lma2 = (LMA ^ (lma_lambda2 - 1)) / lma_lambda2) #Box-Cox transformation applied in a new column

lma_mod_trans2 <- lm(transformed_lma2 ~ type, data = trees)
autoplot(lma_mod_trans2)
shapiro.test(resid(lma_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_lma2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kruskal2 <- trees %>% kruskal_test(LMA ~ type)) #n = 204; df = 3; p = 0.000019 
(lma_effect2 <- trees %>% kruskal_effsize(LMA ~ type)) #effect size = 0.108 ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.108 

#Dunn post-hoc test
(dunn_lma2 <- trees %>% dunn_test(LMA ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly
#invasive/native and invasive/naturalised differ significantly

(lma_boxplot2 <- ggplot(trees, 
                       aes(x = factor(code_two, levels = 
                                        c('Native', 'Naturalised', 'Invasive', 
                                          'C. bullatus')), #reorders the types 
                           y = LMA, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) +
    labs(x = "\n Invasion status", 
         y = expression("LMA (g cm"^-2*")")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  #italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("lma_boxplot2.jpg", lma_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 


#Chlorophyll ----
chl_mod2 <- lm(Chl ~ type, data = trees)
autoplot(chl_mod2)
shapiro.test(resid(chl_mod2)) #residuals not distributed normally
bartlett.test(Chl ~ type, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
chl_boxcox2 <- boxcox(Chl ~ 1, data = trees) #the λ is the highest point on the curve
(chl_lambda2 <- chl_boxcox2$x[which.max(chl_boxcox2$y)]) #λ = -0.06060606
trees <- trees %>% 
  mutate(transformed_chl2 = (Chl ^ (chl_lambda2 - 1)) / chl_lambda2) #Box-Cox transformation applied in a new column

chl_mod_trans2 <- lm(transformed_chl2 ~ type, data = trees)
autoplot(chl_mod_trans2)
shapiro.test(resid(chl_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_chl2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(chl_kruskal2 <- trees %>% kruskal_test(Chl ~ type)) #n = 204; df = 3; p = 0.000014 
(chl_effect2 <- trees %>% kruskal_effsize(Chl ~ type)) #effect size = 0.111; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.111

#Dunn post-hoc test
(dunn_chl2 <- trees %>% dunn_test(Chl ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly
#invasive/native and invasive/naturalised differ significantly

(chl_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 
                                           'Invasive', 'C. bullatus')), #reorders the types 
                            y = Chl, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9"),
                      breaks = c("Native", "Naturalised", 
                                 "Invasive", 
                                 "C. bullatus"),
                      labels = c("Native", "Naturalised", 
                                 "Invasive", expression(italic("C. bullatus")))) + #colours each boxplot this particular colour
    labs(x = "", 
         y = expression("Chlorophyll content (SPAD)"),
         fill = "Invasion status") + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm"),
          legend.position = "none"))

ggsave("chl_boxplot2.jpg", chl_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 


#LDCM ----
ldcm_mod2 <- lm(LDMC ~ type, data = trees)
autoplot(ldcm_mod2)
shapiro.test(resid(ldcm_mod2)) #residuals distributed normally
bartlett.test(LDMC ~ type, data = trees) #homoscedascity
anova(ldcm_mod2) #significant, p = 0.002177; df = 3 

anova_ldmc <- aov(ldcm_mod2)
TukeyHSD(anova_ldmc)
#invasive/alien differ significantly
#native/invasive and naturalised/invasive differ significantly

(ldcm_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 
                                           'Invasive', 'C. bullatus')), #reorders the types 
                            y = LDMC, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("LDMC (g" ~ "g"^-1~")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("ldcm_boxplot2.jpg", ldcm_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 



#Assimilation rate ----
a_mod2 <- lm(A ~ code_two, data = trees_pos)
autoplot(a_mod2)
shapiro.test(resid(a_mod2)) #residuals distributed normally
bartlett.test(A ~ type, data = trees_pos) #homoscedascity
anova(a_mod2) #p = 0.003221; significant; df = 3

anova_a <- aov(a_mod2)
TukeyHSD(anova_a)
#invasive/alien differ significantly
#native/invasive and naturalised/invasive differ significantly

(a_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 
                                         'Invasive', 'C. bullatus')), #reorders the types 
                          y = A, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("A (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("a_boxplot2.jpg", a_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 



#Transpiration rate ----
e_mod2 <- lm(E ~ code_two, data = trees_pos)
autoplot(e_mod2)
shapiro.test(resid(e_mod2)) #residuals distributed normally
bartlett.test(E ~ code_two, data = trees_pos) #homoscedascity
anova(e_mod2) #NS; p-value = 0.8372

(e_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 
                                         'Invasive', 'C. bullatus')), #reorders the types 
                          y = E, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("E (", mu, "mol H"[2]*"O m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("e_boxplot2.jpg", e_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 





#Dark respiration ----
dr_mod2 <- lm(Rleaf ~ code_two, data = trees_pos)
autoplot(dr_mod2)
shapiro.test(resid(dr_mod2)) #residuals not distributed normally
bartlett.test(Rleaf ~ type, data = trees_pos) #heteroscedascity

#Cannot do a Box-Cox transformation because values are negative:
(dr_kruskal2 <- trees_pos %>% kruskal_test(Rleaf ~ type)) #n = 202; df = 3; p = 0.437  
(dr_effect2 <- trees_pos %>% kruskal_effsize(Rleaf ~ type)) #effect size = -0.00142; small magnitude
#report as: small effect size is detected, eta2[H] = -0.00142

y <- expression(paste("R"[leaf], " (", mu, "mol CO"[2], "m"^-2*"s"^-1, ")"))
(dr_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 
                                         'Invasive', 'C. bullatus')), #reorders the types 
                          y = Rleaf, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = y) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("dr_boxplot2.jpg", dr_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 

#GH2O ----
g_mod2 <- lm(g ~ code_two, data = trees_pos)
autoplot(g_mod2)
shapiro.test(resid(g_mod2)) #residuals not distributed normally
bartlett.test(g ~ code_two, data = trees_pos) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
g_boxcox2 <- boxcox(g ~ 1, data = trees_pos) #the λ is the highest point on the curve
(g_lambda2 <- g_boxcox2$x[which.max(g_boxcox2$y)]) #λ = -0.5454545
trees_pos <- trees_pos %>% 
  mutate(transformed_g2 = (g ^ (g_lambda2 - 1)) / g_lambda2) #Box-Cox transformation applied in a new column

g_mod_trans2 <- lm(transformed_g2 ~ type, data = trees_pos)
autoplot(g_mod_trans2)
shapiro.test(resid(g_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_g2 ~ type, data = trees_pos) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(g_kruskal2 <- trees_pos %>% kruskal_test(g ~ code_two)) #n = 202; df = 3; p = 0.000366
(g_effect2 <- trees_pos %>% kruskal_effsize(g ~ code_two)) #effect size = 0.0777; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0777

#Dunn post-hoc test
(dunn_g2 <- trees_pos %>% dunn_test(g ~ code_two, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly from one another
#invasive/naturalised differ significantly from one another

(g_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 
                                         'Invasive', 'C. bullatus')), #reorders the types 
                          y = g, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("g (", mu, "mol H"[2]*"O m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  #italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("g_boxplot2.jpg", g_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 


#C:N ----
cn_mod2 <- lm(CN ~ type, data = cn_trees)
autoplot(cn_mod2)
shapiro.test(resid(cn_mod2)) #residuals not distributed normally
bartlett.test(CN ~ type, data = cn_trees) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox2 <- boxcox(CN ~ 1, data = cn_trees) #the λ is the highest point on the curve
(cn_lambda2 <- cn_boxcox2$x[which.max(cn_boxcox2$y)]) #λ = -0.7474747
cn_trees <- cn_trees %>% 
  mutate(transformed_cn2 = (CN ^ (cn_lambda2 - 1)) / cn_lambda2) #Box-Cox transformation applied in a new column

cn_mod_trans2 <- lm(transformed_cn2 ~ type, data = cn_trees)
autoplot(cn_mod_trans2)
shapiro.test(resid(cn_mod_trans2)) #residuals distributed normally
bartlett.test(transformed_cn2 ~ type, data = cn_trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kruskal2 <- cn_trees %>% kruskal_test(CN ~ type)) #n = 99; df = 3; p = 0.00182 
(cn_effect2 <- cn_trees %>% kruskal_effsize(CN ~ type)) #effect size = 0.126 ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.126 

#Dunn post-hoc test
(dunn_cn2 <- cn_trees %>% dunn_test(CN ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type"))  
#invasive/native differ significantly
#invasive/alien differ significantly

(cn_boxplot2 <- ggplot(cn_trees, 
                      aes(x = factor(code_two, levels = c('Native', 'Naturalised', 
                                                          'Invasive', 'C. bullatus')), #reorders the types 
                          y = CN, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", 
                                 "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", 
                                 "C. bullatus" = "#5EA8D9")) + 
    labs(x = "", 
         y = expression("C/N ratio")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", 
                                              "plain", "italic")),  #italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("cn_boxplot2.jpg", cn_boxplot2, path = "Plots", units = "cm",
       width = 25, height = 13) 


#Box plots grid ----
#These remove the x-axis tick marks, titles, and labels
a_boxplot2 <- a_boxplot2 + theme(axis.title.x = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.ticks.x = element_blank())
dr_boxplot2 <- dr_boxplot2 + theme(axis.title.x = element_blank(), 
                                   axis.text.x = element_blank(), 
                                   axis.ticks.x = element_blank())
e_boxplot2 <- e_boxplot2 + theme(axis.title.x = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.ticks.x = element_blank())
g_boxplot2 <- g_boxplot2 + theme(axis.title.x = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.ticks.x = element_blank())
lma_boxplot2 <- lma_boxplot2 + theme(axis.title.x = element_blank(), 
                                     axis.text.x = element_blank(), 
                                     axis.ticks.x = element_blank())
ldcm_boxplot2 <- ldcm_boxplot2 + theme(axis.title.x = element_blank(), 
                                       axis.text.x = element_blank(), 
                                       axis.ticks.x = element_blank())

(grid1 <- grid.arrange(e_boxplot2, dr_boxplot2, g_boxplot2, a_boxplot2,
                      lma_boxplot2, ldcm_boxplot2, cn_boxplot2, chl_boxplot2, 
                      ncol = 2, widths = c(0.8, 0.8), heights = c(1, 1, 1, 1)))

ggsave("boxplot_grid.jpg", grid1, path = "Plots", units = "cm",
       width = 25, height = 30)

legend <- cowplot::get_legend(chl_boxplot2) #extract legend from previous plot 
#only works if you set the legend.position = "bottom" before running the below code

#run these two also to get no x-axis tick marks:
cn_boxplot2 <- cn_boxplot2 + theme(axis.title.x = element_blank(), 
                                   axis.text.x = element_blank(), 
                                   axis.ticks.x = element_blank())
chl_boxplot2 <- chl_boxplot2 + theme(axis.title.x = element_blank(), 
                                     axis.text.x = element_blank(), 
                                     axis.ticks.x = element_blank())
(grid2 <- grid.arrange(e_boxplot2, a_boxplot2, g_boxplot2, dr_boxplot2,
                      ldcm_boxplot2, lma_boxplot2, cn_boxplot2, chl_boxplot2, 
                      ncol = 2, widths = c(1, 1), heights = c(1, 1, 1, 1)))

(grid_final <- cowplot::plot_grid(
  grid2, legend, 
  nrow = 2, rel_heights = c(1, 0.05)))  #Adjust the relative heights as needed

ggsave("boxplot_final_grid.jpg", grid_final, path = "Plots", units = "cm",
       width = 25, height = 32)



#Step 4 - LMs/GLMs for random effects ----
#in order to select the additional traits, i chose the traits which had the highest r^2 from the NMDS (add the en table in the results)
#or do reading and pick traits that are most interesting for invasion potential (or maybe one of each group)
#LMA LM ----
hist(nns$lma) #looks normal-ish; slightly skewed
shapiro_test(nns$lma) #non-normal

#trying a log-transformation
hist(log(nns$lma)) #normal
shapiro_test(log(nns$lma)) #normal

lma_lm_null <- lm(log(lma) ~ 1, data = nns)
lma_lm1 <- lm(log(lma) ~ type, data = nns)
plot(lma_lm1)
ad.test(residuals(lma_lm1)) #p = 0.1062; normally distributed errors
lma_lm2 <- lm(log(lma) ~ type + age + dbh + ever_dec + canopy_pos, data = nns)
plot(lma_lm2)
ad.test(residuals(lma_lm2)) #p = 0.7818; normally distributed errors
#both of these are ok for the data

#stepwise model selection -> start with a full model and decrease by each trait (least significant)


#remove species because only one sample per species -> might not be representative, also studies that talk ab interspecific variation

#see which model is most parismoius
AIC(lma_lm_null, lma_lm1, lma_lm2) #lma_lm2 is most parsimonious
summary(lma_lm2)
#we can see a significant effect of light and decidousness on LMA

#invasion type is significant, but deciduousness is more important to LMA (based on the above result)
#further research -> discuss the effects of these additional points whether this is more important in native/naturalised
#also how many are deciduous/evergreen from each group for this study

(interspecific_lma <- ggplot(trees, aes(x = type, y = lma, fill = ever_dec)) +
    geom_boxplot() +
    theme_classic() +
    ylab(bquote("LMA (g cm"^{-2}*")")) +
    xlab("Species") +
    theme(legend.position = c(0.95, 0.95),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.size = unit(1.5, "line"),
          legend.title = element_blank()))

(interspecific_lma <- ggplot(trees, aes(x = type, y = lma, fill = canopy_pos)) +
    geom_boxplot() +
    theme_classic() +
    ylab(bquote("LMA (g cm"^{-2}*")")) +
    xlab("Species") +
    theme(legend.position = c(0.95, 0.95),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.size = unit(1.5, "line"),
          legend.title = element_blank()))


ggsave("interspecific_lma_variation.jpg", interspecific_lma, path = "Plots", 
       width = 30, height = 12)

#refer to the appendix when discussing -> just a sentence, not too much

#Chl GLM ----
hist(nns$chl) #looks normal
shapiro_test(nns$chl) #non-normal

#chl values out of 95% interval were removed (or some other number)

#see if removing it is also fine for boxplots + NMDS (but if necessary just include to improve the normality of the chl data)
#if that doesnt work, select the best transformation -> this reduced the p value of the Shapiro test (still significant)

#trying a log-transformation
hist(log(nns$chl)) #looks normal
shapiro_test(log(nns$chl)) #non-normal

#inverse
hist(1/nns$chl) #looks normal
shapiro_test(1/nns$chl) #non-normal

#sqrt
hist(sqrt(nns$chl)) #looks normal
shapiro_test(sqrt(nns$chl)) #non-normal

hist((nns$chl)) #looks normal
shapiro_test(sqrt(nns$chl)) #non-normal


chl_lm_null_reg <- lm(chl ~ 1, data = nns)
chl_lm1_reg <- lm(chl ~ type, data = nns)
plot(chl_lm1_reg)
ad.test(residuals(chl_lm1_reg)) 


chl_lm_null <- lm(log(chl) ~ 1, data = nns)
chl_lm1 <- lm(log(chl) ~ type, data = nns)
plot(chl_lm1)
ad.test(residuals(chl_lm1)) 


chl_lm_null_sqrt <- lm(sqrt(chl) ~ 1, data = nns)
chl_lm1_sqrt <- lm(sqrt(chl) ~ type, data = nns)
plot(chl_lm1_sqrt)
ad.test(residuals(chl_lm1_sqrt)) #p = 0.006; normally distributed errors

chl
chl_lm2 <- lm(chl ~ type + age + ever_dec + canopy_pos, data = nns)
plot(chl_lm2)
ad.test(residuals(chl_lm2)) #p = 0.5144; normally distributed errors

#none of them worked, so glm it is
chl_null_glm <- glm(chl ~ 1, data = trees, family = 'gaussian')
chl_glm <- glm(lma ~ type + age + ever_dec + canopy_pos, data = trees, family = "gaussian") #basic glm
plot(chl_glm)
ad.test(residuals(chl_glm)) #p = 0.0005; non-normally distributed errors

chl_glm2 <- glm(chl ~ type + age + ever_dec + canopy_pos, data = trees, family = Gamma(link = "inverse")) #default gamma distribution link
#used the gamma distribution because the outcome (chl) is continuous but non-normal
plot(chl_glm2) 
ad.test(residuals(chl_glm2)) #p = 0.03; non-normally distributed errors; not to be used

AIC(chl_lm_null_reg, chl_lm1_reg, chl_lm_null, chl_lm1, chl_lm_null_sqrt, 
    chl_lm1_sqrt, chl_lm2, chl_null_glm, chl_glm) 
#chl_lm1 <- lm(log(chl) ~ type, data = trees) is best
#additional variables don't explain the variation in the data


summary(lma_glm2)
#overdispersion test: residual deviance/df = 0.02072451; not overdispersed (bc < 2)
#RP tends to have a higher LMA than natives and naturalised
#This suggests that invasive species may invest more in leaf structural components, resulting in higher LMA.
#tree age has a significant effect on LMA (older trees = lower lma)
#higher dbh = higher lma
#upper leaves and evergreen leaves have higher LMA (boxplot only; model disagrees)

#options:
#1. remove outliers
#2. glm
#3. sqrt transformation







#g model ----
hist(nns$g)
#transformations don't work

boxplot(g ~ type, data = nns)
g_glm <- glm(g ~ type + age + dbh + canopy_pos + ever_dec, data = nns, family = Gamma(link = "inverse")) #default gamma distribution link
g_null_glm <- glm(g ~ 1, data = nns, family = Gamma(link = "inverse"))

ad.test(residuals(g_glm)) 
#not working -> maybe find justification for why they are important

