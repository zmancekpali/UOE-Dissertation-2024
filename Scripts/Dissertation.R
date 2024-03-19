##%#########################################################################%##
#                                                                             #
#                    Dissertation script - Zoja Manček Páli                   #
#                              Started: 30.9.2023                             #
#                                                                             #
##%#########################################################################%##

#WD
setwd("~/") #erases previously set WDs
setwd("Personal repo - zmancekpali/Dissertation") #sets a new one
getwd() #check that it's worked


#Libraries
library(ape)
library(cowplot)
library(dunn.test)
library(e1071)
library(ggfortify)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(lme4)
library(MASS)
library(multcomp)
library(stats)
library(tidyverse)
library(vegan)

#Data
trees <- read.csv("traits_analysis2.csv")
trees <- trees %>% 
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  mutate(code_two = recode(code_two,
                           "CB" = "C. bullatus",
                           "RPS" = "R. pseudoacacia semperflorens")) %>% #recode alien species names
  filter(A >= 0) %>% 
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus', 'R. pseudoacacia semperflorens'))) #rearranges the categories in this order

trees$age <- as.numeric(trees$age)

(trees_counts <- trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))
#1 invasive, 12 naturalised, 20 native, 2 alien


nns <- trees %>% 
  filter(type %in% c('Native', 'Naturalised', "Invasive")) %>% #excluding the alien group for initial analysis
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  arrange(type = factor(type, levels = c('Native', 'Naturalised', 'Invasive'))) %>%  #rearranges the categories in this order
  filter(A >= 0) #removed negative A values (they were dead leaves)

(nns_counts <- nns %>%
  group_by(type) %>%
  summarise(unique_species = n_distinct(code)))
#1 invasive, 12 naturalised, 20 native

traits.palette <- c("#CD6090", "#698B69", "#EEC900")    #defining 3 colours
traits.palette2 <- c("#CD6090", "#698B69", "#EEC900", "#5EA8D9", "#245C82", "#4A3E87", "#5A5DC7")

cn_trees <- read.csv("cn_analysis.csv")
cn_trees <- cn_trees %>% 
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  mutate(code_two = recode(code_two,
                           "CB" = "C. bullatus",
                           "RPS" = "R. pseudoacacia semperflorens")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus',
                                             'R. pseudoacacia semperflorens'))) %>% #rearranges the categories in this order
  mutate(c_n = C/N)

cn_trees$age <- as.numeric(cn_trees$age)

cn_nns <- cn_trees %>% 
  filter(type %in% c('Native', 'Naturalised', "Invasive")) %>% #excluding the alien group for initial analysis
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  arrange(type = factor(type, levels = c('Native', 'Naturalised', 'Invasive'))) %>%   #rearranges the categories in this order
  mutate(c_n = C/N)

(cn_counts <- cn_trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))


#Exploration
head(nns)
str(nns)
head(cn_nns)

#Step 1: see whether NN and I species differ in their traits + respective box plots + post-hoc tests if significant
#LMA ----
lma_mod <- lm(lma ~ type, data = nns)
autoplot(lma_mod)
shapiro.test(resid(lma_mod)) #residuals not distributed normally
bartlett.test(lma ~ type, data = nns) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox <- boxcox(lma ~ 1, data = nns) #the λ is the highest point on the curve
(lma_lambda <- lma_boxcox$x[which.max(lma_boxcox$y)]) #λ = -0.1818182
nns <- nns %>% mutate(transformed_lma = (lma ^ (lma_lambda - 1)) / lma_lambda) #Box-Cox transformation applied in a new column

lma_mod_trans <- lm(transformed_lma ~ type, data = nns)
autoplot(lma_mod_trans)
shapiro.test(resid(lma_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_lma ~ type, data = nns) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kw <- kruskal.test(lma ~ type, data = nns)) #p-value = 0.0005229; significant

(lma_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                           y = lma, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g/cm"^2*")"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("lma_boxplot.jpg", lma_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_lma <- dunn.test(nns$lma, nns$type, method = "bonferroni") #invasives differ significantly from natives yay
#naturalised also differ significantly from natives


#Average chlorophyll ----
chl_mod <- lm(chl ~ type, data = nns)
autoplot(chl_mod)
shapiro.test(resid(chl_mod)) #residuals not distributed normally
bartlett.test(chl ~ type, data = nns) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
chl_boxcox <- boxcox(nns$chl ~ 1)
(chl_lambda <- chl_boxcox$x[which.max(chl_boxcox$y)]) #λ = -0.1010101
nns <- nns %>% mutate(transformed_chl = (chl ^ (chl_lambda - 1)) / chl_lambda) #Box-Cox transformation applied in a new column

chl_mod_trans <- lm(transformed_chl ~ type, data = nns)
autoplot(chl_mod_trans)
shapiro.test(resid(chl_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_chl ~ type, data = nns) #heteroscedascity           

#Transformation did not work, moving on to non-parametric alternative:
kruskal.test(chl ~ type, data = nns) #0.0003969; significant

(chl_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                           y = chl, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("Average chlorophyll (SPAD)"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("chl_boxplot.jpg", chl_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_chl <- dunn.test(nns$chl, nns$type, method = "bonferroni") #invasive differ significantly from natives and naturalised yay
#naturalised and natives show no significant difference

#Assimilation rate ----
a_mod <- lm(A ~ type, data = nns)
autoplot(a_mod)
shapiro.test(resid(a_mod)) #residuals distributed normally
bartlett.test(A ~ type, data = nns) #homoscedascity

anova(a_mod)
#p = 0.0006566

(a_boxplot <- ggplot(nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = A, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Assimilation rate (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none")) 

ggsave("a_boxplot.jpg", a_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_a <- dunn.test(nns$A, nns$type, method = "bonferroni") #invasives differ significantly from natives and naturalised species
#naturalised and natives show no significant difference



#LDCM ----
ldcm_mod <- lm(ldcm ~ type, data = nns)
autoplot(ldcm_mod)
shapiro.test(resid(ldcm_mod)) #residuals distributed normally
bartlett.test(ldcm ~ type, data = nns) #heteroscedascity
anova(ldcm_mod)
#p = 0.001115; significant

(ldmc_boxplot <- ggplot(nns, 
                        aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                            y = ldcm, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Leaf dry matter concentration (g" ~ "g"^-1~")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("ldmc_boxplot.jpg", ldmc_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_ldmc <- dunn.test(nns$ldcm, nns$type, method = "bonferroni") #invasive differ significantly from natives and naturalised species
#naturalised and natives show no significant difference


#Evapotransiration rate ----
e_mod <- lm(E ~ type, data = nns)
autoplot(e_mod)
shapiro.test(resid(e_mod)) #residuals distributed normally
bartlett.test(E ~ type, data = nns) #homoscedascity
anova(e_mod) #NS; p-value = 0.5231

(e_boxplot <- ggplot(nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = E, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Evapotranspiration rate (", mu, "mol CO"[2] ~ "m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("e_boxplot.jpg", e_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Tukey's Honestly Significant Difference post-hoc test
aov_e <- aov(E ~ type, data = nns)
(tukey_e <- TukeyHSD(aov_e))
#no significant differences


#GH20 ----
g_mod <- lm(g ~ type, data = nns)
autoplot(g_mod)
shapiro.test(resid(g_mod)) #residuals not distributed normally
bartlett.test(g ~ type, data = nns) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
g_boxcox <- boxcox(nns$g ~ 1)
(g_lambda <- g_boxcox$x[which.max(g_boxcox$y)]) #λ = -0.6666667
nns <- nns %>% mutate(transformed_g = (g ^ (g_lambda - 1)) / g_lambda) #Box-Cox transformation applied in a new column

g_mod_trans <- lm(transformed_g ~ type, data = nns)
autoplot(g_mod_trans)
shapiro.test(resid(g_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_g ~ type, data = nns) #homoscedascity           

#Transformation did not work, moving on to non-parametric alternative:
kruskal.test(g ~ type, data = nns) #0.02564; significant

(g_boxplot <- ggplot(nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = g, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Stomatal conductance rate (", mu, "mol CO"[2] ~ "m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("g_boxplot.jpg", g_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_g <- dunn.test(nns$g, nns$type, method = "bonferroni") #invasive differ significantly from natives and naturalised yay
#naturalised and natives show no significant difference


#C:N ----
cn_mod <- lm(c_n ~ type, data = cn_nns)
autoplot(cn_mod)
shapiro.test(resid(cn_mod)) #residuals not distributed normally
bartlett.test(c_n ~ type, data = cn_nns) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox <- boxcox(c_n ~ 1, data = cn_nns) #the λ is the highest point on the curve
(cn_lambda <- cn_boxcox$x[which.max(cn_boxcox$y)]) #λ = -0.7070707
cn_nns <- cn_nns %>% mutate(transformed_cn = (c_n ^ (cn_lambda - 1)) / cn_lambda) #Box-Cox transformation applied in a new column

cn_mod_trans <- lm(transformed_cn ~ type, data = cn_nns)
autoplot(cn_mod_trans)
shapiro.test(resid(cn_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_cn ~ type, data = cn_nns) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kw <- kruskal.test(c_n ~ type, data = cn_nns)) #p-value = 0.001083; significant

(cn_boxplot <- ggplot(cn_nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = c_n, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("C:N"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("cn_boxplot.jpg", cn_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_cn <- dunn.test(cn_nns$c_n, cn_nns$type, method = "bonferroni") #invasives differ significantly from natives





#Step 2: compare alien species ----
#LMA ----
lma_mod2 <- lm(lma ~ code_two, data = trees)
autoplot(lma_mod2)
shapiro.test(resid(lma_mod2)) #residuals not distributed normally
bartlett.test(lma ~ code_two, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox2 <- boxcox(lma ~ 1, data = trees) #the λ is the highest point on the curve
(lma_lambda2 <- lma_boxcox2$x[which.max(lma_boxcox2$y)]) #λ = -0.1414141
trees <- trees %>% mutate(transformed_lma2 = (lma ^ (lma_lambda2 - 1)) / lma_lambda2) #Box-Cox transformation applied in a new column

lma_mod_trans2 <- lm(transformed_lma2 ~ type, data = trees)
autoplot(lma_mod_trans2)
shapiro.test(resid(lma_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_lma2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kw2 <- kruskal.test(lma ~ type, data = trees)) #p-value = 6.645e-05; significant

(lma_boxplot2 <- ggplot(trees, 
                       aes(x = factor(code_two, levels = 
                                        c('Native', 'Naturalised', 'Invasive', 
                                          'C. bullatus', 
                                          'R. pseudoacacia semperflorens')), #reorders the types 
                           y = lma, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g/cm"^2*")"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("lma_boxplot2.jpg", lma_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Dunn post-hoc test
dunn_lma_2 <- dunn.test(trees$lma, trees$code_two, method = "bonferroni") 
#C. bullatus differs significantly from native species; RPS does not

#Average chlorophyll ----
chl_mod2 <- lm(chl ~ code_two, data = trees)
autoplot(chl_mod2)
shapiro.test(resid(chl_mod2)) #residuals not distributed normally
bartlett.test(chl ~ code_two, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
chl_boxcox2 <- boxcox(chl ~ 1, data = trees) #the λ is the highest point on the curve
(chl_lambda2 <- chl_boxcox2$x[which.max(chl_boxcox2$y)]) #λ = -0.06060606
trees <- trees %>% mutate(transformed_chl2 = (chl ^ (chl_lambda2 - 1)) / chl_lambda2) #Box-Cox transformation applied in a new column

chl_mod_trans2 <- lm(transformed_chl2 ~ type, data = trees)
autoplot(chl_mod_trans2)
shapiro.test(resid(chl_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_chl2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(chl_kw2 <- kruskal.test(chl ~ type, data = trees)) #p-value = 0.001421; significant


(chl_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 'Invasive', 
                                           'C. bullatus',
                                           'R. pseudoacacia semperflorens')), #reorders the types 
                            y = chl, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop("Average chlorophyll (SPAD)"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("chl_boxplot2.jpg", chl_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Dunn post-hoc test
dunn_chl_2 <- dunn.test(trees$chl, trees$code_two, method = "bonferroni") #invasives differ significantly from natives yay
#CB and RPS differ significantly from natives and naturalised (RPS also from invasive)

#LDCM ----
ldcm_mod2 <- lm(ldcm ~ code_two, data = trees)
autoplot(ldcm_mod2)
shapiro.test(resid(ldcm_mod2)) #residuals distributed normally
bartlett.test(ldcm ~ code_two, data = trees) #homoscedascity
anova(ldcm_mod2)
#significant, p = 0.00016

(ldcm_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 'Invasive', 
                                           'C. bullatus',
                                           'R. pseudoacacia semperflorens')), #reorders the types 
                            y = ldcm, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Leaf dry matter concentration (g" ~ "g"^-1~")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("ldcm_boxplot2.jpg", ldcm_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Dunn post-hoc test
dunn_ldcm_2 <- dunn.test(trees$ldcm, trees$code_two, method = "bonferroni") #invasives differ significantly from natives yay
#neither of the alien species differ significantly from the natives


#Assimilation rate ----
a_mod2 <- lm(A ~ code_two, data = trees)
autoplot(a_mod2)
shapiro.test(resid(a_mod2)) #residuals distributed normally
bartlett.test(A ~ code_two, data = trees) #heteroscedascity (just barely)
anova(a_mod2)
#p = 0.0586

(a_boxplot2 <- ggplot(trees, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus',
                                         'R. pseudoacacia semperflorens')), #reorders the types 
                          y = A, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Assimilation rate (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("a_boxplot2.jpg", a_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Dunn post-hoc test
dunn_a_2 <- dunn.test(trees$A, trees$code_two, method = "bonferroni") #invasives differ significantly from natives yay
#neither of the alien species differ significantly from the natives


#Evapotranspiration rate ----
e_mod2 <- lm(E ~ code_two, data = trees)
autoplot(e_mod2)
shapiro.test(resid(e_mod2)) #residuals distributed normally
bartlett.test(E ~ code_two, data = trees) #homoscedascity
anova(e_mod2) #NS; p-value = 0.8722

(e_boxplot2 <- ggplot(trees, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus',
                                         'R. pseudoacacia semperflorens')), #reorders the types 
                          y = E, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Evapotranspiration rate (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("e_boxplot2.jpg", e_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 





#GH2O ----
g_mod2 <- lm(g ~ code_two, data = trees)
autoplot(g_mod2)
shapiro.test(resid(g_mod2)) #residuals not distributed normally
bartlett.test(g ~ code_two, data = trees) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
g_boxcox2 <- boxcox(g ~ 1, data = trees) #the λ is the highest point on the curve
(g_lambda2 <- g_boxcox2$x[which.max(g_boxcox2$y)]) #λ = -0.7070707
trees <- trees %>% mutate(transformed_g2 = (g ^ (g_lambda2 - 1)) / g_lambda2) #Box-Cox transformation applied in a new column

g_mod_trans2 <- lm(transformed_g2 ~ type, data = trees)
autoplot(g_mod_trans2)
shapiro.test(resid(g_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_g2 ~ type, data = trees) #homoscedascity

#Transformation did not work, moving on to non-parametric alternative:
(g_kw2 <- kruskal.test(g ~ type, data = trees)) #p-value = 0.03221; significant


(g_boxplot2 <- ggplot(trees, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus',
                                         'R. pseudoacacia semperflorens')), #reorders the types 
                          y = g, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                 "R. pseudoacacia semperflorens" = "steelblue3")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Stomatal conductance rate (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("g_boxplot2.jpg", g_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Dunn post-hoc test
dunn_g_2 <- dunn.test(trees$g, trees$code_two, method = "bonferroni") #invasives differ significantly from natives yay
#RPS is sig. different from natives; CB sig. different from native, but invasive is not (p = 0.666)





#C:N ----
cn_mod2 <- lm(c_n ~ code_two, data = cn_trees)
autoplot(cn_mod2)
shapiro.test(resid(cn_mod2)) #residuals not distributed normally
bartlett.test(c_n ~ code_two, data = cn_trees) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox2 <- boxcox(c_n ~ 1, data = cn_trees) #the λ is the highest point on the curve
(cn_lambda2 <- cn_boxcox2$x[which.max(cn_boxcox2$y)]) #λ = 0.5050505
cn_trees <- cn_trees %>% mutate(transformed_cn2 = (c_n ^ (cn_lambda2 - 1)) / cn_lambda2) #Box-Cox transformation applied in a new column

cn_mod_trans2 <- lm(transformed_cn2 ~ code_two, data = cn_trees)
autoplot(cn_mod_trans2)
shapiro.test(resid(cn_mod_trans2)) #residuals distributed normally
bartlett.test(transformed_cn2 ~ code_two, data = cn_trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kw2 <- kruskal.test(c_n ~ code_two, data = cn_trees)) #p-value = 7.488e-05; significant

(cn_boxplot2 <- ggplot(cn_trees, 
                      aes(x = factor(code_two, levels = c('Native', 'Naturalised', 'Invasive', 
                                                          'C. bullatus', 
                                                          'R. pseudoacacia semperflorens')), #reorders the types 
                          y = c_n, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                   "Naturalised" = "#EEC900", "C. bullatus" = "steelblue3",
                                   "R. pseudoacacia semperflorens" = "steelblue3")) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("C:N"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("cn_boxplot2.jpg", cn_boxplot2, path = "Plots", units = "cm", width = 20, height = 15) 

#Dunn post-hoc test
dunn_cn <- dunn.test(cn_trees$c_n, cn_trees$code_two, method = "bonferroni") #invasives differ significantly from natives
#RPS differs significantly from natives, naturalised, and invasive; CB does not (for any)




#Step 3 - Mixed effect models?? ----
#LMA LMER ----
null_lma <- lm(lma ~ 1, data = nns)
model_lma <- lmer(lma ~ type + (1 | ever_dec), data = nns)
model_lma_1 <- lmer(lma ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns) 
model_lma_2 <- lmer(lma ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_lma_3 <- lmer(lma ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_lma, model_lma, model_lma_1, model_lma_2, model_lma_3)
#models 2 and 3 fall within 2 AIC scores; so virtually identical fit to the data
#using model 2 (simplest of the two)
#so for model_lma_2; there is still 9.382 residual std that is not explained by any of these random effects

#Chl LMER ----
null_chl <- lm(chl ~ 1, data = nns)
model_chl_1 <- lmer(chl ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns)
model_chl_2 <- lmer(chl ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_chl_3 <- lmer(chl ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_chl, model_chl_1, model_chl_2, model_chl_3)
#Models 1 and 2 fall within two AIC scores, so i'll use 1 (simplest)
#3.572 residual std that isn't explained by the rest of the variables

#A LMER ----
null_a <- lm(A ~ 1, data = nns)
model_a_1 <- lmer(A ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns)
model_a_2 <- lmer(A ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_a_3 <- lmer(A ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_a, model_a_1, model_a_2, model_a_3)
#models 2 and 3 are the same; will use 2 (simplest of the two)
#1.67752 residual unexplained variance

#LDMC LMER ----
null_ldmc <- lm(ldcm ~ 1, data = nns)
model_ldmc_1 <- lmer(ldcm ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns)
model_ldmc_2 <- lmer(ldcm ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_ldmc_3 <- lmer(ldcm ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_ldmc, model_ldmc_1, model_ldmc_2, model_ldmc_3)
#models 2 and 3 are virtually identical, so will use 2 (simplest)
#0.061836 residual unexplained variance

#E LMER ----
null_e <- lm(E ~ 1, data = nns)
model_e_1 <- lmer(E ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns)
model_e_2 <- lmer(E ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_e_3 <- lmer(E ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_e, model_e_1, model_e_2, model_e_3)
#models 1, 2, and 3 are the same; will use 1 (simplest of the three)
#1.2015 residual unexplained variance

#GH2O LMER ----
null_g <- lm(g ~ 1, data = nns)
model_g_1 <- lmer(g ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = nns)
model_g_2 <- lmer(g ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = nns)
model_g_3 <- lmer(g ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = nns)
AIC(null_g, model_g_1, model_g_2, model_g_3)
#models 1 and 2 are the same; will use 1 (simplest of the two)
#11.171 residual unexplained variance


#C:N LMER ----
null_cn <- lm(c_n ~ 1, data = cn_nns)
model_cn_1 <- lmer(c_n ~ type + (1 | code) + (1 | age) +  (1 | ever_dec), data = cn_nns)
model_cn_2 <- lmer(c_n ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos), data = cn_nns)
model_cn_3 <- lmer(c_n ~ type + (1 | code) + (1 | age) +  (1 | ever_dec) + (1 | canopy_pos) + (1 | dbh), data = cn_nns)
AIC(null_cn, model_cn_1, model_cn_2, model_cn_3)
#all models are virtually the same; will use 1 (simplest of the three)
#2.304 residual unexplained variance





#NMDS----
physiological <- nns %>% select(A, E, g)
morphological <- nns %>% select(lma, ldcm)
chemical <- cn_nns %>% select(c_n)
#Morphological NMDS (nns) ----
numeric_cols_morph <- colnames(morphological)[sapply(morphological, is.numeric)] 
numeric_data_morph <- morphological[, numeric_cols_morph]
numeric_data_morph <- numeric_data_morph %>% select(lma, ldcm)

nmds_morph <- metaMDS(morphological, distance = "euclidean")
nmds_coords_morph <- as.data.frame(scores(nmds_morph, "sites"))
nmds_coords_morph$type <- nns$type

hull.data <- data.frame()
for (i in unique(nmds_coords_morph$type)) {
  temp <- nmds_coords_morph[nmds_coords_morph$type == i, ][chull(nmds_coords_morph[nmds_coords_morph$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_morph <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords_morph, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Morphological Leaf Traits by Invasion Type"))
ggsave("morphological_nmds.jpg", nmds_morph, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix_morph <- vegdist(morphological, method = "bray")
anosim(diss_matrix_morph, nns$type, permutations = 9999) 
#significant, the three types are significantly different in their morphological traits (p = 0.0357);
#however, the R value is close to 0 (0.03814), indicating a slight but significant difference between the groups

#Physiological NMDS (nns) ----
numeric_cols_phys <- colnames(physiological)[sapply(physiological, is.numeric)] 
numeric_data_phys <- physiological[, numeric_cols_phys]
numeric_data_phys <- numeric_data_phys %>% select(A, E, g)

nmds_phys <- metaMDS(physiological, distance = "euclidean")
nmds_coords_phys <- as.data.frame(scores(nmds_phys, "sites"))
nmds_coords_phys$type <- nns$type

hull.data <- data.frame()
for (i in unique(nmds_coords_phys$type)) {
  temp <- nmds_coords_phys[nmds_coords_phys$type == i, ][chull(nmds_coords_phys[nmds_coords_phys$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_phys <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords_phys, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Physiological Leaf Traits by Invasion Type"))
ggsave("physiological_nmds.jpg", nmds_phys, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix_phys <- vegdist(physiological, method = "bray")
anosim(diss_matrix_phys, nns$type, permutations = 9999) 
#not significant, the three types are not significantly different in their physiological traits (p = 0.2215);
#the R is also close to 0 (0.01142), so this relationship is insignificant and a weak difference between the groups


#Chemical NMDS (nns) ----
numeric_cols_chem <- colnames(chemical)[sapply(chemical, is.numeric)] 
numeric_data_chem <- chemical[, numeric_cols_chem]
numeric_data_chem <- numeric_data_chem %>% select(c_n)

nmds_chem <- metaMDS(chemical, distance = "euclidean")
nmds_coords_chem <- as.data.frame(scores(nmds_chem, "sites"))
nmds_coords_chem$type <- cn_nns$type

hull.data <- data.frame()
for (i in unique(nmds_coords_chem$type)) {
  temp <- nmds_coords_chem[nmds_coords_chem$type == i, ][chull(nmds_coords_chem[nmds_coords_chem$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_chem <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords_chem, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Chemical Leaf Traits by Invasion Type"))
ggsave("chemical_nmds.jpg", nmds_chem, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix_chem <- vegdist(chemical, method = "bray")
anosim(diss_matrix_chem, cn_nns$type, permutations = 9999) 
#significant, the three types are significantly different in their chemical traits (p = 8e-04);
#however, the R is close to 0 (0.117), so this relationship is not very strong but still significant







#PCA in progress (badly and sadly) ---- 
phys_subset <- nns %>% select(A, E, g, type)
morph_subset <- nns %>% select(lma, ldcm, type)
chem_subset <- cn_nns %>% select(c_n, type)
morphological_subset <- morph_subset[1:nrow(chemical), ]
physiological_subset <- phys_subset[1:nrow(chemical), ] #subsetting the physiological and morphological
#data to match the number of columns in chemical (because there were 2 samples each in chemical, as opposed to 3 for the other two)
#this is done randomly, just as the chemical traits samples were

combined_tree_data <- cbind(physiological_subset, morphological_subset)
combined_tree_data$type <- as.factor(combined_tree_data$type)
combined_tree_data_no_type <- combined_tree_data[, !names(combined_tree_data) %in% c("type")] #removes type as a variable; unnecesary for the PCA
scaled_data <- scale(combined_tree_data) #scaling; important for a PCA
(pca_result <- prcomp(scaled_data, scale. = TRUE))

#plotting the PCA
pc_scores <- as.data.frame(pca_result$x)
rotation_matrix <- as.data.frame(pca_result$rotation) #extract rotation matrix (loadings)
biplot(pca_result, scale = 0, cex = 0.7)
#Direction of Arrows: The direction of the arrows represents the relationship between the original variables and the principal components. Arrows that point in similar directions indicate positive correlations between the variables and the principal components, while arrows pointing in opposite directions indicate negative correlations.
#Length of Arrows: The length of the arrows represents the importance or weight of each variable in defining the principal components. Longer arrows indicate variables that have a greater influence on the principal components.
#Angle between Arrows: The angle between arrows indicates the correlation (or lack thereof) between the corresponding variables. Arrows that are perpendicular (at a right angle) to each other are uncorrelated with each other. In your case, if the arrows for variables A and E are at a right angle, it suggests that these variables are uncorrelated in the dataset.
#Observation Points: The points in the biplot represent individual observations (samples). Observations that are closer to each other in the biplot space are more similar in terms of their variable values. Observations that are further apart are more dissimilar.

#so i can see that g and E are negatively correlated; as are A and lma (and c_n and lma)
#A and E (and c_n and E) are not correlated (right angle); neither are A and lma; lma and g; and g and A (and g and c_n)
#45 degree angle = uncorrelated also

type_colors <- c("Native" = "#698B69", "Naturalised" = "#EEC900", "Invasive" = "#CD6090")
pca_df <- data.frame(pc_scores[, 1:2], type = combined_tree_data$type)

# Plot the PCA biplot with colored observations using ggplot2
ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = type_colors, breaks = levels(combined_tree_data$type)) +
  theme_minimal() +
  labs(x = "PC1", y = "PC2", color = "Type")




#lucy
pca_data <- nns %>% 
  select(lma, ldcm, A, E, g, chl)

pca <- princomp(pca_data, cor = TRUE, scores = TRUE)
pca$loadings #pca1 explains loads -> add the scree plot and see that the comp 1 explains a lot (plot(pca))
#lma and chl are the most important variables in the pca1; followed by g, E, and A (but the latter in the opp direction)
#add the pca1 and pca2 to the original dataset to plot 
#can go back to the models and use the pca as a response variable
#then decide how to name the groups (i.e. pca1 = morphological and chemical; pca2 = morphological and physiological)

biplot <- biplot(pca)

#workflow
#dimension reduction (pca) -> modelling -> post-hoc tests
#double check that the traits are normally distributed for the PCA (if not, then transform)
hist(nns$lma)
shapiro.test(nns$lma) #non-normal
hist(nns$ldcm)
shapiro.test(nns$ldcm) #normal
hist(nns$A)
shapiro.test(nns$A) #non-normal
hist(nns$E)
shapiro.test(nns$E) #normal
hist(nns$g)
shapiro.test(nns$g) #non-normal
hist(nns$chl)
shapiro.test(nns$chl) #non-normal



