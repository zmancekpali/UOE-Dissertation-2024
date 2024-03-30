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
library(ellipse)
library(e1071)
library(ggfortify)
library(ggpubr)
library(ggrepel)
library(ggsignif)
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
                           "CB" = "C. bullatus")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus'))) #rearranges the categories in this order

trees$age <- as.numeric(trees$age)

(trees_counts <- trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))
#1 invasive, 12 naturalised, 20 native, 1 alien


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

traits.palette <- c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900")    #defining 3 colours
traits.palette2 <- c("#CD6090", "#698B69", "#EEC900", "#5EA8D9", "#245C82", "#4A3E87", "#5A5DC7")

cn_trees <- read.csv("Data/cn_analysis.csv")
cn_trees <- cn_trees %>% 
  dplyr::mutate(canopy_pos = recode(canopy_pos, 
                             "M" = "Lower",
                             "U" = "Upper"), #recode canopy positions
         code_two = recode(code_two,
                           "CB" = "C. bullatus")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus'))) %>% #rearranges the categories in this order
  mutate(c_n = C/N)

cn_trees$age <- as.numeric(cn_trees$age)

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

(lma_kruskal <- nns %>% kruskal_test(lma ~ type)) #n = 196; df = 2
(lma_effect <- nns %>% kruskal_effsize(lma ~ type)) #effect size = 0.0679; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0679
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 6.79% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_lma1 <- nns %>% dunn_test(lma ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 

(res_kruskal_lma1 <- nns %>% kruskal_test(lma ~ type))

my_comparisons <- list(c("Native, Naturalised"), c("Native, Invasive"))

(lma_ggpubr <- ggboxplot(data = nns, x = "type", y = "lma", fill = "type", 
                         bxp.errorbar = TRUE, bxp.errorbar.width = 0.15,
                         legend = "none",
                         xlab = "\n Invasion status") +
                         labs(y = bquote("LMA (g cm"^{-2}*")")) +
    stat_compare_means(method = "kruskal.test", label.y = 160) +        # Add global anova p-value
    scale_fill_manual(values = traits.palette))

(lma_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                           y = lma, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g cm"^-2*")"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none") +
    stat_compare_means(method = "kruskal.test"))


    geom_signif(comparisons = list(c("Native", "Invasive"),
                                   c("Invasive", "Naturalised")),
                map_signif_level = TRUE,
                test = "kruskal",
                y_position = c(160, 140))

ggsave("lma_boxplot1.jpg", lma_boxplot, path = "Plots", units = "cm", width = 20, height = 15) 


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

(chl_kruskal <- nns %>% kruskal_test(chl ~ type)) #n = 196; df = 2
(chl_effect <- nns %>% kruskal_effsize(chl ~ type)) #effect size = 0.0708; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0708
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 7.08% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_chl1 <- nns %>% dunn_test(chl ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 

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
          legend.position = "none") +
    geom_signif(comparisons = list(c("Native", "Invasive"),
                                   c("Native", "Naturalised")),
                map_signif_level = TRUE,
                y_position = c(75, 70)))


(chl_boxplot <- ggplot(nns, 
                       aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                           y = lma, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g cm"^-2*")"))) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none") +
    geom_signif(comparisons = list(c("Native", "Invasive"),
                                   c("Invasive", "Naturalised")),
                map_signif_level = TRUE,
                y_position = c(160, 140)))


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


#Transiration rate ----
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
         y = expression(paste("Evapotranspiration rate (", mu, "mol H"[2]*"O" ~ "m"^-2*"s"^-1, ")"))) +    
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
         y = expression(paste("Evapotranspiration rate (", mu, "mol H"[2]*"O" ~ "m"^-2*"s"^-1, ")"))) +    
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
(lma_kw2 <- kruskal.test(lma ~ type, data = trees)) #p-value = 1.909e-05; significant

(lma_kruskal2 <- trees %>% kruskal_test(lma ~ type)) #n = 202; df = 3
(lma_effect2 <- trees %>% kruskal_effsize(lma ~ type)) #effect size = 0.109; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.109
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 10.9% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_lma2 <- trees %>% dunn_test(lma ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 

(lma_boxplot2 <- ggplot(trees, 
                       aes(x = factor(code_two, levels = 
                                        c('Native', 'Naturalised', 'Invasive', 
                                          'C. bullatus')), #reorders the types 
                           y = lma, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) +
    labs(x = "\n Invasion status", 
         y = expression(atop("LMA (g cm"^-2*")"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("lma_boxplot2.jpg", lma_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


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
(chl_kw2 <- kruskal.test(chl ~ code_two, data = trees)) #p-value = 1.587e-05; significant

(chl_kruskal2 <- trees %>% kruskal_test(chl ~ code_two)) #n = 202; df = 3
(chl_effect2 <- trees %>% kruskal_effsize(chl ~ code_two)) #effect size = 0.111; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.111
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 11.1% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_chl2 <- trees %>% dunn_test(chl ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 


(chl_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 'Invasive', 
                                           'C. bullatus')), #reorders the types 
                            y = chl, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop("Average chlorophyll (SPAD)"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("chl_boxplot2.jpg", chl_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#LDCM ----
ldcm_mod2 <- lm(ldcm ~ code_two, data = trees)
autoplot(ldcm_mod2)
shapiro.test(resid(ldcm_mod2)) #residuals distributed normally
bartlett.test(ldcm ~ code_two, data = trees) #homoscedascity
(anova_ldmc <- aov(ldcm_mod2))
#significant, p = 0.00249; df = 3; 
TukeyHSD(anova_ldmc)
#invasive and native (and invasive and naturalised) differ significantly; invasive and CB differ significantly

(ldcm_boxplot2 <- ggplot(trees, 
                        aes(x = factor(code_two, levels = 
                                         c('Native', 'Naturalised', 'Invasive', 
                                           'C. bullatus')), #reorders the types 
                            y = ldcm, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Leaf dry matter concentration (g" ~ "g"^-1~")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("ldcm_boxplot2.jpg", ldcm_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 



#Assimilation rate ----

#remove dead leaves first (two rows of data; n = 202):
trees_phys <- trees %>% 
  mutate(canopy_pos = recode(canopy_pos, 
                             "L" = "Lower",
                             "U" = "Upper")) %>%  #recode canopy positions from abbreviations
  mutate(code_two = recode(code_two,
                           "CB" = "C. bullatus")) %>% #recode alien species names
  arrange(code_two = factor(type, levels = c('Native', 'Naturalised', 'Invasive', 
                                             'C. bullatus'))) %>% #rearranges the categories in this order
  filter(A >= 0) #removed negative A values (they were dead leaves)


a_mod2 <- lm(A ~ code_two, data = trees_phys)
autoplot(a_mod2)
shapiro.test(resid(a_mod2)) #residuals distributed normally
bartlett.test(A ~ type, data = trees_phys) #homoscedascity
anova_a <- aov(a_mod2)
anova(a_mod2)#p = 0.003221

TukeyHSD(anova_a)
#invasive and CB differ significantly; invasive and native differ significantly; CB and native do not differ
#naturalised and invasive differ significantly; naturalised and native do not

(a_boxplot2 <- ggplot(trees_phys, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus')), #reorders the types 
                          y = A, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Assimilation rate (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("a_boxplot2.jpg", a_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 



#Transpiration rate ----
e_mod2 <- lm(E ~ code_two, data = trees_phys)
autoplot(e_mod2)
shapiro.test(resid(e_mod2)) #residuals distributed normally
bartlett.test(E ~ code_two, data = trees_phys) #homoscedascity
anova(e_mod2) #NS; p-value = 0.8722

(e_boxplot2 <- ggplot(trees_phys, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus',
                                         'R. pseudoacacia semperflorens')), #reorders the types 
                          y = E, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Evapotranspiration rate (", mu, "mol CO"[2]*"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("e_boxplot2.jpg", e_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 





#GH2O ----
g_mod2 <- lm(g ~ code_two, data = trees_phys)
autoplot(g_mod2)
shapiro.test(resid(g_mod2)) #residuals not distributed normally
bartlett.test(g ~ code_two, data = trees_phys) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
g_boxcox2 <- boxcox(g ~ 1, data = trees_phys) #the λ is the highest point on the curve
(g_lambda2 <- g_boxcox2$x[which.max(g_boxcox2$y)]) #λ = -0.5454545
trees_phys <- trees_phys %>% mutate(transformed_g2 = (g ^ (g_lambda2 - 1)) / g_lambda2) #Box-Cox transformation applied in a new column

g_mod_trans2 <- lm(transformed_g2 ~ type, data = trees_phys)
autoplot(g_mod_trans2)
shapiro.test(resid(g_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_g2 ~ type, data = trees_phys) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(g_kw2 <- kruskal.test(g ~ code_two, data = trees_phys)) #p-value = 0.0003655; significant

(g_kruskal2 <- trees_phys %>% kruskal_test(g ~ code_two)) #n = 202; df = 3
(g_effect2 <- trees_phys %>% kruskal_effsize(g ~ code_two)) #effect size = 0.0777; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0777
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 7.77% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_g2 <- trees_phys %>% dunn_test(g ~ code_two, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 


(g_boxplot2 <- ggplot(trees_phys, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus')), #reorders the types 
                          y = g, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(atop(paste("Stomatal conductance rate (", mu, "mol CO"[2]*"m"^-2*~"s"^-1, ")")))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("g_boxplot2.jpg", g_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 




#C:N ----
cn_mod2 <- lm(c_n ~ code_two, data = cn_trees)
autoplot(cn_mod2)
shapiro.test(resid(cn_mod2)) #residuals not distributed normally
bartlett.test(c_n ~ code_two, data = cn_trees) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox2 <- boxcox(c_n ~ 1, data = cn_trees) #the λ is the highest point on the curve
(cn_lambda2 <- cn_boxcox2$x[which.max(cn_boxcox2$y)]) #λ = -0.7474747
cn_trees <- cn_trees %>% mutate(transformed_cn2 = (c_n ^ (cn_lambda2 - 1)) / cn_lambda2) #Box-Cox transformation applied in a new column

cn_mod_trans2 <- lm(transformed_cn2 ~ code_two, data = cn_trees)
autoplot(cn_mod_trans2)
shapiro.test(resid(cn_mod_trans2)) #residuals distributed normally
bartlett.test(transformed_cn2 ~ code_two, data = cn_trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kw2 <- kruskal.test(c_n ~ code_two, data = cn_trees)) #p-value = 0.003251; significant

(cn_kruskal2 <- cn_trees %>% kruskal_test(c_n ~ code_two)) #n = 99; df = 3
(cn_effect2 <- cn_trees %>% kruskal_effsize(c_n ~ code_two)) #effect size = 0.113; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.113
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 11.3% of the variance in LMA
#M. T. Tomczak and Tomczak 2014: Tomczak, Maciej T., and Ewa Tomczak. 2014. “The Need to 
#Report Effect Size Estimates Revisited. an Overview of Some Recommended Measures of Effect Size.” 
#Trends in SportSciences

#Dunn post-hoc test
(dunn_cn2 <- cn_trees %>% dunn_test(c_n ~ code_two, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type"))  
#invasive and native differ significantly, as do CB and invasive

(cn_boxplot2 <- ggplot(cn_trees, 
                      aes(x = factor(code_two, levels = c('Native', 'Naturalised', 'Invasive', 
                                                          'C. bullatus')), #reorders the types 
                          y = c_n, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                   "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + 
    labs(x = "\n Invasion status", 
         y = expression(atop("C:N"))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("cn_boxplot2.jpg", cn_boxplot2, path = "Plots", units = "cm", width = 25, height = 12) 


#Step 3 - GLMs for random effects ----
#LMA LM ----
hist(trees$lma) #looks normal-ish; slightly skewwed
shapiro_test(trees$lma) #non-normal

#trying a log-transformation
hist(log(trees$lma)) #normal
shapiro_test(log(trees$lma)) #normal

lma_lm_null <- lm(log(lma) ~ 1, data = trees)
lma_lm1 <- lm(log(lma) ~ type, data = trees)
plot(lma_lm1)
ad.test(residuals(lma_lm1)) #p = 0.14; normally distributed errors
lma_lm2 <- lm(log(lma) ~ type + age + dbh + ever_dec + code + canopy_pos, data = trees)
plot(lma_lm2)
ad.test(residuals(lma_lm2)) #p = 0.08; normally distributed errors
#both of these are ok for the data

#see which model is most parismoius
AIC(lma_lm_null, lma_lm1, lma_lm2) #lma_lm2 is most parsimonious
summary(lma_lm2)
#we can see a significant effect of age, dbh, light, and decidousness on LMA
#as well as some interspecific variation
(interspecific_lma <- ggplot(trees, aes(x = code, y = lma, fill = type)) +
    geom_boxplot() +
    theme_classic() +
    ylab(bquote("LMA (g cm"^{-2}*")")) +
    xlab("Species") +
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "Alien" = "#5EA8D9",
                                 name = "Invasion type")) +  #colours each boxplot this particular colour)
    theme(legend.position = c(0.95, 0.95),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          legend.key.size = unit(1.5, "line"),
          legend.title = element_blank()))

ggsave("interspecific_lma_variation.jpg", interspecific_lma, path = "Plots", 
       width = 30, height = 12)


#Chl GLM ----
hist(trees$chl) #looks normal
shapiro_test(trees$chl) #non-normal

#trying a log-transformation
hist(log(trees$chl)) #looks normal
shapiro_test(log(trees$chl)) #non-normal

#inverse
hist(1/trees$chl) #looks normal
shapiro_test(1/trees$chl) #non-normal

#sqrt
hist(sqrt(trees$chl)) #looks normal
shapiro_test(sqrt(trees$chl)) #non-normal

chl_lm_null_reg <- lm(chl ~ 1, data = trees)
chl_lm1_reg <- lm(chl ~ type, data = trees)
plot(chl_lm1_reg)
ad.test(residuals(chl_lm1_reg)) #p = 0.0004; normally distributed errors

chl_lm_null <- lm(log(chl) ~ 1, data = trees)
chl_lm1 <- lm(log(chl) ~ type, data = trees)
plot(chl_lm1)
ad.test(residuals(chl_lm1)) #p = 0.008; non-normally distributed errors

chl_lm_null_sqrt <- lm(sqrt(chl) ~ 1, data = trees)
chl_lm1_sqrt <- lm(sqrt(chl) ~ type, data = trees)
plot(chl_lm1_sqrt)
ad.test(residuals(chl_lm1_sqrt)) #p = 0.006; normally distributed errors

chl_lm2 <- lm(chl ~ type + age + ever_dec + canopy_pos, data = trees)
plot(chl_lm2)
ad.test(residuals(chl_lm2)) #p = 0.27; normally distributed errors

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

#A GLM---

#LDMC GLM--

#E GLM--

#GH2O GLM---

#C:N GLM---







#NMDS----
morphological <- nns %>% select(lma, ldcm)
physiological <- nns %>% select(A, E, g)

merged_data <- merge(cn_trees, nns, by = c("code", "canopy_pos"))
chemical <- merged_data %>% select(c_n, chl)
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
#significant, the three types are significantly different in their morphological traits (p = 0.0362);
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
#not significant, the three types are not significantly different in their physiological traits (p = 0.231);
#the R is also close to 0 (0.01142), so this relationship is insignificant and a weak difference between the groups


#Chemical NMDS (nns) ----
numeric_cols_chem <- colnames(chemical)[sapply(chemical, is.numeric)] 
numeric_data_chem <- chemical[, numeric_cols_chem]
numeric_data_chem <- numeric_data_chem %>% select(c_n)

nmds_chem <- metaMDS(chemical, distance = "euclidean")
nmds_coords_chem <- as.data.frame(scores(nmds_chem, "sites"))
nmds_coords_chem$type <- merged_data$code_two.x

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
anosim(diss_matrix_chem, merged_data$code_two.x, permutations = 9999) 
#significant, the three types are significantly different in their chemical traits (p = 4e-04);
#however, the R is close to 0 (0.08338), so this relationship is not very strong but still significant








#PCA in progress (badly and sadly) ---- 
merged_trees <- merge(trees, cn_trees[, c("code", "canopy_pos", "c_n")], by = c("code", "canopy_pos"))
combined_tree_data_numeric <- merged_trees %>%
  select(chl, lma, ldcm, A, E, g, c_n)
combined_tree_data$type_numeric <- as.numeric(combined_tree_data$type)

merged_trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code))


scaled_data <- scale(combined_tree_data_numeric)

(pca_result <- prcomp(scaled_data, scale. = TRUE))
pca_df <- as.data.frame(pca_result$x)
pca_df$type <- merged_trees$type  # Assuming 'type' is a column in combined_tree_data
type_colors <- c("Native" = "#698B69", "Naturalised" = "#EEC900", "Invasive" = "#CD6090", "Alien" = "#5EA8D9")
(plot <- biplot(pca_result, scale = 0, cex = 0.7))

loadings <- pca_result$rotation
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

arrows_df <- data.frame(
  PC1 = loadings[, 1] * sqrt(variance_explained[1]),
  PC2 = loadings[, 2] * sqrt(variance_explained[2])
)

(pca_plot <- ggbiplot(pca_df, obs.scale = 1, var.scale = 1) +
    ggtitle("PCA Biplot") +
    scale_color_manual(values = type_colors) +
    geom_point(aes(color = type), size = 3) +  # Color points by 'type'
    theme_classic() +
    theme(legend.position = c(0.9, 0.2),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.key.size = unit(1.5, "line"),
          legend.title = element_blank()))
  
plot <- ggbiplot(pca_result, choices = c(1, 2), scale = 0, labels = rownames(pca_result$rotation))

pca_result <- prcomp(combined_tree_data_scaled, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$type <- combined_tree_data$type  # Assuming 'type' is a column in combined_tree_data
type_colors <- c("Native" = "#698B69", "Naturalised" = "#EEC900", "Invasive" = "#CD6090")
biplot(pca_result, scale = 0, cex = 0.7, col = type_colors[pca_df$type])


#plotting the PCA
pc_scores <- as.data.frame(pca_result$x)
rotation_matrix <- as.data.frame(pca_result$rotation) #extract rotation matrix (loadings)
biplot(pca_result, scale = 0, cex = 0.7, col = type)
#Direction of Arrows: The direction of the arrows represents the relationship between the original variables and the principal components. Arrows that point in similar directions indicate positive correlations between the variables and the principal components, while arrows pointing in opposite directions indicate negative correlations.
#Length of Arrows: The length of the arrows represents the importance or weight of each variable in defining the principal components. Longer arrows indicate variables that have a greater influence on the principal components.
#Angle between Arrows: The angle between arrows indicates the correlation (or lack thereof) between the corresponding variables. Arrows that are perpendicular (at a right angle) to each other are uncorrelated with each other. In your case, if the arrows for variables A and E are at a right angle, it suggests that these variables are uncorrelated in the dataset.
#Observation Points: The points in the biplot represent individual observations (samples). Observations that are closer to each other in the biplot space are more similar in terms of their variable values. Observations that are further apart are more dissimilar.

#so i can see that g and E are negatively correlated; as are A and lma (and c_n and lma)
#A and E (and c_n and E) are not correlated (right angle); neither are A and lma; lma and g; and g and A (and g and c_n)
#45 degree angle = uncorrelated also

type_colors <- c("Native" = "#698B69", "Naturalised" = "#EEC900", "Invasive" = "#CD6090")
pca_df <- data.frame(pc_scores[, 1:2], type = combined_tree_data$type)
pca_loadings <- pca_result$rotation[, 1:2]
loadings_df <- data.frame(PC1 = rep(0, nrow(pca_loadings)),
                          PC2 = rep(0, nrow(pca_loadings)),
                          xend = pca_loadings[, 1],
                          yend = pca_loadings[, 2])
ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3) +
  geom_segment(data = loadings_df,
               aes(x = PC1, y = PC2, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "black", alpha = 0.5) +
  scale_color_manual(values = type_colors, breaks = levels(combined_tree_data$type)) +
  theme_classic() +
  labs(x = "PC1", y = "PC2", color = "Type")

# Plot the PCA biplot with colored observations and PCA arrows
ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3) +
  geom_segment(x = 0, y = 0, xend = pca_loadings[, 1], yend = pca_loadings[, 2],
               arrow = arrow(length = unit(0.2, "inches")), color = "black", alpha = 0.5) +
  scale_color_manual(values = type_colors, breaks = levels(combined_tree_data$type)) +
  theme_classic() +
  labs(x = "PC1", y = "PC2", color = "Type")

# Plot the PCA biplot with colored observations using ggplot2
ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = type_colors, breaks = levels(combined_tree_data$type)) +
  theme_classic() +
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



