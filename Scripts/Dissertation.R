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
library(plotly)
library(plot3D)
library(rstatix)
library(sjPlot)
library(stats)
library(tidyverse)
library(vegan)
library(vegan3d)

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
                                             'C. bullatus'))) #rearranges the categories in this order

trees_pos <- trees %>% 
  filter(A >= 0) #for physiological measurements

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
  arrange(type = factor(type, levels = c('Native', 'Naturalised', 'Invasive'))) #rearranges the categories in this order
#  filter(A >= 0) #removed negative A values (they were dead leaves)

nns_pos <- nns %>% 
  filter(A >= 0) %>% 
  rename("LMA" = "lma") %>% 
  rename("LDMC" = "ldcm") %>% 
  rename("Chl" = "chl") %>% 

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
                           "CB" = "C. bullatus", 
                           "R. ponticum" = "Invasive")) %>% #recode alien species names
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

cn_nmds_data <- cn_nns %>% 
  rename("CN" = "c_n")

(cn_counts <- cn_trees %>%
    group_by(type) %>%
    summarise(unique_species = n_distinct(code)))

#Mean trait values for each group - for invasion index
(means_trees <- trees_pos %>% 
  group_by(type) %>% 
  summarise(mean_lma = mean(lma),
            mean_ldmc = mean(ldcm), 
            mean_chl = mean(chl),
            mean_A = mean(A),
            mean_E = mean(E),
            mean_g = mean(g)))

(means_cn_trees <- cn_trees %>% 
  group_by(type) %>% 
  summarise(mean_cn = mean(C/N)))

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

lma_mod_log <- lm(log(lma) ~ type, data = nns)
autoplot(lma_mod_log)
shapiro.test(resid(lma_mod_log)) #residuals not distributed normally
bartlett.test(log(lma) ~ type, data = nns) #heteroscedascity

lma_mod_inv <- lm(1/(lma) ~ type, data = nns)
autoplot(lma_mod_inv)
shapiro.test(resid(lma_mod_inv)) #residuals not distributed normally
bartlett.test(1/(lma) ~ type, data = nns) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox <- boxcox(lma ~ 1, data = nns) #the λ is the highest point on the curve
(lma_lambda <- lma_boxcox$x[which.max(lma_boxcox$y)]) #λ = -0.1818182
nns <- nns %>% mutate(transformed_lma = (lma ^ (lma_lambda - 1)) / lma_lambda) #Box-Cox transformation applied in a new column

lma_mod_trans <- lm(transformed_lma ~ type, data = nns)
autoplot(lma_mod_trans)
shapiro.test(resid(lma_mod_trans)) #residuals not distributed normally
bartlett.test(transformed_lma ~ type, data = nns) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kruskal <- nns %>% kruskal_test(lma ~ type)) #n = 198; df = 2; p = 0.000542
(lma_effect <- nns %>% kruskal_effsize(lma ~ type)) #effect size = 0.0669; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0679
#this value indicates the % of variance in the dependent variable (lma) explained by the invasion status
#so this explains 6.69% of the variance in LMA

#Dunn post-hoc test
(dunn_lma1 <- nns %>% dunn_test(lma ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 

my_comparisons <- list(c("Native, Naturalised"), c("Native, Invasive"))

(lma_ggpubr <- ggboxplot(data = nns, x = "type", y = "lma", fill = "type", 
                         bxp.errorbar = TRUE, bxp.errorbar.width = 0.15,
                         legend = "none",
                         xlab = "\n Invasion status") +
                         labs(y = bquote("LMA (g cm"^{-2}*")")) +
    stat_compare_means(method = "kruskal.test", label.y = 160) + #Add global anova p-value
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
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")))

ggsave("lma_boxplot1.jpg", lma_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 


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
(chl_kruskal <- nns %>% kruskal_test(chl ~ type)) #n = 198; df = 2; p = 0.00037 
(chl_effect <- nns %>% kruskal_effsize(chl ~ type)) #effect size = 0.0708; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0708
#this value indicates the % of variance in the dependent variable (chl) explained by the invasion status
#so this explains 7.08% of the variance in chl


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
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
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
ldcm_mod <- lm(ldcm ~ type, data = nns)
autoplot(ldcm_mod)
shapiro.test(resid(ldcm_mod)) #residuals distributed normally
bartlett.test(ldcm ~ type, data = nns) #homoscedascity
anova(ldcm_mod) #p = 0.001115; df = 2

tukey_ldcm <- aov(ldcm_mod)
TukeyHSD(tukey_ldcm, conf.level = 0.95) #native/invasive and naturalised/invasive differ significantly


(ldmc_boxplot <- ggplot(nns, 
                        aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                            y = ldcm, fill = type)) + 
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

#no Tukey test bc no significant differences found overall

(e_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
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
dr_mod <- lm(Dark_resp ~ type, data = nns_pos)
autoplot(dr_mod)
shapiro.test(resid(dr_mod)) #residuals not distributed normally
bartlett.test(Dark_resp ~ type, data = nns_pos) #heteroscedascity

#cannot do a boxcox transformation bc values are negative       

#Transformation did not work, moving on to non-parametric alternative:
(dr_kruskal <- nns_pos %>% kruskal_test(Dark_resp ~ type)) #n = 196; df = 2; p = 0.288  
(dr_effect <- nns_pos %>% kruskal_effsize(Dark_resp ~ type)) #effect size = 0.00253; small magnitude
#report as: small effect size is detected, eta2[H] = 0.00253

(dr_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = Dark_resp, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(paste("Dark respiration (", mu, "mol CO"[2] ~ "m"^-2*"s"^-1, ")"))) +    
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
(g_effect <- nns_pos %>% kruskal_effsize(g ~ type)) #effect size = 0.0276 ; small magnitude
#report as: small effect size is detected, eta2[H] = 0.0276 
#this value indicates the % of variance in the dependent variable (g) explained by the invasion status
#so this explains 7.08% of the variance in g

#Dunn post-hoc test
(dunn_g1 <- nns_pos %>% dunn_test(g ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #invasive/native and invsive/naturalised differ sig.

(g_boxplot <- ggplot(nns_pos, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = g, fill = type)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = traits.palette) + 
    labs(x = "\n Invasion status", 
         y = expression(paste("g (", mu, "mol H"[2]*"O" ~ "m"^-2*"s"^-1, ")"))) +    
    theme_classic() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.position = "none"))

ggsave("g_boxplot.jpg", g_boxplot, path = "Plots", units = "cm", 
       width = 20, height = 15) 


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
(cn_kruskal <- cn_nns %>% kruskal_test(c_n ~ type)) #n = 96; df = 2; p = 0.00108   
(cn_effect <- cn_nns %>% kruskal_effsize(c_n ~ type)) #effect size = 0.125  ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.125  
#this value indicates the % of variance in the dependent variable (chl) explained by the invasion status
#so this explains 7.08% of the variance in chl

#Dunn post-hoc test
(dunn_cn1 <- cn_nns %>% dunn_test(c_n ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) #invasive/native differ sig.

(cn_boxplot <- ggplot(cn_nns, 
                     aes(x = factor(type, levels = c('Native', 'Naturalised', 'Invasive')), #reorders the types 
                         y = c_n, fill = type)) + 
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


#Step 2: compare alien species ----
#LMA ----
lma_mod2 <- lm(lma ~ type, data = trees)
autoplot(lma_mod2)
shapiro.test(resid(lma_mod2)) #residuals not distributed normally
bartlett.test(lma ~ type, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
lma_boxcox2 <- boxcox(lma ~ 1, data = trees) #the λ is the highest point on the curve
(lma_lambda2 <- lma_boxcox2$x[which.max(lma_boxcox2$y)]) #λ = -0.1414141
trees <- trees %>% mutate(transformed_lma2 = (lma ^ (lma_lambda2 - 1)) / lma_lambda2) #Box-Cox transformation applied in a new column

lma_mod_trans2 <- lm(transformed_lma2 ~ type, data = trees)
autoplot(lma_mod_trans2)
shapiro.test(resid(lma_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_lma2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(lma_kruskal2 <- trees %>% kruskal_test(lma ~ type)) #n = 204; df = 3; p = 0.000019 
(lma_effect2 <- trees %>% kruskal_effsize(lma ~ type)) #effect size = 0.108 ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.108 

#Dunn post-hoc test
(dunn_lma2 <- trees %>% dunn_test(lma ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly
#invasive/native and invasive/naturalised differ significantly


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
         y = expression("LMA (g cm"^-2*")")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("lma_boxplot2.jpg", lma_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 


#Average chlorophyll ----
chl_mod2 <- lm(chl ~ type, data = trees)
autoplot(chl_mod2)
shapiro.test(resid(chl_mod2)) #residuals not distributed normally
bartlett.test(chl ~ type, data = trees) #heteroscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
chl_boxcox2 <- boxcox(chl ~ 1, data = trees) #the λ is the highest point on the curve
(chl_lambda2 <- chl_boxcox2$x[which.max(chl_boxcox2$y)]) #λ = -0.06060606
trees <- trees %>% mutate(transformed_chl2 = (chl ^ (chl_lambda2 - 1)) / chl_lambda2) #Box-Cox transformation applied in a new column

chl_mod_trans2 <- lm(transformed_chl2 ~ type, data = trees)
autoplot(chl_mod_trans2)
shapiro.test(resid(chl_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_chl2 ~ type, data = trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(chl_kruskal2 <- trees %>% kruskal_test(chl ~ type)) #n = 204; df = 3; p = 0.000014 
(chl_effect2 <- trees %>% kruskal_effsize(chl ~ type)) #effect size = 0.111; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.111

#Dunn post-hoc test
(dunn_chl2 <- trees %>% dunn_test(chl ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly
#invasive/native and invasive/naturalised differ significantly

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
         y = expression("Average chlorophyll (SPAD)"),
         fill = "Invasion status") + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm"),
          legend.position = "none"))

ggsave("chl_boxplot2.jpg", chl_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 


#LDCM ----
ldcm_mod2 <- lm(ldcm ~ type, data = trees)
autoplot(ldcm_mod2)
shapiro.test(resid(ldcm_mod2)) #residuals distributed normally
bartlett.test(ldcm ~ type, data = trees) #homoscedascity
anova(ldcm_mod2) #significant, p = 0.002177; df = 3; 

anova_ldmc <- aov(ldcm_mod2)
TukeyHSD(anova_ldmc)
#invasive/alien differ significantly
#native/invasive and naturalised/invasive differ significantly

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
         y = expression(paste("LDMC (g" ~ "g"^-1~")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
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
anova(a_mod2) #p = 0.003221

anova_a <- aov(a_mod2)
TukeyHSD(anova_a)
#invasive/alien differ significantly
#native/invasive differ significantly
#naturalised/invasive differ significantly


(a_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus')), #reorders the types 
                          y = A, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("A (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
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
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus',
                                         'R. pseudoacacia semperflorens')), #reorders the types 
                          y = E, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("E (", mu, "mol H"[2]*"O m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("e_boxplot2.jpg", e_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 





#Dark respiration ----
dr_mod2 <- lm(Dark_resp ~ code_two, data = trees_pos)
autoplot(dr_mod2)
shapiro.test(resid(dr_mod2)) #residuals not distributed normally
bartlett.test(Dark_resp ~ type, data = trees_pos) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(dr_kruskal2 <- trees_pos %>% kruskal_test(Dark_resp ~ type)) #n = 202; df = 3; p = 0.437  
(dr_effect2 <- trees_pos %>% kruskal_effsize(Dark_resp ~ type)) #effect size = -0.00142; small magnitude
#report as: small effect size is detected, eta2[H] = -0.00142

(dr_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus')), #reorders the types 
                          y = Dark_resp, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("Dark respiration (", mu, "mol CO"[2]~"m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
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
trees_pos <- trees_pos %>% mutate(transformed_g2 = (g ^ (g_lambda2 - 1)) / g_lambda2) #Box-Cox transformation applied in a new column

g_mod_trans2 <- lm(transformed_g2 ~ type, data = trees_pos)
autoplot(g_mod_trans2)
shapiro.test(resid(g_mod_trans2)) #residuals not distributed normally
bartlett.test(transformed_g2 ~ type, data = trees_pos) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(g_kw2 <- kruskal.test(g ~ code_two, data = trees_pos)) #p-value = 0.0003655; significant

(g_kruskal2 <- trees_pos %>% kruskal_test(g ~ code_two)) #n = 202; df = 3
(g_effect2 <- trees_pos %>% kruskal_effsize(g ~ code_two)) #effect size = 0.0777; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.0777

#Dunn post-hoc test
(dunn_g2 <- trees_pos %>% dunn_test(g ~ code_two, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type")) 
#alien/native and alien/naturalised differ significantly from one another
#invasive/naturalised differ significantly from one another

(g_boxplot2 <- ggplot(trees_pos, 
                      aes(x = factor(code_two, levels = 
                                       c('Native', 'Naturalised', 'Invasive', 
                                         'C. bullatus')), #reorders the types 
                          y = g, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                 "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + #colours each boxplot this particular colour
    labs(x = "\n Invasion status", 
         y = expression(paste("g (", mu, "mol H"[2]*"O m"^-2*~"s"^-1, ")"))) +
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("g_boxplot2.jpg", g_boxplot2, path = "Plots", units = "cm", 
       width = 25, height = 13) 




#C:N ----
cn_mod2 <- lm(c_n ~ type, data = cn_trees)
autoplot(cn_mod2)
shapiro.test(resid(cn_mod2)) #residuals not distributed normally
bartlett.test(c_n ~ type, data = cn_trees) #homoscedascity

#Attempt mathematical transformation first to meet ANOVA assumptions:
cn_boxcox2 <- boxcox(c_n ~ 1, data = cn_trees) #the λ is the highest point on the curve
(cn_lambda2 <- cn_boxcox2$x[which.max(cn_boxcox2$y)]) #λ = -0.7474747
cn_trees <- cn_trees %>% mutate(transformed_cn2 = (c_n ^ (cn_lambda2 - 1)) / cn_lambda2) #Box-Cox transformation applied in a new column

cn_mod_trans2 <- lm(transformed_cn2 ~ type, data = cn_trees)
autoplot(cn_mod_trans2)
shapiro.test(resid(cn_mod_trans2)) #residuals distributed normally
bartlett.test(transformed_cn2 ~ type, data = cn_trees) #heteroscedascity

#Transformation did not work, moving on to non-parametric alternative:
(cn_kruskal2 <- cn_trees %>% kruskal_test(c_n ~ type)) #n = 99; df = 3; p = 0.00325 
(cn_effect2 <- cn_trees %>% kruskal_effsize(c_n ~ type)) #effect size = 0.126 ; moderate magnitude
#report as: moderate effect size is detected, eta2[H] = 0.126 

#Dunn post-hoc test
(dunn_cn2 <- cn_trees %>% dunn_test(c_n ~ type, p.adjust.method = "bonferroni") %>% 
    add_xy_position(x = "type"))  
#invasive/native differ significantly
#invasive/alien differ significantly

(cn_boxplot2 <- ggplot(cn_trees, 
                      aes(x = factor(code_two, levels = c('Native', 'Naturalised', 'Invasive', 
                                                          'C. bullatus')), #reorders the types 
                          y = c_n, fill = code_two)) + 
    geom_boxplot() + #creates the boxplot
    stat_boxplot(geom ='errorbar', width = 0.3) + #adds the whisker ends
    scale_fill_manual(values = c("Invasive" = "#CD6090", "Native" = "#698B69",
                                   "Naturalised" = "#EEC900", "C. bullatus" = "#5EA8D9")) + 
    labs(x = "\n Invasion status", 
         y = expression("C/N ratio")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "italic", "italic")),  # Italicize selected names
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 11), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm")) +
    guides(fill = FALSE))

ggsave("cn_boxplot2.jpg", cn_boxplot2, path = "Plots", units = "cm",
       width = 25, height = 13) 


#Boxplots grid ----
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

(grid <- grid.arrange(e_boxplot2, a_boxplot2, g_boxplot2, dr_boxplot2,
                      ldcm_boxplot2, lma_boxplot2, cn_boxplot2, chl_boxplot2, 
                      ncol = 2, widths = c(1, 1), heights = c(1, 1, 1, 1.3)))

ggsave("boxplot_grid.jpg", grid, path = "Plots", units = "cm",
       width = 20, height = 29.7)

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







#NMDS (nns) ----
merged_trees_nns <- merge(nns_pos, cn_nmds_data[, c("code", "canopy_pos", "CN")], by = c("code", "canopy_pos"))

numeric_cols_nns <- colnames(merged_trees_nns)[sapply(merged_trees_nns, is.numeric)] 
numeric_data_nns <- merged_trees_nns[, numeric_cols_nns]
numeric_data_nns <- numeric_data_nns %>% select(Chl, LMA, LDMC, A, E, g, CN, Dark_resp)

#finding the lowest stress for up to 6 dimensions:
dimcheckMDS(numeric_data_nns,
            distance = "euclidean",
            k = 6) #goeveg package
#generally accepted that stress < 0.2 is a fair fit for ordination, so will use 
#2 dimensions here 

#2-dimensional NMDS (nns)----
nmds_nns <- metaMDS(numeric_data_nns, distance = "euclidean") #2 dimensions, stress = 0.2114154
nmds_coords_nns <- as.data.frame(scores(nmds_nns, "sites"))
nmds_coords_nns$type <- merged_trees_nns$type

plot(nmds_nns, type = "t") #base r NMDS plot
stressplot(nmds_nns) #stressplot; linear R^2 = 0.98; non-linear R^2 = 0.995

gof_nns <- goodness(nmds_nns)
plot(nmds_nns, typ = "t", main = "goodness of fit")
points(nmds_nns, display = "sites", cex = gof_nns * 100)

fit_nns <- envfit(nmds_nns, numeric_data_nns, permu = 999) #environmental variables

plot(nmds_nns, display = "sites")
plot(fit_nns, p.max = 0.05) #The arrow(s) point to the direction of most rapid change in 
#the environmental variable. Often this is called the direction of the gradient.
#The length of the arrow(s) is proportional to the correlation between ordination and environmental variable. 
#often this is called the strength of the gradient.

(stress_nns <- nmds_nns$stress) #0.0731529
goodness(object = nmds_nns)
plot(nmds_nns$diss, nmds_nns$dist)

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
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Leaf Traits by Invasion Type"))
ggsave("nmds_nns_2d.jpg", nmds_nns_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 


diss_matrix_nns <- vegdist(numeric_data_nns, method = "euclidean")
numeric_data_nns$type <- merged_trees_nns$type
anosim(diss_matrix_nns, numeric_data_nns$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.0949), indicating a slight but significant difference between the groups

#drivers of the variation:
nmds_nns <- metaMDS(numeric_data_nns, distance = "euclidean", k = 2) #2 dimensions, stress = 0.2114154
en = envfit(nmds_nns, numeric_data_nns, permutations = 999, na.rm = TRUE)

plot(nmds_nns)
plot(en)
nmds_coords_nns <- as.data.frame(scores(nmds_nns, "sites"))
nmds_coords_nns$type <- merged_trees_nns$type

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
        legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = row.names(en_coord_cont)) +
  labs(colour = "type") +
  labs(title = "NMDS Plot of Leaf Traits by Invasion Type"))


ggsave("drivers_nmds_nns_2d.jpg", drivers_nmds, path = "Plots", units = "cm", 
       width = 20, height = 20)

#importance of each leaf trait:
nmds_nns <- metaMDS(numeric_data_nns, distance = "euclidean", k = 2) #2 dimensions, stress = 0.2114154
(en = envfit(nmds_nns, numeric_data_nns, permutations = 999, na.rm = TRUE))
#lma, A, and E have high r2 values (close to 1), indicating strong correlations with both NMDS1 and NMDS2.
  #This suggests that these traits play a significant role in differentiating between invasive and native species in terms of their leaf characteristics.
  #Additionally, the significant p-values (***) indicate that these correlations are statistically robust.
#chl and c_n also show significant correlations with the NMDS axes, but to a lesser extent compared to lma, A, and E.
#ldcm and g have relatively low r2 values and (ldcm is non-significant), suggesting weaker correlations with the NMDS axes.

#3-dimensional NMDS (nns) ----
nmds_nns3 <- metaMDS(numeric_data_nns, distance = "euclidean", k = 3) #3 dimensions, stress = 0.2114154
nmds_coords_nns3 <- as.data.frame(scores(nmds_nns3, "sites"))
nmds_coords_nns3$type <- merged_trees_nns$type

plot(nmds_nns3, type = "t") #base r NMDS plot
stressplot(nmds_nns3) #stressplot; linear R^2 = 0.893; non-linear R^2 = 0.981

gof_nns3 <- goodness(nmds_nns3)
plot(nmds_nns3, typ = "t", main = "goodness of fit")
points(nmds_nns3, display = "sites", cex = gof_nns3 * 100)

(stress_nns3 <- nmds_nns3$stress) #0.1382597
goodness(object = nmds_nns3)
plot(nmds_nns3$diss, nmds_nns3$dist)

#ggplot NMDS
hull.data3 <- data.frame()
for (i in unique(nmds_coords_nns3$type)) {
  temp <- nmds_coords_nns3[nmds_coords_nns3$type == i, ][chull(nmds_coords_nns3[nmds_coords_nns3$type == i, c("NMDS1", "NMDS2", "NMDS3")]), ]
  hull.data3 <- rbind(hull.data3, temp)
}

(nmds_nns3_plot <- ggplot() +
    geom_polygon(data = hull.data3[hull.data3$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.5) + # polygons for non-invasive types
    geom_polygon(data = hull.data3[hull.data3$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.8) + # polygons for invasive type
    geom_point(data = nmds_coords_nns3, aes(x = NMDS1, y = NMDS2, z = NMDS3, color = type), size = 3) + # points for NMDS1, NMDS2, and NMDS3
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Leaf Traits by Invasion Type"))
ggsave("nmds_nns_3D.jpg", nmds_nns3_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix_nns3 <- vegdist(numeric_data_nns, method = "euclidean", k = 3)
anosim(diss_matrix_nns3, merged_trees_nns$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.0949), indicating a slight but significant difference between the groups
#same as above

(drivers_nmds3 <- ggplot(data = nmds_coords_nns3, aes(x = NMDS1, y = NMDS2, z = NMDS3)) +
    geom_point(aes(colour = type), size = 3) +
    geom_polygon(data = hull.data3[hull.data3$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.5) + # polygons for non-invasive types
    geom_polygon(data = hull.data3[hull.data3$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.8) + # polygons for invasive type
    scale_colour_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) + 
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900")) + 
    theme_classic() +
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.text = element_text(size = 9, colour = "grey30"),
          legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(colour = "type") +
    labs(title = "NMDS Plot of Leaf Traits by Invasion Type (3D)")) #literally identical to the previous graph
  
#geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.5) +
#geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.8) +
#geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2, z = NMDS3, label = row.names(en_coord_cont)), 
#colour = "grey30", fontface = "bold") +


traits.palette <- c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900")    #defining 3 colours

plot_ly(data = nmds_coords_nns3, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
        type = "scatter3d", mode = "markers", color = ~type, colors = traits.palette)

#3d PLOT attempt (you HAVE TO run this in order)
(ef = envfit(nmds_nns3, numeric_data_nns, permutations = 999, na.rm = TRUE, choices = 1:3))

nmds_coords_nns3 <- as.data.frame(scores(nmds_nns3, "sites"))
out <- ordiplot3d(nmds_coords_nns3, col = as.factor(nmds_coords_nns3$type), pch = 18,
                  ax.col = "red", envfit = ef, arr.len = 0.1)
nmds_coords_nns3$type <- merged_trees_nns$type
points(out$xyz.convert(nmds_coords_nns3), pch = 20,
       col = c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900"))
nmds_coords_nns3 <- as.data.frame(scores(nmds_nns3, "sites"))
nmdsD <- scores(nmds_nns3, choices = 1:3, display = c("species"))
x_coords <- out$xyz.convert(nmdsD)$x
y_coords <- out$xyz.convert(nmdsD)$y
z_coords <- out$xyz.convert(nmdsD)$z
text(x = x_coords, y = y_coords, z = z_coords, labels = rownames(nmdsD), pos = 2, col = "blue")



#NMDS (with alien) ----
merged_trees <- merge(trees_pos, cn_trees[, c("code", "canopy_pos", "c_n")], by = c("code", "canopy_pos"))

numeric_cols <- colnames(merged_trees)[sapply(merged_trees, is.numeric)] 
numeric_data <- merged_trees[, numeric_cols]
numeric_data <- numeric_data %>% select(chl, lma, ldcm, A, E, g, c_n, Dark_resp)

#finding the lowest stress for up to 6 dimensions:
dimcheckMDS(numeric_data,
            distance = "euclidean",
            k = 6) #goeveg package


#2-dimensional NMDS (with alien) ----
nmds <- metaMDS(numeric_data, distance = "euclidean", k = 2) #two dimensions, stress = 0.2146863; 3 dimensions, stress = 0.141
nmds_coords <- as.data.frame(scores(nmds, "sites"))
nmds_coords$type <- merged_trees$type

stressplot(nmds) #non-metric R^2 = 0.954, linear R^2 = 0.816

hull.data <- data.frame()
for (i in unique(nmds_coords$type)) {
  temp <- nmds_coords[nmds_coords$type == i, ][chull(nmds_coords[nmds_coords$type == i, c("NMDS1", "NMDS2")]), ]
  hull.data <- rbind(hull.data, temp)
}

(nmds_plot <- ggplot() +
    geom_polygon(data = hull.data[hull.data$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.5) + #add polygons for non-invasive types
    geom_polygon(data = hull.data[hull.data$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, group = type, fill = type), alpha = 0.8) + #add polygons for invasive type
    geom_point(data = nmds_coords, aes(x = NMDS1, y = NMDS2, color = type), size = 3) + # Add points
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900", "Alien" = "#5EA8D9")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900", "Alien" = "#5EA8D9")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Leaf Traits by Invasion Type"))
ggsave("combined_nmds_2d.jpg", nmds_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix <- vegdist(numeric_data, method = "euclidean")
anosim(diss_matrix, merged_trees$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.1405), indicating a slight but significant difference between the groups

#3-dimensional NMDS (with alien) ----
nmds3 <- metaMDS(numeric_data, distance = "euclidean", k = 3) #3 dimensions, stress = 0.1408857 
nmds_coords3 <- as.data.frame(scores(nmds3, "sites"))
nmds_coords3$type <- merged_trees$type

plot(nmds3, type = "t") #base r NMDS plot
stressplot(nmds3) #stressplot; linear R^2 = 0.887; non-linear R^2 = 0.98

gof3 <- goodness(nmds3)
plot(nmds3, typ = "t", main = "goodness of fit")
points(nmds3, display = "sites", cex = gof3 * 100)

plot(nmds3$diss, nmds3$dist)

#ggplot NMDS
hull.data3 <- data.frame()
for (i in unique(nmds_coords3$type)) {
  temp <- nmds_coords3[nmds_coords3$type == i, ][chull(nmds_coords3[nmds_coords3$type == i, c("NMDS1", "NMDS2", "NMDS3")]), ]
  hull.data3 <- rbind(hull.data3, temp)
}

(nmds3_plot <- ggplot() +
    geom_polygon(data = hull.data3[hull.data3$type != "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.5) + # polygons for non-invasive types
    geom_polygon(data = hull.data3[hull.data3$type == "Invasive", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.8) + # polygons for invasive type
    geom_polygon(data = hull.data3[hull.data3$type == "Alien", ], aes(x = NMDS1, y = NMDS2, z = NMDS3, group = type, fill = type), alpha = 0.8) + # polygons for invasive type
    geom_point(data = nmds_coords3, aes(x = NMDS1, y = NMDS2, z = NMDS3, color = type), size = 3) + # points for NMDS1, NMDS2, and NMDS3
    scale_color_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900", "Alien" = "#5EA8D9")) +
    scale_fill_manual(values = c("Native" = "#698B69", "Invasive" = "#CD6090", "Naturalised" = "#EEC900", "Alien" = "#5EA8D9")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9), legend.direction = "vertical", legend.title = element_blank()) +
    labs(title = "NMDS Plot of Leaf Traits by Invasion Type"))
ggsave("nmds_3D.jpg", nmds3_plot, path = "Plots", units = "cm", 
       width = 20, height = 20) 

diss_matrix3 <- vegdist(numeric_data, method = "euclidean", k = 3)
anosim(diss_matrix3, merged_trees$type, permutations = 9999) 
#significant (1e-04), the three types are significantly different in their traits;
#however, the R value is close to 0 (0.0949), indicating a slight but significant difference between the groups
#same as above

traits.palette2 <- c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900", 
                    "Alien" = "#5EA8D9")    #defining 3 colours

plot_ly(data = nmds_coords3, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, 
        type = "scatter3d", mode = "markers", color = ~type, colors = traits.palette2)

#3d PLOT attempt (you HAVE TO run this in order)
(ef2 = envfit(nmds3, numeric_data, permutations = 999, na.rm = TRUE, choices = 1:3))

nmds_coords3 <- as.data.frame(scores(nmds3, "sites"))
out2 <- ordiplot3d(nmds_coords3, col = as.factor(nmds_coords3$type), pch = 18,
                  ax.col = "red", envfit = ef2, arr.len = 0.1)
nmds_coords3$type <- merged_trees$type
points(out$xyz.convert(nmds_coords3), pch = 20,
       col = c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900", 
               "Alien" = "#5EA8D9"))
nmds_coords3 <- as.data.frame(scores(nmds3, "sites"))
nmdsD <- scores(nmds3, choices = 1:3, display = c("species"))
x_coords <- out$xyz.convert(nmdsD)$x
y_coords <- out$xyz.convert(nmdsD)$y
z_coords <- out$xyz.convert(nmdsD)$z
text(x = x_coords, y = y_coords, z = z_coords, labels = rownames(nmdsD), pos = 1, col = "blue")


#PCA in progress (badly and sadly) ---- 
merged_trees <- merge(trees, cn_trees[, c("code", "canopy_pos", "c_n")], by = c("code", "canopy_pos"))
combined_tree_data_numeric <- merged_trees %>%
  select(chl, lma, ldcm, A, E, g, c_n)
combined_tree_data$type_numeric <- as.numeric(combined_tree_data$type)

X <- combined_tree_data_numeric

prin_comp <- prcomp(X, rank. = 3)

components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, merged_trees$type)

tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

tit = 'Total Explained Variance = 99.48'
traits.palette <- c("Invasive" = "#CD6090", "Native" = "#698B69", "Naturalised" = "#EEC900")   
fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~merged_trees$type, colors = traits.palette) %>%
  add_markers(size = 12)


fig <- fig %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#e5ecf6")
  )

fig


##pca old way idk
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
merged_trees <- merge(trees_pos, cn_trees[, c("code", "canopy_pos", "c_n")], by = c("code", "canopy_pos"))

numeric_cols <- colnames(merged_trees)[sapply(merged_trees, is.numeric)] 
numeric_data <- merged_trees[, numeric_cols]
numeric_data <- numeric_data %>% select(chl, lma, ldcm, A, E, g, c_n, Dark_resp)

pca <- princomp(numeric_data, cor = TRUE, scores = TRUE)
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



