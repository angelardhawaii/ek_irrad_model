#script for mixed model analysis of Ek and irradiance

#load the various libraries
library(lme4)
library(lmerTest)
library(effects)
library(car)
library(MuMIn)
library (dplyr)
library(emmeans)
library(DHARMa)
library(performance)
library(patchwork)
library(rstatix)
#for plots and tables
library(ggplot2)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)

ek_irrad_data <- read.csv("/Users/Angela/src/work/limu/irradiance_ek/output/irrad_ek.csv")

# assign run as a factor
ek_irrad_data$Run <- as.factor(ek_irrad_data$Run)

#assign temperature as a factor
ek_irrad_data$Temperature <- as.factor(ek_irrad_data$Temp...C.)

#assigns treatment as characters from integers then to factors
ek_irrad_data$Treatment <- as.factor(as.character(ek_irrad_data$Treatment))

# assign deltaNPQ as a factor
ek_irrad_data$deltaNPQ <- as.factor(ek_irrad_data$deltaNPQ)

#toggle between the species for output. Use Day 9 for final analysis
#recent change: removing the very odd ek.1 value of 559.4 in hypnea dataset
hypnea <- subset(ek_irrad_data, Species == "hm" & RLC.Day == 9 & ek.1 < 226)
hypnea$treatment_graph[hypnea$Treatment == 1] <- "1) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "2) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "4) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "5) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 5] <- "3) 28ppt/53umol"

ulva <- subset(ek_irrad_data, Species == "ul" & RLC.Day == 9)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "4) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "5) 11ppt/80umol"
#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/run5-6_growth_all_042922.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate
gr_ulva <- subset(growth_rate, Species == "Ul")
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & final.weight != 0.3224)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)


#____________________________________________________________
#ULVA 
#run model without interaction between the treatments and temperature - irradiance_over_ek is in minutes
#take RLC.Order our of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
ek_irrad_model_noint_ulva <- lmer(formula = irradiance_over_ek ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = ulva)

#make a histogram and residual plots of the data for ulva
hist(ulva$irradiance_over_ek, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
hist(resid(ek_irrad_model_noint))
plot(resid(ek_irrad_model_noint) ~ fitted(ek_irrad_model_noint))
qqnorm(resid(ek_irrad_model_noint))
qqline(resid(ek_irrad_model_noint))


#check the performance of the model
performance::check_model(ek_irrad_model_noint)

r.squaredGLMM(ek_irrad_model_noint)
summary(ek_irrad_model_noint)

#check for equal variance
bartlett.test(irradiance_over_ek ~ Treatment, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment_ulva <- oneway.test(irradiance_over_ek ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment_ulva
welch_anova_temp_ulva <- oneway.test(irradiance_over_ek ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp_ulva
games_howell_test(ulva, irradiance_over_ek ~ Treatment, conf.level = 0.95, detailed = TRUE)
#run ANOVA and pairwise comparisons
anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
ulva_ek_irrad_model_aov <- aov(irradiance_over_ek ~ Treatment + Temperature, data = ulva)
TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)

plot(allEffects(ek_irrad_model_noint_ulva))

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_irrad_ek_graph <- ggplot(ulva, aes(x=irradiance_over_ek, y=growth_rate)) + geom_point() + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "Ulva lactuca Minutes above Ek vs Growth Rate", x = "Minutes above Ek", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_irrad_ek_graph

ulva %>% ggplot(aes(irradiance_over_ek)) +
        geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

ulva %>% ggplot(aes(treatment_graph, irradiance_over_ek)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Time above Ek (min)", title= "Time above Ek", subtitle = "Ulva lactuca") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        stat_mean() + 
        geom_hline(yintercept=0, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw()



#HYPNEA 
ek_irrad_model_noint_hypnea <- lmer(formula = irradiance_over_ek ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)

#make a histogram and residual plots of the data for ulva
hist(hypnea$irradiance_over_ek, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
hist(resid(ek_irrad_model_noint))
plot(resid(ek_irrad_model_noint) ~ fitted(ek_irrad_model_noint))
qqnorm(resid(ek_irrad_model_noint))
qqline(resid(ek_irrad_model_noint))


#check the performance of the model
performance::check_model(ek_irrad_model_noint)

r.squaredGLMM(ek_irrad_model_noint)
summary(ek_irrad_model_noint)

#check for equal variance
bartlett.test(irradiance_over_ek ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(irradiance_over_ek ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(irradiance_over_ek ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, irradiance_over_ek ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, irradiance_over_ek ~ Temperature, conf.level = 0.95, detailed = TRUE)

plot(allEffects(ek_irrad_model_noint_hypnea))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_irrad_ek_graph <- ggplot(hypnea, aes(x=irradiance_over_ek, y=growth_rate)) + geom_point() + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "Hypnea musciformis Minutes above Ek vs Growth Rate", x = "Minutes above Ek", y = "growth rate (%)") + 
        stat_regline_equation(label.y = 155) + stat_cor(label.y = 145)
hypnea_growth_irrad_ek_graph

hypnea %>% ggplot(aes(irradiance_over_ek)) +
        geom_histogram(binwidth=5, fill = "#a8325e", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

hypnea %>% ggplot(aes(treatment_graph, irradiance_over_ek)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Time above Ek (min)", title= "Time above Ek", subtitle = "Hypnea musciformis") + 
        scale_x_discrete(labels = c("35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        stat_mean() + 
        geom_hline(yintercept=0, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw()

