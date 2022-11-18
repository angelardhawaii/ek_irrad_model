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

ek_irrad_data <- read.csv("input_data/mean_supersaturation_by_period.csv")

# assign run as a factor
ek_irrad_data$Run <- as.factor(ek_irrad_data$Run)

#assign temperature as a factor
ek_irrad_data$Temperature <- as.factor(ek_irrad_data$Temp...C.)

#assigns treatment as characters from integers then to factors
ek_irrad_data$Treatment <- as.factor(as.character(ek_irrad_data$Treatment))

# assign deltaNPQ as a factor
#ek_irrad_data$deltaNPQ <- as.factor(ek_irrad_data$deltaNPQ)

#toggle between the species for output. Use Day 9 for final analysis
#recent change: removing the very odd ek.1 value of 559.4 in hypnea dataset
hypnea <- subset(ek_irrad_data, Species == "hm")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

ulva <- subset(ek_irrad_data, Species == "ul")
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"

#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_102222.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)

gr_ulva <- subset(growth_rate, Species == "Ul")

#for Hypnea remove the same replicate that was removed in the previous dataset (hm6-4 on 11/12 that had no d9 RLC) 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)


#____________________________________________________________
#ULVA 
#run model without interaction between the treatments and temperature - supersat_total is in minutes
#take RLC.Order our of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
ek_irrad_model_ulva <- lmer(formula = supersat_total ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = ulva)
rel_hsat_model_ulva <- lmer(formula = supersat_rel ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = ulva)
#make a histogram and residual plots of the data for ulva
hist(ulva$supersat_rel, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
hist(resid(rel_hsat_model_ulva))
plot(resid(rel_hsat_model_ulva) ~ fitted(rel_hsat_model_ulva))
qqnorm(resid(rel_hsat_model_ulva))
qqline(resid(rel_hsat_model_ulva))
ulva %>% ggplot(aes(supersat_total)) +
        geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#check the performance of the model
performance::check_model(rel_hsat_model_ulva)

r.squaredGLMM(rel_hsat_model_ulva)
summary(rel_hsat_model_ulva)

#check for equal variance
bartlett.test(supersat_rel ~ Treatment, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment_ulva <- oneway.test(supersat_rel ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment_ulva
welch_anova_temp_ulva <- oneway.test(supersat_rel ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp_ulva
games_howell_test(ulva, supersat_rel ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_ek_irrad_model_aov <- aov(supersat_total ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)

plot(allEffects(rel_hsat_model_ulva))

ulva %>% ggplot(aes(treatment_graph, supersat_rel)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = FALSE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Mean Rel. Hsat Time", title= "A", subtitle = "Ulva lactuca") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(0, 1) + stat_mean() + 
        geom_hline(yintercept=0.5, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#summarize the means for rETRmax
ulva %>% group_by(Treatment) %>% summarise_at(vars(supersat_rel), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(supersat_rel), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(day_length_avg), list(mean = mean))

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_irrad_ek_graph <- ggplot(ulva, aes(x=supersat_total, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "A", subtitle = "Ulva lactuca", x = "Mean Hsat Time (min)", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 475, label.y = 225) + 
        stat_cor(label.x = 475, label.y = 215) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
ulva_growth_irrad_ek_graph
ulva_growth_rel_hsat_graph <- ggplot(ulva, aes(x=supersat_rel, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "A", subtitle = "Ulva lactuca", x = "Mean Rel. Hsat", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 0.75, label.y = 200) + 
        stat_cor(label.x = 0.75, label.y = 215) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
ulva_growth_rel_hsat_graph
#HYPNEA
#_____________________________________________________________________________________________________________
ek_irrad_model_hypnea <- lmer(formula = supersat_total ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)
rel_hsat_model_hypnea  <- lmer(formula = supersat_rel ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)
#make a histogram and residual plots of the data for ulva
hist(hypnea$supersat_rel, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
hist(resid(rel_hsat_model_hypnea))
plot(resid(rel_hsat_model_hypnea) ~ fitted(rel_hsat_model_hypnea))
qqnorm(resid(rel_hsat_model_hypnea))
qqline(resid(rel_hsat_model_hypnea))
hypnea %>% ggplot(aes(supersat_total)) +
        geom_histogram(binwidth=5, fill = "#a8325e", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#check the performance of the model
performance::check_model(rel_hsat_model_hypnea)

r.squaredGLMM(rel_hsat_model_hypnea)
summary(rel_hsat_model_hypnea)

#check for equal variance
bartlett.test(supersat_rel ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(supersat_rel ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(supersat_rel ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, supersat_rel ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, supersat_rel ~ Temperature, conf.level = 0.95, detailed = TRUE)

plot(allEffects(rel_hsat_model_hypnea))

hypnea %>% ggplot(aes(treatment_graph, supersat_rel)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Mean Rel. Hsat Time", title= "B", subtitle = "Hypnea musciformis") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(0, 1) + stat_mean() + 
        geom_hline(yintercept=0.5, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw() +
        theme(legend.position = c(0.88,0.88), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_irrad_ek_graph <- ggplot(hypnea, aes(x=supersat_total, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "B", subtitle = "Hypnea musciformis", x = "Mean saturation time (min)", y = "growth rate (%)") + 
        stat_regline_equation(label.x = 450, label.y = 150) + stat_cor(label.x = 450, label.y = 140) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
hypnea_growth_irrad_ek_graph

hypnea_growth_rel_hsat_graph <- ggplot(hypnea, aes(x=supersat_rel, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "B", subtitle = "Hypnea musciformis", x = "Mean Hsat Time (min)", y = "growth rate (%)") + 
        stat_regline_equation(label.x = 0.7, label.y = 100) + stat_cor(label.x = 0.7, label.y = 108) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
hypnea_growth_rel_hsat_graph
