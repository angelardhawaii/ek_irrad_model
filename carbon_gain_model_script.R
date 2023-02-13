#script for mixed model analysis of Diurnal Saturated Photosynthesis Index
#This is a two species analysis and needs to be run in sequence. You can skip
#the first species but will need to run some of the code located in the Ulva section
#in order to run linear regressions with growth.

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

ek_irrad_data$RLC.Order <- as.factor(ek_irrad_data$RLC.Order)

#convert carbon_total from micromol to mol
ek_irrad_data$carbon_mol <- ek_irrad_data$carbon_total / 1e+6


#subset data by species, remove treatment 2.5 for Ulva
ulva <- subset(ek_irrad_data, Species == "ul" & Treatment != 2.5)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"



##ULVA____________________________________________________________
 
#make a histogram and residual plots of the data for ulva
hist(ulva$carbon_mol, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)

ulva %>% ggplot(aes(carbon_mol)) +
        geom_histogram(binwidth=0.25, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#run model without interaction between the treatments and temperature - supersat_avg is in minutes
#take RLC.Order our of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
dspi_model_ulva <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)

#plot residuals
hist(resid(dspi_model_ulva))
plot(resid(dspi_model_ulva) ~ fitted(dspi_model_ulva))
qqnorm(resid(dspi_model_ulva))
qqline(resid(dspi_model_ulva))


#check the performance of the model
performance::check_model(dspi_model_ulva)
r.squaredGLMM(dspi_model_ulva)
summary(dspi_model_ulva)
plot(allEffects(dspi_model_ulva))
tab_model(dspi_model_ulva, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_dspi_treatment_null <- lmer(formula = carbon_mol ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_dspi_model2 <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_dspi_treatment_null, ulva_dspi_model2)
ulva_dspi_temperature_null <- lmer(formula = carbon_mol ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_dspi_model3 <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_dspi_temperature_null, ulva_dspi_model3)


ulva %>% ggplot(aes(treatment_graph, carbon_mol)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment", y= "DSPI (mol m-2)", title= "E", subtitle = "Ulva lactuca") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", 
                                    "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(0, 10) + stat_mean() + 
        geom_hline(yintercept=5.0, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#295102", "#7CB950", "#BDE269")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))


#summarize the means for rETRmax
ulva %>% group_by(Treatment) %>% summarise_at(vars(carbon_mol), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(carbon_mol), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(day_length_avg), list(mean = mean))



#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_011723.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final.weight - growth_rate$Initial.weight) / growth_rate$Initial.weight * 100

#subset growth data
gr_ulva <- subset(growth_rate, Species == "Ul" & treatment != 2.5)
ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)

#plot growth vs dspi
ulva_growth_dspi_graph <- ggplot(ulva, aes(x=carbon_mol, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "A", subtitle = "Ulva lactuca", x = "DSPI (mol)", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 5, label.y = 225) + 
        stat_cor(label.x = 5, label.y = 215) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
ulva_growth_dspi_graph

#plot individual regressions for each treatment vs. growth
ulva_growth_carbon_indiv <- ggplot(ulva, aes(x=carbon_mol, y=growth_rate)) +
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) +
        geom_smooth(method = "lm", col = "black") + theme_bw() +
        labs(title = "A", subtitle = "Ulva lactuca", x = "Mean Daily Carbon Gain (mol m-2)", 
             y = "growth rate (%)") + 
        stat_cor() +
        #stat_regline_equation() + 
        facet_wrap(vars(Treatment), scales = "free")
ulva_growth_carbon_indiv


#HYPNEA____________________________________________________________

hypnea <- subset(ek_irrad_data, Species == "hm" & day1_rlc_time != "11:34:05")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

#make a histogram and residual plots of the data for ulva
hist(hypnea$carbon_mol, main = paste("Hypnea musciformis"), col = "#993333", labels = TRUE)

hypnea %>% ggplot(aes(carbon_mol)) +
        geom_histogram(binwidth=0.25, fill = "#993333", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#run model without interaction between the treatments and temperature - supersat_avg is in minutes
#take RLC.Order our of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
dspi_model_hypnea <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)

hist(resid(dspi_model_hypnea))
plot(resid(dspi_model_hypnea) ~ fitted(dspi_model_hypnea))
qqnorm(resid(dspi_model_hypnea))
qqline(resid(dspi_model_hypnea))

#check the performance of the model
performance::check_model(dspi_model_hypnea)
r.squaredGLMM(dspi_model_hypnea)
summary(dspi_model_hypnea)
plot(allEffects(dspi_model_hypnea))
tab_model(dspi_model_hypnea, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#construct null model to perform likelihood ratio test REML must be FALSE
hypnea_dspi_treatment_null <- lmer(formula = carbon_mol ~ Temperature + (1 | Run) + (1 | Plant.ID) , data = hypnea, REML = FALSE)
hypnea_dspi_model2 <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea, REML = FALSE)
anova(hypnea_dspi_treatment_null, hypnea_dspi_model2)
hypnea_dspi_temperature_null <- lmer(formula = carbon_mol ~ Treatment + (1 | Run) + (1 | Plant.ID), data = hypnea, REML = FALSE)
hypnea_dspi_model3 <- lmer(formula = carbon_mol ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea, REML = FALSE)
anova(hypnea_dspi_temperature_null, hypnea_dspi_model3)

#plots
hypnea %>% ggplot(aes(treatment_graph, carbon_mol)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment", y= "DSPI (mol m-2)", title= "F", subtitle = "Hypnea musciformis") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", 
                                    "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(0, 10) + stat_mean() + 
        geom_hline(yintercept=5.0, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for rETRmax
hypnea %>% group_by(Treatment) %>% summarise_at(vars(carbon_mol), list(mean = mean))
hypnea %>% group_by(Run) %>% summarise_at(vars(carbon_mol), list(mean = mean))
hypnea %>% group_by(Run) %>% summarise_at(vars(day_length_avg), list(mean = mean))

#for Hypnea remove hm6-4 on 11/12 that had no d9 RLC (final weight 0.1017)
# and hm6-4 on 10/29/21 because it was white and also looked dead 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & growth_rate_percent > -87.96837)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_carbon_graph <- ggplot(hypnea, aes(x=carbon_mol, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        xlim(2, 8) + ylim(-100, 190) + 
        labs(title = "B", subtitle = "Hypnea musciformis", x = "DSPI (mol m-2)", 
             y = "growth rate (%)") + 
        stat_regline_equation(aes(size = 14), label.x = 6, label.y = -75) + 
        scale_color_manual(values = c("#57330F", "#22773E", "#9970C2", "#FF3333", "#FF9933", "#99FF33")) +
        stat_cor(aes(size = 14), label.x = 6, label.y = -85) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05, size = 14), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05, size = 14))
hypnea_growth_carbon_graph

#plot individual regressions for each treatment vs. growth
hypnea_growth_carbon_indiv <- ggplot(hypnea, aes(x=carbon_mol, y=growth_rate)) +
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) +
        geom_smooth(method = "lm", col = "black") + theme_bw() +
        labs(title = "B", subtitle = "Hypnea musciformis", x = "Mean Daily Carbon Gain (mol m-2)", 
             y = "growth rate (%)") + 
        stat_cor() +
        #stat_regline_equation() + 
        facet_wrap(vars(Treatment), scales = "free")
hypnea_growth_carbon_indiv





#DON'T USE___________________________________________________________________________
#check for equal variance
bartlett.test(carbon_mol ~ Treatment, data = ulva)


#run Welch's ANOVA if not equal variance
welch_anova_treatment_ulva <- oneway.test(carbon_mol ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment_ulva
welch_anova_temp_ulva <- oneway.test(carbon_mol ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp_ulva
games_howell_test(ulva, carbon_mol ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_ek_irrad_model_aov <- aov(supersat_avg ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)

#check for equal variance
bartlett.test(carbon_mol ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance
welch_anova_treatment_hypnea <- oneway.test(carbon_mol ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment_hypnea
welch_anova_temp_hypnea <- oneway.test(carbon_mol ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp_hypnea
games_howell_test(hypnea, carbon_mol ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_ek_irrad_model_aov <- aov(supersat_avg ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)