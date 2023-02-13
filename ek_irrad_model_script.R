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

# make a new column for delta Ek
ek_irrad_data$deltaEk <- (((ek_irrad_data$day9_ek - ek_irrad_data$ek.1)/ ek_irrad_data$ek.1) * 100)


##ULVA ____________________________________________________________

#toggle between the species for output. 
ulva <- subset(ek_irrad_data, Species == "ul" & Treatment != 2.5)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"

#make a histogram of the data for ulva
hist(ulva$supersat_avg, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)

ulva %>% ggplot(aes(supersat_avg)) +
        geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#run model without interaction between the treatments and temperature - supersat_avg is in minutes
#take RLC.Order our of the random effects because causing problems of singularity. R2 is same with or without (+ (1 | RLC.Order))
hsat_model_ulva <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = ulva)

#check residual plots
hist(resid(hsat_model_ulva))
plot(resid(hsat_model_ulva) ~ fitted(hsat_model_ulva))
qqnorm(resid(hsat_model_ulva))
qqline(resid(hsat_model_ulva))

#check the performance of the model
performance::check_model(hsat_model_ulva)
r.squaredGLMM(hsat_model_ulva)
summary(hsat_model_ulva)
plot(allEffects(hsat_model_ulva))
tab_model(hsat_model_ulva, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_hsat_treatment_null <- lmer(formula = supersat_avg ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_hsat_model2 <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_hsat_treatment_null, ulva_hsat_model2)
ulva_hsat_temperature_null <- lmer(formula = supersat_avg ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_hsat_model3 <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_hsat_temperature_null, ulva_hsat_model3)

ulva %>% ggplot(aes(treatment_graph, supersat_avg)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Hsat Time (min)", title= "A", subtitle = "Ulva lactuca") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(140, 620) + stat_mean() + 
        geom_hline(yintercept=400, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#295102", "#7CB950", "#BDE269")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

ulva %>% ggplot(aes(treatment_graph, deltaEk)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = FALSE) + 
        labs(x="Treatment (salinty/nitrate)", y= "delta Ek (%)", title= "A", subtitle = "Ulva lactuca") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(-80, 120) + stat_mean() + 
        geom_hline(yintercept=0, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#summarize the means for variables of interest
ulva %>% group_by(Treatment) %>% summarise_at(vars(supersat_avg), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(supersat_avg), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(day_length_avg), list(mean = mean))
ulva %>% group_by(Treatment) %>% summarise_at(vars(deltaEk), list(mean = mean))

#add growth rate from other dataset to this one and subset by species for regression
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_011723.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate$growth_rate_percent <- (growth_rate$final.weight - growth_rate$Initial.weight) / growth_rate$Initial.weight * 100

#cannot remove treatment 2.5 at this stage
gr_ulva <- subset(growth_rate, Species == "Ul" & treatment != 2.5)
ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)

#plot a regression between the photosynthetic independent variables of interest and growth rate
#but BE CLEAR, the treatment 2.5 was removed for this analysis with growth rate only

ulva_growth_hsat_graph <- ggplot(ulva, aes(x=supersat_avg, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "Ulva lactuca Hsat vs Growth Rate", x = "Hsat time (minutes)", 
             y = "growth rate (%)") + stat_regline_equation(label.x = 240, label.y = 165) + stat_cor()
ulva_growth_hsat_graph





#HYPNEA
#_____________________________________________________________________________________________________________
hypnea <- subset(ek_irrad_data, Species == "hm" & day1_rlc_time != "11:34:05")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

#make a histogram  for ulva
hist(hypnea$supersat_avg, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)

hypnea %>% ggplot(aes(supersat_avg)) +
        geom_histogram(binwidth=5, fill = "#a8325e", color = "black", size = 0.25, alpha = 0.85) +
        theme_bw()

#Run model for hsat hypnea
hsat_model_hypnea <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)

#and residual plots of the data
hist(resid(hsat_model_hypnea))
plot(resid(hsat_model_hypnea) ~ fitted(hsat_model_hypnea))
qqnorm(resid(hsat_model_hypnea))
qqline(resid(hsat_model_hypnea))

#check the performance of the model
performance::check_model(hsat_model_hypnea)
r.squaredGLMM(hsat_model_hypnea)
summary(hsat_model_hypnea)
tab_model(hsat_model_hypnea, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

plot(allEffects(hsat_model_hypnea))

#construct null model to perform likelihood ratio test REML must be FALSE
hyp_hsat_treatment_null <- lmer(formula = supersat_avg ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hyp_hsat_model2 <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hyp_hsat_treatment_null, hyp_hsat_model2)
hyp_hsat_temperature_null <- lmer(formula = supersat_avg ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hyp_hsat_model3 <- lmer(formula = supersat_avg ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hyp_hsat_temperature_null, hyp_hsat_model3)

hypnea %>% ggplot(aes(treatment_graph, supersat_avg)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
        labs(x="Treatment (salinty/nitrate)", y= "Hsat Time (min)", title= "B", subtitle = "Hypnea musciformis") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(140, 620) + stat_mean() + 
        geom_hline(yintercept=400, color = "red", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
        theme_bw() +
        theme(legend.position = c(0.88,0.88), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#plot deltaEk for Hypnea
hypnea %>% ggplot(aes(treatment_graph, deltaEk)) + 
        geom_boxplot(size=0.5) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = FALSE) + 
        labs(x="Treatment (salinty/nitrate)", y= "delta Ek (%)", title= "B", subtitle = "Hypnea musciformis") + 
        scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
        ylim(-100, 250) + stat_mean() + 
        geom_hline(yintercept=0, color = "orange", size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        theme_bw() +
        theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))


#for Hypnea remove hm6-4 on 11/12 that had no d9 RLC (final weight 0.1017)
# and hm6-4 on 10/09/21 because it was white and also looked dead 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & growth_rate_percent > -87.96837)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_irrad_ek_graph <- ggplot(hypnea, aes(x=supersat_avg, y=growth_rate)) + 
        geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
        geom_smooth(method = "lm", col = "black") + theme_bw() + 
        labs(title = "B", subtitle = "Hypnea musciformis", x = "Hsat Time (min)", y = "growth rate (%)") + 
        stat_regline_equation(label.x = 400, label.y = 150) + stat_cor(label.x = 400, label.y = 140) +
        theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
              plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))
hypnea_growth_irrad_ek_graph

#summarize the means for rETRmax
hypnea %>% group_by(Treatment) %>% summarise_at(vars(supersat_avg), list(mean = mean))
hypnea %>% group_by(Run) %>% summarise_at(vars(supersat_avg), list(mean = mean))
hypnea %>% group_by(Run) %>% summarise_at(vars(day_length_avg), list(mean = mean))
hypnea %>% group_by(Treatment) %>% summarise_at(vars(deltaEk), list(mean = mean))

#no longer use_________________________________________________________________

#check for equal variance
bartlett.test(supersat_avg ~ Treatment, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment_ulva <- oneway.test(supersat_avg ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment_ulva
welch_anova_temp_ulva <- oneway.test(supersat_avg ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp_ulva
games_howell_test(ulva, supersat_avg ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(ek_irrad_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_ek_irrad_model_aov <- aov(supersat_avg ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_ek_irrad_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_ek_irrad_model_aov, "Temperature", ordered = FALSE)

#check for equal variance
bartlett.test(supersat_avg ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(supersat_avg ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(supersat_avg ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, supersat_avg ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, supersat_avg ~ Temperature, conf.level = 0.95, detailed = TRUE)