#this script contains all models, checks, and plots for the components of this
#analysis that include a GLM

library(lme4)
library(tidyverse)
library(lamW)
library(frair)
library(ggeffects)
library(DHARMa)
library(emmeans)
library(patchwork)
source("ggplot_paper_theme.R") #for custom plot themes
#custom colour palettes for plots
pal <- c("#625a94", "#11c2b5")
#pal2 <- c("#3f3a5e", "#73bdb7")
#pal3 <- c("#8a84b3", "#1e7d76")
#pal4 <- c("#fd8d3c","#800026")

# LOADING DATA ----------------------------------------
#the follow csv file is the complete dataset for both Beth's experiment and this 
#one
crab <- read.csv("greencrabforaging.csv") 

#this study only used the varnish clams (VC) and crabs with no sediment so 
#we'll filter out the rest of the data collected at the same time
mycrab <- filter(crab, Clam == "VC", Sed == "no") %>% 
  mutate(clam_alive = Density - Clam_Consumed,
         #standardize response variables
         stand_density = c(scale(Density)),
         stand_cw = c(scale(CW)),
         stand_claw = c(scale(Claw_width)),
         prop_eaten = Clam_Consumed/Density,
         claw_ratio = Claw_width/CW,
         stand_claw_ratio = c(scale(claw_ratio)))

#data on clam sizes from trials and whether each was consumed
female_clam_consumption <- read_csv("clam_size_consumption_female.csv") %>% 
  mutate(Sex = "Female crabs")
clam_consumption <- read_csv("clam_size_consumption_male.csv") %>% 
  #only want varnish clams and no sediment so we'll filter out the rest of the 
  #data collected at the same time
  filter(Species == "VC" & `Sed (1=yes)` == 0) %>% 
  mutate(Sex = "Male crabs") %>% 
  #combine with the female data
  bind_rows(female_clam_consumption) %>% 
  #fix column names
  rename(consumed = `Consumed (1= consumed)`,
         length = `Length (siphon to foot)`,
         width = `Width (hinge up)`,
         sex = Sex) %>% 
  #change consumption to a factor rather than integer
  mutate(consumed = factor(consumed,levels = c("0", "1")),
         consumed_yn = case_when(consumed == "0" ~ "Unconsumed",
                                 TRUE ~ "Consumed")) %>% 
  #remove any clams where we don't have data on if they were consumed or not
  filter(!is.na(consumed))

#crab behaviour data
behaviour <- read.csv("crab_behaviour.csv",  header=T, stringsAsFactors=T) %>% 
  mutate(total_time = time_moving + time_still,
         prop_moving_exact = time_moving/total_time)

crab_sizes <- read_csv("crab_population.csv") 

# FUNCTIONAL RESPONSE MODELS----------------------------------------------------
#males first
mcrab <- filter(mycrab, Sex == "M")

# Test for type II/type III
frair_test(Clam_Consumed ~ Density, data = mcrab)
#evidence for type II

# Frair fit
outII_male <- frair_fit(Clam_Consumed ~ Density, data = mcrab, 
                        response = 'rogersII',
                        start = list(a = 0.2, h = 0.2), fixed = list(T=1))
# A linear fit
outI_male <- frair_fit(Clam_Consumed ~ Density, data = mcrab, 
                       response = 'typeI',
                       start = list(a = 0.2), fixed=list(T=1))

# Visualize fits
plot(outII_male, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(outII_male)
lines(outI_male, lty=3)

# male parameter estimates
male_attack <- outII_male$coefficients[1]
male_handling <- outII_male$coefficients[2]

#extract just the clam densities from the original df
fitsmale <- data.frame(x = mcrab$Density) 

#calculate expected number of clam eatens eaten for each observations, based on 
#frair flexnr equation, using lambert function
fitsmale$Ne <- fitsmale$x - 
  lambertW0(male_attack * male_handling * fitsmale$x * 
              exp(-male_attack * (1 - male_handling * fitsmale$x)))/
                                                  (male_attack * male_handling) 

#now calculate the residuals by examining the predicted and observed values
fitsmale$actual <- mcrab$Clam_Consumed
fitsmale$resid <- fitsmale$Ne - fitsmale$actual

plot(x = fitsmale$x, y = fitsmale$Ne)
plot(x = fitsmale$Ne, y = fitsmale$resid)
abline(h = 0, lty = 'dotted')

# look at original fits returned by mle2
summary(outII_male$fit)
summary(outI_male$fit)

# Compare models using AIC
AIC(outI_male$fit,outII_male$fit) 
#type II is substantially better

# Bootstrap confidence intervals around the predictions
set.seed(309331)
outII_male_boot <- frair_boot(outII_male, start = NULL, strata=mcrab[,6], 
                              nboot=2000,
                              para=TRUE, ncores=NaN, WARN.ONLY=FALSE)
outII_male_boot
confint(outII_male_boot)

# Illustrate bootlines
plot(outII_male_boot, xlim=c(0,16), ylim = c(0, 16), type='n', 
     main='All bootstrapped lines')
lines(outII_male_boot, all_lines=TRUE)
points(outII_male_boot, pch=20)

# Illustrate bootpolys
plot(outII_male_boot, xlim=c(0,16), ylim = c(0, 16), type='n', 
     main='Empirical 95 percent CI')
drawpoly(outII_male_boot, col=hsv(2/6,0.2, 0.8))
points(outII_male_boot, pch=20)
lines(outII_male_boot, all_lines=FALSE)

# nlm to look at asymptote 
mcrab_asymp <- nls(Clam_Consumed ~ SSasymp(Density, Asym, R0, lrc), 
                   data = mcrab)

summary(mcrab_asymp)


#now do the same for females
fcrab <- filter(mycrab, Sex == "F")

# Test for type II/type III
frair_test(Clam_Consumed ~ Density, data = fcrab)

# Frair fit
outII_female <- frair_fit(Clam_Consumed ~ Density, data = fcrab, 
                          response = 'rogersII',
                          start = list(a = 0.2, h = 0.2), fixed = list(T=1))
# A linear fit
outI_female <- frair_fit(Clam_Consumed ~ Density, data = fcrab, 
                         response = 'typeI',
                         start = list(a = 0.2), fixed=list(T=1))

# Visualise fits
plot(outII_male, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(outII_female)
lines(outI_female, lty=3)

# female resid
female_attack <- outII_female$coefficients[1]
female_handling <- outII_female$coefficients[2]

fitsfemale <- data.frame(x = fcrab$Density) 
# Calculate 'a' for each value of x, where x is density of clam

fitsfemale$Ne <- fitsfemale$x - 
  lambertW0(female_attack * female_handling * fitsfemale$x * 
              exp(-female_attack * (1 - female_handling * fitsfemale$x)))/
  (female_attack * female_handling) 
# calculate expected number of clam eaten, based on frair flexnr equation, 
#using lambert function

fitsfemale$actual <- fcrab$Clam_Consumed
fitsfemale$resid <- fitsfemale$Ne - fitsfemale$actual

plot(x = fitsfemale$x, y = fitsfemale$Ne)
plot(x = fitsfemale$Ne, y = fitsfemale$resid)
abline(h = 0, lty = 'dotted')
#pretty nice!

# Have a look at original fits returned by mle2
summary(outII_female$fit)
summary(outI_female$fit)

# Compare models using AIC
AIC(outI_female$fit,outII_female$fit) 
#type II is for sure better

set.seed(309331)
outII_female_boot <- frair_boot(outII_female, start = NULL, strata=fcrab[,6], 
                                nboot=2000,
                                para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_female_boot
confint(outII_female_boot)

# Illustrate bootlines
plot(outII_female_boot, xlim=c(0,16), ylim = c(0, 16), type='n', 
     main='All bootstrapped lines')
lines(outII_female_boot, all_lines=TRUE)
points(outII_female_boot, pch=20)

# Illustrate bootpolys
plot(outII_female_boot, xlim=c(0,16), ylim = c(0, 16), type='n', 
     main='Empirical 95 percent CI')
drawpoly(outII_female_boot, col=hsv(2/6,0.2, 0.8))
points(outII_female_boot, pch=20)
lines(outII_female_boot, all_lines=FALSE)

# nls to look at asymptote 
fcrab_asymp <- nls(Clam_Consumed ~ SSasymp(Density, Asym, R0, lrc), 
                   data = fcrab)
summary(fcrab_asymp)

#add jitter to the points manually, since we can't do it within the frboot.plot
outII_female_boot$x <- jitter(outII_female_boot$x, amount = 0.07)
outII_male_boot$x <- jitter(outII_male_boot$x, amount = 0.07)

#final figure
pdf("Figure1.pdf", width = 9, height = 6)
par(bg = 'white', fg = 'black', mar = c(4.1, 4.1, 3.1, 2.1))
plot(outII_female_boot, xlim=c(0, 16), ylim = c(0, 10), type='n',
     xlab = "",
     ylab = "",
     cex.lab = 1.5,
     cex.axis = 1,
     cex.main = 1.5,
     bty="l")
lines(outII_female_boot, lwd = 3, all_lines=FALSE, col= "#625a94", lty = 2)
lines(outII_male_boot, lwd = 3, all_lines=FALSE, col= "#11c2b5", lty = 1)
drawpoly(outII_female_boot, border = NA, 
         col=adjustcolor("#625a94", alpha.f = 0.4))
drawpoly(outII_male_boot, border = NA, 
         col=adjustcolor("#11c2b5", alpha.f = 0.4))
points(outII_female_boot, pch=24, col="#625a94", 
       bg=adjustcolor("#625a94", alpha.f = 0.4), cex = 1.5)
points(outII_male_boot, pch=21,  col="#196c66", 
       bg=adjustcolor("#11c2b5", alpha.f = 0.4), cex = 1.5)
legend(x = "topleft", legend = c("Male", "Female"), 
       col = c("#11c2b5", "#625a94"), lty = c(1, 2), 
       cex = 1, pt.cex = 1.5, lwd = 1.5, pch = c(20, 17),
       bty="n")
title(xlab = "Initial clam density (number/trial)", mgp = c(2, 1, 0), cex = 1.5)
title(ylab = "Number of clams consumed", mgp = c(2, 1, 0), cex = 1.5)
dev.off()

# BOOTSTRAP NeS ----------------------------------------------------------------
#compare Ne at max experimental density
# Set clam density to the max prey density offered
clam_den <- c(16)

# get the bootstrapped coefficients
a_m <- outII_male_boot$bootcoefs[, 1]
h_m <- outII_male_boot$bootcoefs[, 2]
T_m <- outII_male_boot$bootcoefs[, 3]

# use the coefficients to calculate expected number of clams eaten, based on 
#frair flexnr equation, using the rogersII function
male_boot <- rogersII(16, a_m, h_m, T_m)

# do the same for the females
a_f <- outII_female_boot$bootcoefs[, 1]
h_f <- outII_female_boot$bootcoefs[, 2]
T_f <- outII_female_boot$bootcoefs[, 3]

#can look at density = 16 to see realistic values
female_boot <- rogersII(16, a_f, h_f, T_f)
boot_dist_test <- data.frame(female = female_boot, 
                             male = male_boot,
                             diff = male_boot - female_boot)
#the 2.5% and 97.5% quantiles of a bootstrapped metric (in this case the 
#difference between males and females) represent the upper and lower bounds of
#the 95% confidence interval
lower_ci <- quantile(boot_dist_test$diff, 0.025)
upper_ci <- quantile(boot_dist_test$diff, 0.975)

ggplot(boot_dist_test, aes(diff)) + 
  geom_histogram(binwidth = 0.25) +
  theme_paper_large() +
  geom_vline(xintercept = lower_ci, col = "red") + 
  geom_vline(xintercept = upper_ci, col = "red") + 
  theme(legend.position = "none") +
  labs(x = "Difference between male Ne and female Ne",
       y = "Count")
#ggsave("male_vs_female_Ne.pdf",width = 7, height = 5)

# PERMUTATION TEST FOR NeS AND ASYMPTOTES----------------------------------------
#the CI above can still be used for signifance testing (i.e., if the 95% CI 
#crosses 0, it's not significant), but we can use a permuation test if we want
#a precise p-value
mycrab_perm <- mycrab %>% 
  select(Sex, Density, Clam_Consumed)

#get the mean Nes for males and females from the original models
a_m_mean <- outII_male$coefficients[1]
h_m_mean <- outII_male$coefficients[2]
T_m_mean <- outII_male$coefficients[3]
male_ne_mean <- rogersII(16, a_m_mean, h_m_mean, T_m_mean)
a_f_mean <- outII_female$coefficients[1]
h_f_mean <- outII_female$coefficients[2]
T_f_mean <- outII_female$coefficients[3]
female_ne_mean <- rogersII(16, a_f_mean, h_f_mean, T_f_mean)

#the models predict that males will have a higher Ne than females
mean_ne_diff = male_ne_mean-female_ne_mean

#now permutate the data 2000 times to see the distribution of differences
ne_diffs <- NULL
for(i in 1:2000){
    sampled <- sample(seq(from = 1, to = nrow(mycrab_perm), by = 1), 
                      replace = FALSE)
    temp_df <- tibble(Sex = mycrab_perm$Sex,
                      Density = mycrab_perm$Density[sampled],
                      Clam_Consumed = mycrab_perm$Clam_Consumed[sampled])
    temp_male <- filter(temp_df, Sex == "M")
    outII_male_perm <- frair_fit(Clam_Consumed ~ Density, data = temp_male, 
                            response = 'rogersII',
                            start = list(a = 0.2, h = 0.2), fixed = list(T=1))
    temp_female <- filter(temp_df, Sex == "F")
    outII_female_perm <- frair_fit(Clam_Consumed ~ Density, data = temp_female, 
                                 response = 'rogersII',
                                 start = list(a = 0.2, h = 0.2), 
                                 fixed = list(T=1))
    # get the coefficients
    a_m_perm <- outII_male_perm$coefficients[1]
    h_m_perm <- outII_male_perm$coefficients[2]
    T_m_perm <- outII_male_perm$coefficients[3]
    male_ne <- rogersII(16, a_m_perm, h_m_perm, T_m_perm)
    a_f_perm <- outII_female_perm$coefficients[1]
    h_f_perm <- outII_female_perm$coefficients[2]
    T_f_perm <- outII_female_perm$coefficients[3]
    female_ne <- rogersII(16, a_f_perm, h_f_perm, T_f_perm)
    ne_diffs[i] <- male_ne-female_ne
}
hist(ne_diffs)

#now see what proportion of the time the permutated difference is greater than 
#or equal to the mean difference
props <- ne_diffs %>% 
  as.data.frame() %>% 
  filter(`.` >= mean_ne_diff) %>% 
  summarize(n = n())

#multiply by 2 for a two-sided p-value
props[1,1]/2000 * 2  

# ATTACK RATE, HANDLING TIME, MAX CONSUMPTION COMPARISON -----------
#use bootcoefs to calculate CI around FRR
female_boots <- as.data.frame(outII_female_boot$bootcoefs) %>% 
  transmute(frr_female = a/h,
            a_female = a,
            h_female = h,
            max_rate_female = 1/h)
male_boots <- as.data.frame(outII_male_boot$bootcoefs) %>% 
  transmute(frr_male = a/h,
         a_male = a,
         h_male = h,
         max_rate_male = 1/h)
compare_boots <- bind_cols(female_boots, male_boots) %>% 
  mutate(frr_diff = frr_male - frr_female,
         frr_ratio = frr_male/frr_female,
         a_diff = a_male - a_female,
         h_diff = h_male - h_female,
         max_rate_diff = max_rate_male - max_rate_female)
compare_summary <- compare_boots %>% 
  summarize(a_diff_lower = quantile(a_diff, 0.025),
            a_diff_mean = mean(a_diff),
            a_diff_upper = quantile(a_diff, 0.975),
            h_diff_lower = quantile(h_diff, 0.025),
            h_diff_mean = mean(h_diff),
            h_diff_upper = quantile(h_diff, 0.975),
            max_rate_diff_lower = quantile(max_rate_diff, 0.025),
            max_rate_diff_mean = mean(max_rate_diff),
            max_rate_diff_upper = quantile(max_rate_diff, 0.975),
            frr_diff_lower = quantile(frr_diff, 0.025),
            frr_diff_mean = mean(frr_diff),
            frr_diff_upper = quantile(frr_diff, 0.975),
            frr_ratio_lower = quantile(frr_ratio, 0.025),
            frr_ratio_mean = mean(frr_ratio),
            frr_ratio_upper = quantile(frr_ratio, 0.975),
            frr_male_lower = quantile(frr_male, 0.025),
            frr_male_mean = mean(frr_male),
            frr_male_upper = quantile(frr_male, 0.975),
            frr_female_lower = quantile(frr_female, 0.025),
            frr_female_mean = mean(frr_female),
            frr_female_upper = quantile(frr_female, 0.975))

ggplot(compare_boots, aes(frr_diff)) + 
  geom_histogram()
ggplot(compare_boots, aes(a_diff)) + 
  geom_histogram()
ggplot(compare_boots, aes(h_diff)) + 
  geom_histogram()

compare_lower <- quantile(compare_boots$frr_diff, 0.025)
compare_upper <- quantile(compare_boots$frr_diff, 0.975)

female_frr <- mean(female_boots$frr)
female_lower_frr <- quantile(female_boots$frr_female, 0.025)
female_upper_frr <- quantile(female_boots$frr_female, 0.975)
male_frr <- mean(male_boots$frr)
male_lower_frr <- quantile(male_boots$frr_male, 0.025)
male_upper_frr <- quantile(male_boots$frr_male, 0.975)

# use frair to compare a and h between males and females
frair_compare(outII_female, outII_male, start = NULL)

female_params <- as.data.frame(coef(outII_female_boot)) %>% 
  rownames_to_column() %>% 
  mutate(lower_bca = NA,
         upper_bca = NA) %>% 
  filter(rowname != "T") %>% 
  mutate(sex = "female") %>% 
  rename(value = `coef(outII_female_boot)`)

female_params[1,3] <- confint(outII_female_boot)[[1]]$bca$lower
female_params[1,4] <- confint(outII_female_boot)[[1]]$bca$upper
female_params[2,3] <- confint(outII_female_boot)[[2]]$bca$lower
female_params[2,4] <- confint(outII_female_boot)[[2]]$bca$upper


male_params <- as.data.frame(coef(outII_male_boot)) %>% 
  rownames_to_column() %>% 
  mutate(lower_bca = NA,
         upper_bca = NA) %>% 
  filter(rowname != "T") %>% 
  mutate(sex = "male") %>% 
  rename(value = `coef(outII_male_boot)`)

male_params[1,3] <- confint(outII_male_boot)[[1]]$bca$lower
male_params[1,4] <- confint(outII_male_boot)[[1]]$bca$upper
male_params[2,3] <- confint(outII_male_boot)[[2]]$bca$lower
male_params[2,4] <- confint(outII_male_boot)[[2]]$bca$upper

both_params <- bind_rows(female_params, male_params)

a_param <- both_params %>% 
  filter(rowname == "a")

a <- ggplot(a_param, aes(sex, value, colour = sex, shape = sex)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower_bca, ymax = upper_bca), width = 0.2) +
  theme_paper_large() +
  labs(y = "Attack rate",
       x = "") +
  theme(legend.position = "none",
        axis.text.x = element_blank())+
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(17,19)) +
  scale_x_discrete(labels = c("Female", "Male"))

h_param <- both_params %>% 
  filter(rowname == "h")

h <- ggplot(h_param, aes(sex, value, colour = sex, shape = sex)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower_bca, ymax = upper_bca), width = 0.2) +
  theme_paper_large() +
  labs(y = "Handling time",
       x = "") +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(17,19)) +
  scale_x_discrete(labels = c("Female", "Male"))

ne_param <- boot_dist_test %>% 
  summarize(sex_female_mean = mean(female),
            sex_male_mean = mean(male),
            sex_female_upper = quantile(female, 0.975),
            sex_female_lower = quantile(female, 0.025),
            sex_male_upper = quantile(male, 0.975),
            sex_male_lower = quantile(male, 0.025)) 
#janky fix to this df
ne_param_plotting <- tibble(sex = c("female","male"),
                            mean = c(ne_param[1,1], ne_param[1,2]),
                            upper = c(ne_param[1,3], ne_param[1,5]),
                            lower = c(ne_param[1,4], ne_param[1,6]))

ne <- ggplot(ne_param_plotting, aes(sex, mean, colour = sex, shape = sex)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_paper_large() +
  labs(y = "Ne at density = 16",
       x = "Sex") +
  theme(legend.position = "none") +
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(17,19)) +
  scale_x_discrete(labels = c("Female", "Male"))


#can also bootstrap 1/ht for each sex
max_rate_male <- tibble(value = 1/h_m*T_m) %>% 
  summarize(mean = mean(value),
            upper_ci = quantile(value,0.025),
            lower_ci = quantile(value,0.975)) %>% 
  mutate(sex = "Male")
max_rates <- tibble(value = 1/h_f*T_f) %>% 
  summarize(mean = mean(value),
            upper_ci = quantile(value,0.025),
            lower_ci = quantile(value,0.975)) %>% 
  mutate(sex = "Female") %>% 
  bind_rows(max_rate_male)

max_rates_plot <- ggplot(max_rates, aes(sex, mean, colour = sex, shape = sex)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  theme_paper_large() +
  labs(y = "Maximum consumption rate",
       x = "Sex") +
  theme(legend.position = "none") +
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(17,19)) +
  scale_x_discrete(labels = c("Female", "Male"))
max_rates_plot

a + h + max_rates_plot +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a") +
  theme(plot.margin = margin(r = 0, l = 0)) & 
  theme(plot.tag.position = c(0.17, 0.97),
        plot.tag = element_text(size = 16, hjust = 0, vjust = 0))
#ggsave("Figure2.pdf", width = 5, height = 9)
#ggsave("Figure2.png", width = 5, height = 9)
# COMPARISON OF CRAB SIZE IN FUNCTIONAL RESPONSE EXPERIMENT ------------
cara_compare <- lm(CW ~ Sex, data = mycrab)
summary(cara_compare) #carapace width is almost significantly different
simulateResiduals(cara_compare, plot = TRUE)
#residuals look good

claw_compare <- lm(Claw_width ~ Sex, data = mycrab)
summary(claw_compare) # claw size is significantly different 
simulateResiduals(claw_compare, plot = TRUE)

# COMPARISON OF CRAB SIZE IN WILD ---------------------------------------------
cara_compare_wild <- lm(cw ~ sex, data = crab_sizes)
summary(cara_compare_wild) #carapace width is almost significantly different
simulateResiduals(cara_compare_wild, plot = TRUE)
#residuals look good

claw_compare_wild <- lm(claw ~ sex, data = crab_sizes)
summary(claw_compare_wild) # claw size is significantly different 
simulateResiduals(claw_compare_wild, plot = TRUE)
#not perfect but reasonable for the simple comparison we're making
# GLM MODEL USING SIZE AS PREDICTORS OF PROPORTION EATEN --------------------
# see what size measure is a better predictor of proportion eaten
prop_eaten_model_cara <- glm(prop_eaten ~ stand_density + stand_cw + Sex + 
                          Sex:stand_cw, 
                        data = mycrab, 
                        weights = Density, 
                        family = binomial(link = logit))
summary(prop_eaten_model_cara)

prop_eaten_model_claw <- glm(prop_eaten ~ stand_density + stand_claw + Sex + 
                          Sex:stand_claw, 
                        data = mycrab, 
                        weights = Density, 
                        family = binomial(link = logit))
summary(prop_eaten_model_claw)

#compare with AIC corrected for small sample size
MuMIn::AICc(prop_eaten_model_cara, prop_eaten_model_claw)
#ok they are virtually identical in how well they explain the data and they have
#slightly different explanations, so we'll present both. We'll cover claw size 
#in the main text since we have a stronger biological rationale for using it
#as a predictor, and cover carapace width in the supplement

#and a quick multicollinearity check without the interaction:
coll_test1 <-glm(prop_eaten ~ stand_density + 
                                             stand_cw + Sex, 
                                           data = mycrab, 
                                           weights = Density, 
                                           family = binomial(link = logit))
car::vif(coll_test1)
coll_test2 <- glm(prop_eaten ~ stand_density + 
                                             stand_claw + Sex, 
                                          data = mycrab, 
                                          weights = Density, 
                                          family = binomial(link = logit))
car::vif(coll_test2)
#and both are fine
# PLOT OF PROP_EATEN_MODEL ----------------------
#CARAPACE WIDTH
df_predict <- ggpredict(prop_eaten_model_cara, 
                        terms = c("stand_cw[n=1000]", "Sex"),
                        condition = c(stand_density = 0)) %>% 
  rename(stand_cw = x, 
         Sex = group) %>%
  mutate(CW = stand_cw*(sd(mycrab$CW)) + mean(mycrab$CW)) %>% 
  #limit predictions to values in the data
  filter((Sex == 'M' & (stand_cw >= min(filter(mycrab, Sex == "M")$stand_cw) & 
                          stand_cw <= max(filter(mycrab, Sex == "M")$stand_cw))) |
           (Sex == 'F' & (stand_cw >= min(filter(mycrab, Sex == "F")$stand_cw) & 
                            stand_cw <= max(filter(mycrab, Sex == "F")$stand_cw))))

ggplot() + 
  geom_ribbon(data = df_predict, aes(x = CW, ymin = conf.low, ymax = conf.high, 
                                     fill = Sex),
              alpha = 0.3) +
  geom_line(data = df_predict, aes(x = CW, y = predicted, 
                                   colour = Sex, 
                                   lty = Sex)) +
  geom_point(data = mycrab, aes(x = CW, y = prop_eaten, 
                                color = Sex, pch = Sex)) +
  labs(x = "Carapace width (mm)", y = "Proportion consumed") +
  scale_color_manual(name = "Sex", values = pal, labels = c("Female", "Male")) +
  scale_fill_manual(values = pal, labels = c("Female", "Male"), guide = NULL) +
  scale_linetype_manual(name = "Sex", values = c("dashed","solid"),
                        labels = c("Female", "Male")) +
  scale_shape_manual(values = c(17, 19),labels = c("Female", "Male")) +
  ylim(0,1) +
  theme_paper_large() +
  theme(legend.position = c(0.08,0.95),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size=unit(2,"lines")) 
#ggsave("FigureS2.png",width = 9, height = 6)
#ggsave("FigureS2.pdf",width = 9, height = 6)

#CLAW WIDTH
claw_prop_eaten_model_predict <- ggpredict(prop_eaten_model_claw, 
                                           terms = c("stand_claw", "Sex"),
                                           condition = c(stand_density = 0)) %>% 
  rename(stand_claw = x, Sex = group) %>%
  mutate(Claw_width = stand_claw*(sd(mycrab$Claw_width)) +
           mean(mycrab$Claw_width)) %>% 
  filter((Sex == 'M' & (stand_claw >= min(filter(mycrab, Sex == "M")$stand_claw) & 
                          stand_claw <= max(filter(mycrab, Sex == "M")$stand_claw))) |
           (Sex == 'F' & (stand_claw >= min(filter(mycrab, Sex == "F")$stand_claw) & 
                            stand_claw <= max(filter(mycrab, Sex == "F")$stand_claw))))

ggplot() + 
  geom_ribbon(data = claw_prop_eaten_model_predict, 
              aes(x = Claw_width, ymin = conf.low, ymax =  conf.high, 
                  fill = Sex), alpha = 0.3) +
  geom_line(data = claw_prop_eaten_model_predict, 
            aes(x = Claw_width, y = predicted, colour = Sex, lty = Sex)) +
  geom_point(data = mycrab, aes(x = Claw_width, y = prop_eaten, 
                                color = Sex, shape = Sex)) +
  labs(x = "Crusher claw height (mm)", y = "Proportion consumed") +
  scale_color_manual(name = "Sex", values = pal, labels = c("Female", "Male")) +
  scale_fill_manual(values = pal, labels = c("Female", "Male"), guide = NULL) +
  scale_linetype_manual(name = "Sex", values = c("dashed","solid"),
                        labels = c("Female", "Male")) +
  scale_shape_manual(values = c(17, 19),labels = c("Female", "Male")) +
  ylim(0,1) +
  theme_paper_large() +
  theme(legend.position = c(0.08,0.95),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size=unit(2,"lines")) 
#ggsave("Figure3.pdf",width = 9, height = 6)
#ggsave("Figure3.png",width = 9, height = 6)

# CONSUMPTION OF CLAM SIZE -------------------------------
# clam length
lm_clam_consumption <- lm(length ~ sex * consumed, data = clam_consumption)
summary(lm_clam_consumption)
simulateResiduals(lm_clam_consumption, plot = TRUE)
emmeans(lm_clam_consumption, pairwise ~ sex | consumed)$contrasts
emmeans(lm_clam_consumption, pairwise ~ consumed | sex)$contrasts

lm_clam_consumption_fit <- ggpredict(lm_clam_consumption, 
                                     terms = c("sex", "consumed")) %>% 
  mutate(sex = x, 
         length = predicted,
         consumed = case_when(group == 0 ~ "Unconsumed",
                              TRUE ~ "Consumed"))
ggplot() + 
  geom_jitter(data = clam_consumption,
              aes(x = consumed_yn, y = length, colour = sex, shape = sex), 
              alpha = 0.2, height=0) +
  geom_point(data = lm_clam_consumption_fit, 
             aes(x = consumed, y = length, colour = sex, shape = sex), 
             size = 3) +
  geom_errorbar(data = lm_clam_consumption_fit, 
                aes(x = consumed, ymin = conf.low, ymax = conf.high, 
                    colour = sex), 
                width = 0.5) +
  labs(y="Clam length (mm)", x = "") +
  facet_wrap(~sex)+
  theme_paper_large() +
  scale_color_manual(values = pal, labels = c("Female", "Male")) +
  scale_shape_manual(values = c(17,19)) +
  theme(legend.position = "none")
#ggsave("Figure4.pdf",width = 9, height = 6)
#ggsave("Figure4.png",width = 9, height = 6)

# clam width
#lm_clam_consumption_w <- lm(width ~ Sex * consumed, data = clam_consumption)
#summary(lm_clam_consumption_w)
#simulateResiduals(lm_clam_consumption_w, plot = TRUE)
#
#lm_clam_consumption_fit_w <- ggpredict(lm_clam_consumption_w, 
#                                     terms = c("Sex", "consumed")) %>% 
#  mutate(sex = x, 
#         width = predicted,
#         consumed = group)
#
#ggplot() + 
#  geom_jitter(data = clam_consumption,
#              aes(x = consumed, y = width), 
#              alpha = 0.2, height=0) +
#  geom_point(data = lm_clam_consumption_fit_w, 
#             aes(x = consumed, y = width), 
#             size = 3) +
#  geom_errorbar(data = lm_clam_consumption_fit_w, 
#                aes(x = consumed, ymin = conf.low, ymax = conf.high), 
#                width = 0.5) +
#  labs (y="Width of clams consumed (mm)", x = "Consumed") +
#  facet_wrap(~sex)+
#  theme_classic()

# EXPLORATORY BEHAVIOUR --------------------------------
#TESTING DIFFERENCES IN PATH LENGTH 
lm_path <- lm(path_length ~ sex*cw, data = behaviour)
summary(lm_path)
simulateResiduals(lm_path, plot = T)

lm_path_fit <- ggpredict(lm_path, terms = "sex") %>% 
  mutate(sex = x, path_length = predicted)

# TESTING DIFFERENCES IN PROPORTION OF TIME MOVING
lm_prop_moving <- lm(prop_moving_exact ~ sex*cw, data = behaviour)
summary(lm_prop_moving)
simulateResiduals(lm_prop_moving, plot = T)
#looks good - note that the lm was a way better fit than a logistic

lm_prop_moving_fit <- ggpredict(lm_prop_moving, terms = "sex") %>% 
  mutate(sex = x, prop_moving = predicted)

# TESTING DIFFERENCES IN SPEED
lm_speed <- lm(speed ~ sex*cw, data = behaviour)
summary(lm_speed)
simulateResiduals(lm_speed, plot = T)

lm_speed_fit <- ggpredict(lm_speed, terms = "sex") %>% 
  mutate(sex = x, speed = predicted)

path_plot <- ggplot() + 
  geom_jitter(data = behaviour, 
              aes(x = sex, y = path_length, colour = sex, shape = sex), 
              alpha = 0.4, height=0) +
  geom_point(data = lm_path_fit, 
             aes(x = sex, y = path_length, colour = sex, shape = sex), 
             size = 5) +
  geom_errorbar(data = lm_path_fit, 
                aes(x = sex, ymin = conf.low, ymax = conf.high, colour = sex), 
                width = 0.2) +
  labs (y="Path length (cm)", x = "") +
  ylim(200,1200) +
  theme_paper_large() +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(17,19)) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) + 
  scale_x_discrete(labels = c("Female", "Male"))
path_plot

moving_plot <- ggplot() + 
  geom_jitter(data = behaviour, 
              aes(x = sex, y = prop_moving, colour = sex, shape = sex), 
              alpha = 0.4, height=0) +
  geom_point(data = lm_prop_moving_fit, 
             aes(x = sex, y = prop_moving, colour = sex, shape = sex), 
             size = 5) +
  geom_errorbar(data = lm_prop_moving_fit, 
                aes(x = sex, ymin = conf.low, ymax = conf.high, colour = sex), 
                width = 0.2) +
  labs (y="Proportion time moving", x = "") +
  ylim(0,1) +
  theme_paper_large()+
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(17,19)) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) + 
  scale_x_discrete(labels = c("Female", "Male"))
moving_plot 


speed_plot <- ggplot() + 
  geom_jitter(data = behaviour, 
              aes(x = sex, y = speed, colour = sex, shape = sex), 
              alpha = 0.4, height=0) +
  geom_point(data = lm_speed_fit, 
             aes(x = sex, y = speed, colour = sex, shape = sex), 
             size = 5) +
  geom_errorbar(data = lm_speed_fit, 
                aes(x = sex, ymin = conf.low, ymax = conf.high, colour = sex), 
                width = 0.2) +
  labs (y="Speed (cm/s)", x = "Sex") +
  theme_paper_large() +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(17,19)) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("Female", "Male"))

speed_plot 

# Behaviour plot

path_plot + moving_plot + speed_plot + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a") +
  theme(plot.margin = margin(r = 0, l = 0)) & 
  theme(plot.tag.position = c(0.17, 0.97),
        plot.tag = element_text(size = 16, hjust = 0, vjust = 0))
#ggsave("Figure5.pdf", width = 5, height = 9)
#ggsave("Figure5.png", width = 5, height = 9)

# SUPPLEMENTAL FIGURE  outliers------------------------------------------------
# Plot with the the two "outliers"/lowest female values removes
mycrab_noout <- mycrab %>% 
  filter(CW > 55) %>% 
  mutate(stand_density = c(scale(Density)),
         stand_cw = c(scale(CW)))

prop_eaten_model_noout <- glm(prop_eaten ~ stand_density + stand_cw + Sex + 
                                Sex:stand_cw, 
                              data = mycrab_noout, 
                              weights = Density,
                              family = binomial(link = logit))

summary(prop_eaten_model_noout)

df_predict_noout <- ggpredict(prop_eaten_model_noout, 
                              terms = c("stand_cw", "Sex")) %>% 
  rename(stand_cw = x, 
         Sex = group) %>%
  mutate(CW = stand_cw*(sd(mycrab_noout$CW)) + mean(mycrab_noout$CW))

pal <- c("#625a94", "#11c2b5")

ggplot() + 
  geom_ribbon(data = df_predict_noout, aes(x = CW, ymin = conf.low, 
                                           ymax = conf.high, fill = Sex), 
              alpha = 0.3) +
  geom_line(data = df_predict_noout, aes(x = CW, y = predicted, colour = Sex)) +
  geom_point(data = mycrab_noout, aes(x = CW, y = prop_eaten, color = Sex)) +
  labs(x = "Carapace width (mm)", y = "Proportion consumed", 
       size = "Clam\nDensity") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ylim(0,1) +
  theme_classic()

# SUPPLEMENTAL FIGURE ratio---------------------------------------------
prop_eaten_model_ratio <- glm(prop_eaten ~ stand_density + stand_claw_ratio + 
                                Sex + Sex:stand_claw_ratio, 
                              data = mycrab, 
                              weights = Density, 
                              family = binomial(link = logit))
summary(prop_eaten_model_ratio)

df_predict_ratio <- ggpredict(prop_eaten_model_ratio, 
                              terms = c("stand_claw_ratio[n=100]", "Sex")) %>% 
  rename(stand_claw_ratio = x, 
         Sex = group) %>%
  mutate(claw_ratio = stand_claw_ratio*(sd(mycrab$claw_ratio)) + 
           mean(mycrab$claw_ratio)) %>% 
  filter((Sex == "M" & claw_ratio > 0.206)| (Sex == "F" & claw_ratio < 0.274))

pal <- c("#625a94", "#11c2b5")

ggplot() + 
  geom_ribbon(data = df_predict_ratio, aes(x = claw_ratio, ymin = conf.low, 
                                           ymax = conf.high, fill = Sex), 
              alpha = 0.3) +
  geom_line(data = df_predict_ratio, aes(x = claw_ratio, y = predicted, 
                                         colour = Sex)) +
  geom_point(data = mycrab, aes(x = claw_ratio, y = prop_eaten, color = Sex)) +
  labs(x = "Ratio claw to carapace width (mm)", y = "Proportion consumed") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ylim(0,1) +
  theme_classic()
